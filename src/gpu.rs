//! GPU (wgpu compute-shader) implementation of the MM pairwise kernel.
//!
//! The kernel evaluates the same Coulomb + Lennard-Jones pair sums as the CPU
//! path in `mmpbsa::calc_mm`.  It requires a device with the
//! `Features::SHADER_F64` feature so the results stay consistent with the CPU
//! double-precision sums; adapters without f64 shader support are treated as
//! "GPU not supported" and the caller falls back to the CPU path.  The
//! screening `exp` is evaluated in plain arithmetic because driver f64
//! transcendentals are low-precision or broken on common stacks (see the
//! comment in the shader).

use std::sync::OnceLock;

use bytemuck::{Pod, Zeroable};
use ndarray::Array1;
use wgpu::util::{BufferInitDescriptor, DeviceExt};

/// Mirror of the WGSL `Atom` struct (40 bytes, align 8).
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuAtom {
    charge: f64,
    x: f64,
    y: f64,
    z: f64,
    type_id: u32,
    _pad: u32,
}

/// Mirror of the WGSL `Params` struct (32 bytes, align 8).
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuParams {
    n_a: u32,
    n_b: u32,
    use_screen: u32,
    n_types: u32,
    r_cutoff_sq: f64,
    lambda_d: f64,
}

const WGSL_SHADER: &str = r#"
// f64 types require the SHADER_F64 device feature; there is no WGSL
// `enable f64;` directive in naga (wgpu's WGSL front end).

struct Atom {
    charge: f64,
    x: f64,
    y: f64,
    z: f64,
    type_id: u32,
    _pad: u32,
}

struct Params {
    n_a: u32,
    n_b: u32,
    use_screen: u32,
    n_types: u32,
    r_cutoff_sq: f64,
    lambda_d: f64,
}

@group(0) @binding(0) var<storage> params: Params;
@group(0) @binding(1) var<storage> atoms_a: array<Atom>;
@group(0) @binding(2) var<storage> atoms_b: array<Atom>;
@group(0) @binding(3) var<storage> c6: array<f64>;
@group(0) @binding(4) var<storage> c12: array<f64>;
// Per-atom (electrostatic, van-der-Waals) sums: [0..n_a) receptor, [n_a..n_a+n_b) ligand.
@group(0) @binding(5) var<storage, read_write> out: array<vec2<f64>>;
// Bonded exclusion CSR (main_self only): excl_offsets[n + 1] slice ends into
// the sorted per-atom exclusion lists in excl_list.
@group(0) @binding(6) var<storage> excl_offsets: array<u32>;
@group(0) @binding(7) var<storage> excl_list: array<u32>;

// Driver f64 transcendentals (exp/exp2) are implemented with reduced
// precision or are outright wrong for some arguments on common GPU stacks,
// so exp for the screening factor is evaluated here in plain arithmetic:
// exp(x) = 2^(x*log2 e) = 2^k * exp(f*ln2) with k = floor(t), f = t - k,
// the last factor expanded in a Taylor series.  Only multiplies, divides
// and floor are used; 2^k is exact for k <= 0 (products of powers of two).
const INV_LOG2E: f64 = 1.4426950408889634;
const LN2: f64 = 0.6931471805599453;

// 2^k for -255 <= k <= 0, exact: multiplications by powers of two only.
// Bit i of k contributes 2^(2^i); the factors follow by exact squaring.
fn pow2_nonpos(k: i32) -> f64 {
    var p = f64(1.0);
    var factor = f64(2.0);
    let m = u32(-k);
    for (var i = 0u; i < 8u; i++) {
        if ((m >> i) & 1u) == 1u {
            p = p / factor;
        }
        factor = factor * factor;
    }
    return p;
}

// exp(x) for x <= 0 to ~1e-14 relative accuracy (x < -176 underflows to 0,
// like the libm exp the CPU path uses would for these magnitudes).
fn exp_neg(x: f64) -> f64 {
    if (x < f64(-176.0)) {
        return f64(0.0);
    }
    let t = x * INV_LOG2E;
    let k = i32(floor(t));
    let g = (t - floor(t)) * LN2;
    // e^g via the nested Taylor form 1 + g/1 (1 + g/2 (1 + ... (1 + g/17))).
    var h = f64(1.0);
    for (var n = 17; n >= 1; n--) {
        h = f64(1.0) + g * h / f64(n);
    }
    return pow2_nonpos(k) * h;
}

fn pair_energy(a: Atom, b: Atom) -> vec2<f64> {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    let r_sq = dx * dx + dy * dy + dz * dz;
    if (r_sq > params.r_cutoff_sq) {
        return vec2<f64>(0.0, 0.0);
    }
    let r = sqrt(r_sq) * 0.1;  // A -> nm
    let r_inv = 1.0 / r;
    var screen = f64(1.0);
    if (params.use_screen == 1u) {
        screen = exp_neg(-r / params.lambda_d);
    }
    let e_elec = a.charge * b.charge * r_inv * screen;
    let r6_inv = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
    let idx = a.type_id * params.n_types + b.type_id;
    let e_vdw = (c12[idx] * r6_inv - c6[idx]) * r6_inv;
    return vec2<f64>(e_elec, e_vdw);
}

// Dispatch A: one thread per receptor atom, summing over all ligand atoms.
@compute @workgroup_size(64)
fn main_a(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n_a) {
        return;
    }
    var acc = vec2<f64>(0.0, 0.0);
    for (var j = 0u; j < params.n_b; j += 1u) {
        acc = acc + pair_energy(atoms_a[i], atoms_b[j]);
    }
    out[i] = acc;
}

// Dispatch B: one thread per ligand atom, summing over all receptor atoms.
// pair_energy keeps the (receptor, ligand) argument order of dispatch A
// because the c6/c12 lookup is not guaranteed symmetric.
@compute @workgroup_size(64)
fn main_b(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n_b) {
        return;
    }
    var acc = vec2<f64>(0.0, 0.0);
    for (var j = 0u; j < params.n_a; j += 1u) {
        acc = acc + pair_energy(atoms_a[j], atoms_b[i]);
    }
    out[params.n_a + i] = acc;
}

// Dispatch S: one thread per atom, summing over all OTHER atoms of the same
// group (the zero-distance self pair would poison the sum with infinities).
// Bonded exclusions (1-2/1-3 and the explicit [ pairs ] handled elsewhere)
// arrive as a CSR list per atom, sorted ascending; the inner loop walks a
// cursor through it while j advances, so every excluded partner is skipped
// exactly once.  The thread's atom stays the first pair_energy argument, so
// the c6/c12 lookup is [own type, other type].
@compute @workgroup_size(64)
fn main_self(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let i = gid.x;
    if (i >= params.n_a) {
        return;
    }
    _ = lid;
    var acc = vec2<f64>(0.0, 0.0);
    var ex = excl_offsets[i];
    let ex_end = excl_offsets[i + 1u];
    for (var j = 0u; j < params.n_b; j += 1u) {
        if (j == i) {
            continue;
        }
        if (ex < ex_end && excl_list[ex] == j) {
            ex += 1u;
            continue;
        }
        acc = acc + pair_energy(atoms_a[i], atoms_b[j]);
    }
    out[i] = acc;
}
"#;

const WORKGROUP_SIZE: u32 = 64;

pub struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline_a: wgpu::ComputePipeline,
    pipeline_b: wgpu::ComputePipeline,
    pipeline_self: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    adapter_name: String,
}

static GPU_CONTEXT: OnceLock<Option<GpuContext>> = OnceLock::new();

/// Returns the shared GPU context, initializing it on first use.  `None` means
/// GPU is not supported on this machine (no adapter, or no f64 shader support);
/// the reason is printed once and the caller must fall back to the CPU path.
pub fn mm_context() -> Option<&'static GpuContext> {
    GPU_CONTEXT
        .get_or_init(|| match GpuContext::new() {
            Ok(ctx) => {
                println!(
                    "GPU enabled for MM calculations on: {}.",
                    ctx.adapter_name
                );
                Some(ctx)
            }
            Err(e) => {
                println!("GPU requested but not supported ({}); falling back to CPU.", e);
                None
            }
        })
        .as_ref()
}

/// PCI vendor ID of AMD/ATI.
const VENDOR_AMD: u32 = 0x1002;

/// Whether this adapter is an AMD GPU.  AMD drivers break on the s_mmpbsa
/// f64 compute pipelines (on a Radeon 780M Windows Vulkan stack the pipeline
/// compilation intermittently lost the device, the builtin f64 exp produced
/// wrong results, and multigrid solves wedged mid-run), so AMD adapters
/// always take the CPU path regardless of what the driver advertises.
fn is_amd_adapter(info: &wgpu::AdapterInfo) -> bool {
    (info.vendor & 0xffff) == VENDOR_AMD
        || info.name.contains("AMD")
        || info.name.contains("Radeon")
        || info.name.contains("Advanced Micro Devices")
}

impl GpuContext {
    fn new() -> Result<GpuContext, String> {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .map_err(|e| format!("no GPU adapter found: {}", e))?;
        let adapter_info = adapter.get_info();
        if is_amd_adapter(&adapter_info) {
            return Err(format!(
                "AMD GPU '{}' detected, the GPU path is disabled for AMD drivers",
                adapter_info.name
            ));
        }
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err(format!(
                "adapter '{}' lacks 64-bit float shader support",
                adapter.get_info().name
            ));
        }
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("s_mmpbsa MM device"),
            required_features: wgpu::Features::SHADER_F64,
            required_limits: wgpu::Limits::default(),
            memory_hints: Default::default(),
            trace: Default::default(),
        }))
        .map_err(|e| format!("cannot open GPU device: {}", e))?;

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("s_mmpbsa MM shader"),
            source: wgpu::ShaderSource::Wgsl(std::borrow::Cow::Borrowed(WGSL_SHADER)),
        });
        let bind_group_layout = Self::create_bind_group_layout(&device);
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("s_mmpbsa MM layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });
        let make_pipeline = |entry_point: &'static str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry_point),
                layout: Some(&layout),
                module: &module,
                entry_point: Some(entry_point),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let pipeline_a = make_pipeline("main_a");
        let pipeline_b = make_pipeline("main_b");
        let pipeline_self = make_pipeline("main_self");
        let adapter_name = format!("{} ({:?})", adapter_info.name, adapter_info.backend);
        Ok(GpuContext {
            device,
            queue,
            pipeline_a,
            pipeline_b,
            pipeline_self,
            bind_group_layout,
            adapter_name,
        })
    }

    fn create_bind_group_layout(device: &wgpu::Device) -> wgpu::BindGroupLayout {
        let storage_read = |binding: u32| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: true },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let storage_rw = |binding: u32| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: false },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("s_mmpbsa MM bind group layout"),
            entries: &[
                storage_read(0), // params
                storage_read(1), // receptor atoms
                storage_read(2), // ligand atoms
                storage_read(3), // c6 matrix
                storage_read(4), // c12 matrix
                storage_rw(5),   // per-atom (elec, vdw) sums
                storage_read(6), // exclusion CSR offsets (main_self)
                storage_read(7), // exclusion CSR lists (main_self)
            ],
        })
    }
}

/// A flat `(atom index, charge, type id, x, y, z)` record collected for one
/// molecule group, mirroring the tuples built by `mmpbsa::calc_mm`.
pub type AtomEntry = (usize, f64, usize, f64, f64, f64);

/// Runs the MM pairwise kernel on the GPU and returns the unscaled per-atom
/// electrostatic and van-der-Waals sums, in the same (unscaled) units the CPU
/// path produces before applying `elec_scale` / `vdw_scale`.
///
/// `c6` / `c12` are the Lennard-Jones matrices indexed by `[type_id, type_id]`.
pub fn calc_mm_gpu(
    ctx: &GpuContext,
    n_atoms: usize,
    rec_props: &[AtomEntry],
    lig_props: &[AtomEntry],
    c6: &ndarray::Array2<f64>,
    c12: &ndarray::Array2<f64>,
    r_cutoff_sq: f64,
    use_screen: bool,
    lambda_d: f64,
) -> Result<(Array1<f64>, Array1<f64>), String> {
    let (n_a, n_b) = (rec_props.len(), lig_props.len());
    if n_a == 0 || n_b == 0 || n_atoms == 0 {
        return Ok((Array1::zeros(n_atoms), Array1::zeros(n_atoms)));
    }
    let n_types = c6.dim().0;
    if c6.dim() != c12.dim() || c6.dim().1 != n_types {
        return Err("LJ matrices are not square or differ in shape".to_string());
    }

    let atoms_a: Vec<GpuAtom> = rec_props
        .iter()
        .map(|&(_, q, t, x, y, z)| GpuAtom {
            charge: q,
            x,
            y,
            z,
            type_id: t as u32,
            _pad: 0,
        })
        .collect();
    let atoms_b: Vec<GpuAtom> = lig_props
        .iter()
        .map(|&(_, q, t, x, y, z)| GpuAtom {
            charge: q,
            x,
            y,
            z,
            type_id: t as u32,
            _pad: 0,
        })
        .collect();
    let params = GpuParams {
        n_a: n_a as u32,
        n_b: n_b as u32,
        use_screen: use_screen as u32,
        n_types: n_types as u32,
        r_cutoff_sq,
        lambda_d,
    };
    let c6_flat = flatten_row_major(c6);
    let c12_flat = flatten_row_major(c12);

    // Output layout: [0..n_a) receptor sums, [n_a..n_a+n_b) ligand sums,
    // each entry an interleaved (elec, vdw) pair.
    let n_out = n_a + n_b;
    let params_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let atoms_a_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM atoms A"),
            contents: bytemuck::cast_slice(&atoms_a),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let atoms_b_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM atoms B"),
            contents: bytemuck::cast_slice(&atoms_b),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let c6_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM c6"),
            contents: bytemuck::cast_slice(&c6_flat),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let c12_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM c12"),
            contents: bytemuck::cast_slice(&c12_flat),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let out_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("MM out"),
        size: (n_out * std::mem::size_of::<[f64; 2]>()) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let read_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("MM readback"),
        size: out_buf.size(),
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("MM bind group"),
        layout: &ctx.bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: params_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: atoms_a_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: atoms_b_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: c6_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 4,
                resource: c12_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 5,
                resource: out_buf.as_entire_binding(),
            },
            // main_a / main_b never read the exclusion CSR; it is only bound
            // to satisfy the shared layout (main_self builds its own bind
            // group with the real buffers).
            wgpu::BindGroupEntry {
                binding: 6,
                resource: params_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 7,
                resource: params_buf.as_entire_binding(),
            },
        ],
    });

    let mut encoder = ctx
        .device
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("MM encoder"),
        });
    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("MM pass"),
            timestamp_writes: None,
        });
        pass.set_bind_group(0, &bind_group, &[]);
        let groups_a = (n_a as u32).div_ceil(WORKGROUP_SIZE);
        if groups_a > 0 {
            pass.set_pipeline(&ctx.pipeline_a);
            pass.dispatch_workgroups(groups_a, 1, 1);
        }
        let groups_b = (n_b as u32).div_ceil(WORKGROUP_SIZE);
        if groups_b > 0 {
            pass.set_pipeline(&ctx.pipeline_b);
            pass.dispatch_workgroups(groups_b, 1, 1);
        }
    }
    encoder.copy_buffer_to_buffer(&out_buf, 0, &read_buf, 0, out_buf.size());
    ctx.queue.submit(Some(encoder.finish()));
    let raw = map_read_f64(&ctx.device, &read_buf)?;

    // Scatter the group-local sums back onto their global atom indices.
    let mut de_elec = Array1::zeros(n_atoms);
    let mut de_vdw = Array1::zeros(n_atoms);
    for (local, &(gid, _, _, _, _, _)) in rec_props.iter().enumerate() {
        de_elec[gid] = raw[local][0];
        de_vdw[gid] = raw[local][1];
    }
    for (local, &(gid, _, _, _, _, _)) in lig_props.iter().enumerate() {
        de_elec[gid] = raw[n_a + local][0];
        de_vdw[gid] = raw[n_a + local][1];
    }
    Ok((de_elec, de_vdw))
}

/// Runs the MM pairwise kernel on the GPU for the self-interaction of one
/// atom group: every atom against all other atoms of the same set, the self
/// pair and the bonded exclusions (1-2/1-3, plus the explicit `[ pairs ]`
/// whose energy the caller adds separately) skipped.  Returns the unscaled
/// per-atom electrostatic and van-der-Waals sums, in the same (unscaled)
/// units the CPU path produces before applying `elec_scale` / `vdw_scale`.
///
/// `excl_offsets` / `excl_list` hold the per-atom exclusion lists of the
/// group in CSR form with **group-local, ascending** indices; `c6` / `c12`
/// are indexed by `[type_id, type_id]` with the summing atom first, matching
/// the WGSL `pair_energy` lookup.
pub fn calc_mm_self_gpu(
    ctx: &GpuContext,
    n_atoms: usize,
    self_props: &[AtomEntry],
    c6: &ndarray::Array2<f64>,
    c12: &ndarray::Array2<f64>,
    r_cutoff_sq: f64,
    use_screen: bool,
    lambda_d: f64,
    excl_offsets: &[u32],
    excl_list: &[u32],
) -> Result<(Array1<f64>, Array1<f64>), String> {
    let n = self_props.len();
    if n == 0 || n_atoms == 0 {
        return Ok((Array1::zeros(n_atoms), Array1::zeros(n_atoms)));
    }
    if excl_offsets.len() != n + 1 {
        return Err("exclusion CSR offsets do not match the group size".to_string());
    }
    let n_types = c6.dim().0;
    if c6.dim() != c12.dim() || c6.dim().1 != n_types {
        return Err("LJ matrices are not square or differ in shape".to_string());
    }

    let atoms: Vec<GpuAtom> = self_props
        .iter()
        .map(|&(_, q, t, x, y, z)| GpuAtom {
            charge: q,
            x,
            y,
            z,
            type_id: t as u32,
            _pad: 0,
        })
        .collect();
    let params = GpuParams {
        n_a: n as u32,
        n_b: n as u32,
        use_screen: use_screen as u32,
        n_types: n_types as u32,
        r_cutoff_sq,
        lambda_d,
    };
    let c6_flat = flatten_row_major(c6);
    let c12_flat = flatten_row_major(c12);

    let params_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM self params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let atoms_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM self atoms"),
            contents: bytemuck::cast_slice(&atoms),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let c6_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM self c6"),
            contents: bytemuck::cast_slice(&c6_flat),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let c12_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM self c12"),
            contents: bytemuck::cast_slice(&c12_flat),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let excl_offsets_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM self excl offsets"),
            contents: bytemuck::cast_slice(excl_offsets),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let excl_list_buf = ctx
        .device
        .create_buffer_init(&BufferInitDescriptor {
            label: Some("MM self excl list"),
            contents: bytemuck::cast_slice(excl_list),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
    let out_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("MM self out"),
        size: (n * std::mem::size_of::<[f64; 2]>()) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let read_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("MM self readback"),
        size: out_buf.size(),
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    // The same read-only atom buffer feeds both the "a" and "b" slots of the
    // shader.
    let bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("MM self bind group"),
        layout: &ctx.bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: params_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: atoms_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: atoms_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: c6_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 4,
                resource: c12_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 5,
                resource: out_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 6,
                resource: excl_offsets_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 7,
                resource: excl_list_buf.as_entire_binding(),
            },
        ],
    });

    let mut encoder = ctx
        .device
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("MM self encoder"),
        });
    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("MM self pass"),
            timestamp_writes: None,
        });
        pass.set_bind_group(0, &bind_group, &[]);
        let groups = (n as u32).div_ceil(WORKGROUP_SIZE);
        if groups > 0 {
            pass.set_pipeline(&ctx.pipeline_self);
            pass.dispatch_workgroups(groups, 1, 1);
        }
    }
    encoder.copy_buffer_to_buffer(&out_buf, 0, &read_buf, 0, out_buf.size());
    ctx.queue.submit(Some(encoder.finish()));
    let raw = map_read_f64(&ctx.device, &read_buf)?;

    // Scatter the group-local sums back onto their global atom indices.
    let mut de_elec = Array1::zeros(n_atoms);
    let mut de_vdw = Array1::zeros(n_atoms);
    for (local, &(gid, _, _, _, _, _)) in self_props.iter().enumerate() {
        de_elec[gid] = raw[local][0];
        de_vdw[gid] = raw[local][1];
    }
    Ok((de_elec, de_vdw))
}

fn map_read_f64(device: &wgpu::Device, buffer: &wgpu::Buffer) -> Result<Vec<[f64; 2]>, String> {
    let (tx, rx) = std::sync::mpsc::channel();
    buffer
        .slice(..)
        .map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
    device
        .poll(wgpu::PollType::Wait)
        .map_err(|e| format!("GPU poll failed: {}", e))?;
    rx.recv()
        .map_err(|_| "GPU map callback was dropped".to_string())?
        .map_err(|e| format!("GPU output mapping failed: {}", e))?;
    let view = buffer.slice(..).get_mapped_range();
    let mut data: Vec<[f64; 2]> = Vec::with_capacity(view.len() / 16);
    for chunk in view.chunks_exact(16) {
        data.push([
            f64::from_ne_bytes(chunk[0..8].try_into().unwrap()),
            f64::from_ne_bytes(chunk[8..16].try_into().unwrap()),
        ]);
    }
    drop(view);
    buffer.unmap();
    Ok(data)
}

/// Row-major flattening that works for any ndarray memory layout.
fn flatten_row_major(a: &ndarray::Array2<f64>) -> Vec<f64> {
    let mut out = Vec::with_capacity(a.len());
    for row in a.rows() {
        for &v in row {
            out.push(v);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngExt;

    /// Reference CPU implementation of the unscaled MM pair sums, mirroring
    /// `mmpbsa::calc_mm_cpu`.
    fn reference(
        n_atoms: usize,
        rec: &[AtomEntry],
        lig: &[AtomEntry],
        c6: &ndarray::Array2<f64>,
        c12: &ndarray::Array2<f64>,
        r_cutoff_sq: f64,
        use_screen: bool,
        lambda_d: f64,
    ) -> (Array1<f64>, Array1<f64>) {
        let mut elec = Array1::zeros(n_atoms);
        let mut vdw = Array1::zeros(n_atoms);
        for &(i, qi, ci, xi, yi, zi) in rec {
            for &(j, qj, cj, xj, yj, zj) in lig {
                let (dx, dy, dz) = (xi - xj, yi - yj, zi - zj);
                let r_sq = dx * dx + dy * dy + dz * dz;
                if r_sq <= r_cutoff_sq {
                    let r = r_sq.sqrt() * 0.1;
                    let r_inv = 1.0 / r;
                    let screen = if use_screen { (-r / lambda_d).exp() } else { 1.0 };
                    let e_elec = qi * qj * r_inv * screen;
                    // same left-to-right multiply chain as the WGSL kernel
                    let r6_inv = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                    let e_vdw = (c12[[ci, cj]] * r6_inv - c6[[ci, cj]]) * r6_inv;
                    elec[i] += e_elec;
                    elec[j] += e_elec;
                    vdw[i] += e_vdw;
                    vdw[j] += e_vdw;
                }
            }
        }
        (elec, vdw)
    }

    /// Atoms sit on a jittered grid so pair distances stay physical (>~1 A);
    /// with uniformly random coordinates some pairs land so close that their
    /// huge LJ terms cancel in the sum and any 1-ulp difference in the r^6
    /// evaluation shows up as a visible difference in the totals.
    fn random_system(n_rec: usize, n_lig: usize, n_types: usize) -> (
        Vec<AtomEntry>, Vec<AtomEntry>, ndarray::Array2<f64>, ndarray::Array2<f64>,
    ) {
        let mut rng = rand::rng();
        let mut place = |ids: std::ops::Range<usize>, offset: f64| {
            ids.map(|k| {
                let cell = k as f64;
                (
                    k, rng.random_range(-1.0..1.0), rng.random_range(0..n_types),
                    offset + 2.0 * (cell - 5.0 * (cell / 5.0).floor()) + rng.random_range(-0.3..0.3),
                    2.0 * (cell / 5.0).floor() + rng.random_range(-0.3..0.3),
                    offset + rng.random_range(-0.3..0.3),
                )
            })
            .collect::<Vec<AtomEntry>>()
        };
        let rec = place(0..n_rec, 0.0);
        let lig = place(n_rec..n_rec + n_lig, 1.0);
        let c6 = ndarray::Array2::from_shape_fn((n_types, n_types), |_| rng.random_range(0.0..0.01));
        let c12 = ndarray::Array2::from_shape_fn((n_types, n_types), |_| rng.random_range(0.0..1e-5));
        (rec, lig, c6, c12)
    }

    fn assert_close(a: &Array1<f64>, b: &Array1<f64>, what: &str) {
        assert_eq!(a.len(), b.len(), "{}: length mismatch", what);
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            let scale = x.abs().max(y.abs()).max(1e-12);
            assert!(
                (x - y).abs() / scale < 1e-12,
                "{} mismatch at atom {}: gpu={} cpu={} (rel diff {:.3e})",
                what, i, x, y, (x - y).abs() / scale
            );
        }
    }

    /// GPU vs CPU agreement.  Skipped (not failed) on machines where no
    /// f64-capable GPU is available, so CPU-only environments still pass.
    #[test]
    fn gpu_matches_cpu_on_random_systems() {
        let Some(ctx) = mm_context() else {
            println!("No f64-capable GPU found; skipping the GPU/CPU comparison test.");
            return;
        };
        for n_types in [1usize, 4] {
            let (rec, lig, c6, c12) = random_system(37, 11, n_types);
            for (r_cutoff_sq, use_screen) in [(f64::INFINITY, true), (100.0, false), (2500.0, true)] {
                let (gpu_e, gpu_v) = calc_mm_gpu(
                    ctx, 48, &rec, &lig, &c6, &c12, r_cutoff_sq, use_screen, 0.8,
                ).expect("GPU MM run failed");
                let (cpu_e, cpu_v) = reference(48, &rec, &lig, &c6, &c12, r_cutoff_sq, use_screen, 0.8);
                assert_close(&gpu_e, &cpu_e, "elec");
                assert_close(&gpu_v, &cpu_v, "vdw");
            }
        }
    }

    /// Atoms outside both groups must stay zero on the GPU path, and empty
    /// groups must produce all-zero results.
    /// The self-interaction kernel must skip the bonded exclusions given as
    /// a CSR list and otherwise match a per-atom CPU reference.
    #[test]
    fn gpu_self_skips_exclusions_and_matches_cpu() {
        let Some(ctx) = mm_context() else {
            println!("No f64-capable GPU found; skipping the GPU self test.");
            return;
        };
        let (rec, _lig, c6, c12) = random_system(23, 5, 4);
        let n_atoms = 23;

        // Exclude a handful of pairs in group-local indices.
        let excluded: Vec<(usize, usize)> = vec![
            (0, 1), (0, 5), (2, 9), (7, 11), (4, 17), (13, 21), (6, 18),
        ];
        let mut excl_by_atom: Vec<Vec<u32>> = vec![Vec::new(); n_atoms];
        let mut excluded_set = std::collections::HashSet::new();
        for &(a, b) in &excluded {
            excl_by_atom[a].push(b as u32);
            excl_by_atom[b].push(a as u32);
            // Both directions, so the reference loop skips whichever order it
            // meets the pair in.
            excluded_set.insert((a, b));
            excluded_set.insert((b, a));
        }
        for list in &mut excl_by_atom {
            list.sort_unstable();
        }
        let mut excl_offsets = vec![0u32; n_atoms + 1];
        let mut excl_list: Vec<u32> = Vec::new();
        for (li, list) in excl_by_atom.iter().enumerate() {
            excl_list.extend_from_slice(list);
            excl_offsets[li + 1] = excl_list.len() as u32;
        }

        for (r_cutoff_sq, use_screen) in [(f64::INFINITY, true), (100.0, false)] {
            let (gpu_e, gpu_v) = calc_mm_self_gpu(
                ctx, n_atoms, &rec, &c6, &c12, r_cutoff_sq, use_screen, 0.8,
                &excl_offsets, &excl_list,
            ).expect("GPU self MM run failed");

            let mut cpu_e = Array1::zeros(n_atoms);
            let mut cpu_v = Array1::zeros(n_atoms);
            for &(i, qi, ci, xi, yi, zi) in &rec {
                for &(j, qj, cj, xj, yj, zj) in &rec {
                    if j == i {
                        continue;
                    }
                    if excluded_set.contains(&(i, j)) {
                        continue;
                    }
                    let (dx, dy, dz) = (xi - xj, yi - yj, zi - zj);
                    let r_sq = dx * dx + dy * dy + dz * dz;
                    if r_sq <= r_cutoff_sq {
                        let r = r_sq.sqrt() * 0.1;
                        let r_inv = 1.0 / r;
                        let screen = if use_screen { (-r / 0.8).exp() } else { 1.0 };
                        let e_elec = qi * qj * r_inv * screen;
                        let r6_inv = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                        let e_vdw = (c12[[ci, cj]] * r6_inv - c6[[ci, cj]]) * r6_inv;
                        cpu_e[i] += e_elec;
                        cpu_v[i] += e_vdw;
                    }
                }
            }
            assert_close(&gpu_e, &cpu_e, "self+excl elec");
            assert_close(&gpu_v, &cpu_v, "self+excl vdw");
        }
    }

    /// Minimal deterministic probe: three atoms, the pair (0, 1) excluded.
    #[test]
    fn gpu_self_minimal_exclusion_probe() {
        let Some(ctx) = mm_context() else {
            println!("No f64-capable GPU found; skipping the minimal probe.");
            return;
        };
        let rec: Vec<AtomEntry> = vec![
            (0, 0.4, 0, 0.0, 0.0, 0.0),
            (1, 0.4, 0, 2.5, 0.0, 0.0),
            (2, 0.4, 0, 0.0, 3.0, 0.0),
        ];
        let c6 = ndarray::Array2::from_elem((1, 1), 0.001);
        let c12 = ndarray::Array2::from_elem((1, 1), 1e-6);

        // Exclude the pair (0, 1): atom 1 keeps only its pair with atom 2.
        let excl_offsets: Vec<u32> = vec![0, 1, 2, 2];
        let excl_list: Vec<u32> = vec![1, 0];

        let (gpu_e, _) = calc_mm_self_gpu(
            ctx, 3, &rec, &c6, &c12, f64::INFINITY, false, 0.8,
            &excl_offsets, &excl_list,
        ).expect("GPU self MM run failed");

        let pair = |ri: &[f64; 3], rj: &[f64; 3]| -> f64 {
            let (dx, dy, dz) = (ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]);
            let r = (dx * dx + dy * dy + dz * dz).sqrt() * 0.1;
            0.4 * 0.4 / r
        };
        let c = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [0.0, 3.0, 0.0]];
        // Atom 0: only pair with atom 2.  Atom 1: only pair with atom 2.
        // Atom 2: pairs with 0 and 1.
        let want = [pair(&c[0], &c[2]), pair(&c[1], &c[2]), pair(&c[0], &c[2]) + pair(&c[1], &c[2])];
        for i in 0..3 {
            println!("[mini] atom {i}: gpu={} want={}", gpu_e[i], want[i]);
            assert!((gpu_e[i] - want[i]).abs() < 1e-12, "atom {i}: gpu {} vs want {}", gpu_e[i], want[i]);
        }
    }
}
