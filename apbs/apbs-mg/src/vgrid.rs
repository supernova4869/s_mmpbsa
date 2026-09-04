// APBS vgrid.rs - Grid data structure and I/O
// Port of src/mg/vgrid.h / vgrid.c

use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use apbs_generic::error::{ApbsError, ApbsResult};

const APBS_PACKAGE_STRING: &str = concat!("APBS Rust ", env!("CARGO_PKG_VERSION"));

/// Number of decimal places for comparisons
const _VGRID_DIGITS: f64 = 1.0e6;

/// Uniform Cartesian grid for storing potentials, dielectric maps, etc.
#[derive(Debug, Clone)]
pub struct Vgrid {
    /// Grid points in x
    pub nx: usize,
    /// Grid points in y
    pub ny: usize,
    /// Grid points in z
    pub nz: usize,
    /// Grid spacing in x (A)
    pub hx: f64,
    /// Grid spacing in y (A)
    pub hy: f64,
    /// Grid spacing in z (A)
    pub hzed: f64,
    /// Lower grid corner x coordinate
    pub xmin: f64,
    /// Lower grid corner y coordinate
    pub ymin: f64,
    /// Lower grid corner z coordinate
    pub zmin: f64,
    /// Upper grid corner x coordinate
    pub xmax: f64,
    /// Upper grid corner y coordinate
    pub ymax: f64,
    /// Upper grid corner z coordinate
    pub zmax: f64,
    /// nx*ny*nz array of data
    pub data: Vec<f64>,
    /// Data was read from file
    pub readdata: bool,
    /// Data was included at construction
    pub ctordata: bool,
}

impl Vgrid {
    /// Create a new grid
    pub fn new(
        nx: usize, ny: usize, nz: usize,
        hx: f64, hy: f64, hzed: f64,
        xmin: f64, ymin: f64, zmin: f64,
    ) -> Self {
        let xmax = xmin + (nx - 1) as f64 * hx;
        let ymax = ymin + (ny - 1) as f64 * hy;
        let zmax = zmin + (nz - 1) as f64 * hzed;

        Self {
            nx, ny, nz,
            hx, hy, hzed,
            xmin, ymin, zmin,
            xmax, ymax, zmax,
            data: vec![0.0; nx * ny * nz],
            readdata: false,
            ctordata: false,
        }
    }

    /// Create a grid with provided data
    pub fn with_data(
        nx: usize, ny: usize, nz: usize,
        hx: f64, hy: f64, hzed: f64,
        xmin: f64, ymin: f64, zmin: f64,
        data: Vec<f64>,
    ) -> Self {
        assert_eq!(data.len(), nx * ny * nz);
        let xmax = xmin + (nx - 1) as f64 * hx;
        let ymax = ymin + (ny - 1) as f64 * hy;
        let zmax = zmin + (nz - 1) as f64 * hzed;

        Self {
            nx, ny, nz,
            hx, hy, hzed,
            xmin, ymin, zmin,
            xmax, ymax, zmax,
            data,
            readdata: false,
            ctordata: true,
        }
    }

    /// 3D index to flat index
    #[inline]
    pub fn ijk(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.nx * self.ny + j * self.nx + i
    }

    /// Interpolate value at a point using trilinear interpolation
    pub fn value(&self, pt: [f64; 3]) -> ApbsResult<f64> {
        // Check bounds
        if pt[0] < self.xmin - self.hx || pt[0] > self.xmax + self.hx ||
           pt[1] < self.ymin - self.hy || pt[1] > self.ymax + self.hy ||
           pt[2] < self.zmin - self.hzed || pt[2] > self.zmax + self.hzed {
            return Err(ApbsError::Grid(format!(
                "Point [{}, {}, {}] is outside grid bounds",
                pt[0], pt[1], pt[2]
            )));
        }

        // Compute fractional indices
        let fx = (pt[0] - self.xmin) / self.hx;
        let fy = (pt[1] - self.ymin) / self.hy;
        let fz = (pt[2] - self.zmin) / self.hzed;

        let ix = fx.floor() as i32;
        let iy = fy.floor() as i32;
        let iz = fz.floor() as i32;

        let dx = fx - ix as f64;
        let dy = fy - iy as f64;
        let dz = fz - iz as f64;

        // Clamp to grid
        let ix = ix.max(0).min(self.nx as i32 - 2) as usize;
        let iy = iy.max(0).min(self.ny as i32 - 2) as usize;
        let iz = iz.max(0).min(self.nz as i32 - 2) as usize;

        // Trilinear interpolation
        let c000 = self.data[self.ijk(ix, iy, iz)];
        let c100 = self.data[self.ijk(ix + 1, iy, iz)];
        let c010 = self.data[self.ijk(ix, iy + 1, iz)];
        let c110 = self.data[self.ijk(ix + 1, iy + 1, iz)];
        let c001 = self.data[self.ijk(ix, iy, iz + 1)];
        let c101 = self.data[self.ijk(ix + 1, iy, iz + 1)];
        let c011 = self.data[self.ijk(ix, iy + 1, iz + 1)];
        let c111 = self.data[self.ijk(ix + 1, iy + 1, iz + 1)];

        let c00 = c000 * (1.0 - dx) + c100 * dx;
        let c10 = c010 * (1.0 - dx) + c110 * dx;
        let c01 = c001 * (1.0 - dx) + c101 * dx;
        let c11 = c011 * (1.0 - dx) + c111 * dx;

        let c0 = c00 * (1.0 - dy) + c10 * dy;
        let c1 = c01 * (1.0 - dy) + c11 * dy;

        Ok(c0 * (1.0 - dz) + c1 * dz)
    }

    /// Compute gradient at a point using central differences
    pub fn gradient(&self, pt: [f64; 3]) -> ApbsResult<[f64; 3]> {
        let eps = 1.0e-6;
        let mut grad = [0.0; 3];

        for d in 0..3 {
            let mut pt_plus = pt;
            let mut pt_minus = pt;
            pt_plus[d] += eps;
            pt_minus[d] -= eps;

            let v_plus = self.value(pt_plus)?;
            let v_minus = self.value(pt_minus)?;
            grad[d] = (v_plus - v_minus) / (2.0 * eps);
        }

        Ok(grad)
    }

    /// Compute L2 norm of grid data
    pub fn norm_l2(&self) -> f64 {
        let mut sum = 0.0;
        for &v in &self.data {
            sum += v * v;
        }
        sum.sqrt()
    }

    /// Compute L1 norm of grid data
    pub fn norm_l1(&self) -> f64 {
        self.data.iter().map(|v| v.abs()).sum()
    }

    /// Compute L-inf norm of grid data
    pub fn norm_linf(&self) -> f64 {
        self.data.iter().map(|v| v.abs()).fold(0.0f64, f64::max)
    }

    /// Integrate grid data (sum * cell volume)
    pub fn integrate(&self) -> f64 {
        let vol = self.hx * self.hy * self.hzed;
        self.data.iter().sum::<f64>() * vol
    }

    /// Read OpenDX format file
    /// Port of Vgrid_readDX from vgrid.c
    pub fn read_dx(&mut self, filename: &str) -> ApbsResult<()> {
        let file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let reader = BufReader::new(file);
        let mut data = Vec::new();
        let mut header_done = false;
        let mut delta_count = 0;

        for line in reader.lines() {
            let line = line.map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            if trimmed.starts_with("object") && !header_done {
                // Parse: object 1 class gridpositions counts NX NY NZ
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() >= 8 && fields[1] == "1" && fields[2] == "class" && fields[3] == "gridpositions" {
                    self.nx = fields[5].parse().unwrap_or(self.nx);
                    self.ny = fields[6].parse().unwrap_or(self.ny);
                    self.nz = fields[7].parse().unwrap_or(self.nz);
                }
                continue;
            }

            if trimmed.starts_with("origin") && !header_done {
                // Parse: origin X Y Z
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() >= 4 {
                    self.xmin = fields[1].parse().unwrap_or(self.xmin);
                    self.ymin = fields[2].parse().unwrap_or(self.ymin);
                    self.zmin = fields[3].parse().unwrap_or(self.zmin);
                }
                continue;
            }

            if trimmed.starts_with("delta") && !header_done {
                // Parse: delta D1 D2 D3
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() >= 4 {
                    let d1: f64 = fields[1].parse().unwrap_or(0.0);
                    let d2: f64 = fields[2].parse().unwrap_or(0.0);
                    let d3: f64 = fields[3].parse().unwrap_or(0.0);
                    match delta_count {
                        0 => { self.hx = d1; }
                        1 => { self.hy = d2; }
                        2 => { self.hzed = d3; }
                        _ => {}
                    }
                    delta_count += 1;
                }
                continue;
            }

            // After header is done, collect data values
            if trimmed.starts_with("object") && header_done {
                // Another object line (e.g., field definition) - skip
                continue;
            }
            if trimmed.contains("data follows") || trimmed.contains("binary") {
                header_done = true;
                continue;
            }

            // If we haven't seen gridpositions yet, skip
            if self.nx == 0 || self.ny == 0 || self.nz == 0 {
                continue;
            }

            // Parse data values (may be multiple per line)
            if header_done || delta_count >= 3 {
                for token in trimmed.split_whitespace() {
                    if let Ok(v) = token.parse::<f64>() {
                        data.push(v);
                        if data.len() >= self.nx * self.ny * self.nz {
                            break;
                        }
                    }
                }
                if data.len() >= self.nx * self.ny * self.nz {
                    break;
                }
            }
        }

        // The data from readDX is in column-major order (i, j, k)
        // Convert to row-major: u = k*nx*ny + j*nx + i
        let len = self.nx * self.ny * self.nz;
        if data.len() == len {
            self.data = vec![0.0; len];
            let mut incr = 0;
            for i in 0..self.nx {
                for j in 0..self.ny {
                    for k in 0..self.nz {
                        if incr < data.len() {
                            let u = k * self.nx * self.ny + j * self.nx + i;
                            self.data[u] = data[incr];
                            incr += 1;
                        }
                    }
                }
            }
        } else if !data.is_empty() {
            // Data might already be in row-major or different format
            self.data = data;
        }

        self.xmax = self.xmin + (self.nx - 1) as f64 * self.hx;
        self.ymax = self.ymin + (self.ny - 1) as f64 * self.hy;
        self.zmax = self.zmin + (self.nz - 1) as f64 * self.hzed;
        self.readdata = true;

        Ok(())
    }

    /// Write OpenDX format file
    pub fn write_dx(&self, filename: &str, title: &str) -> ApbsResult<()> {
        let mut file = File::create(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        self.write_dx_text(&mut file, title)
    }

    /// Write gzipped DX file
    pub fn write_gz(&self, filename: &str, title: &str) -> ApbsResult<()> {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let file = File::create(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let mut encoder = GzEncoder::new(file, Compression::default());

        self.write_dx_text(&mut encoder, title)?;
        encoder.finish().map_err(|e| ApbsError::Io(e.to_string()))?;
        Ok(())
    }

    fn write_dx_text<W: Write>(&self, writer: &mut W, title: &str) -> ApbsResult<()> {
        writeln!(writer, "# Data from {}", APBS_PACKAGE_STRING).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "# ").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "# {}", title).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "# ").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "object 1 class gridpositions counts {} {} {}", self.nx, self.ny, self.nz)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "origin {:12.6e} {:12.6e} {:12.6e}", self.xmin, self.ymin, self.zmin)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "delta {:12.6e} {:12.6e} {:12.6e}", self.hx, 0.0, 0.0)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "delta {:12.6e} {:12.6e} {:12.6e}", 0.0, self.hy, 0.0)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "delta {:12.6e} {:12.6e} {:12.6e}", 0.0, 0.0, self.hzed)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "object 2 class gridconnections counts {} {} {}", self.nx, self.ny, self.nz)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "object 3 class array type double rank 0 items {} data follows", self.data.len())
            .map_err(|e| ApbsError::Io(e.to_string()))?;

        let mut icol = 0usize;
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    let idx = self.ijk(i, j, k);
                    write!(writer, "{:12.6e} ", self.data[idx])
                        .map_err(|e| ApbsError::Io(e.to_string()))?;
                    icol += 1;
                    if icol == 3 {
                        icol = 0;
                        writeln!(writer).map_err(|e| ApbsError::Io(e.to_string()))?;
                    }
                }
            }
        }
        if icol != 0 {
            writeln!(writer).map_err(|e| ApbsError::Io(e.to_string()))?;
        }

        writeln!(writer, "attribute \"dep\" string \"positions\"").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "object \"regular positions regular connections\" class field")
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "component \"positions\" value 1").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "component \"connections\" value 2").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(writer, "component \"data\" value 3").map_err(|e| ApbsError::Io(e.to_string()))?;

        Ok(())
    }

    /// Write APBS/C-compatible flat scalar format.
    pub fn write_flat_values(&self, filename: &str, title: &str, count: usize) -> ApbsResult<()> {
        let mut file = File::create(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        writeln!(file, "# Data from {}", APBS_PACKAGE_STRING).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "# ").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "# {}", title).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "# ").map_err(|e| ApbsError::Io(e.to_string()))?;
        for value in self.data.iter().take(count.min(self.data.len())) {
            writeln!(file, "{:12.6e}", value).map_err(|e| ApbsError::Io(e.to_string()))?;
        }
        Ok(())
    }

    /// Write flat format file
    pub fn write_flat(&self, filename: &str) -> ApbsResult<()> {
        let mut file = File::create(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let idx = self.ijk(i, j, k);
                    writeln!(file, "{} {} {} {:.6e}",
                        self.xmin + i as f64 * self.hx,
                        self.ymin + j as f64 * self.hy,
                        self.zmin + k as f64 * self.hzed,
                        self.data[idx]
                    ).map_err(|e| ApbsError::Io(e.to_string()))?;
                }
            }
        }

        Ok(())
    }

    /// Read binary OpenDX format file.
    /// Port of Vgrid_readDXBIN from vgrid.c line 810.
    pub fn read_dxbin(&mut self, filename: &str) -> ApbsResult<()> {
        use std::io::Read;

        let file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        // Read text header
        let mut buf_reader = BufReader::new(&file);
        let mut line = String::new();

        // Skip comment lines
        loop {
            line.clear();
            buf_reader.read_line(&mut line)
                .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
            if !line.starts_with('#') && !line.trim().is_empty() {
                break;
            }
        }

        // Parse: object 1 class gridpositions counts NX NY NZ
        {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() >= 7 && fields[1] == "1" && fields[2] == "class" && fields[3] == "gridpositions" {
                self.nx = fields[5].parse().map_err(|_| ApbsError::Io(format!("Failed to parse nx from: {}", line)))?;
                self.ny = fields[6].parse().map_err(|_| ApbsError::Io(format!("Failed to parse ny from: {}", line)))?;
                // nz might be on same line or next token
                if fields.len() > 7 {
                    self.nz = fields[7].parse().map_err(|_| ApbsError::Io(format!("Failed to parse nz from: {}", line)))?;
                }
            } else {
                return Err(ApbsError::Io(format!("Failed to read dimensions from: {}", line)));
            }
        }

        // If nz wasn't in the counts line, read next line
        // (typically the counts line has all three, but handle edge case)

        // Read origin line
        line.clear();
        buf_reader.read_line(&mut line)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        {
            let trimmed = line.trim();
            if trimmed.starts_with("origin") {
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() >= 4 {
                    self.xmin = fields[1].parse().map_err(|_| ApbsError::Io(format!("Failed to parse xmin")))?;
                    self.ymin = fields[2].parse().map_err(|_| ApbsError::Io(format!("Failed to parse ymin")))?;
                    self.zmin = fields[3].parse().map_err(|_| ApbsError::Io(format!("Failed to parse zmin")))?;
                }
            }
        }

        // Read 3 delta lines
        for _ in 0..3 {
            line.clear();
            buf_reader.read_line(&mut line)
                .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
            let trimmed = line.trim();
            if trimmed.starts_with("delta") {
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() >= 4 {
                    let d1: f64 = fields[1].parse().unwrap_or(0.0);
                    let d2: f64 = fields[2].parse().unwrap_or(0.0);
                    let d3: f64 = fields[3].parse().unwrap_or(0.0);
                    self.hx += d1;
                    self.hy += d2;
                    self.hzed += d3;
                }
            }
        }

        // Skip to "binary data follows" line
        loop {
            line.clear();
            buf_reader.read_line(&mut line)
                .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
            if line.contains("binary") || line.contains("data follows") {
                break;
            }
            if line.trim().is_empty() {
                continue;
            }
        }

        // Now switch to raw binary reading from the original file handle
        drop(buf_reader);

        // Re-open for binary reading from current position
        // We need to track how many bytes the text header consumed
        // Actually, we need to reopen since BufReader consumed the buffer
        let file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        // We need to skip past the text header to get to the binary data
        // Re-read and count bytes consumed
        let mut header_bytes = 0usize;
        let mut header_lines = 0;
        {
            let mut reader = BufReader::new(&file);
            loop {
                let mut l = String::new();
                let n = reader.read_line(&mut l)
                    .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
                header_bytes += n;
                header_lines += 1;
                if l.contains("data follows") {
                    break;
                }
                if header_lines > 100 {
                    return Err(ApbsError::Io(format!("Header too long in {}", filename)));
                }
            }
        }

        // Re-open and seek past header
        let mut file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        use std::io::Seek;
        file.seek(std::io::SeekFrom::Start(header_bytes as u64))
            .map_err(|e| ApbsError::Io(format!("Seek error: {}", e)))?;

        let tot = self.nx * self.ny * self.nz;
        self.data = vec![0.0; tot];

        // Read binary doubles in column-major order (i, j, k) and store at u = k*nx*ny + j*nx + i
        let mut buf = [0u8; 8];
        let mut counter = 0;
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    file.read_exact(&mut buf)
                        .map_err(|e| ApbsError::Io(format!("Read error at {}: {}", filename, e)))?;
                    let val = f64::from_le_bytes(buf);
                    let u = k * self.nx * self.ny + j * self.nx + i;
                    self.data[u] = val;
                    counter += 1;
                }
            }
        }

        if counter != tot {
            return Err(ApbsError::Io(format!(
                "Read {} doubles, expected {}", counter, tot
            )));
        }

        // Compute grid maxima
        self.xmax = self.xmin + (self.nx - 1) as f64 * self.hx;
        self.ymax = self.ymin + (self.ny - 1) as f64 * self.hy;
        self.zmax = self.zmin + (self.nz - 1) as f64 * self.hzed;
        self.readdata = true;

        Ok(())
    }

    /// Write binary OpenDX format file.
    /// Port of Vgrid_writeDXBIN from vgrid.c line 1458.
    pub fn write_dxbin(&self, filename: &str, title: &str) -> ApbsResult<()> {
        use std::io::Write;

        let mut file = File::create(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        // Write text header
        writeln!(file, "# Data from {}", APBS_PACKAGE_STRING).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "# ").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "# {}", title).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "# ").map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "object 1 class gridpositions counts {} {} {}", self.nx, self.ny, self.nz)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "origin {:.6e} {:.6e} {:.6e}", self.xmin, self.ymin, self.zmin)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "delta {:.6e} 0 0", self.hx)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "delta 0 {:.6e} 0", self.hy)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "delta 0 0 {:.6e}", self.hzed)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "object 2 class gridconnections counts {} {} {}", self.nx, self.ny, self.nz)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "object 3 class array type double rank 0 items {} binary data follows",
            self.data.len())
            .map_err(|e| ApbsError::Io(e.to_string()))?;

        // Write binary data in column-major order (i, j, k) -> u = k*nx*ny + j*nx + i
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    let u = k * self.nx * self.ny + j * self.nx + i;
                    let bytes = self.data[u].to_le_bytes();
                    file.write_all(&bytes)
                        .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
                }
            }
        }

        // Write field definition
        writeln!(file).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "attribute \"dep\" string \"positions\"")
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "object \"regular positions regular connections\" class field")
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "component \"positions\" value 1")
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "component \"connections\" value 2")
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "component \"data\" value 3")
            .map_err(|e| ApbsError::Io(e.to_string()))?;

        Ok(())
    }

    /// Write UHBD grid format.
    pub fn write_uhbd(&self, filename: &str, title: &str) -> ApbsResult<()> {
        if self.hx != self.hy || self.hy != self.hzed {
            return Err(ApbsError::UnsupportedFormat(
                "UHBD output requires uniform grid spacing".to_string(),
            ));
        }

        let mut file = File::create(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        writeln!(file, "{:>72}", title).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(
            file,
            "{:12.5e}{:12.5e}{:7}{:7}{:7}{:7}{:7}",
            1.0, 0.0, -1, 0, self.nz, 1, self.nz
        ).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(
            file,
            "{:7}{:7}{:7}{:12.5e}{:12.5e}{:12.5e}{:12.5e}",
            self.nx, self.ny, self.nz, self.hx, self.xmin - self.hx, self.ymin - self.hx, self.zmin - self.hx
        ).map_err(|e| ApbsError::Io(e.to_string()))?;
        writeln!(file, "{:12.5e}{:12.5e}{:12.5e}{:12.5e}", 0.0, 0.0, 0.0, 0.0)
            .map_err(|e| ApbsError::Io(e.to_string()))?;
        write!(file, "{:12.5e}{:12.5e}{:7}{:7}", 0.0, 0.0, 0, 0)
            .map_err(|e| ApbsError::Io(e.to_string()))?;

        for k in 0..self.nz {
            writeln!(file, "\n{:7}{:7}{:7}", k + 1, self.nx, self.ny)
                .map_err(|e| ApbsError::Io(e.to_string()))?;
            let mut icol = 0usize;
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let idx = self.ijk(i, j, k);
                    write!(file, " {:12.5e}", self.data[idx])
                        .map_err(|e| ApbsError::Io(e.to_string()))?;
                    icol += 1;
                    if icol == 6 {
                        icol = 0;
                        writeln!(file).map_err(|e| ApbsError::Io(e.to_string()))?;
                    }
                }
            }
            if icol != 0 {
                writeln!(file).map_err(|e| ApbsError::Io(e.to_string()))?;
            }
        }

        Ok(())
    }

    /// Read gzip-compressed OpenDX format file.
    /// Port of Vgrid_readGZ from vgrid.c line 462.
    pub fn read_gz(&mut self, filename: &str) -> ApbsResult<()> {
        use flate2::read::GzDecoder;
        use std::io::Read;

        let file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let mut decoder = GzDecoder::new(file);
        let mut contents = Vec::new();
        decoder.read_to_end(&mut contents)
            .map_err(|e| ApbsError::Io(format!("GzDecoder read error: {}", e)))?;

        // Parse as text
        let text = String::from_utf8(contents)
            .map_err(|e| ApbsError::Io(format!("UTF-8 error: {}", e)))?;

        let mut lines_iter = text.lines();
        let mut header = 0;

        self.hx = 0.0;
        self.hy = 0.0;
        self.hzed = 0.0;

        // Parse header (7 non-comment, non-blank lines)
        while header < 7 {
            let line = match lines_iter.next() {
                Some(l) => l,
                None => break,
            };
            let trimmed = line.trim();
            if trimmed.starts_with('#') || trimmed.is_empty() {
                continue;
            }

            match header {
                0 => {
                    // object 1 class gridpositions counts NX NY NZ
                    if let Some(rest) = trimmed.strip_prefix("object 1 class gridpositions counts ") {
                        let fields: Vec<&str> = rest.split_whitespace().collect();
                        if fields.len() >= 3 {
                            self.nx = fields[0].parse().unwrap_or(0);
                            self.ny = fields[1].parse().unwrap_or(0);
                            self.nz = fields[2].parse().unwrap_or(0);
                        }
                    }
                }
                1 => {
                    // origin X Y Z
                    if let Some(rest) = trimmed.strip_prefix("origin ") {
                        let fields: Vec<&str> = rest.split_whitespace().collect();
                        if fields.len() >= 3 {
                            self.xmin = fields[0].parse().unwrap_or(0.0);
                            self.ymin = fields[1].parse().unwrap_or(0.0);
                            self.zmin = fields[2].parse().unwrap_or(0.0);
                        }
                    }
                }
                2 | 3 | 4 => {
                    // delta lines - accumulate spacings
                    if let Some(rest) = trimmed.strip_prefix("delta ") {
                        let fields: Vec<&str> = rest.split_whitespace().collect();
                        if fields.len() >= 3 {
                            self.hx += fields[0].parse::<f64>().unwrap_or(0.0);
                            self.hy += fields[1].parse::<f64>().unwrap_or(0.0);
                            self.hzed += fields[2].parse::<f64>().unwrap_or(0.0);
                        }
                    }
                }
                _ => {} // lines 5-6: skip (gridconnections, etc.)
            }
            header += 1;
        }

        // Parse data values from remaining text (3 doubles per line, column-major order)
        let mut temp = Vec::with_capacity(self.nx * self.ny * self.nz);
        for line in lines_iter {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.contains("object") {
                continue;
            }
            for token in trimmed.split_whitespace() {
                if let Ok(v) = token.parse::<f64>() {
                    temp.push(v);
                }
            }
        }

        // Convert from column-major to row-major
        let len = self.nx * self.ny * self.nz;
        self.data = vec![0.0; len];
        let mut incr = 0;
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    if incr < temp.len() {
                        let u = k * self.nx * self.ny + j * self.nx + i;
                        self.data[u] = temp[incr];
                        incr += 1;
                    }
                }
            }
        }

        // Compute grid maxima
        self.xmax = self.xmin + (self.nx - 1) as f64 * self.hx;
        self.ymax = self.ymin + (self.ny - 1) as f64 * self.hy;
        self.zmax = self.zmin + (self.nz - 1) as f64 * self.hzed;
        self.readdata = true;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vgrid_create() {
        let grid = Vgrid::new(5, 5, 5, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        assert_eq!(grid.nx, 5);
        assert_eq!(grid.data.len(), 125);
    }

    #[test]
    fn test_vgrid_interpolation() {
        let mut grid = Vgrid::new(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        // Set center point
        let idx = grid.ijk(1, 1, 1);
        grid.data[idx] = 1.0;

        // Value at center should be 1.0
        let v = grid.value([1.0, 1.0, 1.0]).unwrap();
        assert!((v - 1.0).abs() < 1e-10);

        // Value at corner should be 0.0
        let v = grid.value([0.0, 0.0, 0.0]).unwrap();
        assert!(v.abs() < 1e-10);
    }

    #[test]
    fn test_vgrid_write_read_dx() {
        let mut grid = Vgrid::new(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        let idx = grid.ijk(1, 1, 1);
        grid.data[idx] = 42.0;

        let tmpfile = std::env::temp_dir().join("test_vgrid.dx");
        grid.write_dx(tmpfile.to_str().unwrap(), "test").unwrap();

        let mut grid2 = Vgrid::new(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        grid2.read_dx(tmpfile.to_str().unwrap()).unwrap();

        std::fs::remove_file(tmpfile).ok();
    }

    #[test]
    fn test_vgrid_write_dx_uses_c_file_order() {
        let grid = Vgrid::with_data(
            2, 2, 2,
            1.0, 1.0, 1.0,
            0.0, 0.0, 0.0,
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        );

        let tmpfile = std::env::temp_dir().join("test_vgrid_order.dx");
        grid.write_dx(tmpfile.to_str().unwrap(), "order").unwrap();
        let text = std::fs::read_to_string(&tmpfile).unwrap();
        std::fs::remove_file(tmpfile).ok();

        let mut vals = Vec::new();
        for line in text.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if trimmed.starts_with("object")
                || trimmed.starts_with("origin")
                || trimmed.starts_with("delta")
                || trimmed.starts_with("attribute")
                || trimmed.starts_with("component")
            {
                continue;
            }
            for field in trimmed.split_whitespace() {
                if let Ok(v) = field.parse::<f64>() {
                    vals.push(v);
                }
            }
        }

        assert_eq!(vals, vec![0.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0]);
    }
}
