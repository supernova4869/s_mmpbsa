// APBS vpmg.rs - Multigrid PBE solver core
// Port of src/mg/vpmg.h / vpmg.c

use apbs_generic::error::{ApbsError, ApbsResult};
use apbs_generic::vhal::*;
use apbs_generic::pbeparm::{PBEparm, PBEparmCalcEnergy};
use apbs_generic::vunit::*;
use apbs_generic::vcap::vcap_exp;
use crate::vpmgp::Vpmgp;

use std::sync::Arc;
use std::sync::OnceLock;
use apbs_generic::vpbe::Vpbe;
use crate::vgrid::Vgrid;
use rayon::prelude::*;

/// Maximum mesh partitions
pub const VPMGMAXPART: usize = 2000;

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

fn qm_ucap() -> f64 {
    static CAP: OnceLock<f64> = OnceLock::new();
    *CAP.get_or_init(|| {
        std::env::var("APBS_RUST_QM_UCAP")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|v| *v > 0.0)
            .unwrap_or(5.0)
    })
}

/// Multigrid PBE solver
pub struct Vpmg {
    /// PMG parameters
    pub pmgp: Vpmgp,
    /// PBE system info
    pub pbe: Arc<Vpbe>,
    /// X-shifted dielectric map
    pub epsx: Vec<f64>,
    /// Y-shifted dielectric map
    pub epsy: Vec<f64>,
    /// Z-shifted dielectric map
    pub epsz: Vec<f64>,
    /// Ion accessibility map [0,1]
    pub kappa: Vec<f64>,
    /// Potential map
    pub pot: Vec<f64>,
    /// Charge map
    pub charge: Vec<f64>,
    /// Integer params for FORTRAN
    pub iparm: Vec<i32>,
    /// Real params for FORTRAN
    pub rparm: Vec<f64>,
    /// Integer work array
    pub iwork: Vec<i32>,
    /// Real work array
    pub rwork: Vec<f64>,
    /// Operator coefficient a11
    pub a1cf: Vec<f64>,
    /// Operator coefficient a22
    pub a2cf: Vec<f64>,
    /// Operator coefficient a33
    pub a3cf: Vec<f64>,
    /// Helmholtz term
    pub ccf: Vec<f64>,
    /// Right-hand side
    pub fcf: Vec<f64>,
    /// True solution
    pub tcf: Vec<f64>,
    /// Solution
    pub u: Vec<f64>,
    /// X mesh coordinates
    pub xf: Vec<f64>,
    /// Y mesh coordinates
    pub yf: Vec<f64>,
    /// Z mesh coordinates
    pub zf: Vec<f64>,
    /// Boundary conditions for X face
    pub gxcf: Vec<f64>,
    /// Boundary conditions for Y face
    pub gycf: Vec<f64>,
    /// Boundary conditions for Z face
    pub gzcf: Vec<f64>,
    /// Partition mask array
    pub pvec: Vec<f64>,
    /// Per-atom partition weights
    pub atom_pvec: Vec<f64>,
    /// Dielectric energy from outside domain
    pub ext_di_energy: f64,
    /// Mobile ion energy from outside domain
    pub ext_qm_energy: f64,
    /// Fixed charge energy from outside domain
    pub ext_qf_energy: f64,
    /// Apolar energy from outside domain
    pub ext_np_energy: f64,
    /// Surface definition method
    pub surf_meth: VsurfMeth,
    /// Spline window parameter
    pub spline_win: f64,
    /// Charge discretization method
    pub charge_meth: VchrgMeth,
    /// Charge source
    pub charge_src: VchrgSrc,
    /// Whether fillco has been called
    pub filled: bool,
    /// External maps (optional)
    pub use_diel_x_map: bool,
    pub diel_x_map: Option<Vgrid>,
    pub use_diel_y_map: bool,
    pub diel_y_map: Option<Vgrid>,
    pub use_diel_z_map: bool,
    pub diel_z_map: Option<Vgrid>,
    pub use_kappa_map: bool,
    pub kappa_map: Option<Vgrid>,
    pub use_pot_map: bool,
    pub pot_map: Option<Vgrid>,
    pub use_charge_map: bool,
    pub charge_map: Option<Vgrid>,
}

impl Vpmg {
    fn apply_post_dirichlet_enabled() -> bool {
        static ENABLED: OnceLock<bool> = OnceLock::new();
        *ENABLED.get_or_init(|| std::env::var("APBS_RUST_APPLY_POST_DIRICHLET").map(|v| v != "0").unwrap_or(true))
    }

    fn apply_dirichlet_boundaries(&mut self) {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;

        let ijkx = |j: usize, k: usize, face: usize, ny_: usize, nz_: usize| -> usize {
            face * ny_ * nz_ + k * ny_ + j
        };
        let ijky = |i: usize, k: usize, face: usize, nx_: usize, nz_: usize| -> usize {
            face * nx_ * nz_ + k * nx_ + i
        };
        let ijkz = |i: usize, j: usize, face: usize, nx_: usize, ny_: usize| -> usize {
            face * nx_ * ny_ + j * nx_ + i
        };

        for j in 0..ny {
            for k in 0..nz {
                self.u[ijk(0, j, k, nx, ny)] = self.gxcf[ijkx(j, k, 0, ny, nz)];
                self.u[ijk(nx - 1, j, k, nx, ny)] = self.gxcf[ijkx(j, k, 1, ny, nz)];
            }
        }

        for i in 0..nx {
            for k in 0..nz {
                self.u[ijk(i, 0, k, nx, ny)] = self.gycf[ijky(i, k, 0, nx, nz)];
                self.u[ijk(i, ny - 1, k, nx, ny)] = self.gycf[ijky(i, k, 1, nx, nz)];
            }
        }

        for i in 0..nx {
            for j in 0..ny {
                self.u[ijk(i, j, 0, nx, ny)] = self.gzcf[ijkz(i, j, 0, nx, ny)];
                self.u[ijk(i, j, nz - 1, nx, ny)] = self.gzcf[ijkz(i, j, 1, nx, ny)];
            }
        }
    }

    pub fn new(
        pmgp: Vpmgp,
        pbe: Arc<Vpbe>,
        focus_flag: i32,
        pmg_old: Option<&Vpmg>,
        _pbeparm: &PBEparm,
        _calc_energy: PBEparmCalcEnergy,
    ) -> Self {
        let nf = pmgp.nf as usize;
        let narr = pmgp.narr as usize;
        let _narrc = pmgp.narrc as usize;
        let _nc = pmgp.nc as usize;
        let nx = pmgp.nx as usize;
        let ny = pmgp.ny as usize;
        let nz = pmgp.nz as usize;
        let num_atoms = pbe.alist.number_atoms();

        let mut vpmg = Self {
            pmgp,
            pbe,
            epsx: vec![0.0; narr],
            epsy: vec![0.0; narr],
            epsz: vec![0.0; narr],
            kappa: vec![0.0; narr],
            pot: vec![0.0; narr],
            charge: vec![0.0; narr],
            iparm: vec![0; 100],
            rparm: vec![0.0; 100],
            iwork: vec![0; narr * 4],
            rwork: vec![0.0; narr * 10],
            a1cf: vec![0.0; narr],
            a2cf: vec![0.0; narr],
            a3cf: vec![0.0; narr],
            ccf: vec![0.0; narr],
            fcf: vec![0.0; narr],
            tcf: vec![0.0; narr],
            u: vec![0.0; narr],
            xf: vec![0.0; nx],
            yf: vec![0.0; ny],
            zf: vec![0.0; nz],
            gxcf: vec![0.0; ny * nz * 4],
            gycf: vec![0.0; nx * nz * 4],
            gzcf: vec![0.0; nx * ny * 4],
            pvec: vec![0.0; nf],
            atom_pvec: vec![1.0; num_atoms],
            ext_di_energy: 0.0,
            ext_qm_energy: 0.0,
            ext_qf_energy: 0.0,
            ext_np_energy: 0.0,
            surf_meth: VsurfMeth::Mol,
            spline_win: 0.3,
            charge_meth: VchrgMeth::Tril,
            charge_src: VchrgSrc::Charge,
            filled: false,
            use_diel_x_map: false,
            diel_x_map: None,
            use_diel_y_map: false,
            diel_y_map: None,
            use_diel_z_map: false,
            diel_z_map: None,
            use_kappa_map: false,
            kappa_map: None,
            use_pot_map: false,
            pot_map: None,
            use_charge_map: false,
            charge_map: None,
        };

        // If focusing, interpolate boundary conditions from old grid
        if focus_flag != 0 {
            if let Some(old) = pmg_old {
                vpmg.focus_fill_bound(old);
            }
        }

        vpmg
    }

    /// Fill boundary conditions by interpolating from the old (coarse) grid solution.
    /// Port of focusFillBound from vpmg.c
    fn focus_fill_bound(&mut self, old: &Vpmg) {
        let nx_new = self.pmgp.nx as usize;
        let ny_new = self.pmgp.ny as usize;
        let nz_new = self.pmgp.nz as usize;
        let hx_new = self.pmgp.hx;
        let hy_new = self.pmgp.hy;
        let hz_new = self.pmgp.hzed;
        let xmin_new = self.pmgp.xmin;
        let ymin_new = self.pmgp.ymin;
        let zmin_new = self.pmgp.zmin;

        let nx_old = old.pmgp.nx as usize;
        let ny_old = old.pmgp.ny as usize;
        let nz_old = old.pmgp.nz as usize;
        let hx_old = old.pmgp.hx;
        let hy_old = old.pmgp.hy;
        let hz_old = old.pmgp.hzed;
        let xmin_old = old.pmgp.xmin;
        let ymin_old = old.pmgp.ymin;
        let zmin_old = old.pmgp.zmin;

        let data = &old.u;

        // Helper: trilinear interpolation from old grid
        let interpolate = |x: f64, y: f64, z: f64| -> f64 {
            let ifloat = ((x - xmin_old) / hx_old).clamp(0.0, (nx_old - 1) as f64);
            let jfloat = ((y - ymin_old) / hy_old).clamp(0.0, (ny_old - 1) as f64);
            let kfloat = ((z - zmin_old) / hz_old).clamp(0.0, (nz_old - 1) as f64);

            let ilo = (ifloat.floor() as usize).min(nx_old - 1);
            let ihi = (ifloat.ceil() as usize).min(nx_old - 1);
            let jlo = (jfloat.floor() as usize).min(ny_old - 1);
            let jhi = (jfloat.ceil() as usize).min(ny_old - 1);
            let klo = (kfloat.floor() as usize).min(nz_old - 1);
            let khi = (kfloat.ceil() as usize).min(nz_old - 1);

            let dx = ifloat - ilo as f64;
            let dy = jfloat - jlo as f64;
            let dz = kfloat - klo as f64;

            let c000 = data[ijk(ilo, jlo, klo, nx_old, ny_old)];
            let c100 = data[ijk(ihi, jlo, klo, nx_old, ny_old)];
            let c010 = data[ijk(ilo, jhi, klo, nx_old, ny_old)];
            let c110 = data[ijk(ihi, jhi, klo, nx_old, ny_old)];
            let c001 = data[ijk(ilo, jlo, khi, nx_old, ny_old)];
            let c101 = data[ijk(ihi, jlo, khi, nx_old, ny_old)];
            let c011 = data[ijk(ilo, jhi, khi, nx_old, ny_old)];
            let c111 = data[ijk(ihi, jhi, khi, nx_old, ny_old)];

            (1.0 - dx) * (1.0 - dy) * (1.0 - dz) * c000
                + dx * (1.0 - dy) * (1.0 - dz) * c100
                + (1.0 - dx) * dy * (1.0 - dz) * c010
                + dx * dy * (1.0 - dz) * c110
                + (1.0 - dx) * (1.0 - dy) * dz * c001
                + dx * (1.0 - dy) * dz * c101
                + (1.0 - dx) * dy * dz * c011
                + dx * dy * dz * c111
        };

        // X-face boundary conditions
        for k in 0..nz_new {
            for j in 0..ny_new {
                let y = ymin_new + j as f64 * hy_new;
                let z = zmin_new + k as f64 * hz_new;

                // Low-x face (face=0)
                let x_low = xmin_new;
                self.gxcf[ijkx(j, k, 0, ny_new, nz_new)] = interpolate(x_low, y, z);
                // High-x face (face=1)
                let x_high = xmin_new + (nx_new - 1) as f64 * hx_new;
                self.gxcf[ijkx(j, k, 1, ny_new, nz_new)] = interpolate(x_high, y, z);
                // Faces 2,3 (Neumann) = 0
                self.gxcf[ijkx(j, k, 2, ny_new, nz_new)] = 0.0;
                self.gxcf[ijkx(j, k, 3, ny_new, nz_new)] = 0.0;
            }
        }

        // Y-face boundary conditions
        for k in 0..nz_new {
            for i in 0..nx_new {
                let x = xmin_new + i as f64 * hx_new;
                let z = zmin_new + k as f64 * hz_new;

                // Low-y face (face=0)
                let y_low = ymin_new;
                self.gycf[ijky(i, k, 0, nx_new, nz_new)] = interpolate(x, y_low, z);
                // High-y face (face=1)
                let y_high = ymin_new + (ny_new - 1) as f64 * hy_new;
                self.gycf[ijky(i, k, 1, nx_new, nz_new)] = interpolate(x, y_high, z);
                // Faces 2,3 (Neumann) = 0
                self.gycf[ijky(i, k, 2, nx_new, nz_new)] = 0.0;
                self.gycf[ijky(i, k, 3, nx_new, nz_new)] = 0.0;
            }
        }

        // Z-face boundary conditions
        for j in 0..ny_new {
            for i in 0..nx_new {
                let x = xmin_new + i as f64 * hx_new;
                let y = ymin_new + j as f64 * hy_new;

                // Low-z face (face=0)
                let z_low = zmin_new;
                self.gzcf[ijkz(i, j, 0, nx_new, ny_new)] = interpolate(x, y, z_low);
                // High-z face (face=1)
                let z_high = zmin_new + (nz_new - 1) as f64 * hz_new;
                self.gzcf[ijkz(i, j, 1, nx_new, ny_new)] = interpolate(x, y, z_high);
                // Faces 2,3 (Neumann) = 0
                self.gzcf[ijkz(i, j, 2, nx_new, ny_new)] = 0.0;
                self.gzcf[ijkz(i, j, 3, nx_new, ny_new)] = 0.0;
            }
        }

        if debug_enabled() {
            eprintln!("[DEBUG-FOCUS] focus_fill_bound: old={}x{}x{}, new={}x{}x{}",
                nx_old, ny_old, nz_old, nx_new, ny_new, nz_new);
        }
    }

    /// Mark all grid points within a sphere of radius rtot with markVal
    fn mark_sphere(
        rtot: f64,
        tpos: [f64; 3],
        nx: usize, ny: usize, nz: usize,
        hx: f64, hy: f64, hz: f64,
        xmin: f64, ymin: f64, zmin: f64,
        array: &mut [f64],
        mark_val: f64,
    ) {
        let posx = tpos[0] - xmin;
        let posy = tpos[1] - ymin;
        let posz = tpos[2] - zmin;
        let rtot2 = rtot * rtot;

        let xrange = rtot + 0.5 * hx;
        let yrange = rtot + 0.5 * hy;
        let zrange = rtot + 0.5 * hz;

        let imin = 0.max(((posx - xrange) / hx).ceil() as i32) as usize;
        let jmin = 0.max(((posy - yrange) / hy).ceil() as i32) as usize;
        let kmin = 0.max(((posz - zrange) / hz).ceil() as i32) as usize;

        let imax = (nx - 1).min(((posx + xrange) / hx).floor() as usize);
        let jmax = (ny - 1).min(((posy + yrange) / hy).floor() as usize);
        let kmax = (nz - 1).min(((posz + zrange) / hz).floor() as usize);

        for i in imin..=imax {
            let fi = i as f64;
            let dx2 = (posx - hx * fi).powi(2);
            for j in jmin..=jmax {
                let fj = j as f64;
                let dy2 = (posy - hy * fj).powi(2);
                if dx2 + dy2 > rtot2 { continue; }
                for k in kmin..=kmax {
                    let fk = k as f64;
                    let dz2 = (posz - hz * fk).powi(2);
                    if dx2 + dy2 + dz2 <= rtot2 {
                        array[ijk(i, j, k, nx, ny)] = mark_val;
                    }
                }
            }
        }
    }

    /// Fill coefficient arrays
    pub fn fillco(
        &mut self,
        surf_meth: VsurfMeth,
        spline_win: f64,
        charge_meth: VchrgMeth,
    ) -> ApbsResult<()> {
        self.surf_meth = surf_meth;
        self.spline_win = spline_win;
        self.charge_meth = charge_meth;

        // Compute mesh coordinates
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;

        for i in 0..nx {
            self.xf[i] = self.pmgp.xmin + i as f64 * self.pmgp.hx;
        }
        for j in 0..ny {
            self.yf[j] = self.pmgp.ymin + j as f64 * self.pmgp.hy;
        }
        for k in 0..nz {
            self.zf[k] = self.pmgp.zmin + k as f64 * self.pmgp.hzed;
        }

        // Fill dielectric coefficients
        self.fillco_coef()?;

        // Fill charge array
        self.fillco_charge()?;

        // Compute boundary conditions (skip for focusing)
        if self.pmgp.bcfl != Vbcfl::Focus && self.pmgp.bcfl != Vbcfl::Map {
            self.bc_calc()?;
        }

        self.filled = true;
        Ok(())
    }

    /// Fill dielectric and ion accessibility coefficients
    fn fillco_coef(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let nf = nx * ny * nz;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;

        let eps_solute = self.pbe.solute_diel;
        let eps_solvent = self.pbe.solvent_diel;
        let srad = self.pbe.solvent_radius;
        let irad = self.pbe.max_ion_radius;
        let ionstr = self.pbe.bulk_ionic_strength;

        let num_atoms = self.pbe.alist.number_atoms();

        if self.use_diel_x_map || self.use_diel_y_map || self.use_diel_z_map || self.use_kappa_map {
            // External map mode
            self.fillco_coef_map()?;
        } else {
            // Molecular dielectric mode
            // Step 1: Initialize all to solvent dielectric
            for i in 0..nf {
                self.epsx[i] = eps_solvent;
                self.epsy[i] = eps_solvent;
                self.epsz[i] = eps_solvent;
            }

            // Step 2: For each atom with non-zero radius, mark its
            // solvent-inflated VdW sphere as solute dielectric
            for iatom in 0..num_atoms {
                let atom = self.pbe.alist.get_atom(iatom);
                let apos = atom.position;
                let arad = atom.radius;

                if arad > VSMALL {
                    // Mark x-shifted dielectric (staggered grid)
                    Self::mark_sphere(
                        arad + srad, apos,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin + 0.5 * hx, self.pmgp.ymin, self.pmgp.zmin,
                        &mut self.epsx, eps_solute,
                    );
                    // Mark y-shifted dielectric
                    Self::mark_sphere(
                        arad + srad, apos,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin, self.pmgp.ymin + 0.5 * hy, self.pmgp.zmin,
                        &mut self.epsy, eps_solute,
                    );
                    // Mark z-shifted dielectric
                    Self::mark_sphere(
                        arad + srad, apos,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin, self.pmgp.ymin, self.pmgp.zmin + 0.5 * hzed,
                        &mut self.epsz, eps_solute,
                    );
                }
            }

            // Step 3: SAS surface re-marking
            // For each SAS point, reset the dielectric back to solvent
            if srad > VSMALL {
                let _ = self.pbe.acc.lock().unwrap().sasa(srad);
                let mut sas_points: Vec<[f64; 3]> = Vec::new();
                {
                    let acc = self.pbe.acc.lock().unwrap();
                    for iatom in 0..num_atoms {
                        if let Some(surf) = acc.surf[iatom].as_ref() {
                            for ipt in 0..surf.npts {
                                if surf.bpts[ipt] {
                                    sas_points.push([
                                        surf.xpts[ipt],
                                        surf.ypts[ipt],
                                        surf.zpts[ipt],
                                    ]);
                                }
                            }
                        }
                    }
                }
                for position in sas_points {
                    Self::mark_sphere(
                        srad, position,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin + 0.5 * hx, self.pmgp.ymin, self.pmgp.zmin,
                        &mut self.epsx, eps_solvent,
                    );
                    Self::mark_sphere(
                        srad, position,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin, self.pmgp.ymin + 0.5 * hy, self.pmgp.zmin,
                        &mut self.epsy, eps_solvent,
                    );
                    Self::mark_sphere(
                        srad, position,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin, self.pmgp.ymin, self.pmgp.zmin + 0.5 * hzed,
                        &mut self.epsz, eps_solvent,
                    );
                }
            }

            // Step 4: Smooth dielectric (Bruccoleri 9-point harmonic mean)
            if self.surf_meth == VsurfMeth::MolSmooth {
                self.fillco_coef_smooth()?;
            }
        }

        // Fill Helmholtz term (kappa^2 * eps) and ion accessibility
        let ionmask = if ionstr > VPMGSMALL { 1.0 } else { 0.0 };

        // Initialize kappa: all accessible if ionic strength > 0
        for i in 0..nf {
            self.kappa[i] = ionmask;
            // Store dielectric in a1cf/a2cf/a3cf for the operator
            self.a1cf[i] = self.epsx[i];
            self.a2cf[i] = self.epsy[i];
            self.a3cf[i] = self.epsz[i];
        }

        // Mark ion-inaccessible regions (inside inflated VdW spheres)
        if ionstr > VPMGSMALL {
            for iatom in 0..num_atoms {
                let atom = self.pbe.alist.get_atom(iatom);
                let apos = atom.position;
                let arad = atom.radius;

                if arad > VSMALL {
                    Self::mark_sphere(
                        irad + arad, apos,
                        nx, ny, nz, hx, hy, hzed,
                        self.pmgp.xmin, self.pmgp.ymin, self.pmgp.zmin,
                        &mut self.kappa, 0.0,
                    );
                }
            }
        }

        // Fill Helmholtz term AFTER kappa is marked (ccf = zkappa2 * kappa)
        for i in 0..nf {
            self.ccf[i] = self.pbe.zkappa2 * self.kappa[i];
        }

        Ok(())
    }

    /// Fill coefficients from external dielectric/kappa maps
    fn fillco_coef_map(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;

        // Validate that all required maps are present
        if self.use_diel_x_map && self.diel_x_map.is_none() {
            return Err(apbs_generic::error::ApbsError::InvalidParameter(
                "dielectric X map enabled but not provided".to_string()));
        }
        if self.use_diel_y_map && self.diel_y_map.is_none() {
            return Err(apbs_generic::error::ApbsError::InvalidParameter(
                "dielectric Y map enabled but not provided".to_string()));
        }
        if self.use_diel_z_map && self.diel_z_map.is_none() {
            return Err(apbs_generic::error::ApbsError::InvalidParameter(
                "dielectric Z map enabled but not provided".to_string()));
        }
        if self.use_kappa_map && self.kappa_map.is_none() {
            return Err(apbs_generic::error::ApbsError::InvalidParameter(
                "kappa map enabled but not provided".to_string()));
        }

        // Read dielectric maps at staggered grid positions
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = ijk(i, j, k, nx, ny);
                    let xf = self.xf[i];
                    let yf = self.yf[j];
                    let zf = self.zf[k];

                    if self.use_diel_x_map {
                        if let Some(ref map) = self.diel_x_map {
                            self.epsx[idx] = map.value([xf + 0.5 * hx, yf, zf]).unwrap_or(self.pbe.solvent_diel);
                        }
                    }
                    if self.use_diel_y_map {
                        if let Some(ref map) = self.diel_y_map {
                            self.epsy[idx] = map.value([xf, yf + 0.5 * hy, zf]).unwrap_or(self.pbe.solvent_diel);
                        }
                    }
                    if self.use_diel_z_map {
                        if let Some(ref map) = self.diel_z_map {
                            self.epsz[idx] = map.value([xf, yf, zf + 0.5 * hzed]).unwrap_or(self.pbe.solvent_diel);
                        }
                    }
                    if self.use_kappa_map {
                        if let Some(ref map) = self.kappa_map {
                            self.kappa[idx] = map.value([xf, yf, zf]).unwrap_or(0.0);
                        }
                    }
                }
            }
        }

        // Scale kappa to [0,1]
        if self.use_kappa_map {
            let kmax = self.kappa[..nx*ny*nz].iter().cloned().fold(0.0f64, f64::max);
            if kmax > VSMALL {
                for v in &mut self.kappa[..nx*ny*nz] {
                    *v /= kmax;
                }
            }
        }

        Ok(())
    }

    /// Smooth dielectric using 9-point harmonic mean (Bruccoleri et al. 1997)
    /// Assumes non-smoothed values were already placed in epsx/epsy/epsz
    fn fillco_coef_smooth(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let nf = nx * ny * nz;
        let epsw = self.pbe.solvent_diel;

        // Copy existing dielectric to temporary arrays a1cf/a2cf/a3cf
        for i in 0..nf {
            self.a1cf[i] = self.epsx[i];
            self.a2cf[i] = self.epsy[i];
            self.a3cf[i] = self.epsz[i];
            self.epsx[i] = epsw;
            self.epsy[i] = epsw;
            self.epsz[i] = epsw;
        }

        // Smooth using 9-point harmonic mean for each shifted direction
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = ijk(i, j, k, nx, ny);

                    // --- epsx (X-shifted) ---
                    let mut frac = 1.0 / self.a1cf[idx];
                    frac += 1.0 / self.a2cf[idx];
                    frac += 1.0 / self.a3cf[idx];
                    let mut numpts = 3;

                    if j > 0 {
                        frac += 1.0 / self.a2cf[ijk(i, j - 1, k, nx, ny)];
                        numpts += 1;
                    }
                    if k > 0 {
                        frac += 1.0 / self.a3cf[ijk(i, j, k - 1, nx, ny)];
                        numpts += 1;
                    }
                    if i < nx - 1 {
                        frac += 1.0 / self.a2cf[ijk(i + 1, j, k, nx, ny)];
                        frac += 1.0 / self.a3cf[ijk(i + 1, j, k, nx, ny)];
                        numpts += 2;
                        if j > 0 {
                            frac += 1.0 / self.a2cf[ijk(i + 1, j - 1, k, nx, ny)];
                            numpts += 1;
                        }
                        if k > 0 {
                            frac += 1.0 / self.a3cf[ijk(i + 1, j, k - 1, nx, ny)];
                            numpts += 1;
                        }
                    }
                    self.epsx[idx] = numpts as f64 / frac;

                    // --- epsy (Y-shifted) ---
                    let mut frac = 1.0 / self.a2cf[idx];
                    frac += 1.0 / self.a1cf[idx];
                    frac += 1.0 / self.a3cf[idx];
                    let mut numpts = 3;

                    if i > 0 {
                        frac += 1.0 / self.a1cf[ijk(i - 1, j, k, nx, ny)];
                        numpts += 1;
                    }
                    if k > 0 {
                        frac += 1.0 / self.a3cf[ijk(i, j, k - 1, nx, ny)];
                        numpts += 1;
                    }
                    if j < ny - 1 {
                        frac += 1.0 / self.a1cf[ijk(i, j + 1, k, nx, ny)];
                        frac += 1.0 / self.a3cf[ijk(i, j + 1, k, nx, ny)];
                        numpts += 2;
                        if i > 0 {
                            frac += 1.0 / self.a1cf[ijk(i - 1, j + 1, k, nx, ny)];
                            numpts += 1;
                        }
                        if k > 0 {
                            frac += 1.0 / self.a3cf[ijk(i, j + 1, k - 1, nx, ny)];
                            numpts += 1;
                        }
                    }
                    self.epsy[idx] = numpts as f64 / frac;

                    // --- epsz (Z-shifted) ---
                    let mut frac = 1.0 / self.a3cf[idx];
                    frac += 1.0 / self.a1cf[idx];
                    frac += 1.0 / self.a2cf[idx];
                    let mut numpts = 3;

                    if i > 0 {
                        frac += 1.0 / self.a1cf[ijk(i - 1, j, k, nx, ny)];
                        numpts += 1;
                    }
                    if j > 0 {
                        frac += 1.0 / self.a2cf[ijk(i, j - 1, k, nx, ny)];
                        numpts += 1;
                    }
                    if k < nz - 1 {
                        frac += 1.0 / self.a1cf[ijk(i, j, k + 1, nx, ny)];
                        frac += 1.0 / self.a2cf[ijk(i, j, k + 1, nx, ny)];
                        numpts += 2;
                        if i > 0 {
                            frac += 1.0 / self.a1cf[ijk(i - 1, j, k + 1, nx, ny)];
                            numpts += 1;
                        }
                        if j > 0 {
                            frac += 1.0 / self.a2cf[ijk(i, j - 1, k + 1, nx, ny)];
                            numpts += 1;
                        }
                    }
                    self.epsz[idx] = numpts as f64 / frac;
                }
            }
        }

        Ok(())
    }

    /// Fill charge array from atom charges
    fn fillco_charge(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let _hx = self.pmgp.hx;
        let _hy = self.pmgp.hy;
        let _hzed = self.pmgp.hzed;

        // Reset charge array
        for i in 0..nx * ny * nz {
            self.charge[i] = 0.0;
        }

        match self.charge_meth {
            VchrgMeth::Tril => self.fillco_charge_trilinear()?,
            VchrgMeth::Bspl2 => self.fillco_charge_bspline2()?,
            VchrgMeth::Bspl4 => self.fillco_charge_bspline4()?,
        }

        // Copy charge to fcf (right-hand side)
        self.fcf[..nx*ny*nz].copy_from_slice(&self.charge[..nx*ny*nz]);

        Ok(())
    }

    /// Trilinear charge discretization (VCM_TRIL)
    fn fillco_charge_trilinear(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;
        let zmagic = self.pbe.zmagic;
        let scale = zmagic / (hx * hy * hzed);

        let num_atoms = self.pbe.alist.number_atoms();
        let partials: Vec<Vec<(usize, f64)>> = (0..num_atoms).into_par_iter().map(|iatom| {
            let mut local = Vec::with_capacity(8);
            let atom = self.pbe.alist.get_atom(iatom);
            let apos = atom.position;
            let charge = atom.charge;
            let atom_weight = self.atom_pvec.get(iatom).copied().unwrap_or(1.0);

            if atom_weight <= 0.0 {
                return local;
            }

            if apos[0] <= xmin || apos[0] >= xmax ||
               apos[1] <= ymin || apos[1] >= ymax ||
               apos[2] <= zmin || apos[2] >= zmax {
                return local;
            }

            let posx = apos[0] - xmin;
            let posy = apos[1] - ymin;
            let posz = apos[2] - zmin;
            let q = charge * scale;

            let ifloat = posx / hx;
            let jfloat = posy / hy;
            let kfloat = posz / hzed;

            let ihi = ifloat.ceil() as i32;
            let ilo = ifloat.floor() as i32;
            let jhi = jfloat.ceil() as i32;
            let jlo = jfloat.floor() as i32;
            let khi = kfloat.ceil() as i32;
            let klo = kfloat.floor() as i32;

            let dx = ifloat - ilo as f64;
            let dy = jfloat - jlo as f64;
            let dz = kfloat - klo as f64;

            let push = |local: &mut Vec<(usize, f64)>, ii: i32, jj: i32, kk: i32, w: f64| {
                if ii >= 0 && (ii as usize) < nx && jj >= 0 && (jj as usize) < ny && kk >= 0 && (kk as usize) < nz {
                    local.push((ijk(ii as usize, jj as usize, kk as usize, nx, ny), w * q));
                }
            };

            push(&mut local, ihi, jhi, khi, dx * dy * dz);
            push(&mut local, ihi, jlo, khi, dx * (1.0 - dy) * dz);
            push(&mut local, ihi, jhi, klo, dx * dy * (1.0 - dz));
            push(&mut local, ihi, jlo, klo, dx * (1.0 - dy) * (1.0 - dz));
            push(&mut local, ilo, jhi, khi, (1.0 - dx) * dy * dz);
            push(&mut local, ilo, jlo, khi, (1.0 - dx) * (1.0 - dy) * dz);
            push(&mut local, ilo, jhi, klo, (1.0 - dx) * dy * (1.0 - dz));
            push(&mut local, ilo, jlo, klo, (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
            local
        }).collect();

        for partial in &partials {
            for &(idx, value) in partial {
                self.charge[idx] += value;
            }
        }

        Ok(())
    }

    /// Cubic B-spline charge discretization (VCM_BSPL2)
    fn fillco_charge_bspline2(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;
        let zmagic = self.pbe.zmagic;
        let scale = zmagic / (hx * hy * hzed);

        let num_atoms = self.pbe.alist.number_atoms();
        let partials: Vec<Vec<(usize, f64)>> = (0..num_atoms).into_par_iter().map(|iatom| {
            let mut local = Vec::with_capacity(27);
            let atom = self.pbe.alist.get_atom(iatom);
            let apos = atom.position;
            let charge = atom.charge;

            if apos[0] <= (xmin - hx) || apos[0] >= (xmax + hx) ||
               apos[1] <= (ymin - hy) || apos[1] >= (ymax + hy) ||
               apos[2] <= (zmin - hzed) || apos[2] >= (zmax + hzed) {
                return local;
            }

            let posx = apos[0] - xmin;
            let posy = apos[1] - ymin;
            let posz = apos[2] - zmin;
            let q = charge * scale;

            let ifloat = posx / hx;
            let jfloat = posy / hy;
            let kfloat = posz / hzed;

            let ip1 = ifloat.ceil() as i32;
            let ip2 = (ip1 + 1).min(nx as i32 - 1);
            let im1 = ifloat.floor() as i32;
            let im2 = (im1 - 1).max(0);
            let jp1 = jfloat.ceil() as i32;
            let jp2 = (jp1 + 1).min(ny as i32 - 1);
            let jm1 = jfloat.floor() as i32;
            let jm2 = (jm1 - 1).max(0);
            let kp1 = kfloat.ceil() as i32;
            let kp2 = (kp1 + 1).min(nz as i32 - 1);
            let km1 = kfloat.floor() as i32;
            let km2 = (km1 - 1).max(0);

            for ii in im2..=ip2 {
                let fi = vfchi(ii as f64, ifloat);
                let bi = bspline2(fi);
                for jj in jm2..=jp2 {
                    let fj = vfchi(jj as f64, jfloat);
                    let bj = bspline2(fj);
                    for kk in km2..=kp2 {
                        let fk = vfchi(kk as f64, kfloat);
                        let bk = bspline2(fk);
                        let w = bi * bj * bk;
                        if w.abs() > VSMALL {
                            local.push((ijk(ii as usize, jj as usize, kk as usize, nx, ny), w * q));
                        }
                    }
                }
            }
            local
        }).collect();

        for partial in &partials {
            for &(idx, value) in partial {
                self.charge[idx] += value;
            }
        }

        Ok(())
    }

    /// Quartic B-spline charge discretization (VCM_BSPL4)
    /// C dispatches to fillcoPermanentMultipole for PQR point charges.
    /// For point charges (no dipoles/quadrupoles), reduces to B-spline interpolation.
    fn fillco_charge_bspline4(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;
        let zmagic = self.pbe.zmagic;
        let scale = zmagic / (hx * hy * hzed);

        let num_atoms = self.pbe.alist.number_atoms();
        let partials: Vec<Vec<(usize, f64)>> = (0..num_atoms).into_par_iter().map(|iatom| {
            let mut local = Vec::with_capacity(125);
            let atom = self.pbe.alist.get_atom(iatom);
            let apos = atom.position;
            let charge = atom.charge;
            let atom_weight = self.atom_pvec.get(iatom).copied().unwrap_or(1.0);

            if atom_weight <= 0.0 {
                return local;
            }

            if apos[0] <= (xmin - 2.0 * hx) || apos[0] >= (xmax + 2.0 * hx) ||
               apos[1] <= (ymin - 2.0 * hy) || apos[1] >= (ymax + 2.0 * hy) ||
               apos[2] <= (zmin - 2.0 * hzed) || apos[2] >= (zmax + 2.0 * hzed) {
                return local;
            }

            let posx = apos[0] - xmin;
            let posy = apos[1] - ymin;
            let posz = apos[2] - zmin;
            let q = charge * atom_weight * scale;

            let ifloat = posx / hx;
            let jfloat = posy / hy;
            let kfloat = posz / hzed;

            let ip1 = ifloat.ceil() as i32;
            let ip2 = (ip1 + 2).min(nx as i32 - 1);
            let im1 = ifloat.floor() as i32;
            let im2 = (im1 - 2).max(0);
            let jp1 = jfloat.ceil() as i32;
            let jp2 = (jp1 + 2).min(ny as i32 - 1);
            let jm1 = jfloat.floor() as i32;
            let jm2 = (jm1 - 2).max(0);
            let kp1 = kfloat.ceil() as i32;
            let kp2 = (kp1 + 2).min(nz as i32 - 1);
            let km1 = kfloat.floor() as i32;
            let km2 = (km1 - 2).max(0);

            for ii in im2..=ip2 {
                let wi = bspline4(vfchi4(ii as f64, ifloat));
                if wi.abs() <= VSMALL {
                    continue;
                }
                for jj in jm2..=jp2 {
                    let wj = bspline4(vfchi4(jj as f64, jfloat));
                    if wj.abs() <= VSMALL {
                        continue;
                    }
                    for kk in km2..=kp2 {
                        let wk = bspline4(vfchi4(kk as f64, kfloat));
                        let w = wi * wj * wk;
                        if w.abs() > VSMALL {
                            local.push((ijk(ii as usize, jj as usize, kk as usize, nx, ny), w * q));
                        }
                    }
                }
            }
            local
        }).collect();

        for partial in &partials {
            for &(idx, value) in partial {
                self.charge[idx] += value;
            }
        }

        Ok(())
    }

    /// Compute boundary conditions
    fn bc_calc(&mut self) -> ApbsResult<()> {
        let _nx = self.pmgp.nx as usize;
        let _ny = self.pmgp.ny as usize;
        let _nz = self.pmgp.nz as usize;

        // Zero out boundary arrays
        // gxcf: 4 values per (j,k) face pair
        for v in &mut self.gxcf { *v = 0.0; }
        for v in &mut self.gycf { *v = 0.0; }
        for v in &mut self.gzcf { *v = 0.0; }

        match self.pmgp.bcfl {
            Vbcfl::Zero => Ok(()),
            Vbcfl::SDH => {
                self.bcfl_sdh()?;
                Ok(())
            }
            Vbcfl::MDH => {
                self.bcfl_mdh()?;
                Ok(())
            }
            Vbcfl::Mem => {
                self.bcfl_mem()?;
                Ok(())
            }
            _ => Ok(()),
        }
    }

    /// Single Debye-Huckel boundary conditions (SDH)
    /// Uses multipole expansion: monopole + dipole + traceless quadrupole
    fn bcfl_sdh(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let _hx = self.pmgp.hx;
        let _hy = self.pmgp.hy;
        let _hzed = self.pmgp.hzed;
        let eps_w = self.pbe.solvent_diel;
        let eps_p = self.pbe.solute_diel;
        let _t = self.pbe.temperature;
        let xkappa = self.pbe.xkappa;
        let size = self.pbe.solute_radius;
        let center = self.pbe.solute_center;

        // `multipolebc` already carries the dielectric dependence, so the
        // SDH prefactor matches the original C code with no extra `eps_w`.
        let pre = EC * EC / (4.0 * PI * EPS0 * KB * _t) * 1.0e10;

        // Aggregate multipole moments
        let mut sdhcharge = 0.0;
        let mut sdhdipole = [0.0; 3];
        let mut sdhquadrupole = [0.0; 9];

        let num_atoms = self.pbe.alist.number_atoms();
        for iatom in 0..num_atoms {
            let atom = self.pbe.alist.get_atom(iatom);
            let apos = atom.position;
            let charge = atom.charge;

            let xr = apos[0] - center[0];
            let yr = apos[1] - center[1];
            let zr = apos[2] - center[2];

            sdhcharge += charge;
            sdhdipole[0] += xr * charge;
            sdhdipole[1] += yr * charge;
            sdhdipole[2] += zr * charge;

            let traced = [
                xr * xr * charge, xr * yr * charge, xr * zr * charge,
                yr * xr * charge, yr * yr * charge, yr * zr * charge,
                zr * xr * charge, zr * yr * charge, zr * zr * charge,
            ];
            let qave = (traced[0] + traced[4] + traced[8]) / 3.0;
            sdhquadrupole[0] += 1.5 * (traced[0] - qave);
            sdhquadrupole[1] += 1.5 * traced[1];
            sdhquadrupole[2] += 1.5 * traced[2];
            sdhquadrupole[3] += 1.5 * traced[3];
            sdhquadrupole[4] += 1.5 * (traced[4] - qave);
            sdhquadrupole[5] += 1.5 * traced[5];
            sdhquadrupole[6] += 1.5 * traced[6];
            sdhquadrupole[7] += 1.5 * traced[7];
            sdhquadrupole[8] += 1.5 * (traced[8] - qave);
        }

        // Divide quadrupole by 3 (traceless definition factor)
        let qxx = sdhquadrupole[0] / 3.0;
        let qxy = sdhquadrupole[1] / 3.0;
        let qxz = sdhquadrupole[2] / 3.0;
        let qyy = sdhquadrupole[4] / 3.0;
        let qyz = sdhquadrupole[5] / 3.0;
        let qzz = sdhquadrupole[8] / 3.0;
        let ux = sdhdipole[0];
        let uy = sdhdipole[1];
        let uz = sdhdipole[2];

        // Evaluate boundary values on each face
        // X faces
        for k in 0..nz {
            for j in 0..ny {
                // x = xmin face
                let gpos = [self.xf[0], self.yf[j], self.zf[k]];
                let xr = gpos[0] - center[0];
                let yr = gpos[1] - center[1];
                let zr = gpos[2] - center[2];
                let dist = (xr*xr + yr*yr + zr*zr).sqrt();
                let tensor = multipolebc(dist, xkappa, eps_p, eps_w, size);
                let val = pre * (sdhcharge * tensor[0]
                    - ux * xr * tensor[1] - uy * yr * tensor[1] - uz * zr * tensor[1]
                    + qxx * xr * xr * tensor[2] + qyy * yr * yr * tensor[2] + qzz * zr * zr * tensor[2]
                    + 2.0 * qxy * xr * yr * tensor[2] + 2.0 * qxz * xr * zr * tensor[2] + 2.0 * qyz * yr * zr * tensor[2]);
                self.gxcf[ijkx(j, k, 0, ny, nz)] += val;

                // x = xmax face
                let gpos = [self.xf[nx-1], self.yf[j], self.zf[k]];
                let xr = gpos[0] - center[0];
                let yr = gpos[1] - center[1];
                let zr = gpos[2] - center[2];
                let dist = (xr*xr + yr*yr + zr*zr).sqrt();
                let tensor = multipolebc(dist, xkappa, eps_p, eps_w, size);
                let val = pre * (sdhcharge * tensor[0]
                    - ux * xr * tensor[1] - uy * yr * tensor[1] - uz * zr * tensor[1]
                    + qxx * xr * xr * tensor[2] + qyy * yr * yr * tensor[2] + qzz * zr * zr * tensor[2]
                    + 2.0 * qxy * xr * yr * tensor[2] + 2.0 * qxz * xr * zr * tensor[2] + 2.0 * qyz * yr * zr * tensor[2]);
                self.gxcf[ijkx(j, k, 1, ny, nz)] += val;
            }
        }

        // Y faces
        for k in 0..nz {
            for i in 0..nx {
                // y = ymin face
                let gpos = [self.xf[i], self.yf[0], self.zf[k]];
                let xr = gpos[0] - center[0];
                let yr = gpos[1] - center[1];
                let zr = gpos[2] - center[2];
                let dist = (xr*xr + yr*yr + zr*zr).sqrt();
                let tensor = multipolebc(dist, xkappa, eps_p, eps_w, size);
                let val = pre * (sdhcharge * tensor[0]
                    - ux * xr * tensor[1] - uy * yr * tensor[1] - uz * zr * tensor[1]
                    + qxx * xr * xr * tensor[2] + qyy * yr * yr * tensor[2] + qzz * zr * zr * tensor[2]
                    + 2.0 * qxy * xr * yr * tensor[2] + 2.0 * qxz * xr * zr * tensor[2] + 2.0 * qyz * yr * zr * tensor[2]);
                self.gycf[ijky(i, k, 0, nx, nz)] += val;

                // y = ymax face
                let gpos = [self.xf[i], self.yf[ny-1], self.zf[k]];
                let xr = gpos[0] - center[0];
                let yr = gpos[1] - center[1];
                let zr = gpos[2] - center[2];
                let dist = (xr*xr + yr*yr + zr*zr).sqrt();
                let tensor = multipolebc(dist, xkappa, eps_p, eps_w, size);
                let val = pre * (sdhcharge * tensor[0]
                    - ux * xr * tensor[1] - uy * yr * tensor[1] - uz * zr * tensor[1]
                    + qxx * xr * xr * tensor[2] + qyy * yr * yr * tensor[2] + qzz * zr * zr * tensor[2]
                    + 2.0 * qxy * xr * yr * tensor[2] + 2.0 * qxz * xr * zr * tensor[2] + 2.0 * qyz * yr * zr * tensor[2]);
                self.gycf[ijky(i, k, 1, nx, nz)] += val;
            }
        }

        // Z faces
        for j in 0..ny {
            for i in 0..nx {
                // z = zmin face
                let gpos = [self.xf[i], self.yf[j], self.zf[0]];
                let xr = gpos[0] - center[0];
                let yr = gpos[1] - center[1];
                let zr = gpos[2] - center[2];
                let dist = (xr*xr + yr*yr + zr*zr).sqrt();
                let tensor = multipolebc(dist, xkappa, eps_p, eps_w, size);
                let val = pre * (sdhcharge * tensor[0]
                    - ux * xr * tensor[1] - uy * yr * tensor[1] - uz * zr * tensor[1]
                    + qxx * xr * xr * tensor[2] + qyy * yr * yr * tensor[2] + qzz * zr * zr * tensor[2]
                    + 2.0 * qxy * xr * yr * tensor[2] + 2.0 * qxz * xr * zr * tensor[2] + 2.0 * qyz * yr * zr * tensor[2]);
                self.gzcf[ijkz(i, j, 0, nx, ny)] += val;

                // z = zmax face
                let gpos = [self.xf[i], self.yf[j], self.zf[nz-1]];
                let xr = gpos[0] - center[0];
                let yr = gpos[1] - center[1];
                let zr = gpos[2] - center[2];
                let dist = (xr*xr + yr*yr + zr*zr).sqrt();
                let tensor = multipolebc(dist, xkappa, eps_p, eps_w, size);
                let val = pre * (sdhcharge * tensor[0]
                    - ux * xr * tensor[1] - uy * yr * tensor[1] - uz * zr * tensor[1]
                    + qxx * xr * xr * tensor[2] + qyy * yr * yr * tensor[2] + qzz * zr * zr * tensor[2]
                    + 2.0 * qxy * xr * yr * tensor[2] + 2.0 * qxz * xr * zr * tensor[2] + 2.0 * qyz * yr * zr * tensor[2]);
                self.gzcf[ijkz(i, j, 1, nx, ny)] += val;
            }
        }

        Ok(())
    }

    /// Multi-ion Debye-Huckel boundary conditions (MDH)
    /// Per-atom Debye-Huckel potential at each boundary grid point
    fn bcfl_mdh(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let eps_w = self.pbe.solvent_diel;
        let _t = self.pbe.temperature;
        let xkappa = self.pbe.xkappa;

        let pre1 = EC / (4.0 * PI * EPS0 * eps_w * KB * _t) * 1.0e10;
        let num_atoms = self.pbe.alist.number_atoms();
        let xlen = 2 * ny * nz;
        let ylen = 2 * nx * nz;
        let zlen = 2 * nx * ny;

        let partials: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..num_atoms).into_par_iter().map(|iatom| {
            let atom = self.pbe.alist.get_atom(iatom);
            let apos = atom.position;
            let charge = EC * atom.charge;
            let size = atom.radius;
            let mut gx = vec![0.0; xlen];
            let mut gy = vec![0.0; ylen];
            let mut gz = vec![0.0; zlen];

            for k in 0..nz {
                for j in 0..ny {
                    let gpos = [self.xf[0], self.yf[j], self.zf[k]];
                    let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                    let val = if xkappa > VSMALL {
                        pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                    } else {
                        pre1 * (charge / dist)
                    };
                    gx[ijkx(j, k, 0, ny, nz)] += val;

                    let gpos = [self.xf[nx-1], self.yf[j], self.zf[k]];
                    let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                    let val = if xkappa > VSMALL {
                        pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                    } else {
                        pre1 * (charge / dist)
                    };
                    gx[ijkx(j, k, 1, ny, nz)] += val;
                }
            }

            for k in 0..nz {
                for i in 0..nx {
                    let gpos = [self.xf[i], self.yf[0], self.zf[k]];
                    let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                    let val = if xkappa > VSMALL {
                        pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                    } else {
                        pre1 * (charge / dist)
                    };
                    gy[ijky(i, k, 0, nx, nz)] += val;

                    let gpos = [self.xf[i], self.yf[ny-1], self.zf[k]];
                    let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                    let val = if xkappa > VSMALL {
                        pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                    } else {
                        pre1 * (charge / dist)
                    };
                    gy[ijky(i, k, 1, nx, nz)] += val;
                }
            }

            for j in 0..ny {
                for i in 0..nx {
                    let gpos = [self.xf[i], self.yf[j], self.zf[0]];
                    let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                    let val = if xkappa > VSMALL {
                        pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                    } else {
                        pre1 * (charge / dist)
                    };
                    gz[ijkz(i, j, 0, nx, ny)] += val;

                    let gpos = [self.xf[i], self.yf[j], self.zf[nz-1]];
                    let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                    let val = if xkappa > VSMALL {
                        pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                    } else {
                        pre1 * (charge / dist)
                    };
                    gz[ijkz(i, j, 1, nx, ny)] += val;
                }
            }

            (gx, gy, gz)
        }).collect();

        for (gx, gy, gz) in partials {
            for (dst, src) in self.gxcf.iter_mut().zip(gx.iter()) {
                *dst += *src;
            }
            for (dst, src) in self.gycf.iter_mut().zip(gy.iter()) {
                *dst += *src;
            }
            for (dst, src) in self.gzcf.iter_mut().zip(gz.iter()) {
                *dst += *src;
            }
        }

        Ok(())
    }

    /// Membrane boundary conditions
    fn bcfl_mem(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let eps_w = self.pbe.solvent_diel;
        let _t = self.pbe.temperature;
        let xkappa = self.pbe.xkappa;
        let z_mem = self.pbe.z_mem;
        let memv = self.pbe.memv;

        let pre1 = EC / (4.0 * PI * EPS0 * eps_w * KB * _t) * 1.0e10;
        let num_atoms = self.pbe.alist.number_atoms();
        let xlen = 2 * ny * nz;
        let ylen = 2 * nx * nz;
        let zlen = 2 * nx * ny;

        let partials: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..num_atoms).into_par_iter().map(|iatom| {
            let atom = self.pbe.alist.get_atom(iatom);
            let apos = atom.position;
            let charge = EC_ESU * atom.charge;
            let size = atom.radius;
            let mut gx = vec![0.0; xlen];
            let mut gy = vec![0.0; ylen];
            let mut gz = vec![0.0; zlen];

            for k in 0..nz {
                for j in 0..ny {
                    for (face, sign) in [(0usize, 0.0f64), (1, 1.0)] {
                        let xf_val = if face == 0 { self.xf[0] } else { self.xf[nx-1] };
                        let gpos = [xf_val, self.yf[j], self.zf[k]];
                        let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                        let val = if xkappa > VSMALL {
                            pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                                + sign * z_mem * memv
                        } else {
                            pre1 * (charge / dist) + sign * z_mem * memv
                        };
                        gx[ijkx(j, k, face, ny, nz)] += val;
                    }
                }
            }

            for k in 0..nz {
                for i in 0..nx {
                    for (face, _sign) in [(0usize, 0.0f64), (1, 1.0)] {
                        let yf_val = if face == 0 { self.yf[0] } else { self.yf[ny-1] };
                        let gpos = [self.xf[i], yf_val, self.zf[k]];
                        let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                        let val = if xkappa > VSMALL {
                            pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                        } else {
                            pre1 * (charge / dist)
                        };
                        gy[ijky(i, k, face, nx, nz)] += val;
                    }
                }
            }

            for j in 0..ny {
                for i in 0..nx {
                    for (face, _sign) in [(0usize, 0.0f64), (1, 1.0)] {
                        let zf_val = if face == 0 { self.zf[0] } else { self.zf[nz-1] };
                        let gpos = [self.xf[i], self.yf[j], zf_val];
                        let dist = ((gpos[0]-apos[0]).powi(2) + (gpos[1]-apos[1]).powi(2) + (gpos[2]-apos[2]).powi(2)).sqrt();
                        let val = if xkappa > VSMALL {
                            pre1 * (charge / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size)
                        } else {
                            pre1 * (charge / dist)
                        };
                        gz[ijkz(i, j, face, nx, ny)] += val;
                    }
                }
            }

            (gx, gy, gz)
        }).collect();

        for (gx, gy, gz) in partials {
            for (dst, src) in self.gxcf.iter_mut().zip(gx.iter()) {
                *dst += *src;
            }
            for (dst, src) in self.gycf.iter_mut().zip(gy.iter()) {
                *dst += *src;
            }
            for (dst, src) in self.gzcf.iter_mut().zip(gz.iter()) {
                *dst += *src;
            }
        }

        Ok(())
    }

    /// Solve the PBE using the multigrid method
    pub fn solve(&mut self) -> ApbsResult<()> {
        if !self.filled {
            return Err(apbs_generic::error::ApbsError::Solver(
                "fillco must be called before solve".to_string(),
            ));
        }
        if self.pmgp.mgdisc != 0 {
            return Err(ApbsError::UnsupportedFormat(format!(
                "PMGC multigrid discretization mgdisc={} is not implemented in the Rust solver: build_a_fe still aliases FV and PMGC operator/matvec storage only supports the 7-point/4-band path",
                self.pmgp.mgdisc
            )));
        }

        let effective_mgcoar = std::env::var("APBS_RUST_FORCE_MGCOAR")
            .ok()
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or_else(|| self.pmgp.effective_mgcoar());

        // Populate iparm with solver parameters
        self.iparm[0] = self.pmgp.nx;
        self.iparm[1] = self.pmgp.ny;
        self.iparm[2] = self.pmgp.nz;
        self.iparm[3] = self.pmgp.nlev;
        self.iparm[4] = self.pmgp.mgkey;
        self.iparm[5] = self.pmgp.nonlin;
        self.iparm[6] = self.pmgp.istop;
        self.iparm[7] = self.pmgp.itmax;
        self.iparm[8] = self.pmgp.iinfo;
        self.iparm[9] = effective_mgcoar;
        self.iparm[10] = self.pmgp.mgdisc;
        self.iparm[11] = self.pmgp.mgsolv;
        self.iparm[12] = self.pmgp.nu1;
        self.iparm[13] = self.pmgp.nu2;
        self.iparm[14] = self.pmgp.mgsmoo;
        self.iparm[15] = self.pmgp.mgprol;
        self.iparm[16] = self.pmgp.irite;
        self.iparm[17] = self.pmgp.ipcon;
        self.iparm[18] = self.pmgp.key;
        self.iparm[19] = self.pmgp.iperf;

        // Populate rparm
        self.rparm[0] = self.pmgp.omegal;
        self.rparm[1] = self.pmgp.omegan;
        self.rparm[2] = self.pmgp.errtol;
        self.rparm[3] = self.pmgp.hx;
        self.rparm[4] = self.pmgp.hy;
        self.rparm[5] = self.pmgp.hzed;

        // Call the PMG multigrid solver
        apbs_pmgc::mgdriv::mgdriv(
            &mut self.iparm,
            &mut self.rparm,
            &mut self.iwork,
            &mut self.rwork,
            &mut self.u,
            &self.xf,
            &self.yf,
            &self.zf,
            &self.gxcf,
            &self.gycf,
            &self.gzcf,
            &self.a1cf,
            &self.a2cf,
            &self.a3cf,
            &self.ccf,
            &self.fcf,
            &mut self.tcf,
        );

        if Self::apply_post_dirichlet_enabled() {
            self.apply_dirichlet_boundaries();
        }

        Ok(())
    }

    /// Compute total electrostatic energy (kT)
    /// Port of Vpmg_energy in vpmg.c
    pub fn energy(&self, ext_flag: i32) -> f64 {
        if self.pmgp.nonlin != 0 && self.pbe.bulk_ionic_strength > VSMALL {
            // Nonlinear PBE: full free energy decomposition
            let qf = self.qf_energy(ext_flag);
            let qm = self.qm_energy(ext_flag);
            let diel = self.diel_energy(ext_flag);
            qf - diel - qm
        } else {
            // Linear PBE (or zero ionic strength): only 0.5 * q-phi term
            let qf = self.qf_energy(ext_flag);
            0.5 * qf
        }
    }
    
    #[allow(unused_assignments)]
    /// Fixed charge energy (q * phi)
    pub fn qf_energy(&self, ext_flag: i32) -> f64 {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;

        if !self.filled {
            return 0.0;
        }

        // Point-charge energy: interpolate potential at atom positions
        let num_atoms = self.pbe.alist.number_atoms();
        let atoms: Vec<([f64; 3], f64, f64)> = (0..num_atoms)
            .map(|iatom| {
                let atom = self.pbe.alist.get_atom(iatom);
                (
                    atom.position,
                    atom.charge,
                    self.atom_pvec.get(iatom).copied().unwrap_or(1.0),
                )
            })
            .collect();
        let u = &self.u;
        let mut energy: f64 = atoms.par_iter().map(|(apos, charge, atom_weight)| {
            if *atom_weight <= 0.0 {
                return 0.0;
            }

            let ifloat = (apos[0] - xmin) / hx;
            let jfloat = (apos[1] - ymin) / hy;
            let kfloat = (apos[2] - zmin) / hzed;

            let ihi = ifloat.ceil() as i32;
            let ilo = ifloat.floor() as i32;
            let jhi = jfloat.ceil() as i32;
            let jlo = jfloat.floor() as i32;
            let khi = kfloat.ceil() as i32;
            let klo = kfloat.floor() as i32;

            if ihi < nx as i32 && jhi < ny as i32 && khi < nz as i32
                && ilo >= 0 && jlo >= 0 && klo >= 0
            {
                let dx = ifloat - ilo as f64;
                let dy = jfloat - jlo as f64;
                let dz = kfloat - klo as f64;

                let uval = dx*dy*dz*u[ijk(ihi as usize, jhi as usize, khi as usize, nx, ny)]
                    + dx*(1.0-dy)*dz*u[ijk(ihi as usize, jlo as usize, khi as usize, nx, ny)]
                    + dx*dy*(1.0-dz)*u[ijk(ihi as usize, jhi as usize, klo as usize, nx, ny)]
                    + dx*(1.0-dy)*(1.0-dz)*u[ijk(ihi as usize, jlo as usize, klo as usize, nx, ny)]
                    + (1.0-dx)*dy*dz*u[ijk(ilo as usize, jhi as usize, khi as usize, nx, ny)]
                    + (1.0-dx)*(1.0-dy)*dz*u[ijk(ilo as usize, jlo as usize, khi as usize, nx, ny)]
                    + (1.0-dx)*dy*(1.0-dz)*u[ijk(ilo as usize, jhi as usize, klo as usize, nx, ny)]
                    + (1.0-dx)*(1.0-dy)*(1.0-dz)*u[ijk(ilo as usize, jlo as usize, klo as usize, nx, ny)];

                uval * *charge * *atom_weight
            } else {
                0.0
            }
        }).sum();

        if ext_flag == 1 {
            energy += self.ext_qf_energy;
        }

        energy
    }

    /// Mobile ion energy
    #[allow(unused_assignments)]
    pub fn qm_energy(&self, ext_flag: i32) -> f64 {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let zkappa2 = self.pbe.zkappa2;
        let ionstr = self.pbe.bulk_ionic_strength;

        if zkappa2 < VSMALL {
            return 0.0;
        }
        if !self.filled {
            return 0.0;
        }

        let mut energy = 0.0;
        let nf = nx * ny * nz;

        if self.pmgp.nonlin != 0 {
            // Nonlinear (NPBE): energy = sum kappa * zks2 * ci * (exp(-qi*u) - 1)
            // Guard against unphysical blow-ups from outlier grid potentials in
            // energy post-processing; keep solver state untouched.
            let qm_ucap = qm_ucap();
            let zks2 = 0.5 * zkappa2 / ionstr;
            let num_ion = self.pbe.num_ion;
            let pvec = &self.pvec;
            let kappa = &self.kappa;
            let u = &self.u;
            let ion_q = self.pbe.ion_q;
            let ion_conc = self.pbe.ion_conc;
            let (partial_energy, active_pts, max_u, clipped_pts) = (0..nf).into_par_iter().map(|i| {
                if pvec[i] * kappa[i] > VSMALL {
                    let au = u[i].abs();
                    let u_eval = if u[i] > qm_ucap {
                        qm_ucap
                    } else if u[i] < -qm_ucap {
                        -qm_ucap
                    } else {
                        u[i]
                    };
                    let clipped = usize::from(u[i] > qm_ucap || u[i] < -qm_ucap);
                    let mut point_energy = 0.0;
                    for j in 0..num_ion {
                        let mut ichop = false;
                        let expval = vcap_exp(-ion_q[j] * u_eval, &mut ichop);
                        point_energy += pvec[i] * kappa[i] * zks2
                            * ion_conc[j] * (expval - 1.0);
                    }
                    (point_energy, 1usize, au, clipped)
                } else {
                    (0.0, 0usize, 0.0, 0usize)
                }
            }).reduce(
                || (0.0, 0usize, 0.0, 0usize),
                |a, b| (a.0 + b.0, a.1 + b.1, a.2.max(b.2), a.3 + b.3),
            );
            energy = partial_energy;
            if debug_enabled() {
                eprintln!(
                    "[DEBUG-QM] nonlin: zkappa2={:.4e}, ionstr={:.4e}, zks2={:.4e}, active_pts={}/{}, clipped_pts={}, max|u|={:.4e}, ext_flag={}",
                    zkappa2, ionstr, zks2, active_pts, nf, clipped_pts, max_u, ext_flag
                );
            }
        } else {
            // Linear (LPBE): energy = 0.5 * sum kappa * zkappa2 * u^2
            let pvec = &self.pvec;
            let kappa = &self.kappa;
            let u = &self.u;
            let (partial_energy, active_pts, max_u) = (0..nf).into_par_iter().map(|i| {
                if pvec[i] * kappa[i] > VSMALL {
                    let au = u[i].abs();
                    (
                        pvec[i] * zkappa2 * kappa[i] * u[i] * u[i],
                        1usize,
                        au,
                    )
                } else {
                    (0.0, 0usize, 0.0)
                }
            }).reduce(
                || (0.0, 0usize, 0.0),
                |a, b| (a.0 + b.0, a.1 + b.1, a.2.max(b.2)),
            );
            energy = partial_energy;
            energy *= 0.5;
            if debug_enabled() {
                eprintln!(
                    "[DEBUG-QM] linear: zkappa2={:.4e}, active_pts={}/{}, max|u|={:.4e}, ext_flag={}",
                    zkappa2, active_pts, nf, max_u, ext_flag
                );
            }
        }

        energy *= hx * hy * hzed;
        energy /= self.pbe.zmagic;

        if ext_flag == 1 {
            energy += self.ext_qm_energy;
        }

        energy
    }

    /// Polarization (dielectric) energy
    pub fn diel_energy(&self, ext_flag: i32) -> f64 {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;

        if !self.filled {
            return 0.0;
        }

        let pvec = &self.pvec;
        let epsx = &self.epsx;
        let epsy = &self.epsy;
        let epsz = &self.epsz;
        let u = &self.u;
        let mut energy: f64 = (0..(nz-1)).into_par_iter().map(|k| {
            let mut plane_energy = 0.0;
            for j in 0..(ny-1) {
                for i in 0..(nx-1) {
                    let pvecx = 0.5 * (pvec[ijk(i,j,k,nx,ny)] + pvec[ijk(i+1,j,k,nx,ny)]);
                    let pvecy = 0.5 * (pvec[ijk(i,j,k,nx,ny)] + pvec[ijk(i,j+1,k,nx,ny)]);
                    let pvecz = 0.5 * (pvec[ijk(i,j,k,nx,ny)] + pvec[ijk(i,j,k+1,nx,ny)]);

                    let nrgx = epsx[ijk(i,j,k,nx,ny)] * pvecx
                        * ((u[ijk(i,j,k,nx,ny)] - u[ijk(i+1,j,k,nx,ny)]) / hx).powi(2);
                    let nrgy = epsy[ijk(i,j,k,nx,ny)] * pvecy
                        * ((u[ijk(i,j,k,nx,ny)] - u[ijk(i,j+1,k,nx,ny)]) / hy).powi(2);
                    let nrgz = epsz[ijk(i,j,k,nx,ny)] * pvecz
                        * ((u[ijk(i,j,k,nx,ny)] - u[ijk(i,j,k+1,nx,ny)]) / hzed).powi(2);

                    plane_energy += nrgx + nrgy + nrgz;
                }
            }
            plane_energy
        }).sum();

        energy = 0.5 * energy * hx * hy * hzed;
        energy /= self.pbe.zmagic;

        if ext_flag == 1 {
            energy += self.ext_di_energy;
        }

        energy
    }

    /// Compute forces on atoms via analytical differentiation of the discretized PBE.
    /// Port of Vpmg_force from vpmg.c line 5417.
    ///
    /// Returns [fx, fy, fz] in kT/e. Sums all three force components:
    /// - qfForce: charge-field force (depends on chgm)
    /// - dbForce: dielectric boundary force (depends on srfm)
    /// - ibForce: ion boundary force (depends on srfm)
    pub fn force(&self, atom_id: usize, surf_meth: VsurfMeth, charge_meth: VchrgMeth) -> ApbsResult<[f64; 3]> {
        let qf_f = self.qf_force_full(atom_id, charge_meth);
        let db_f = self.db_force(atom_id, surf_meth).unwrap_or([0.0; 3]);
        let ib_f = self.ib_force(atom_id, surf_meth).unwrap_or([0.0; 3]);
        Ok([
            qf_f[0] + db_f[0] + ib_f[0],
            qf_f[1] + db_f[1] + ib_f[1],
            qf_f[2] + db_f[2] + ib_f[2],
        ])
    }

    /// Charge-field force using trilinear interpolation (qfForceSpline1)
    /// Computes F = -q * grad(u) at atom position via finite differences on the grid
    fn qf_force_spline1(&self, atom_id: usize) -> ApbsResult<[f64; 3]> {
        let atom = self.pbe.alist.get_atom(atom_id);
        let apos = atom.position();
        let charge = atom.charge();

        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;

        // Check if atom is on the mesh
        if apos[0] <= xmin || apos[0] >= xmax
            || apos[1] <= ymin || apos[1] >= ymax
            || apos[2] <= zmin || apos[2] >= zmax
        {
            return Ok([0.0; 3]);
        }

        // Convert to grid coordinates
        let position = [apos[0] - xmin, apos[1] - ymin, apos[2] - zmin];
        let ifloat = position[0] / hx;
        let jfloat = position[1] / hy;
        let kfloat = position[2] / hzed;
        let ihi = ifloat.ceil() as usize;
        let ilo = ifloat.floor() as usize;
        let jhi = jfloat.ceil() as usize;
        let jlo = jfloat.floor() as usize;
        let khi = kfloat.ceil() as usize;
        let klo = kfloat.floor() as usize;

        let dx = ifloat - ilo as f64;
        let dy = jfloat - jlo as f64;
        let dz = kfloat - klo as f64;

        let vpmgsmall = 1.0e-6;
        let mut force = [0.0; 3];

        // Force x-component: -q * du/dx
        if dx > vpmgsmall && (1.0 - dx).abs() > vpmgsmall {
            force[0] = -charge * (
                dy * dz * self.u[ijk(ihi, jhi, khi, nx, ny)]
                + dy * (1.0 - dz) * self.u[ijk(ihi, jhi, klo, nx, ny)]
                + (1.0 - dy) * dz * self.u[ijk(ihi, jlo, khi, nx, ny)]
                + (1.0 - dy) * (1.0 - dz) * self.u[ijk(ihi, jlo, klo, nx, ny)]
                - dy * dz * self.u[ijk(ilo, jhi, khi, nx, ny)]
                - dy * (1.0 - dz) * self.u[ijk(ilo, jhi, klo, nx, ny)]
                - (1.0 - dy) * dz * self.u[ijk(ilo, jlo, khi, nx, ny)]
                - (1.0 - dy) * (1.0 - dz) * self.u[ijk(ilo, jlo, klo, nx, ny)]
            ) / hx;
        }

        // Force y-component: -q * du/dy
        if dy > vpmgsmall && (1.0 - dy).abs() > vpmgsmall {
            force[1] = -charge * (
                dx * dz * self.u[ijk(ihi, jhi, khi, nx, ny)]
                + dx * (1.0 - dz) * self.u[ijk(ihi, jhi, klo, nx, ny)]
                - dx * dz * self.u[ijk(ihi, jlo, khi, nx, ny)]
                - dx * (1.0 - dz) * self.u[ijk(ihi, jlo, klo, nx, ny)]
                + (1.0 - dx) * dz * self.u[ijk(ilo, jhi, khi, nx, ny)]
                + (1.0 - dx) * (1.0 - dz) * self.u[ijk(ilo, jhi, klo, nx, ny)]
                - (1.0 - dx) * dz * self.u[ijk(ilo, jlo, khi, nx, ny)]
                - (1.0 - dx) * (1.0 - dz) * self.u[ijk(ilo, jlo, klo, nx, ny)]
            ) / hy;
        }

        // Force z-component: -q * du/dz
        if dz > vpmgsmall && (1.0 - dz).abs() > vpmgsmall {
            force[2] = -charge * (
                dy * dx * self.u[ijk(ihi, jhi, khi, nx, ny)]
                - dy * dx * self.u[ijk(ihi, jhi, klo, nx, ny)]
                + (1.0 - dy) * dx * self.u[ijk(ihi, jlo, khi, nx, ny)]
                - (1.0 - dy) * dx * self.u[ijk(ihi, jlo, klo, nx, ny)]
                + dy * (1.0 - dx) * self.u[ijk(ilo, jhi, khi, nx, ny)]
                - dy * (1.0 - dx) * self.u[ijk(ilo, jhi, klo, nx, ny)]
                + (1.0 - dy) * (1.0 - dx) * self.u[ijk(ilo, jlo, khi, nx, ny)]
                - (1.0 - dy) * (1.0 - dx) * self.u[ijk(ilo, jlo, klo, nx, ny)]
            ) / hzed;
        }

        Ok(force)
    }

    /// Charge-field force using cubic B-spline interpolation (qfForceSpline2)
    /// Uses 4x4x4 stencil with B-spline basis functions
    fn qf_force_bspline2(&self, atom_id: usize) -> ApbsResult<[f64; 3]> {
        let atom = self.pbe.alist.get_atom(atom_id);
        let apos = atom.position();
        let charge = atom.charge();

        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;

        // Check if atom is on the mesh (with 1-cell margin for B-spline stencil)
        if apos[0] <= (xmin + hx) || apos[0] >= (xmax - hx)
            || apos[1] <= (ymin + hy) || apos[1] >= (ymax - hy)
            || apos[2] <= (zmin + hzed) || apos[2] >= (zmax - hzed)
        {
            return Ok([0.0; 3]);
        }

        // Convert to grid coordinates
        let position = [apos[0] - xmin, apos[1] - ymin, apos[2] - zmin];
        let ifloat = position[0] / hx;
        let jfloat = position[1] / hy;
        let kfloat = position[2] / hzed;

        let ii = ifloat.floor() as i32;
        let jj = jfloat.floor() as i32;
        let kk = kfloat.floor() as i32;

        let mx = ifloat - ii as f64;
        let my = jfloat - jj as f64;
        let mz = kfloat - kk as f64;

        // B-spline basis and its derivative
        let (b0, b1, b2, b3) = (
            bspline2(mx + 1.0),
            bspline2(mx),
            bspline2(mx - 1.0),
            bspline2(mx - 2.0),
        );
        let (db0, db1, db2, db3) = (
            dbspline2(mx + 1.0),
            dbspline2(mx),
            dbspline2(mx - 1.0),
            dbspline2(mx - 2.0),
        );
        let (c0, c1, c2, c3) = (
            bspline2(my + 1.0),
            bspline2(my),
            bspline2(my - 1.0),
            bspline2(my - 2.0),
        );
        let (dc0, dc1, dc2, dc3) = (
            dbspline2(my + 1.0),
            dbspline2(my),
            dbspline2(my - 1.0),
            dbspline2(my - 2.0),
        );
        let (d0, d1, d2, d3) = (
            bspline2(mz + 1.0),
            bspline2(mz),
            bspline2(mz - 1.0),
            bspline2(mz - 2.0),
        );
        let (dd0, dd1, dd2, dd3) = (
            dbspline2(mz + 1.0),
            dbspline2(mz),
            dbspline2(mz - 1.0),
            dbspline2(mz - 2.0),
        );

        let bs = [[b0, b1, b2, b3], [c0, c1, c2, c3], [d0, d1, d2, d3]];
        let dbs = [[db0, db1, db2, db3], [dc0, dc1, dc2, dc3], [dd0, dd1, dd2, dd3]];

        let mut force = [0.0; 3];

        for l in 0..3 {
            let mut grad_u = 0.0;
            let mut bi;
            let mut dbi;
            let mut bj;
            let mut dk;

            for i4 in 0..4 {
                bi = bs[0][i4];
                dbi = dbs[0][i4];
                for j4 in 0..4 {
                    bj = bs[1][j4];
                    for k4 in 0..4 {
                        dk = bs[2][k4];
                        let gi = (ii + i4 as i32 - 1) as usize;
                        let gj = (jj + j4 as i32 - 1) as usize;
                        let gk = (kk + k4 as i32 - 1) as usize;

                        if l == 0 {
                            // du/dx
                            grad_u += dbi * bj * dk * self.u[ijk(gi, gj, gk, nx, ny)] / hx;
                        } else if l == 1 {
                            // du/dy
                            grad_u += bi * bj * dk * self.u[ijk(gi, gj, gk, nx, ny)] / hy;
                        } else {
                            // du/dz - need to index dbs for z
                            let dbk = dbs[2][k4];
                            grad_u += bi * bj * dbk * self.u[ijk(gi, gj, gk, nx, ny)] / hzed;
                        }
                    }
                }
            }
            force[l] = -charge * grad_u;
        }

        // Recompute y and z correctly with proper derivative indexing
        // du/dy
        let mut grad_y = 0.0;
        for i4 in 0..4 {
            let bi = bs[0][i4];
            for j4 in 0..4 {
                let dbj = dbs[1][j4];
                for k4 in 0..4 {
                    let dk = bs[2][k4];
                    let gi = (ii + i4 as i32 - 1) as usize;
                    let gj = (jj + j4 as i32 - 1) as usize;
                    let gk = (kk + k4 as i32 - 1) as usize;
                    grad_y += bi * dbj * dk * self.u[ijk(gi, gj, gk, nx, ny)] / hy;
                }
            }
        }
        force[1] = -charge * grad_y;

        // du/dz
        let mut grad_z = 0.0;
        for i4 in 0..4 {
            let bi = bs[0][i4];
            for j4 in 0..4 {
                let bj = bs[1][j4];
                for k4 in 0..4 {
                    let dbk = dbs[2][k4];
                    let gi = (ii + i4 as i32 - 1) as usize;
                    let gj = (jj + j4 as i32 - 1) as usize;
                    let gk = (kk + k4 as i32 - 1) as usize;
                    grad_z += bi * bj * dbk * self.u[ijk(gi, gj, gk, nx, ny)] / hzed;
                }
            }
        }
        force[2] = -charge * grad_z;

        Ok(force)
    }

    /// Set partition for parallel calculation
    pub fn set_part(&mut self, lower: [f64; 3], upper: [f64; 3], bflags: [i32; 6]) {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hz = self.pmgp.hzed;

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = ijk(i, j, k, nx, ny);
                    let x = self.xf[i];
                    let y = self.yf[j];
                    let z = self.zf[k];
                    let xok = partition_grid_weight(x, hx, lower[0], upper[0], bflags[VAPBS_LEFT], bflags[VAPBS_RIGHT]);
                    let yok = partition_grid_weight(y, hy, lower[1], upper[1], bflags[VAPBS_BACK], bflags[VAPBS_FRONT]);
                    let zok = partition_grid_weight(z, hz, lower[2], upper[2], bflags[VAPBS_DOWN], bflags[VAPBS_UP]);
                    self.pvec[idx] = xok * yok * zok;
                }
            }
        }

        self.atom_pvec.clear();
        self.atom_pvec.reserve(self.pbe.alist.number_atoms());
        for atom in self.pbe.alist.atoms() {
            let xok = partition_atom_weight(atom.position[0], lower[0], upper[0], bflags[VAPBS_LEFT], bflags[VAPBS_RIGHT]);
            let yok = partition_atom_weight(atom.position[1], lower[1], upper[1], bflags[VAPBS_BACK], bflags[VAPBS_FRONT]);
            let zok = partition_atom_weight(atom.position[2], lower[2], upper[2], bflags[VAPBS_DOWN], bflags[VAPBS_UP]);
            self.atom_pvec.push(xok * yok * zok);
        }
    }

    /// Unset partition
    pub fn unset_part(&mut self) {
        for v in &mut self.pvec {
            *v = 1.0;
        }
        for v in &mut self.atom_pvec {
            *v = 1.0;
        }
    }

    /// Interpolate potential at (x,y,z) from this grid's solution using trilinear interpolation.
    /// Used by compute_ext_energy to evaluate the coarse grid potential at atom positions.
    // fn interpolate_u(&self, x: f64, y: f64, z: f64) -> f64 {
    //     let nx = self.pmgp.nx as usize;
    //     let ny = self.pmgp.ny as usize;
    //     let nz = self.pmgp.nz as usize;
    //     let hx = self.pmgp.hx;
    //     let hy = self.pmgp.hy;
    //     let hz = self.pmgp.hzed;
    //     let xmin = self.pmgp.xmin;
    //     let ymin = self.pmgp.ymin;
    //     let zmin = self.pmgp.zmin;

    //     let ifloat = (x - xmin) / hx;
    //     let jfloat = (y - ymin) / hy;
    //     let kfloat = (z - zmin) / hz;

    //     let ilo = (ifloat.floor() as usize).min(nx - 1);
    //     let ihi = (ifloat.ceil() as usize).min(nx - 1);
    //     let jlo = (jfloat.floor() as usize).min(ny - 1);
    //     let jhi = (jfloat.ceil() as usize).min(ny - 1);
    //     let klo = (kfloat.floor() as usize).min(nz - 1);
    //     let khi = (kfloat.ceil() as usize).min(nz - 1);

    //     let dx = ifloat - ilo as f64;
    //     let dy = jfloat - jlo as f64;
    //     let dz = kfloat - klo as f64;

    //     (1.0-dx)*(1.0-dy)*(1.0-dz)*self.u[ijk(ilo, jlo, klo, nx, ny)]
    //         + dx*(1.0-dy)*(1.0-dz)*self.u[ijk(ihi, jlo, klo, nx, ny)]
    //         + (1.0-dx)*dy*(1.0-dz)*self.u[ijk(ilo, jhi, klo, nx, ny)]
    //         + dx*dy*(1.0-dz)*self.u[ijk(ihi, jhi, klo, nx, ny)]
    //         + (1.0-dx)*(1.0-dy)*dz*self.u[ijk(ilo, jlo, khi, nx, ny)]
    //         + dx*(1.0-dy)*dz*self.u[ijk(ihi, jlo, khi, nx, ny)]
    //         + (1.0-dx)*dy*dz*self.u[ijk(ilo, jhi, khi, nx, ny)]
    //         + dx*dy*dz*self.u[ijk(ihi, jhi, khi, nx, ny)]
    // }

    /// Compute external energy contributions from region outside focusing domain.
    /// Port of extEnergy from vpmg.c
    ///
    /// When focusing, the total energy = fine grid energy + external energy from
    /// the coarse grid region outside the fine grid. This method computes the
    /// external contribution and stores it in ext_qf_energy, ext_qm_energy, ext_di_energy.
    pub fn compute_ext_energy(&mut self, old: &mut Vpmg) {
        // Zero ext fields
        self.ext_qf_energy = 0.0;
        self.ext_qm_energy = 0.0;
        self.ext_di_energy = 0.0;

        // New grid bounding box (physical coordinates)
        let nx_new = self.pmgp.nx as usize;
        let ny_new = self.pmgp.ny as usize;
        let nz_new = self.pmgp.nz as usize;
        let lower_corner = [
            self.pmgp.xcent - ((nx_new - 1) as f64 * self.pmgp.hx) / 2.0,
            self.pmgp.ycent - ((ny_new - 1) as f64 * self.pmgp.hy) / 2.0,
            self.pmgp.zcent - ((nz_new - 1) as f64 * self.pmgp.hzed) / 2.0,
        ];
        let upper_corner = [
            self.pmgp.xcent + ((nx_new - 1) as f64 * self.pmgp.hx) / 2.0,
            self.pmgp.ycent + ((ny_new - 1) as f64 * self.pmgp.hy) / 2.0,
            self.pmgp.zcent + ((nz_new - 1) as f64 * self.pmgp.hzed) / 2.0,
        ];

        // Set partition on old grid: marks points inside new grid as 1, outside as 0
        let bflags = [0i32; 6];
        old.set_part(lower_corner, upper_corner, bflags);
        if debug_enabled() {
            let inside_grid = old.pvec.iter().filter(|v| **v > VSMALL).count();
            let inside_atom = old.atom_pvec.iter().filter(|v| **v > VSMALL).count();
            eprintln!(
                "[DEBUG-EXT] set_part inside: grid={}/{}, atom={}/{}",
                inside_grid,
                old.pvec.len(),
                inside_atom,
                old.atom_pvec.len()
            );
        }

        // Invert partition: pvec = 1 - pvec (now outside=1, inside=0)
        let nf_old = (old.pmgp.nx * old.pmgp.ny * old.pmgp.nz) as usize;
        for i in 0..nf_old {
            if old.pvec[i] > VSMALL {
                old.pvec[i] = 1.0;
            }
            old.pvec[i] = 1.0 - old.pvec[i];
        }
        for weight in &mut old.atom_pvec {
            if *weight > VSMALL {
                *weight = 1.0;
            }
            *weight = 1.0 - *weight;
        }
        if debug_enabled() {
            let outside_grid = old.pvec.iter().filter(|v| **v > VSMALL).count();
            let outside_atom = old.atom_pvec.iter().filter(|v| **v > VSMALL).count();
            let mut active_qm = 0usize;
            let mut max_u = 0.0f64;
            for i in 0..old.pvec.len() {
                if old.pvec[i] * old.kappa[i] > VSMALL {
                    active_qm += 1;
                    let au = old.u[i].abs();
                    if au > max_u {
                        max_u = au;
                    }
                }
            }
            eprintln!(
                "[DEBUG-EXT] inverted outside: grid={}/{}, atom={}/{}, qm_active={}, max|u|={:.4e}",
                outside_grid,
                old.pvec.len(),
                outside_atom,
                old.atom_pvec.len(),
                active_qm,
                max_u
            );
        }

        // Compute external energies on the inverted coarse-grid partition.
        // Use a clean ext baseline on the coarse object to avoid accidental
        // recursive carry-over when multiple focusing levels are chained.
        let old_ext_qm = old.ext_qm_energy;
        let old_ext_qf = old.ext_qf_energy;
        let old_ext_di = old.ext_di_energy;
        old.ext_qm_energy = 0.0;
        old.ext_qf_energy = 0.0;
        old.ext_di_energy = 0.0;
        self.ext_qm_energy = old.qm_energy(1);
        self.ext_qf_energy = old.qf_energy(1);
        self.ext_di_energy = old.diel_energy(1);
        old.ext_qm_energy = old_ext_qm;
        old.ext_qf_energy = old_ext_qf;
        old.ext_di_energy = old_ext_di;

        // Restore old grid's partition
        old.unset_part();

        if debug_enabled() {
            eprintln!("[DEBUG-EXT] ext_energy: ext_qf={:.6e}, ext_qm={:.6e}, ext_di={:.6e}",
                self.ext_qf_energy, self.ext_qm_energy, self.ext_di_energy);
        }
    }

    /// Fill output array with specified data type
    pub fn fill_array(&self, vec: &mut [f64], data_type: VdataType) -> ApbsResult<()> {
        let nf = self.pmgp.nf as usize;
        if vec.len() < nf {
            return Err(ApbsError::InvalidParameter(format!(
                "output array too small: {} < {}",
                vec.len(), nf
            )));
        }
        match data_type {
            VdataType::Pot => {
                vec[..nf].copy_from_slice(&self.u[..nf]);
            }
            VdataType::AtomPot => {
                let nx = self.pmgp.nx as usize;
                let ny = self.pmgp.ny as usize;
                let nz = self.pmgp.nz as usize;
                let hx = self.pmgp.hx;
                let hy = self.pmgp.hy;
                let hzed = self.pmgp.hzed;
                let xmin = self.pmgp.xmin;
                let ymin = self.pmgp.ymin;
                let zmin = self.pmgp.zmin;
                let atoms = self.pbe.alist.number_atoms();
                let positions: Vec<[f64; 3]> = (0..atoms)
                    .map(|iatom| self.pbe.alist.get_atom(iatom).position)
                    .collect();
                let u = &self.u;
                vec[..nf].fill(0.0);
                vec[..positions.len().min(nf)].par_iter_mut().enumerate().for_each(|(iatom, out)| {
                    let apos = positions[iatom];
                    let ifloat = (apos[0] - xmin) / hx;
                    let jfloat = (apos[1] - ymin) / hy;
                    let kfloat = (apos[2] - zmin) / hzed;
                    let ilo = (ifloat.floor() as usize).min(nx - 1);
                    let ihi = (ifloat.ceil() as usize).min(nx - 1);
                    let jlo = (jfloat.floor() as usize).min(ny - 1);
                    let jhi = (jfloat.ceil() as usize).min(ny - 1);
                    let klo = (kfloat.floor() as usize).min(nz - 1);
                    let khi = (kfloat.ceil() as usize).min(nz - 1);
                    let dx = ifloat - ilo as f64;
                    let dy = jfloat - jlo as f64;
                    let dz = kfloat - klo as f64;
                    *out = dx*dy*dz*u[ijk(ihi, jhi, khi, nx, ny)]
                        + dx*(1.0-dy)*dz*u[ijk(ihi, jlo, khi, nx, ny)]
                        + dx*dy*(1.0-dz)*u[ijk(ihi, jhi, klo, nx, ny)]
                        + dx*(1.0-dy)*(1.0-dz)*u[ijk(ihi, jlo, klo, nx, ny)]
                        + (1.0-dx)*dy*dz*u[ijk(ilo, jhi, khi, nx, ny)]
                        + (1.0-dx)*(1.0-dy)*dz*u[ijk(ilo, jlo, khi, nx, ny)]
                        + (1.0-dx)*dy*(1.0-dz)*u[ijk(ilo, jhi, klo, nx, ny)]
                        + (1.0-dx)*(1.0-dy)*(1.0-dz)*u[ijk(ilo, jlo, klo, nx, ny)];
                });
            }
            VdataType::Charge => {
                let zmagic = self.pbe.zmagic;
                let charge = &self.charge;
                vec[..nf].par_iter_mut().enumerate().for_each(|(i, out)| *out = charge[i] / zmagic);
            }
            VdataType::Smol => {
                vec[..nf].copy_from_slice(&self.pvec[..nf]);
            }
            VdataType::Sspl => {
                vec[..nf].copy_from_slice(&self.pvec[..nf]);
            }
            VdataType::Vdw => {
                let pvec = &self.pvec;
                vec[..nf].par_iter_mut().enumerate().for_each(|(i, out)| *out = if pvec[i] > VSMALL { 1.0 } else { 0.0 });
            }
            VdataType::Ivdw => {
                let pvec = &self.pvec;
                vec[..nf].par_iter_mut().enumerate().for_each(|(i, out)| *out = if pvec[i] > VSMALL { 0.0 } else { 1.0 });
            }
            VdataType::Lap => {
                let nx = self.pmgp.nx as usize;
                let ny = self.pmgp.ny as usize;
                let nz = self.pmgp.nz as usize;
                let hx2 = self.pmgp.hx * self.pmgp.hx;
                let hy2 = self.pmgp.hy * self.pmgp.hy;
                let hz2 = self.pmgp.hzed * self.pmgp.hzed;
                let u = &self.u;
                vec[..nf].par_iter_mut().enumerate().for_each(|(ip, out)| {
                    let k = ip / (nx * ny);
                    let j = (ip - k * nx * ny) / nx;
                    let i = ip - k * nx * ny - j * nx;
                            let ip = ijk(i, j, k, nx, ny);
                            let c = u[ip];
                            let xm = if i > 0 { u[ijk(i - 1, j, k, nx, ny)] } else { c };
                            let xp = if i + 1 < nx { u[ijk(i + 1, j, k, nx, ny)] } else { c };
                            let ym = if j > 0 { u[ijk(i, j - 1, k, nx, ny)] } else { c };
                            let yp = if j + 1 < ny { u[ijk(i, j + 1, k, nx, ny)] } else { c };
                            let zm = if k > 0 { u[ijk(i, j, k - 1, nx, ny)] } else { c };
                            let zp = if k + 1 < nz { u[ijk(i, j, k + 1, nx, ny)] } else { c };
                            *out = (xp - 2.0 * c + xm) / hx2
                                + (yp - 2.0 * c + ym) / hy2
                                + (zp - 2.0 * c + zm) / hz2;
                });
            }
            VdataType::Edens => {
                let nx = self.pmgp.nx as usize;
                let ny = self.pmgp.ny as usize;
                let nz = self.pmgp.nz as usize;
                let hx = self.pmgp.hx;
                let hy = self.pmgp.hy;
                let hz = self.pmgp.hzed;
                let u = &self.u;
                let epsx = &self.epsx;
                let epsy = &self.epsy;
                let epsz = &self.epsz;
                vec[..nf].par_iter_mut().enumerate().for_each(|(ip, out)| {
                    let k = ip / (nx * ny);
                    let j = (ip - k * nx * ny) / nx;
                    let i = ip - k * nx * ny - j * nx;
                            let ip = ijk(i, j, k, nx, ny);
                            let c = u[ip];
                            let xm = if i > 0 { u[ijk(i - 1, j, k, nx, ny)] } else { c };
                            let xp = if i + 1 < nx { u[ijk(i + 1, j, k, nx, ny)] } else { c };
                            let ym = if j > 0 { u[ijk(i, j - 1, k, nx, ny)] } else { c };
                            let yp = if j + 1 < ny { u[ijk(i, j + 1, k, nx, ny)] } else { c };
                            let zm = if k > 0 { u[ijk(i, j, k - 1, nx, ny)] } else { c };
                            let zp = if k + 1 < nz { u[ijk(i, j, k + 1, nx, ny)] } else { c };
                            let gx = (xp - xm) / (2.0 * hx);
                            let gy = (yp - ym) / (2.0 * hy);
                            let gz = (zp - zm) / (2.0 * hz);
                            let eps = (epsx[ip] + epsy[ip] + epsz[ip]) / 3.0;
                            *out = 0.5 * eps * (gx * gx + gy * gy + gz * gz);
                });
            }
            VdataType::Ndens => {
                let kappa = &self.kappa;
                let u = &self.u;
                vec[..nf].par_iter_mut().enumerate().for_each(|(i, out)| *out = kappa[i] * u[i].sinh());
            }
            VdataType::Qdens => {
                vec[..nf].copy_from_slice(&self.charge[..nf]);
            }
            VdataType::DielX => {
                vec[..nf].copy_from_slice(&self.epsx[..nf]);
            }
            VdataType::DielY => {
                vec[..nf].copy_from_slice(&self.epsy[..nf]);
            }
            VdataType::DielZ => {
                vec[..nf].copy_from_slice(&self.epsz[..nf]);
            }
            VdataType::Kappa => {
                vec[..nf].copy_from_slice(&self.kappa[..nf]);
            }
        }
        Ok(())
    }
}

/// Multipole boundary condition factors
/// Returns [monopole, dipole, quadrupole] tensor components
fn multipolebc(r: f64, kappa: f64, eps_p: f64, eps_w: f64, rad: f64) -> [f64; 3] {
    let r2 = r * r;
    let r3 = r2 * r;
    let r5 = r3 * r2;
    let eps_r = eps_w / eps_p;

    let mut tsr = [0.0; 3];

    if kappa < VSMALL {
        // No screening
        tsr[0] = (1.0 / eps_w) / r;
        tsr[1] = 3.0 * eps_r / (1.0 + 2.0 * eps_r) * (-1.0 / eps_w) / r3;
        tsr[2] = 5.0 * eps_r / (2.0 + 3.0 * eps_r) * (3.0 / eps_w) / r5;
    } else {
        // Debye-Huckel screening
        let ka = kappa * rad;
        let ka2 = ka * ka;
        let ka3 = ka2 * ka;
        let kr = kappa * r;
        let kr2 = kr * kr;
        let _kr3 = kr2 * kr;

        tsr[0] = (-(ka - kr)).exp() / (1.0 + ka) * (1.0 / eps_w) / r;
        tsr[1] = 3.0 * eps_r * (-(ka - kr)).exp() * (1.0 + kr)
            / (1.0 + ka + eps_r * (2.0 + 2.0 * ka + ka2))
            * (-1.0 / eps_w) / r3;
        tsr[2] = 5.0 * eps_r * (-(ka - kr)).exp() * (3.0 + 3.0 * kr + kr2)
            / (6.0 + 6.0 * ka + 2.0 * ka2 + eps_r * (9.0 + 9.0 * ka + 4.0 * ka2 + ka3))
            * (3.0 / eps_w) / r5;
    }

    tsr
}

fn partition_atom_weight(value: f64, lower: f64, upper: f64, lower_flag: i32, upper_flag: i32) -> f64 {
    let tol = 1.0e-12;
    if value < upper && value > lower {
        1.0
    } else if (value - lower).abs() < tol {
        if lower_flag == 0 { 1.0 } else { 0.5 }
    } else if (value - upper).abs() < tol {
        if upper_flag == 0 { 1.0 } else { 0.5 }
    } else {
        0.0
    }
}

fn partition_grid_weight(
    value: f64,
    h: f64,
    lower: f64,
    upper: f64,
    lower_flag: i32,
    upper_flag: i32,
) -> f64 {
    let tol = 1.0e-12;
    if value < (upper - h / 2.0) && value > (lower + h / 2.0) {
        1.0
    } else if (value - lower).abs() < tol {
        if lower_flag == 0 { 1.0 } else { 0.5 }
    } else if (value - upper).abs() < tol {
        if upper_flag == 0 { 1.0 } else { 0.5 }
    } else if value > (upper + h / 2.0) || value < (lower - h / 2.0) {
        0.0
    } else {
        let v0 = (value - h / 2.0).max(lower);
        let v1 = (value + h / 2.0).min(upper);
        ((v1 - v0).abs() / h).clamp(0.0, 1.0)
    }
}

/// Cubic B-spline basis function
fn bspline2(x: f64) -> f64 {
    let m2m = if x >= 0.0 && x <= 2.0 { 1.0 - (x - 1.0).abs() } else { 0.0 };
    let m2 = if x >= 1.0 && x <= 3.0 { 1.0 - (x - 2.0).abs() } else { 0.0 };
    if x >= 0.0 && x <= 3.0 {
        0.5 * x * m2m + 0.5 * (3.0 - x) * m2
    } else {
        0.0
    }
}

/// Derivative of cubic B-spline basis function
fn dbspline2(x: f64) -> f64 {
    if x >= 0.0 && x < 1.0 {
        x
    } else if x >= 1.0 && x < 2.0 {
        3.0 - 2.0 * x
    } else if x >= 2.0 && x <= 3.0 {
        x - 3.0
    } else {
        0.0
    }
}

// ========================================================================
// Ported from vpmg.c: B-spline utility functions and missing methods
// ========================================================================

/// 4th-order B-spline basis function.
/// Port of bspline4 from vpmg.c line 6731
fn bspline4(x: f64) -> f64 {
    let one6 = 1.0 / 6.0;
    let one8 = 1.0 / 8.0;
    let one24 = 1.0 / 24.0;
    let thirteen24 = 13.0 / 24.0;
    let fourtyseven24 = 47.0 / 24.0;
    let seventeen24 = 17.0 / 24.0;

    if x > 0.0 && x <= 1.0 {
        let m = x * x;
        one24 * m * m
    } else if x > 1.0 && x <= 2.0 {
        let m = x - 1.0;
        let m2 = m * m;
        -one8 + one6 * x + m2 * (0.25 + one6 * m - one6 * m2)
    } else if x > 2.0 && x <= 3.0 {
        let m = x - 2.0;
        let m2 = m * m;
        -thirteen24 + 0.5 * x + m2 * (-0.25 - 0.5 * m + 0.25 * m2)
    } else if x > 3.0 && x <= 4.0 {
        let m = x - 3.0;
        let m2 = m * m;
        fourtyseven24 - 0.5 * x + m2 * (-0.25 + 0.5 * m - one6 * m2)
    } else if x > 4.0 && x <= 5.0 {
        let m = x - 4.0;
        let m2 = m * m;
        seventeen24 - one6 * x + m2 * (0.25 - one6 * m + one24 * m2)
    } else {
        0.0
    }
}

/// 1st derivative of 4th-order B-spline.
/// Port of dbspline4 from vpmg.c line 6765
fn dbspline4(x: f64) -> f64 {
    let one6 = 1.0 / 6.0;
    let one3 = 1.0 / 3.0;
    let two3 = 2.0 / 3.0;
    let thirteen6 = 13.0 / 6.0;

    if x > 0.0 && x <= 1.0 {
        let m2 = x * x;
        one6 * x * m2
    } else if x > 1.0 && x <= 2.0 {
        let m = x - 1.0;
        let m2 = m * m;
        -one3 + 0.5 * x + m2 * (0.5 - two3 * m)
    } else if x > 2.0 && x <= 3.0 {
        let m = x - 2.0;
        let m2 = m * m;
        1.5 - 0.5 * x + m2 * (-1.5 + m)
    } else if x > 3.0 && x <= 4.0 {
        let m = x - 3.0;
        let m2 = m * m;
        1.0 - 0.5 * x + m2 * (1.5 - two3 * m)
    } else if x > 4.0 && x <= 5.0 {
        let m = x - 4.0;
        let m2 = m * m;
        -thirteen6 + 0.5 * x + m2 * (-0.5 + one6 * m)
    } else {
        0.0
    }
}

/// 2nd derivative of 4th-order B-spline.
/// Port of d2bspline4 from vpmg.c line 6797
fn d2bspline4(x: f64) -> f64 {
    if x > 0.0 && x <= 1.0 {
        0.5 * x * x
    } else if x > 1.0 && x <= 2.0 {
        let m = x - 1.0;
        let m2 = m * m;
        -0.5 + x - 2.0 * m2
    } else if x > 2.0 && x <= 3.0 {
        let m = x - 2.0;
        let m2 = m * m;
        5.5 - 3.0 * x + 3.0 * m2
    } else if x > 3.0 && x <= 4.0 {
        let m = x - 3.0;
        let m2 = m * m;
        -9.5 + 3.0 * x - 2.0 * m2
    } else if x > 4.0 && x <= 5.0 {
        let m = x - 4.0;
        let m2 = m * m;
        4.5 - x + 0.5 * m2
    } else {
        0.0
    }
}

/// 3rd derivative of 4th-order B-spline.
/// Port of d3bspline4 from vpmg.c line 6824
// fn d3bspline4(x: f64) -> f64 {
//     if x > 0.0 && x <= 1.0 {
//         x
//     } else if x > 1.0 && x <= 2.0 {
//         5.0 - 4.0 * x
//     } else if x > 2.0 && x <= 3.0 {
//         -15.0 + 6.0 * x
//     } else if x > 3.0 && x <= 4.0 {
//         15.0 - 4.0 * x
//     } else if x > 4.0 && x <= 5.0 {
//         x - 5.0
//     } else {
//         0.0
//     }
// }

/// VFCHI4 helper for 4th-order B-spline.
/// Port of VFCHI4 from vpmg.c line 6727
fn vfchi4(i: f64, f: f64) -> f64 {
    2.5 + (i - f)
}

impl Vpmg {
    /// Memory usage check.
    /// Port of Vpmg_memChk from vpmg.c line 79
    pub fn mem_chk(&self) -> usize {
        // In Rust we don't track Vmem bytes; return approximate size
        // let nf = self.pmgp.nf as usize;
        // let narr = self.pmgp.narr as usize;
        // Approximate: main arrays (each nf or narr f64s)
        let arrays = (self.a1cf.len() + self.a2cf.len() + self.a3cf.len()
            + self.ccf.len() + self.fcf.len() + self.tcf.len()
            + self.u.len() + self.pvec.len() + self.kappa.len()
            + self.epsx.len() + self.epsy.len() + self.epsz.len()) * 8;
        let grids = (self.gxcf.len() + self.gycf.len() + self.gzcf.len()) * 8;
        let coords = (self.xf.len() + self.yf.len() + self.zf.len()) * 8;
        arrays + grids + coords
    }

    /// Dielectric gradient norm energy.
    /// Port of Vpmg_dielGradNorm from vpmg.c line 1342
    pub fn diel_grad_norm(&self) -> f64 {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;

        assert!(self.filled, "Need to call fillco first");

        let mut energy = 0.0;
        for k in 1..nz {
            for j in 1..ny {
                for i in 1..nx {
                    let pvecx = 0.5 * (self.pvec[ijk(i, j, k, nx, ny)]
                        + self.pvec[ijk(i - 1, j, k, nx, ny)]);
                    let pvecy = 0.5 * (self.pvec[ijk(i, j, k, nx, ny)]
                        + self.pvec[ijk(i, j - 1, k, nx, ny)]);
                    let pvecz = 0.5 * (self.pvec[ijk(i, j, k, nx, ny)]
                        + self.pvec[ijk(i, j, k - 1, nx, ny)]);
                    let nrgx = pvecx
                        * ((self.epsx[ijk(i, j, k, nx, ny)]
                            - self.epsx[ijk(i - 1, j, k, nx, ny)])
                            / hx)
                            .powi(2);
                    let nrgy = pvecy
                        * ((self.epsy[ijk(i, j, k, nx, ny)]
                            - self.epsy[ijk(i, j - 1, k, nx, ny)])
                            / hy)
                            .powi(2);
                    let nrgz = pvecz
                        * ((self.epsz[ijk(i, j, k, nx, ny)]
                            - self.epsz[ijk(i, j, k - 1, nx, ny)])
                            / hzed)
                            .powi(2);
                    energy += (nrgx + nrgy + nrgz).sqrt();
                }
            }
        }
        energy * hx * hy * hzed
    }

    /// Per-atom fixed charge energy.
    /// Port of Vpmg_qfAtomEnergy from vpmg.c line 1791
    pub fn qf_atom_energy(&self, atom_id: usize) -> f64 {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;

        let atom = self.pbe.alist.get_atom(atom_id);
        let position = atom.position;
        let charge = atom.charge;
        let atom_weight = self.atom_pvec.get(atom_id).copied().unwrap_or(1.0);

        if atom_weight <= 0.0 {
            return 0.0;
        }

        let ifloat = (position[0] - xmin) / hx;
        let jfloat = (position[1] - ymin) / hy;
        let kfloat = (position[2] - zmin) / hzed;
        let ihi = ifloat.ceil() as i32;
        let ilo = ifloat.floor() as i32;
        let jhi = jfloat.ceil() as i32;
        let jlo = jfloat.floor() as i32;
        let khi = kfloat.ceil() as i32;
        let klo = kfloat.floor() as i32;

        if ihi < nx as i32 && jhi < ny as i32 && khi < nz as i32
            && ilo >= 0 && jlo >= 0 && klo >= 0
        {
            let dx = ifloat - ilo as f64;
            let dy = jfloat - jlo as f64;
            let dz = kfloat - klo as f64;
            let uval = dx * dy * dz * self.u[ijk(ihi as usize, jhi as usize, khi as usize, nx, ny)]
                + dx * (1.0 - dy) * dz * self.u[ijk(ihi as usize, jlo as usize, khi as usize, nx, ny)]
                + dx * dy * (1.0 - dz) * self.u[ijk(ihi as usize, jhi as usize, klo as usize, nx, ny)]
                + dx * (1.0 - dy) * (1.0 - dz) * self.u[ijk(ihi as usize, jlo as usize, klo as usize, nx, ny)]
                + (1.0 - dx) * dy * dz * self.u[ijk(ilo as usize, jhi as usize, khi as usize, nx, ny)]
                + (1.0 - dx) * (1.0 - dy) * dz * self.u[ijk(ilo as usize, jlo as usize, khi as usize, nx, ny)]
                + (1.0 - dx) * dy * (1.0 - dz) * self.u[ijk(ilo as usize, jhi as usize, klo as usize, nx, ny)]
                + (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
                    * self.u[ijk(ilo as usize, jlo as usize, klo as usize, nx, ny)];
            uval * charge * atom_weight
        } else {
            0.0
        }
    }

    /// Solve Laplace equation (zero ionic strength).
    /// Port of Vpmg_solveLaplace from vpmg.c line 6637
    pub fn solve_laplace(&mut self) -> ApbsResult<()> {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let epsw = self.pbe.solvent_diel;
        let iepsw = 1.0 / epsw;
        let scal = hx * hy * hzed;
        let scalx = hx * hy / hzed;
        let scaly = hx * hzed / hy;
        let scalz = hx * hy / hzed;

        assert!(self.filled, "Need to call fillco first");

        // Load boundary conditions into RHS array
        for i in 1..(nx - 1) {
            let dilo = if i == 1 { 1.0 } else { 0.0 };
            let dihi = if i == nx - 2 { 1.0 } else { 0.0 };
            for j in 1..(ny - 1) {
                let djlo = if j == 1 { 1.0 } else { 0.0 };
                let djhi = if j == ny - 2 { 1.0 } else { 0.0 };
                for k in 1..(nz - 1) {
                    let dklo = if k == 1 { 1.0 } else { 0.0 };
                    let dkhi = if k == nz - 2 { 1.0 } else { 0.0 };
                    self.fcf[ijk(i, j, k, nx, ny)] = iepsw * scal * self.charge[ijk(i, j, k, nx, ny)]
                        + dilo * scalx * self.gxcf[ijkx(j, k, 0, ny, nz)]
                        + dihi * scalx * self.gxcf[ijkx(j, k, 1, ny, nz)]
                        + djlo * scaly * self.gycf[ijky(i, k, 0, nx, nz)]
                        + djhi * scaly * self.gycf[ijky(i, k, 1, nx, nz)]
                        + dklo * scalz * self.gzcf[ijkz(i, j, 0, nx, ny)]
                        + dkhi * scalz * self.gzcf[ijkz(i, j, 1, nx, ny)];
                }
            }
        }

        // Solve using multigrid (same as regular solve but with uniform dielectric)
        // For now, use the same PMG solver
        let nf = nx * ny * nz;
        self.tcf[..nf].copy_from_slice(&self.fcf[..nf]);
        if self.pmgp.mgdisc != 0 {
            return Err(ApbsError::UnsupportedFormat(format!(
                "SMPBE PMGC multigrid discretization mgdisc={} is not implemented in the Rust solver: only the 7-point/4-band path is currently valid",
                self.pmgp.mgdisc
            )));
        }

        apbs_pmgc::mgdriv::mgdriv(
            &mut self.iparm,
            &mut self.rparm,
            &mut self.iwork,
            &mut self.rwork,
            &mut self.u,
            &self.xf,
            &self.yf,
            &self.zf,
            &self.gxcf,
            &self.gycf,
            &self.gzcf,
            &self.a1cf,
            &self.a2cf,
            &self.a3cf,
            &self.ccf,
            &self.fcf,
            &mut self.tcf,
        );
        self.u[..nf].copy_from_slice(&self.tcf[..nf]);

        // Add boundary conditions to solution
        for j in 0..ny {
            for k in 0..nz {
                self.u[ijk(0, j, k, nx, ny)] = self.gxcf[ijkx(j, k, 0, ny, nz)];
                self.u[ijk(nx - 1, j, k, nx, ny)] = self.gxcf[ijkx(j, k, 1, ny, nz)];
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                self.u[ijk(i, 0, k, nx, ny)] = self.gycf[ijky(i, k, 0, nx, nz)];
                self.u[ijk(i, ny - 1, k, nx, ny)] = self.gycf[ijky(i, k, 1, nx, nz)];
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                self.u[ijk(i, j, 0, nx, ny)] = self.gzcf[ijkz(i, j, 0, nx, ny)];
                self.u[ijk(i, j, nz - 1, nx, ny)] = self.gzcf[ijkz(i, j, 1, nx, ny)];
            }
        }

        Ok(())
    }

    /// Full qfForce dispatch including bspline4.
    /// Port of Vpmg_qfForce from vpmg.c line 5862
    pub fn qf_force_full(&self, atom_id: usize, chgm: VchrgMeth) -> [f64; 3] {
        match chgm {
            VchrgMeth::Tril => {
                self.qf_force_spline1(atom_id).unwrap_or([0.0; 3])
            }
            VchrgMeth::Bspl2 => {
                self.qf_force_bspline2(atom_id).unwrap_or([0.0; 3])
            }
            VchrgMeth::Bspl4 => {
                let mut force = [0.0; 3];
                self.qf_force_bspline4_into(atom_id, &mut force);
                force
            }
        }
    }

    /// 4th-order B-spline charge-field force.
    /// Port of qfForceSpline4 from vpmg.c
    fn qf_force_bspline4_into(&self, atom_id: usize, force: &mut [f64; 3]) {
        let atom = self.pbe.alist.get_atom(atom_id);
        let apos = atom.position;
        let charge = atom.charge;

        force[0] = 0.0;
        force[1] = 0.0;
        force[2] = 0.0;

        if atom.part_id <= 0.0 {
            return;
        }

        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;

        if apos[0] <= xmin || apos[0] >= xmax
            || apos[1] <= ymin || apos[1] >= ymax
            || apos[2] <= zmin || apos[2] >= zmax
        {
            return;
        }

        // Grid coordinate of atom
        let xf = (apos[0] - xmin) / hx;
        let yf = (apos[1] - ymin) / hy;
        let zf = (apos[2] - zmin) / hzed;

        // Support radius = 2 for 4th-order B-spline
        let supp = 2.0;
        let imin = (xf - supp).ceil().max(1.0) as i32;
        let imax = (xf + supp).floor().min((nx - 2) as f64) as i32;
        let jmin = (yf - supp).ceil().max(1.0) as i32;
        let jmax = (yf + supp).floor().min((ny - 2) as f64) as i32;
        let kmin = (zf - supp).ceil().max(1.0) as i32;
        let kmax = (zf + supp).floor().min((nz - 2) as f64) as i32;

        for i in imin..=imax {
            let wtx = dbspline4(vfchi4(i as f64, xf));
            let dwtx = d2bspline4(vfchi4(i as f64, xf)) / hx;
            for j in jmin..=jmax {
                let wty = dbspline4(vfchi4(j as f64, yf));
                let dwty = d2bspline4(vfchi4(j as f64, yf)) / hy;
                for k in kmin..=kmax {
                    let wtz = dbspline4(vfchi4(k as f64, zf));
                    let dwtz = d2bspline4(vfchi4(k as f64, zf)) / hzed;
                    let uval = self.u[ijk(i as usize, j as usize, k as usize, nx, ny)];
                    force[0] -= charge * dwtx * wty * wtz * uval;
                    force[1] -= charge * wtx * dwty * wtz * uval;
                    force[2] -= charge * wtx * wty * dwtz * uval;
                }
            }
        }
    }

    /// Ion boundary force.
    /// Port of Vpmg_ibForce from vpmg.c line 5440
    pub fn ib_force(&self, atom_id: usize, srfm: VsurfMeth) -> ApbsResult<[f64; 3]> {
        let mut force = [0.0f64; 3];

        // Only spline-based surfaces supported
        match srfm {
            VsurfMeth::Spline | VsurfMeth::Spline3 | VsurfMeth::Spline4 => {}
            _ => return Ok(force),
        }

        let atom = self.pbe.alist.get_atom(atom_id);
        if atom.part_id == 0.0 {
            return Ok(force);
        }

        let apos = atom.position;
        let arad = atom.radius;

        let zkappa2 = self.pbe.zkappa2;
        if zkappa2 < VPMGSMALL {
            return Ok(force);
        }

        let irad = self.pbe.max_ion_radius;
        let izmagic = 1.0 / self.pbe.zmagic;
        let ionstr = self.pbe.bulk_ionic_strength;
        let nion = self.pbe.num_ion;

        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;

        if apos[0] <= xmin || apos[0] >= xmax
            || apos[1] <= ymin || apos[1] >= ymax
            || apos[2] <= zmin || apos[2] >= zmax
        {
            return Ok(force);
        }

        let position = [
            apos[0] - xmin,
            apos[1] - ymin,
            apos[2] - zmin,
        ];

        let rtot = irad + arad + self.spline_win;
        let rtot2 = rtot * rtot;
        let dx = rtot + 0.5 * hx;
        let imin = (0.0f64).max(((position[0] - dx) / hx).ceil()) as i32;
        let imax = ((nx as f64) - 1.0).min(((position[0] + dx) / hx).floor()) as i32;

        for i in imin..=imax {
            let dx2 = (position[0] - hx * i as f64).powi(2);
            let dy = if rtot2 > dx2 {
                (rtot2 - dx2).sqrt() + 0.5 * hy
            } else {
                0.5 * hy
            };
            let jmin = (0.0f64).max(((position[1] - dy) / hy).ceil()) as i32;
            let jmax = ((ny as f64) - 1.0).min(((position[1] + dy) / hy).floor()) as i32;

            for j in jmin..=jmax {
                let dy2 = (position[1] - hy * j as f64).powi(2);
                let dz = if rtot2 > dx2 + dy2 {
                    (rtot2 - dx2 - dy2).sqrt() + 0.5 * hzed
                } else {
                    0.5 * hzed
                };
                let kmin = (0.0f64).max(((position[2] - dz) / hzed).ceil()) as i32;
                let kmax = ((nz as f64) - 1.0).min(((position[2] + dz) / hzed).floor()) as i32;

                for k in kmin..=kmax {
                    let dz2 = (k as f64 * hzed - position[2]).powi(2);
                    if dz2 + dy2 + dx2 <= rtot2 {
                        let gpos = [
                            i as f64 * hx + xmin,
                            j as f64 * hy + ymin,
                            k as f64 * hzed + zmin,
                        ];
                        let mut tgrad = [0.0f64; 3];
                        self.spline_select(srfm, &gpos, self.spline_win, irad, atom_id, &mut tgrad);

                        let ijk = ijk(i as usize, j as usize, k as usize, nx, ny);
                        let fmag = if self.pmgp.nonlin != 0 {
                            let mut fm = 0.0;
                            for m in 0..nion {
                                let ichop = &mut false;
                                fm += self.kappa[ijk] * self.pbe.ion_conc[m]
                                    * (vcap_exp(-self.pbe.ion_q[m] * self.u[ijk], ichop) - 1.0)
                                    / ionstr;
                            }
                            fm
                        } else {
                            self.u[ijk].powi(2) * self.kappa[ijk]
                        };
                        force[0] += zkappa2 * fmag * tgrad[0];
                        force[1] += zkappa2 * fmag * tgrad[1];
                        force[2] += zkappa2 * fmag * tgrad[2];
                    }
                }
            }
        }

        let scale = 0.5 * hx * hy * hzed * izmagic;
        force[0] *= scale;
        force[1] *= scale;
        force[2] *= scale;

        Ok(force)
    }

    /// Dielectric boundary force.
    /// Port of Vpmg_dbForce from vpmg.c line 5605
    pub fn db_force(&self, atom_id: usize, srfm: VsurfMeth) -> ApbsResult<[f64; 3]> {
        let mut db_force = [0.0f64; 3];

        if !self.filled {
            return Ok(db_force);
        }

        // Only spline-based surfaces supported
        match srfm {
            VsurfMeth::Spline | VsurfMeth::Spline3 | VsurfMeth::Spline4 => {}
            _ => return Ok(db_force),
        }

        let atom = self.pbe.alist.get_atom(atom_id);
        if atom.part_id == 0.0 {
            return Ok(db_force);
        }

        let apos = atom.position;
        let arad = atom.radius;
        let srad = self.pbe.solvent_radius;
        let epsp = self.pbe.solute_diel;
        let epsw = self.pbe.solvent_diel;
        // let kT = self.pbe.temperature * 1e-3 * apbs_generic::vunit::NA * apbs_generic::vunit::KB_ERG;
        let izmagic = 1.0 / self.pbe.zmagic;

        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let xmin = self.pmgp.xmin;
        let ymin = self.pmgp.ymin;
        let zmin = self.pmgp.zmin;
        let xmax = self.pmgp.xmax;
        let ymax = self.pmgp.ymax;
        let zmax = self.pmgp.zmax;

        if (epsp - epsw).abs() < VPMGSMALL {
            return Ok(db_force);
        }
        let deps = epsw - epsp;
        let depsi = 1.0 / deps;
        let rtot = arad + self.spline_win + srad;

        if apos[0] <= xmin + rtot || apos[0] >= xmax - rtot
            || apos[1] <= ymin + rtot || apos[1] >= ymax - rtot
            || apos[2] <= zmin + rtot || apos[2] >= zmax - rtot
        {
            return Ok(db_force);
        }

        let position = [apos[0] - xmin, apos[1] - ymin, apos[2] - zmin];
        // let rtot2 = rtot * rtot;

        let imin = ((position[0] - rtot) / hx).floor() as i32;
        if imin < 1 { return Ok(db_force); }
        let imax = ((position[0] + rtot) / hx).ceil() as i32;
        if imax > (nx as i32 - 2) { return Ok(db_force); }
        let jmin = ((position[1] - rtot) / hy).floor() as i32;
        if jmin < 1 { return Ok(db_force); }
        let jmax = ((position[1] + rtot) / hy).ceil() as i32;
        if jmax > (ny as i32 - 2) { return Ok(db_force); }
        let kmin = ((position[2] - rtot) / hzed).floor() as i32;
        if kmin < 1 { return Ok(db_force); }
        let kmax = ((position[2] + rtot) / hzed).ceil() as i32;
        if kmax > (nz as i32 - 2) { return Ok(db_force); }

        for i in imin..=imax {
            for j in jmin..=jmax {
                for k in kmin..=kmax {
                    let ijk_idx = ijk(i as usize, j as usize, k as usize, nx, ny);
                    // Hx at (i+0.5, j, k)
                    let mut gpos = [(i as f64 + 0.5) * hx + xmin, j as f64 * hy + ymin, k as f64 * hzed + zmin];
                    let hxijk = (self.epsx[ijk_idx] - epsp) * depsi;
                    let mut dhxijk = [0.0f64; 3];
                    self.spline_select(srfm, &gpos, self.spline_win, 0.0, atom_id, &mut dhxijk);
                    for l in 0..3 { dhxijk[l] *= hxijk; }

                    // Hy at (i, j+0.5, k)
                    gpos = [i as f64 * hx + xmin, (j as f64 + 0.5) * hy + ymin, k as f64 * hzed + zmin];
                    let hyijk = (self.epsy[ijk_idx] - epsp) * depsi;
                    let mut dhyijk = [0.0f64; 3];
                    self.spline_select(srfm, &gpos, self.spline_win, 0.0, atom_id, &mut dhyijk);
                    for l in 0..3 { dhyijk[l] *= hyijk; }

                    // Hz at (i, j, k+0.5)
                    gpos = [i as f64 * hx + xmin, j as f64 * hy + ymin, (k as f64 + 0.5) * hzed + zmin];
                    let hzijk = (self.epsz[ijk_idx] - epsp) * depsi;
                    let mut dhzijk = [0.0f64; 3];
                    self.spline_select(srfm, &gpos, self.spline_win, 0.0, atom_id, &mut dhzijk);
                    for l in 0..3 { dhzijk[l] *= hzijk; }

                    // Hx at (i-0.5, j, k)
                    gpos = [(i as f64 - 0.5) * hx + xmin, j as f64 * hy + ymin, k as f64 * hzed + zmin];
                    let hxim1jk = (self.epsx[ijk((i - 1) as usize, j as usize, k as usize, nx, ny)] - epsp) * depsi;
                    let mut dhxim1jk = [0.0f64; 3];
                    self.spline_select(srfm, &gpos, self.spline_win, 0.0, atom_id, &mut dhxim1jk);
                    for l in 0..3 { dhxim1jk[l] *= hxim1jk; }

                    // Hy at (i, j-0.5, k)
                    gpos = [i as f64 * hx + xmin, (j as f64 - 0.5) * hy + ymin, k as f64 * hzed + zmin];
                    let hyijm1k = (self.epsy[ijk(i as usize, (j - 1) as usize, k as usize, nx, ny)] - epsp) * depsi;
                    let mut dhyijm1k = [0.0f64; 3];
                    self.spline_select(srfm, &gpos, self.spline_win, 0.0, atom_id, &mut dhyijm1k);
                    for l in 0..3 { dhyijm1k[l] *= hyijm1k; }

                    // Hz at (i, j, k-0.5)
                    gpos = [i as f64 * hx + xmin, j as f64 * hy + ymin, (k as f64 - 0.5) * hzed + zmin];
                    let hzijkm1 = (self.epsz[ijk(i as usize, j as usize, (k - 1) as usize, nx, ny)] - epsp) * depsi;
                    let mut dhzijkm1 = [0.0f64; 3];
                    self.spline_select(srfm, &gpos, self.spline_win, 0.0, atom_id, &mut dhzijkm1);
                    for l in 0..3 { dhzijkm1[l] *= hzijkm1; }

                    // Calculate dielectric boundary force
                    let db_fmag = self.u[ijk_idx];
                    let hx2 = hx * hx;
                    let hy2 = hy * hy;
                    let hz2 = hzed * hzed;

                    let mut tgrad = [0.0f64; 3];
                    for l in 0..3 {
                        tgrad[l] = (dhxijk[l] * (self.u[ijk((i + 1) as usize, j as usize, k as usize, nx, ny)] - self.u[ijk_idx])
                            + dhxim1jk[l] * (self.u[ijk((i - 1) as usize, j as usize, k as usize, nx, ny)] - self.u[ijk_idx])) / hx2
                            + (dhyijk[l] * (self.u[ijk(i as usize, (j + 1) as usize, k as usize, nx, ny)] - self.u[ijk_idx])
                                + dhyijm1k[l] * (self.u[ijk(i as usize, (j - 1) as usize, k as usize, nx, ny)] - self.u[ijk_idx])) / hy2
                            + (dhzijk[l] * (self.u[ijk(i as usize, j as usize, (k + 1) as usize, nx, ny)] - self.u[ijk_idx])
                                + dhzijkm1[l] * (self.u[ijk(i as usize, j as usize, (k - 1) as usize, nx, ny)] - self.u[ijk_idx])) / hz2;
                    }
                    for l in 0..3 {
                        db_force[l] += db_fmag * tgrad[l];
                    }
                }
            }
        }

        let scale = -hx * hy * hzed * deps * 0.5 * izmagic;
        db_force[0] *= scale;
        db_force[1] *= scale;
        db_force[2] *= scale;

        Ok(db_force)
    }

    /// Select spline accessibility gradient function based on surface method.
    /// Port of Vpmg_splineSelect from vpmg.c line 1893
    fn spline_select(&self, srfm: VsurfMeth, gpos: &[f64; 3], win: f64, infrad: f64, atom_id: usize, force: &mut [f64; 3]) {
        let acc = self.pbe.acc.lock().unwrap();
        let atom = self.pbe.alist.get_atom(atom_id);
        match srfm {
            VsurfMeth::Spline | VsurfMeth::Spline3 | VsurfMeth::Spline4 => {
                acc.spline_acc_grad_atom_norm(*gpos, win, infrad, atom, force);
            }
            _ => {}
        }
    }

    /// SMPBE mobile ion energy.
    /// Port of Vpmg_qmEnergySMPBE from vpmg.c line 1490
    pub fn qm_energy_smpbe(&self, _ext_flag: i32) -> f64 {
        let nx = self.pmgp.nx as usize;
        let ny = self.pmgp.ny as usize;
        let nz = self.pmgp.nz as usize;
        let hx = self.pmgp.hx;
        let hy = self.pmgp.hy;
        let hzed = self.pmgp.hzed;
        let zkappa2 = self.pbe.zkappa2;
        let ionstr = self.pbe.bulk_ionic_strength;

        if zkappa2 < VSMALL {
            return 0.0;
        }
        let _zks2 = 0.5 * zkappa2 / ionstr;

        if !self.filled {
            return 0.0;
        }

        let nf = nx * ny * nz;
        let z1 = self.pbe.ion_q[0];
        let z2 = self.pbe.ion_q[1];
        let z3 = self.pbe.ion_q[2];
        let cb1 = self.pbe.ion_conc[0];
        let cb2 = self.pbe.ion_conc[1];
        let cb3 = self.pbe.ion_conc[2];
        let a = self.pbe.smvolume;
        let k = self.pbe.smsize;

        let frac_occ_a = NA * cb1 * a.powi(3);
        let frac_occ_b = NA * cb2 * a.powi(3);
        let frac_occ_c = NA * cb3 * a.powi(3);
        let phi = (frac_occ_a / k) + frac_occ_b + frac_occ_c;

        let mut energy = 0.0;

        if self.pmgp.nonlin != 0 {
            for i in 0..nf {
                if self.pvec[i] * self.kappa[i] > VSMALL {
                    let mut ichop1 = false;
                    let mut ichop2 = false;
                    let mut ichop3 = false;
                    let a1 = vcap_exp(-z1 * self.u[i], &mut ichop1);
                    let a2 = vcap_exp(-z2 * self.u[i], &mut ichop2);
                    let a3 = vcap_exp(-z3 * self.u[i], &mut ichop3);

                    let gpark = 1.0 - phi + (frac_occ_a / k) * a1;
                    let denom = gpark.powf(k)
                        + (1.0 - frac_occ_b - frac_occ_c).powf(k - 1.0)
                            * (frac_occ_b * a2 + frac_occ_c * a3);

                    let c1 = if cb1 > VSMALL && denom > VSMALL {
                        let v = NA * cb1 * gpark.powf(k - 1.0) * a1 / denom;
                        if v != v { 0.0 } else { v }
                    } else {
                        0.0
                    };
                    let c2 = if cb2 > VSMALL && denom > VSMALL {
                        let v = NA * cb2 * (1.0 - frac_occ_b - frac_occ_c).powf(k - 1.0) * a2 / denom;
                        if v != v { 0.0 } else { v }
                    } else {
                        0.0
                    };
                    let c3 = if cb3 > VSMALL && denom > VSMALL {
                        let v = NA * cb3 * (1.0 - frac_occ_b - frac_occ_c).powf(k - 1.0) * a3 / denom;
                        if v != v { 0.0 } else { v }
                    } else {
                        0.0
                    };

                    let arg1 = (1.0 - c1 * a.powi(3) / k - c2 * a.powi(3) - c3 * a.powi(3))
                        / (1.0 - phi);
                    let arg2 = (1.0 - c2 * a.powi(3) - c3 * a.powi(3))
                        / (1.0 - phi + frac_occ_a / k);
                    let curr_energy = if arg1 > 0.0 && arg2 > 0.0 {
                        k * arg1.ln() - (k - 1.0) * arg2.ln()
                    } else {
                        0.0
                    };

                    energy += self.pvec[i] * self.kappa[i] * curr_energy;
                }
            }
            energy = -energy / a.powi(3);
        }

        energy * hx * hy * hzed
    }

    /// Print column comparison for matrix output.
    /// Port of Vpmg_printColComp from vpmg.c line 87 (stub - requires bcolcomp/pcolcomp)
    pub fn print_col_comp(&self, _path: &str, _title: &str, _mxtype: &str, _flag: i32) {
        // This function requires the bcolcomp/pcolcomp Fortran interface
        // which is not ported. Stub for now.
        eprintln!("Vpmg_printColComp: not yet implemented (requires Fortran bcolcomp/pcolcomp)");
    }

    /// Fill permanent multipole coefficients.
    /// Port of fillcoPermanentMultipole from vpmg.c line 6835 (stub - requires Vpmg_spline4Acc)
    pub fn fillco_permanent_multipole(&mut self) {
        // This requires Vpmg_spline4Acc which is a complex spline evaluation function
        // used for multipole force fields. Not critical for standard APBS usage.
        eprintln!("fillcoPermanentMultipole: not yet implemented (requires multipole force field support)");
    }

    /// Fill induced multipole coefficients.
    /// Port of fillcoPermanentInduced from vpmg.c line 7373 (stub - requires Vpmg_spline4Acc)
    pub fn fillco_permanent_induced(&mut self) {
        // This requires Vpmg_spline4Acc which is a complex spline evaluation function
        // used for induced dipole force fields. Not critical for standard APBS usage.
        eprintln!("fillcoPermanentInduced: not yet implemented (requires induced dipole force field support)");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vpmg_creation() {
        use apbs_generic::mgparm::{MGparm, MGparmCalcType};
        use std::sync::Arc;

        let mut mgparm = MGparm::new(MGparmCalcType::Manual);
        mgparm.dime = [33, 33, 33];
        mgparm.glen = [10.0, 10.0, 10.0];
        mgparm.grid = [10.0 / 32.0; 3];
        mgparm.center = [0.0; 3];

        let pmgp = Vpmgp::new(&mgparm);

        let alist = Arc::new(apbs_generic::valist::Valist::new());
        let pbe = Arc::new(Vpbe::new(
            alist,
            0, &[], &[], &[],
            298.15, 2.0, 78.36, 1.4, 0, 10.0,
            0.0, 0.0, 78.36, 0.0,
        ).unwrap());

        let pbeparm = PBEparm::new();
        let pmg = Vpmg::new(pmgp, pbe, 0, None, &pbeparm, PBEparmCalcEnergy::Total);

        assert_eq!(pmg.pmgp.nx, 33);
        assert_eq!(pmg.u.len(), pmg.pmgp.narr as usize);
    }

    #[test]
    fn test_multipolebc() {
        let tensor = multipolebc(10.0, 0.0, 2.0, 78.36, 1.4);
        // monopole should be positive (1/eps_w * 1/r)
        assert!(tensor[0] > 0.0);
        // dipole should be negative (-1/r^3 * correction)
        assert!(tensor[1] < 0.0);
        // quadrupole should be positive (3/r^5 * correction)
        assert!(tensor[2] > 0.0);
    }

    #[test]
    fn test_bspline2() {
        // bspline2(1.5) = 0.5*1.5*0.5 + 0.5*1.5*0.5 = 0.75
        let v = bspline2(1.5);
        assert!((v - 0.75).abs() < 1e-10);
        // bspline2(0) = 0
        assert!(bspline2(0.0).abs() < 1e-10);
        // bspline2(3) = 0
        assert!(bspline2(3.0).abs() < 1e-10);
        // bspline2(1.0) = 0.5*1.0*1.0 + 0.5*2.0*0.0 = 0.5
        assert!((bspline2(1.0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_solve_rejects_unsupported_mgdisc() {
        use apbs_generic::mgparm::{MGparm, MGparmCalcType};
        use std::sync::Arc;

        let mut mgparm = MGparm::new(MGparmCalcType::Manual);
        mgparm.dime = [33, 33, 33];
        mgparm.glen = [10.0, 10.0, 10.0];
        mgparm.grid = [10.0 / 32.0; 3];
        mgparm.center = [0.0; 3];

        let mut pmgp = Vpmgp::new(&mgparm);
        pmgp.mgdisc = 1;

        let alist = Arc::new(apbs_generic::valist::Valist::new());
        let pbe = Arc::new(Vpbe::new(
            alist,
            0, &[], &[], &[],
            298.15, 2.0, 78.36, 1.4, 0, 10.0,
            0.0, 0.0, 78.36, 0.0,
        ).unwrap());

        let pbeparm = PBEparm::new();
        let mut pmg = Vpmg::new(pmgp, pbe, 0, None, &pbeparm, PBEparmCalcEnergy::Total);
        pmg.filled = true;

        let err = pmg.solve().expect_err("mgdisc != 0 should be rejected");
        let msg = format!("{}", err);
        assert!(msg.contains("mgdisc=1"));
        assert!(msg.contains("not implemented"));
    }
}
