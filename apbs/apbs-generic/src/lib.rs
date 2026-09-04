// APBS Generic - Core data structures and utilities
// Port of src/generic/

pub mod error;
pub mod vhal;
pub mod vunit;
pub mod vmatrix;
pub mod vstring;
pub mod vcap;
pub mod vatom;
pub mod valist;
pub mod vclist;
pub mod vparam;
pub mod pbeparm;
pub mod mgparm;
pub mod apolparm;
pub mod nosh;

// Stub modules for larger implementations
pub mod vpbe {
    use crate::error::ApbsResult;
    use crate::valist::Valist;
    use crate::vclist::Vclist;
    use crate::vacc::Vacc;
    
    

    /// Poisson-Boltzmann Equation master object
    pub struct Vpbe {
        pub alist: std::sync::Arc<Valist>,
        pub clist: std::sync::Arc<Vclist>,
        pub acc: std::sync::Mutex<Vacc>,
        pub temperature: f64,
        pub solute_diel: f64,
        pub solvent_diel: f64,
        pub solvent_radius: f64,
        pub bulk_ionic_strength: f64,
        pub max_ion_radius: f64,
        pub num_ion: usize,
        pub ion_conc: [f64; crate::vhal::MAXION],
        pub ion_radii: [f64; crate::vhal::MAXION],
        pub ion_q: [f64; crate::vhal::MAXION],
        pub xkappa: f64,
        pub deblen: f64,
        pub zkappa2: f64,
        pub zmagic: f64,
        pub solute_center: [f64; 3],
        pub solute_radius: f64,
        pub solute_xlen: f64,
        pub solute_ylen: f64,
        pub solute_zlen: f64,
        pub solute_charge: f64,
        pub smvolume: f64,
        pub smsize: f64,
        pub ipkey: i32,
        pub z_mem: f64,
        pub l_mem: f64,
        pub membrane_diel: f64,
        pub memv: f64,
    }

    impl Vpbe {
        pub fn new(
            alist: std::sync::Arc<Valist>,
            ion_num: usize,
            ion_conc: &[f64],
            ion_radii: &[f64],
            ion_q: &[f64],
            temperature: f64,
            solute_diel: f64,
            solvent_diel: f64,
            solvent_radius: f64,
            _focus_flag: i32,
            sdens: f64,
            z_mem: f64,
            l_mem: f64,
            membrane_diel: f64,
            memv: f64,
        ) -> ApbsResult<Self> {
            // Compute solute geometry
            let mut center = [0.0f64; 3];
            let mut mincrd = [f64::MAX; 3];
            let mut maxcrd = [f64::MIN; 3];
            let mut max_rad = 0.0f64;
            let mut charge = 0.0f64;
            let num_atoms = alist.number_atoms();

            for atom in alist.atoms() {
                for d in 0..3 {
                    if atom.position[d] < mincrd[d] { mincrd[d] = atom.position[d]; }
                    if atom.position[d] > maxcrd[d] { maxcrd[d] = atom.position[d]; }
                }
                if atom.radius > max_rad { max_rad = atom.radius; }
                charge += atom.charge;
            }

            for d in 0..3 {
                if num_atoms == 0 {
                    mincrd[d] = -1.0;
                    maxcrd[d] = 1.0;
                }
                center[d] = (mincrd[d] + maxcrd[d]) / 2.0;
            }

            let solute_xlen = maxcrd[0] - mincrd[0];
            let solute_ylen = maxcrd[1] - mincrd[1];
            let solute_zlen = maxcrd[2] - mincrd[2];
            let solute_radius = (solute_xlen * solute_xlen + solute_ylen * solute_ylen + solute_zlen * solute_zlen).sqrt() / 2.0;

            // Compute bulk ionic strength
            let mut bulk_ionic_strength = 0.0;
            let mut max_ion_radius = 0.0;
            for i in 0..ion_num {
                bulk_ionic_strength += 0.5 * ion_conc[i] * ion_q[i] * ion_q[i];
                if ion_radii[i] > max_ion_radius {
                    max_ion_radius = ion_radii[i];
                }
            }

            // Compute Debye-Huckel parameters
            // kappa^2 = (8*pi*N_A*e_c^2*I_s)/(1000*eps_w*k_B*T)
            use crate::vunit::{NA_ESU, EC_ESU, KB_ERG};
            let eps_w = solvent_diel;
            // Match legacy APBS unit conversion: extra 1e-16 factor in kappa^2.
            let kappa2 = (8.0 * std::f64::consts::PI * NA_ESU * EC_ESU * EC_ESU * bulk_ionic_strength * 1.0e-16)
                / (1000.0 * eps_w * KB_ERG * temperature);
            let xkappa = kappa2.sqrt();
            let deblen = if xkappa > 0.0 { 1.0 / xkappa } else { f64::INFINITY };
            let zkappa2 = eps_w * kappa2;
            let zmagic = 4.0 * std::f64::consts::PI * EC_ESU * EC_ESU / (KB_ERG * temperature) * 1.0e8;

            // Create cell list
            // C code: radius = max(maxIonRadius, solventRadius) + MAX_SPLINE_WINDOW
            let max_spline_window = 0.5;
            let clist_max_radius = max_ion_radius.max(solvent_radius) + max_spline_window;
            let solute_dim = solute_xlen.max(solute_ylen).max(solute_zlen);
            let hash_dim = ((solute_dim / 0.5) as usize).max(3).min(crate::vhal::MAX_HASH_DIM);
            let npts = [hash_dim, hash_dim, hash_dim];
            let clist = std::sync::Arc::new(Vclist::new_auto(&alist, clist_max_radius, npts)?);

            // Create accessibility object
            let acc = std::sync::Mutex::new(Vacc::new(&alist, &clist, sdens));

            let mut ion_conc_arr = [0.0; crate::vhal::MAXION];
            let mut ion_radii_arr = [0.0; crate::vhal::MAXION];
            let mut ion_q_arr = [0.0; crate::vhal::MAXION];
            for i in 0..ion_num.min(crate::vhal::MAXION) {
                ion_conc_arr[i] = ion_conc[i];
                ion_radii_arr[i] = ion_radii[i];
                ion_q_arr[i] = ion_q[i];
            }

            Ok(Self {
                alist,
                clist,
                acc,
                temperature,
                solute_diel,
                solvent_diel,
                solvent_radius,
                bulk_ionic_strength,
                max_ion_radius,
                num_ion: ion_num,
                ion_conc: ion_conc_arr,
                ion_radii: ion_radii_arr,
                ion_q: ion_q_arr,
                xkappa,
                deblen,
                zkappa2,
                zmagic,
                solute_center: center,
                solute_radius,
                solute_xlen,
                solute_ylen,
                solute_zlen,
                solute_charge: charge,
                smvolume: 0.0,
                smsize: 0.0,
                ipkey: 0,
                z_mem,
                l_mem,
                membrane_diel,
                memv,
            })
        }

        pub fn get_coulomb_energy1(&self) -> f64 {
            let eps = self.solvent_diel;
            let kt = crate::vunit::KB_ERG * self.temperature;
            let coef = crate::vunit::EC_ESU * crate::vunit::EC_ESU / (4.0 * std::f64::consts::PI * crate::vunit::EPS0 * eps * 1.0e-10) / kt;

            let natoms = self.alist.number_atoms();
            let mut energy = 0.0;
            for i in 0..natoms {
                let ai = self.alist.get_atom(i);
                for j in (i + 1)..natoms {
                    let aj = self.alist.get_atom(j);
                    let dx = ai.position[0] - aj.position[0];
                    let dy = ai.position[1] - aj.position[1];
                    let dz = ai.position[2] - aj.position[2];
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    if r > 0.0 {
                        energy += ai.charge * aj.charge / r;
                    }
                }
            }
            coef * energy
        }

        /// Compute SASA (solvent-accessible surface area)
        pub fn sasa(&self, radius: f64) -> f64 {
            self.acc.lock().unwrap().sasa(radius)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::Vacc;
        use crate::valist::Valist;
        use crate::vatom::Vatom;
        use crate::vclist::Vclist;
        use std::sync::Arc;

        fn single_atom_vacc(radius: f64) -> Vacc {
            let mut atom = Vatom::new();
            atom.set_atom_id(0);
            atom.set_position([0.0, 0.0, 0.0]);
            atom.set_radius(radius);

            let mut alist = Valist::new();
            alist.atoms.push(atom);
            alist.get_statistics();

            let alist = Arc::new(alist);
            let clist = Arc::new(Vclist::new_auto(&alist, 2.0, [10, 10, 10]).unwrap());
            Vacc::new(&alist, &clist, 10.0)
        }

        #[test]
        fn ivdw_accessibility_matches_c_semantics() {
            let vacc = single_atom_vacc(1.5);
            assert_eq!(vacc.ivdw_acc([5.0, 0.0, 0.0], 1.4), 1.0);
            assert_eq!(vacc.ivdw_acc([0.1, 0.0, 0.0], 1.4), 0.0);
        }

        #[test]
        fn isolated_atom_surface_keeps_reference_points() {
            let mut vacc = single_atom_vacc(1.5);
            let ref_npts = vacc.ref_sphere.npts;
            let surf = vacc.atom_sas_points(0, 1.4).unwrap();
            assert!(surf.npts > 0);
            assert_eq!(surf.npts, ref_npts);
        }
    }
}

pub mod vgreen {
    /// Green's function computation (stub)
    pub struct Vgreen;

    impl Vgreen {
        pub fn new() -> Self {
            Self
        }
    }
}

pub mod geoflowparm {
    use crate::error::ApbsResult;

    /// Geometric flow parameters
    #[derive(Debug, Clone, Default)]
    pub struct GEOFLOWparm {
        pub parsed: bool,
    }

    impl GEOFLOWparm {
        pub fn new() -> Self { Self::default() }
        pub fn check(&self) -> ApbsResult<()> { Ok(()) }
    }
}

pub mod pbamparm {
    use crate::error::ApbsResult;

    /// PBAM parameters
    #[derive(Debug, Clone, Default)]
    pub struct PBAMparm {
        pub parsed: bool,
    }

    impl PBAMparm {
        pub fn new() -> Self { Self::default() }
        pub fn check(&self) -> ApbsResult<()> { Ok(()) }
    }
}

pub mod pbsamparm {
    use crate::error::ApbsResult;

    /// PBSAM parameters
    #[derive(Debug, Clone, Default)]
    pub struct PBSAMparm {
        pub parsed: bool,
    }

    impl PBSAMparm {
        pub fn new() -> Self { Self::default() }
        pub fn check(&self) -> ApbsResult<()> { Ok(()) }
    }
}

pub mod vacc {
    use crate::valist::Valist;
    use crate::vclist::Vclist;
    use crate::vatom::Vatom;
    use crate::vhal::VSMALL;
    use crate::apolparm::APOLparm;
    use crate::error::ApbsError;
    use rayon::prelude::*;

    /// Per-atom surface point set
    #[derive(Debug, Clone)]
    pub struct VaccSurf {
        pub xpts: Vec<f64>,
        pub ypts: Vec<f64>,
        pub zpts: Vec<f64>,
        pub bpts: Vec<bool>,
        pub area: f64,
        pub npts: usize,
        pub probe_radius: f64,
    }

    impl VaccSurf {
        pub fn new(npts: usize, probe_radius: f64) -> Self {
            Self {
                xpts: vec![0.0; npts],
                ypts: vec![0.0; npts],
                zpts: vec![0.0; npts],
                bpts: vec![false; npts],
                area: 0.0,
                npts,
                probe_radius,
            }
        }

        /// Generate a reference sphere with approximately npts uniformly distributed points
        pub fn ref_sphere(npts: usize) -> Self {
            // Port of VaccSurf_refSphere from vacc.c line 924
            let mut surf = Self::new(npts, 1.0);
            let frac = npts as f64 / 4.0;
            let ntheta = (std::f64::consts::PI * frac).sqrt().round() as i32;
            let dtheta = std::f64::consts::PI / ntheta as f64;
            let nphimax = 2 * ntheta;
            let mut count = 0;

            for itheta in 0..ntheta {
                let theta = dtheta * itheta as f64;
                let sintheta = theta.sin();
                let nphi = (sintheta * nphimax as f64).round() as i32;
                if nphi > 0 {
                    let dphi = 2.0 * std::f64::consts::PI / nphi as f64;
                    for iphi in 0..nphi {
                        if count >= npts as i32 {
                            break;
                        }
                        let phi = dphi * iphi as f64;
                        let sinphi = phi.sin();
                        let cosphi = phi.cos();
                        surf.xpts[count as usize] = cosphi * sintheta;
                        surf.ypts[count as usize] = sinphi * sintheta;
                        surf.zpts[count as usize] = theta.cos();
                        surf.bpts[count as usize] = true;
                        count += 1;
                    }
                }
            }

            surf.npts = count as usize;
            surf.area = 4.0 * std::f64::consts::PI;
            surf
        }
    }

    /// Solvent/ion accessibility oracle
    pub struct Vacc {
        pub alist: std::sync::Arc<Valist>,
        pub clist: std::sync::Arc<Vclist>,
        pub atom_flags: Vec<bool>,
        pub ref_sphere: VaccSurf,
        pub surf: Vec<Option<VaccSurf>>,
        pub surf_density: f64,
    }

    impl Vacc {
        pub fn new(alist: &std::sync::Arc<Valist>, clist: &std::sync::Arc<Vclist>, surf_density: f64) -> Self {
            let num_atoms = alist.number_atoms();

            // Compute reference sphere size: nsphere = ceil(4*PI*maxrad^2 * surf_density)
            // where maxrad = max atom radius + Vclist max_radius
            let mut maxrad = 0.0f64;
            for i in 0..num_atoms {
                let rad = alist.get_atom(i).radius;
                if rad > maxrad { maxrad = rad; }
            }
            maxrad += clist.max_radius();
            let maxarea = 4.0 * std::f64::consts::PI * maxrad * maxrad;
            let nsphere = (maxarea * surf_density).ceil() as usize;
            let nsphere = nsphere.max(100); // minimum 100 points

            let ref_sphere = VaccSurf::ref_sphere(nsphere);

            Self {
                alist: std::sync::Arc::clone(alist),
                clist: std::sync::Arc::clone(clist),
                atom_flags: vec![false; num_atoms],
                ref_sphere,
                surf: vec![None; num_atoms],
                surf_density,
            }
        }

        /// Van der Waals accessibility test
        pub fn vdw_acc(&self, center: [f64; 3]) -> f64 {
            if let Some(cell) = self.clist.get_cell(center) {
                for &atom_idx in &cell.atoms {
                    let atom = self.alist.get_atom(atom_idx);
                    let dx = center[0] - atom.position[0];
                    let dy = center[1] - atom.position[1];
                    let dz = center[2] - atom.position[2];
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    let rad = atom.radius;
                    if dist2 < rad * rad {
                        return 0.0;
                    }
                }
            }
            1.0
        }

        /// Inflated van der Waals accessibility test
        pub fn ivdw_acc(&self, center: [f64; 3], radius: f64) -> f64 {
            if self.ivdw_acc_exclus(center, radius, -1) { 1.0 } else { 0.0 }
        }

        fn ivdw_acc_exclus(&self, center: [f64; 3], radius: f64, exclude_atom: i32) -> bool {
            Self::ivdw_acc_exclus_from(&self.clist, &self.alist, center, radius, exclude_atom)
        }

        fn ivdw_acc_exclus_from(
            clist: &std::sync::Arc<Vclist>,
            alist: &std::sync::Arc<Valist>,
            center: [f64; 3],
            radius: f64,
            exclude_atom: i32,
        ) -> bool {
            if let Some(cell) = clist.get_cell(center) {
                for &atom_idx in &cell.atoms {
                    if atom_idx as i32 == exclude_atom {
                        continue;
                    }
                    let atom = alist.get_atom(atom_idx);
                    let dx = center[0] - atom.position[0];
                    let dy = center[1] - atom.position[1];
                    let dz = center[2] - atom.position[2];
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    let rad = radius + atom.radius;
                    if dist2 < rad * rad {
                        return false;
                    }
                }
                return true;
            }
            true
        }

        /// Compute total SASA
        pub fn sasa(&mut self, radius: f64) -> f64 {
            // Lazy initialization of per-atom surfaces
            if self.surf.iter().any(|s| s.is_none()) {
                self.build_surfaces(radius);
            }

            let mut total = 0.0;
            for s in &self.surf {
                if let Some(surf) = s {
                    total += surf.area;
                }
            }
            total
        }

        fn build_surfaces(&mut self, prad: f64) {
            let num_atoms = self.alist.number_atoms();
            if num_atoms == 0 {
                return;
            }

            let alist = std::sync::Arc::clone(&self.alist);
            let clist = std::sync::Arc::clone(&self.clist);
            let ref_sphere = self.ref_sphere.clone();
            // Each per-atom surface is computed from the atom list and the
            // neighbor cell list only, so atoms can be processed in parallel.
            let surfaces: Vec<VaccSurf> = (0..num_atoms)
                .into_par_iter()
                .map(|i| {
                    let atom = alist.get_atom(i);
                    Self::atom_surf_from(&alist, &clist, atom, &ref_sphere, prad)
                })
                .collect();
            for (i, surface) in surfaces.into_iter().enumerate() {
                self.surf[i] = Some(surface);
            }
        }

        fn atom_surf(&self, atom: &Vatom, ref_sphere: &VaccSurf, prad: f64) -> VaccSurf {
            Self::atom_surf_from(&self.alist, &self.clist, atom, ref_sphere, prad)
        }

        fn atom_surf_from(
            alist: &std::sync::Arc<Valist>,
            clist: &std::sync::Arc<Vclist>,
            atom: &Vatom,
            ref_sphere: &VaccSurf,
            prad: f64,
        ) -> VaccSurf {
            let arad = atom.radius;
            let rad = arad + prad;
            let mut count = 0;

            let mut result = VaccSurf::new(ref_sphere.npts, prad);

            for ipoint in 0..ref_sphere.npts {
                let x = atom.position[0] + ref_sphere.xpts[ipoint] * rad;
                let y = atom.position[1] + ref_sphere.ypts[ipoint] * rad;
                let z = atom.position[2] + ref_sphere.zpts[ipoint] * rad;

                if Self::ivdw_acc_exclus_from(clist, alist, [x, y, z], prad, atom.atom_id()) {
                    result.xpts[count] = x;
                    result.ypts[count] = y;
                    result.zpts[count] = z;
                    result.bpts[count] = true;
                    count += 1;
                }
            }

            result.npts = count;
            result.area = 4.0 * std::f64::consts::PI * rad * rad * (count as f64 / ref_sphere.npts.max(1) as f64);
            result
        }

        /// Get SAS surface points for an atom (builds surfaces if needed)
        pub fn atom_sas_points(&mut self, atom_idx: usize, prad: f64) -> Option<&VaccSurf> {
            if self.surf[atom_idx].is_none() {
                let atom = self.alist.get_atom(atom_idx);
                let ref_clone = self.ref_sphere.clone();
                self.surf[atom_idx] = Some(self.atom_surf(atom, &ref_clone, prad));
            }
            self.surf[atom_idx].as_ref()
        }

        /// Gradient of spline accessibility for a single atom (normalized by chi).
        /// Port of Vacc_splineAccGradAtomNorm from vacc.c line 314
        pub fn spline_acc_grad_atom_norm(&self, center: [f64; 3], win: f64, infrad: f64, atom: &Vatom, grad: &mut [f64; 3]) {
            grad[0] = 0.0;
            grad[1] = 0.0;
            grad[2] = 0.0;

            if atom.radius <= 0.0 {
                return;
            }

            let apos = atom.position;
            let arad = atom.radius + infrad;
            let dx = center[0] - apos[0];
            let dy = center[1] - apos[1];
            let dz = center[2] - apos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            if dist < arad - win {
                return;
            }
            if dist > arad + win {
                return;
            }
            if (dist - (arad - win)).abs() < VSMALL || (dist - (arad + win)).abs() < VSMALL {
                return;
            }

            let w2i = 1.0 / (win * win);
            let w3i = 1.0 / (win * win * win);
            let sm = dist - arad + win;
            let sm2 = sm * sm;
            let mychi = 0.75 * sm2 * w2i - 0.25 * sm * sm2 * w3i;
            let mygrad = 1.5 * sm * w2i - 0.75 * sm2 * w3i;

            if mychi > 0.0 {
                let scale = -(mygrad / mychi) / dist;
                grad[0] = scale * dx;
                grad[1] = scale * dy;
                grad[2] = scale * dz;
            }
        }

        /// Spline-based accessibility for a single atom
        pub fn spline_acc_atom(&self, center: [f64; 3], win: f64, infrad: f64, atom: &Vatom) -> f64 {
            let dx = center[0] - atom.position[0];
            let dy = center[1] - atom.position[1];
            let dz = center[2] - atom.position[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let arad = atom.radius + infrad;

            if dist < arad - win {
                0.0
            } else if dist > arad + win {
                1.0
            } else {
                let sm = dist - arad + win;
                0.75 * sm * sm / (win * win) - 0.25 * sm * sm * sm / (win * win * win)
            }
        }

        /// Private helper: multi-atom spline accessibility subroutine.
        /// Port of splineAcc from vacc.c line 491.
        /// Loops through atoms in the cell, using atom_flags to avoid double-counting,
        /// and multiplies individual spline_acc_atom values.
        fn spline_acc(&mut self, center: [f64; 3], win: f64, infrad: f64, cell: &crate::vclist::VclistCell) -> f64 {
            let mut value = 1.0f64;

            for &atom_idx in &cell.atoms {
                let atom = self.alist.get_atom(atom_idx);
                let atom_id = atom.atom_id() as usize;

                // Check to see if we've counted this atom already
                if !self.atom_flags[atom_id] {
                    self.atom_flags[atom_id] = true;
                    value *= self.spline_acc_atom(center, win, infrad, atom);

                    if value < VSMALL {
                        return value;
                    }
                }
            }

            value
        }

        /// Multi-atom spline accessibility.
        /// Port of Vacc_splineAcc from vacc.c line 526.
        /// Uses the cell list to find nearby atoms, multiplies individual
        /// spline_acc_atom values together.
        pub fn spline_acc_multi(&mut self, center: [f64; 3], win: f64, infrad: f64) -> f64 {
            // Check that the cell list has sufficient radius
            assert!(
                self.clist.max_radius() >= (win + infrad),
                "Vacc_splineAcc: Vclist has max_radius={}; insufficient for win={}, infrad={}",
                self.clist.max_radius(), win, infrad
            );

            // Get a cell; if no cell, return 1.0 (no nearby atoms)
            let cell = match self.clist.get_cell(center) {
                Some(c) => c.clone(),
                None => return 1.0,
            };

            // Reset the list of atom flags for atoms in this cell
            for &atom_idx in &cell.atoms {
                let atom = self.alist.get_atom(atom_idx);
                let atom_id = atom.atom_id() as usize;
                self.atom_flags[atom_id] = false;
            }

            self.spline_acc(center, win, infrad, &cell)
        }

        /// Gradient of multi-atom spline accessibility.
        /// Port of Vacc_splineAccGrad from vacc.c line 559.
        pub fn spline_acc_grad(&mut self, center: [f64; 3], win: f64, infrad: f64, grad: &mut [f64; 3]) {
            // Check that the cell list has sufficient radius
            assert!(
                self.clist.max_radius() >= (win + infrad),
                "Vacc_splineAccGrad: Vclist max_radius={}; insufficient for win={}, infrad={}",
                self.clist.max_radius(), win, infrad
            );

            // Reset the gradient
            for g in grad.iter_mut() {
                *g = 0.0;
            }

            // Get the cell; check for nullity
            let cell = match self.clist.get_cell(center) {
                Some(c) => c.clone(),
                None => return,
            };

            // Reset the list of atom flags
            for &atom_idx in &cell.atoms {
                let atom = self.alist.get_atom(atom_idx);
                let atom_id = atom.atom_id() as usize;
                self.atom_flags[atom_id] = false;
            }

            // Get the local accessibility
            let acc = self.spline_acc(center, win, infrad, &cell);

            // Accumulate the gradient of all local atoms
            if acc > VSMALL {
                let mut tgrad = [0.0f64; 3];
                for &atom_idx in &cell.atoms {
                    let atom = self.alist.get_atom(atom_idx);
                    self.spline_acc_grad_atom_norm(center, win, infrad, atom, &mut tgrad);
                }
                for i in 0..3 {
                    grad[i] += tgrad[i];
                }
            }
            for g in grad.iter_mut() {
                *g *= -acc;
            }
        }

        /// Get per-atom SASA (solvent-accessible surface area).
        /// Port of Vacc_atomSASA from vacc.c line 778.
        /// Returns the surface area for the specified atom, building
        /// surfaces if needed.
        pub fn atom_sasa(&mut self, atom_idx: usize, radius: f64) -> f64 {
            // Build all surfaces if none have been built yet
            if self.surf.iter().all(|s| s.is_none()) {
                self.build_surfaces(radius);
            }

            // See if this surface needs to be built or rebuilt
            let needs_build = match &self.surf[atom_idx] {
                None => true,
                Some(s) => (s.probe_radius - radius).abs() > VSMALL,
            };

            if needs_build {
                let atom = self.alist.get_atom(atom_idx);
                let ref_clone = self.ref_sphere.clone();
                self.surf[atom_idx] = Some(self.atom_surf(atom, &ref_clone, radius));
            }

            self.surf[atom_idx].as_ref().unwrap().area
        }

        /// Total solvent accessible volume.
        /// Port of Vacc_totalSAV from vacc.c line 1501.
        /// Numerical integration using trapezoidal rule over the accessibility
        /// function across the domain.
        pub fn total_sav(&self, radius: f64, apolparm: &APOLparm) -> f64 {
            let mut sav = 0.0f64;
            let vol_density = 2.0f64;

            let lower_corner = self.clist.lower_corner;
            let upper_corner = self.clist.upper_corner;

            // Compute grid spacings
            let mut spacs = [0.0f64; 3];
            let mut npts = [0usize; 3];
            for i in 0..3 {
                let len = upper_corner[i] - lower_corner[i];
                let fn_val = len * vol_density + 1.0;
                npts[i] = fn_val.ceil() as usize;
                spacs[i] = if npts[i] > 1 {
                    len / (npts[i] as f64 - 1.0)
                } else {
                    len
                };
                if apolparm.setgrid {
                    if apolparm.grid[i] > spacs[i] {
                        // Warning: grid value larger than recommended
                    }
                    spacs[i] = apolparm.grid[i];
                }
            }

            // Trapezoidal integration over the 3D grid
            let mut x = lower_corner[0];
            while x <= upper_corner[0] + VSMALL {
                let wx = if (x - lower_corner[0]).abs() < VSMALL
                    || (x - upper_corner[0]).abs() < VSMALL
                {
                    0.5
                } else {
                    1.0
                };

                let mut y = lower_corner[1];
                while y <= upper_corner[1] + VSMALL {
                    let wy = if (y - lower_corner[1]).abs() < VSMALL
                        || (y - upper_corner[1]).abs() < VSMALL
                    {
                        0.5
                    } else {
                        1.0
                    };

                    let mut z = lower_corner[2];
                    while z <= upper_corner[2] + VSMALL {
                        let wz = if (z - lower_corner[2]).abs() < VSMALL
                            || (z - upper_corner[2]).abs() < VSMALL
                        {
                            0.5
                        } else {
                            1.0
                        };

                        let w = wx * wy * wz;
                        sav += w * (1.0 - self.ivdw_acc([x, y, z], radius));

                        z += spacs[2];
                    } // z loop

                    y += spacs[1];
                } // y loop

                x += spacs[0];
            } // x loop

            let w = spacs[0] * spacs[1] * spacs[2];
            sav *= w;

            sav
        }

        /// Per-atom WCA (Weeks-Chandler-Andersen) energy.
        /// Port of Vacc_wcaEnergyAtom from vacc.c line 1578.
        /// Computes the WCA energy for a single atom by numerical integration
        /// over a local grid patch around the atom.
        pub fn wca_energy_atom(&self, apolparm: &APOLparm, iatom: usize) -> Result<f64, ApbsError> {
            let pad = 14.0f64;
            // let vol_density = 2.0f64;

            // let lower_corner = self.clist.lower_corner;
            // let upper_corner = self.clist.upper_corner;

            let atom = self.alist.get_atom(iatom);
            let pos = atom.position;

            let srad = apolparm.srad;
            let rho = apolparm.bconc;
            let watsigma = apolparm.watsigma;
            let watepsilon = apolparm.watepsilon;
            let psig = atom.radius;
            let epsilon_atom = atom.epsilon;
            let sigma = psig + watsigma;
            let epsilon = (epsilon_atom * watepsilon).sqrt();

            // LJ parameters
            let sigma6 = sigma.powi(6);
            let sigma12 = sigma.powi(12);

            // Local grid bounds around atom
            let xmin = pos[0] - pad;
            let xmax = pos[0] + pad;
            let ymin = pos[1] - pad;
            let ymax = pos[1] + pad;
            let zmin = pos[2] - pad;
            let zmax = pos[2] + pad;

            // Compute grid spacings using extended domain
            let mut spacs = [0.5f64; 3];
            for i in 0..3 {
                // let len = (upper_corner[i] + pad) - (lower_corner[i] - pad);
                if apolparm.setgrid {
                    if apolparm.grid[i] > spacs[i] {
                        // Warning: grid value larger than recommended
                    }
                    spacs[i] = apolparm.grid[i];
                }
            }

            let mut energy = 0.0f64;

            // Numerical integration over local grid
            let mut x = xmin;
            while x <= xmax + VSMALL {
                let wx = if (x - xmin).abs() < VSMALL || (x - xmax).abs() < VSMALL {
                    0.5
                } else {
                    1.0
                };

                let mut y = ymin;
                while y <= ymax + VSMALL {
                    let wy = if (y - ymin).abs() < VSMALL || (y - ymax).abs() < VSMALL {
                        0.5
                    } else {
                        1.0
                    };

                    let mut z = zmin;
                    while z <= zmax + VSMALL {
                        let wz = if (z - zmin).abs() < VSMALL || (z - zmax).abs() < VSMALL {
                            0.5
                        } else {
                            1.0
                        };

                        let w = wx * wy * wz;
                        let vec = [x, y, z];

                        let chi = self.ivdw_acc(vec, srad);

                        let mut eni = 0.0f64;
                        if chi.abs() > VSMALL {
                            let x2 = (vec[0] - pos[0]).powi(2);
                            let y2 = (vec[1] - pos[1]).powi(2);
                            let z2 = (vec[2] - pos[2]).powi(2);
                            let r = (x2 + y2 + z2).sqrt();

                            if r <= 14.0 && r >= sigma {
                                eni = chi * rho * epsilon
                                    * (-2.0 * sigma6 / r.powi(6) + sigma12 / r.powi(12));
                            } else if r <= 14.0 {
                                eni = -1.0 * epsilon * chi * rho;
                            } else {
                                eni = 0.0;
                            }
                        }

                        energy += eni * w;

                        z += spacs[2];
                    } // z loop

                    y += spacs[1];
                } // y loop

                x += spacs[0];
            } // x loop

            let w = spacs[0] * spacs[1] * spacs[2];
            energy *= w;

            Ok(energy)
        }

        /// Total WCA (Weeks-Chandler-Andersen) energy.
        /// Port of Vacc_wcaEnergy from vacc.c line 1719.
        /// Sums per-atom WCA energies. Stores result in apolparm.wca_energy.
        pub fn wca_energy(&self, apolparm: &mut APOLparm) -> Result<(), ApbsError> {
            let rho = apolparm.bconc;

            // Sanity check: watsigma and watepsilon must be set
            // In the Rust APOLparm, these always have default values, so we
            // check they are non-zero as a proxy for "set"

            if rho.abs() < VSMALL {
                apolparm.wca_energy = 0.0;
                return Ok(());
            }

            let mut tenergy = 0.0f64;
            let num_atoms = self.alist.number_atoms();

            for iatom in 0..num_atoms {
                let energy = self.wca_energy_atom(apolparm, iatom)?;
                tenergy += energy;
            }

            apolparm.wca_energy = tenergy;

            Ok(())
        }

        /// Helper: compute SASA for an atom at a given position (without modifying the stored atom).
        /// Used by atomd_sasa for finite-difference derivatives.
        fn atom_sasa_at_pos(&self, atom_idx: usize, pos: [f64; 3], prad: f64) -> f64 {
            let atom = self.alist.get_atom(atom_idx);
            let mut tmp = atom.clone();
            tmp.position = pos;
            let surf = self.atom_surf(&tmp, &self.ref_sphere, prad);
            surf.area
        }

        /// Derivative of SAV (solvent accessible volume) with respect to atomic displacement.
        /// Port of Vacc_atomdSAV from vacc.c line 1198.
        pub fn atomd_sav(&self, srad: f64, atom_idx: usize, dsa: &mut [f64; 3]) {
            let ref_sphere = &self.ref_sphere;
            let atom = self.alist.get_atom(atom_idx);
            let iatom = atom.atom_id();
            let t_pos = atom.position;
            let t_rad = atom.radius;

            dsa[0] = 0.0;
            dsa[1] = 0.0;
            dsa[2] = 0.0;

            if t_rad == 0.0 {
                return;
            }

            let area = 4.0 * std::f64::consts::PI * (t_rad + srad) * (t_rad + srad)
                / (ref_sphere.npts as f64);

            let mut dx = 0.0;
            let mut dy = 0.0;
            let mut dz = 0.0;

            for ipt in 0..ref_sphere.npts {
                let vec = [
                    (t_rad + srad) * ref_sphere.xpts[ipt] + t_pos[0],
                    (t_rad + srad) * ref_sphere.ypts[ipt] + t_pos[1],
                    (t_rad + srad) * ref_sphere.zpts[ipt] + t_pos[2],
                ];
                if self.ivdw_acc_exclus(vec, srad, iatom) {
                    dx += vec[0] - t_pos[0];
                    dy += vec[1] - t_pos[1];
                    dz += vec[2] - t_pos[2];
                }
            }

            let denom = t_rad + srad;
            if denom != 0.0 {
                dsa[0] = dx * area / denom;
                dsa[1] = dy * area / denom;
                dsa[2] = dz * area / denom;
            }
        }

        /// Derivative of SASA (solvent accessible surface area) via finite differences.
        /// Port of Vacc_atomdSASA from vacc.c line 1318.
        pub fn atomd_sasa(&self, dpos: f64, srad: f64, atom_idx: usize, dsa: &mut [f64; 3]) {
            let atom = self.alist.get_atom(atom_idx);
            let t_pos = atom.position;

            // Shift by -/+ on x
            let axb1 = self.atom_sasa_at_pos(atom_idx, [t_pos[0] - dpos, t_pos[1], t_pos[2]], srad);
            let axt1 = self.atom_sasa_at_pos(atom_idx, [t_pos[0] + dpos, t_pos[1], t_pos[2]], srad);

            // Shift by -/+ on y
            let ayb1 = self.atom_sasa_at_pos(atom_idx, [t_pos[0], t_pos[1] - dpos, t_pos[2]], srad);
            let ayt1 = self.atom_sasa_at_pos(atom_idx, [t_pos[0], t_pos[1] + dpos, t_pos[2]], srad);

            // Shift by -/+ on z
            let azb1 = self.atom_sasa_at_pos(atom_idx, [t_pos[0], t_pos[1], t_pos[2] - dpos], srad);
            let azt1 = self.atom_sasa_at_pos(atom_idx, [t_pos[0], t_pos[1], t_pos[2] + dpos], srad);

            // Finite difference
            dsa[0] = (axt1 - axb1) / (2.0 * dpos);
            dsa[1] = (ayt1 - ayb1) / (2.0 * dpos);
            dsa[2] = (azt1 - azb1) / (2.0 * dpos);
        }

        /// WCA (Weeks-Chandler-Andersen) force on a single atom.
        /// Port of Vacc_wcaForceAtom from vacc.c line 1754.
        pub fn wca_force_atom(&self, apolparm: &APOLparm, iatom: usize, force: &mut [f64; 3]) -> Result<(), ApbsError> {
            let pad = 14.0f64;
            // let vol_density = 2.0f64;

            let lower_corner = self.clist.lower_corner;
            let upper_corner = self.clist.upper_corner;

            let atom = self.alist.get_atom(iatom);
            let pos = atom.position;

            let rho = apolparm.bconc;
            let watsigma = apolparm.watsigma;
            let watepsilon = apolparm.watepsilon;
            let psig = atom.radius;
            let epsilon_atom = atom.epsilon;
            let sigma = psig + watsigma;
            let epsilon = (epsilon_atom * watepsilon).sqrt();

            // LJ parameters
            let sigma6 = sigma.powi(6);
            let sigma12 = sigma.powi(12);

            // Grid spacings
            let mut spacs = [0.5f64; 3];
            for i in 0..3 {
                let _len = (upper_corner[i] + pad) - (lower_corner[i] - pad);
                if apolparm.setgrid {
                    if apolparm.grid[i] > spacs[i] {
                        // Warning: grid value larger than recommended
                    }
                    spacs[i] = apolparm.grid[i];
                }
            }

            // Local grid bounds around atom
            let xmin = pos[0] - pad;
            let xmax = pos[0] + pad;
            let ymin = pos[1] - pad;
            let ymax = pos[1] + pad;
            let zmin = pos[2] - pad;
            let zmax = pos[2] + pad;

            force[0] = 0.0;
            force[1] = 0.0;
            force[2] = 0.0;

            // Numerical integration over local grid
            let mut x = xmin;
            while x <= xmax + VSMALL {
                let wx = if (x - xmin).abs() < VSMALL || (x - xmax).abs() < VSMALL {
                    0.5
                } else {
                    1.0
                };

                let mut y = ymin;
                while y <= ymax + VSMALL {
                    let wy = if (y - ymin).abs() < VSMALL || (y - ymax).abs() < VSMALL {
                        0.5
                    } else {
                        1.0
                    };

                    let mut z = zmin;
                    while z <= zmax + VSMALL {
                        let wz = if (z - zmin).abs() < VSMALL || (z - zmax).abs() < VSMALL {
                            0.5
                        } else {
                            1.0
                        };

                        let w = wx * wy * wz;
                        let vec = [x, y, z];

                        let chi = self.ivdw_acc(vec, apolparm.srad);

                        let mut fpt = [0.0f64; 3];
                        if chi.abs() > VSMALL {
                            let x2 = (vec[0] - pos[0]).powi(2);
                            let y2 = (vec[1] - pos[1]).powi(2);
                            let z2 = (vec[2] - pos[2]).powi(2);
                            let r = (x2 + y2 + z2).sqrt();

                            if r <= 14.0 && r >= sigma {
                                let fo = 12.0 * chi * rho * epsilon
                                    * (sigma6 / r.powi(7) - sigma12 / r.powi(13));
                                fpt[0] = -1.0 * (pos[0] - vec[0]) * fo / r;
                                fpt[1] = -1.0 * (pos[1] - vec[1]) * fo / r;
                                fpt[2] = -1.0 * (pos[2] - vec[2]) * fo / r;
                            }
                        }

                        force[0] += w * fpt[0];
                        force[1] += w * fpt[1];
                        force[2] += w * fpt[2];

                        z += spacs[2];
                    } // z loop

                    y += spacs[1];
                } // y loop

                x += spacs[0];
            } // x loop

            let w = spacs[0] * spacs[1] * spacs[2];
            force[0] *= w;
            force[1] *= w;
            force[2] *= w;

            Ok(())
        }
    }
}
