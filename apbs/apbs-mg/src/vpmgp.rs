// APBS vpmgp.rs - Multigrid solver parameters
// Port of src/mg/vpmgp.h / vpmgp.c

use apbs_generic::mgparm::MGparm;
use apbs_generic::vhal::Vbcfl;

/// Multigrid solver parameters
#[derive(Debug, Clone)]
pub struct Vpmgp {
    // User-specified
    pub nx: i32,
    pub ny: i32,
    pub nz: i32,
    pub nlev: i32,
    pub hx: f64,
    pub hy: f64,
    pub hzed: f64,
    pub nonlin: i32,

    // Derived
    pub nxc: i32,
    pub nyc: i32,
    pub nzc: i32,
    pub nf: i32,
    pub nc: i32,
    pub narrc: i32,
    pub n_rpc: i32,
    pub n_iz: i32,
    pub n_ipc: i32,
    pub nrwk: usize,
    pub niwk: i32,
    pub narr: i32,
    pub ipkey: i32,

    // Parameters with defaults
    pub xcent: f64,
    pub ycent: f64,
    pub zcent: f64,
    pub errtol: f64,
    pub itmax: i32,
    pub istop: i32,
    pub iinfo: i32,
    pub bcfl: Vbcfl,
    pub key: i32,
    pub iperf: i32,
    pub meth: i32,
    pub mgkey: i32,
    pub nu1: i32,
    pub nu2: i32,
    pub mgsmoo: i32,
    pub mgprol: i32,
    pub mgcoar: i32,
    pub mgsolv: i32,
    pub mgdisc: i32,
    pub omegal: f64,
    pub omegan: f64,
    pub irite: i32,
    pub ipcon: i32,

    // Domain
    pub xlen: f64,
    pub ylen: f64,
    pub zlen: f64,
    pub xmin: f64,
    pub ymin: f64,
    pub zmin: f64,
    pub xmax: f64,
    pub ymax: f64,
    pub zmax: f64,
}

impl Vpmgp {
    pub fn effective_mgcoar(&self) -> i32 {
        self.mgcoar
    }

    pub fn new(mgparm: &MGparm) -> Self {
        let nx = mgparm.dime[0];
        let ny = mgparm.dime[1];
        let nz = mgparm.dime[2];

        let hx = mgparm.grid[0];
        let hy = mgparm.grid[1];
        let hzed = mgparm.grid[2];

        // Match APBS C: derive actual domain lengths from spacing and counts
        // rather than trusting the user-provided `glen` literals.
        let xlen = (nx - 1) as f64 * hx;
        let ylen = (ny - 1) as f64 * hy;
        let zlen = (nz - 1) as f64 * hzed;

        let xcent = mgparm.center[0];
        let ycent = mgparm.center[1];
        let zcent = mgparm.center[2];

        let xmin = xcent - xlen / 2.0;
        let ymin = ycent - ylen / 2.0;
        let zmin = zcent - zlen / 2.0;
        let xmax = xcent + xlen / 2.0;
        let ymax = ycent + ylen / 2.0;
        let zmax = zcent + zlen / 2.0;

        let nf = nx * ny * nz;
        let max_nlev = Self::compute_nlev(nx, ny, nz);
        let nlev = if mgparm.nlev > 0 {
            mgparm.nlev.min(max_nlev)
        } else {
            max_nlev
        };

        let mut p = Self {
            nx, ny, nz,
            nlev,
            hx, hy, hzed,
            nonlin: 0,
            nxc: 0, nyc: 0, nzc: 0,
            nf,
            nc: 0,
            narrc: 0,
            n_rpc: 0, n_iz: 0, n_ipc: 0,
            nrwk: 0,
            niwk: 0,
            narr: nf,
            ipkey: 0,
            xcent, ycent, zcent,
            errtol: 1.0e-6,
            itmax: 200,
            istop: 1,
            iinfo: 1,
            bcfl: Vbcfl::SDH,
            key: 0,
            iperf: 0,
            meth: mgparm.method,
            mgkey: 0,
            nu1: 2,
            nu2: 2,
            mgsmoo: 1,
            mgprol: 0,
            mgcoar: 2,
            mgsolv: 1,
            mgdisc: 0,
            omegal: 1.0,
            omegan: 0.9,
            irite: 8,
            ipcon: 3,
            xlen, ylen, zlen,
            xmin, ymin, zmin,
            xmax, ymax, zmax,
        };

        p.size();
        p
    }

    /// Compute number of multigrid levels
    fn compute_nlev(nx: i32, ny: i32, nz: i32) -> i32 {
        let mut lev = 0;
        loop {
            lev += 1;
            let iden = 1_i32 << (lev - 1);

            let nxc = (nx - 1) / iden + 1;
            let nyc = (ny - 1) / iden + 1;
            let nzc = (nz - 1) / iden + 1;

            let done =
                ((nxc - 1) * iden != (nx - 1)) || (nxc <= 2) ||
                ((nyc - 1) * iden != (ny - 1)) || (nyc <= 2) ||
                ((nzc - 1) * iden != (nz - 1)) || (nzc <= 2);

            if done {
                return (lev - 1).max(1);
            }
        }
    }

    /// Compute array sizes
    pub fn size(&mut self) {
        // Compute coarse grid dimensions
        let mut nxc = self.nx;
        let mut nyc = self.ny;
        let mut nzc = self.nz;
        for _ in 1..self.nlev {
            nxc = nxc / 2 + 1;
            nyc = nyc / 2 + 1;
            nzc = nzc / 2 + 1;
        }
        self.nxc = nxc;
        self.nyc = nyc;
        self.nzc = nzc;

        // Compute total unknowns across all levels
        let mut narr = 0;
        let mut nx_l = self.nx;
        let mut ny_l = self.ny;
        let mut nz_l = self.nz;
        for _ in 0..self.nlev {
            narr += nx_l * ny_l * nz_l;
            nx_l = nx_l / 2 + 1;
            ny_l = ny_l / 2 + 1;
            nz_l = nz_l / 2 + 1;
        }
        self.narr = narr;
        self.narrc = narr - self.nf;
        self.nc = nxc * nyc * nzc;

        // Compute work array sizes (simplified)
        let numdia = 4i32; // 7-point stencil
        let narrv = narr;
        let narrc = self.narrc;

        self.n_ipc = self.nlev * 20;
        self.n_rpc = self.nlev * 20;
        self.n_iz = self.nlev * 50;

        // Operator storage: numdia diags for each level
        let nac = if self.effective_mgcoar() == 2 {
            14 * narrc + numdia * self.nf
        } else {
            numdia * narr
        };

        // Prolongation: 27 * narrc
        let npc = 27 * narrc;

        // Work arrays
        self.nrwk = (narrv + nac + npc + self.nc + self.nf) as usize;
        self.niwk = self.n_iz + self.n_ipc + narr + nac + npc + self.nc + self.nf;

        // Banded solver storage
        if self.mgsolv == 1 {
            let nb = (self.nx - 2) * (self.ny - 2);
            let nband = if self.mgdisc == 0 { nb + 1 } else { nb * 3 + 1 };
            self.nrwk += (nband * self.nc) as usize;
        }
    }

    /// Compute coarse grid dimensions
    pub fn make_coarse(_num_level: i32, nx_old: i32, ny_old: i32, nz_old: i32) -> (i32, i32, i32) {
        let nx_new = nx_old / 2 + 1;
        let ny_new = ny_old / 2 + 1;
        let nz_new = nz_old / 2 + 1;
        (nx_new, ny_new, nz_new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use apbs_generic::mgparm::MGparm;

    #[test]
    fn test_vpmgp_new() {
        let mut mgparm = MGparm::new(apbs_generic::mgparm::MGparmCalcType::Manual);
        mgparm.dime = [33, 33, 33];
        mgparm.glen = [10.0, 10.0, 10.0];
        mgparm.grid = [10.0 / 32.0; 3];
        mgparm.center = [0.0; 3];

        let pmgp = Vpmgp::new(&mgparm);
        assert_eq!(pmgp.nx, 33);
        assert_eq!(pmgp.ny, 33);
        assert_eq!(pmgp.nz, 33);
        assert!(pmgp.nlev >= 1);
    }

    #[test]
    fn test_effective_mgcoar_preserves_user_setting() {
        let mut mgparm = MGparm::new(apbs_generic::mgparm::MGparmCalcType::Manual);
        mgparm.dime = [33, 33, 33];
        mgparm.glen = [10.0, 10.0, 10.0];
        mgparm.grid = [10.0 / 32.0; 3];
        mgparm.center = [0.0; 3];

        let mut pmgp = Vpmgp::new(&mgparm);
        pmgp.mgcoar = 2;
        pmgp.mgdisc = 0;
        assert_eq!(pmgp.effective_mgcoar(), 2);

        pmgp.mgdisc = 1;
        assert_eq!(pmgp.effective_mgcoar(), 2);
    }

    #[test]
    fn test_compute_nlev_matches_apbs_c_vmaxlev_examples() {
        assert_eq!(Vpmgp::compute_nlev(97, 129, 129), 6);
        assert_eq!(Vpmgp::compute_nlev(129, 129, 129), 7);
        assert_eq!(Vpmgp::compute_nlev(33, 33, 33), 5);
    }
}
