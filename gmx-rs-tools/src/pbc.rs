//! Periodic boundary condition helpers, mirroring the parts of
//! `gromacs/pbcutil/pbc.cpp` and `gromacs/pbcutil/pbcmethods.cpp` used by
//! `gmx trjconv`.

use crate::frame::{Matrix, PbcType, Rvec};

/// Box centre selection, mirroring the `ecenter*` values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ECenter {
    Tric,
    Rect,
    Zero,
}

impl ECenter {
    pub fn from_name(name: &str) -> ECenter {
        match name {
            "rect" => ECenter::Rect,
            "zero" => ECenter::Zero,
            _ => ECenter::Tric,
        }
    }
}

pub fn calc_box_center(ecenter: ECenter, boxm: &Matrix) -> Rvec {
    let mut c = [0.0f32; 3];
    match ecenter {
        ECenter::Tric => {
            for m in 0..3 {
                for d in 0..3 {
                    c[d] += 0.5 * boxm[m][d];
                }
            }
        }
        ECenter::Rect => {
            for d in 0..3 {
                c[d] = 0.5 * boxm[d][d];
            }
        }
        ECenter::Zero => {}
    }
    c
}

pub fn npbc_dims(pbc_type: PbcType) -> usize {
    if pbc_type == PbcType::XY {
        2
    } else {
        3
    }
}

/// `put_atoms_in_box`: puts all atoms inside the unit cell (rect/tric aware).
pub fn put_atoms_in_box(pbc_type: PbcType, boxm: &Matrix, x: &mut [Rvec]) {
    let npbcdim = npbc_dims(pbc_type);
    let mut inv = [0.0f32; 3];
    for m in 0..npbcdim {
        if boxm[m][m] != 0.0 {
            inv[m] = 1.0 / boxm[m][m];
        }
    }
    let triclinic = crate::frame::is_triclinic(boxm);
    for xi in x.iter_mut() {
        if triclinic {
            for m in (0..npbcdim).rev() {
                let shift = (xi[m] * inv[m]).floor();
                for d in 0..=m {
                    xi[d] -= shift * boxm[m][d];
                }
            }
        } else {
            for d in 0..npbcdim {
                let shift = (xi[d] * inv[d]).floor();
                xi[d] -= shift * boxm[d][d];
            }
        }
    }
}

/// `put_atoms_in_triclinic_unitcell()`.
pub fn put_atoms_in_triclinic_unitcell(ecenter: ECenter, boxm: &Matrix, x: &mut [Rvec]) {
    let box_center = calc_box_center(ecenter, boxm);
    let shm01 = boxm[1][0] / boxm[1][1];
    let shm02 = (boxm[1][1] * boxm[2][0] - boxm[2][1] * boxm[1][0]) / (boxm[1][1] * boxm[2][2]);
    let shm12 = boxm[2][1] / boxm[2][2];

    let mut shift_center = [0.0f32; 3];
    for d in 0..3 {
        for m in 0..3 {
            shift_center[d] += boxm[m][d];
        }
    }
    for d in 0..3 {
        shift_center[d] *= 0.5;
    }
    for d in 0..3 {
        shift_center[d] = box_center[d] - shift_center[d];
    }

    shift_center[0] = shm01 * shift_center[1] + shm02 * shift_center[2];
    shift_center[1] = shm12 * shift_center[2];
    shift_center[2] = 0.0;

    for xi in x.iter_mut() {
        for m in (0..3).rev() {
            let mut shift = shift_center[m];
            if m == 0 {
                shift += shm01 * xi[1] + shm02 * xi[2];
            } else if m == 1 {
                shift += shm12 * xi[2];
            }
            while xi[m] - shift < 0.0 {
                for d in 0..=m {
                    xi[d] += boxm[m][d];
                }
            }
            while xi[m] - shift >= boxm[m][m] {
                for d in 0..=m {
                    xi[d] -= boxm[m][d];
                }
            }
        }
    }
}

/// Minimum image distance `a - b`.
pub fn pbc_dx(pbc_type: PbcType, boxm: &Matrix, a: &Rvec, b: &Rvec) -> Rvec {
    let mut dx = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    if pbc_type == PbcType::No || pbc_type == PbcType::Unset {
        return dx;
    }
    let npbcdim = npbc_dims(pbc_type);
    if crate::frame::is_triclinic(boxm) {
        for m in (0..npbcdim).rev() {
            if boxm[m][m] != 0.0 {
                let sh = (dx[m] / boxm[m][m] + 0.5).floor();
                for d in 0..=m {
                    dx[d] -= sh * boxm[m][d];
                }
            }
        }
    } else {
        for d in 0..npbcdim {
            if boxm[d][d] != 0.0 {
                let sh = (dx[d] / boxm[d][d] + 0.5).floor();
                dx[d] -= sh * boxm[d][d];
            }
        }
    }
    dx
}

/// `put_atoms_in_compact_unitcell()`.
pub fn put_atoms_in_compact_unitcell(pbc_type: PbcType, ecenter: ECenter, boxm: &Matrix, x: &mut [Rvec]) {
    let box_center = calc_box_center(ecenter, boxm);
    // set_pbc() precomputes the extra triclinic trial vectors.
    let tric_vec = triclinic_shift_vectors(pbc_type, boxm);
    for xi in x.iter_mut() {
        let dx = pbc_dx_gmx(pbc_type, boxm, xi, &box_center, &tric_vec);
        for d in 0..3 {
            xi[d] = box_center[d] + dx[d];
        }
    }
}

/// Distance vector to the *closest* periodic image.
/// `sc_skewnessMargin` from `pbc.cpp`.
const SKEWNESS_MARGIN: f32 = 1.001;
/// `MAX_NTRICVEC` from `pbc.h`.
const MAX_NTRIC_VEC: usize = 12;

/// `max_cutoff2()`: the largest squared distance for which the plain minimum
/// image reduction is guaranteed to give the shortest vector.
pub fn max_cutoff2(pbc_type: PbcType, boxm: &Matrix) -> f32 {
    let norm2 = |v: &Rvec| v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    let mut min_hv2 = 0.25 * norm2(&boxm[0]).min(norm2(&boxm[1]));
    if pbc_type != PbcType::XY {
        min_hv2 = min_hv2.min(0.25 * norm2(&boxm[2]));
    }
    let min_ss = if pbc_type == PbcType::XY {
        boxm[0][0].min(boxm[1][1])
    } else {
        boxm[0][0].min(boxm[1][1] - boxm[2][1].abs()).min(boxm[2][2])
    };
    min_hv2.min(min_ss * min_ss)
}

/// Extra trial shift vectors computed by `low_set_pbc()` for triclinic boxes.
///
/// These make it possible to find the shortest vector in a truncated
/// octahedron, where the plain minimum image reduction is not sufficient.
pub fn triclinic_shift_vectors(pbc_type: PbcType, boxm: &Matrix) -> Vec<Rvec> {
    let mut out = Vec::new();
    if !crate::frame::is_triclinic(boxm) || pbc_type == PbcType::No || pbc_type == PbcType::Unset
    {
        return out;
    }
    let mut b_pbc = [1i32; 3];
    if pbc_type == PbcType::XY {
        b_pbc[2] = 0;
    }
    let npbcdim = b_pbc.iter().filter(|v| **v != 0).count();
    if npbcdim < 2 {
        return out;
    }
    // `order[]` is {0, -1, 1} in the C code, iterated in that order.
    let order = [0i32, -1, 1];
    let hbox = [0.5 * boxm[0][0], 0.5 * boxm[1][1], 0.5 * boxm[2][2]];
    'outer: for &k in order.iter() {
        if b_pbc[2] == 0 && k != 0 {
            continue;
        }
        for &j in order.iter() {
            if b_pbc[1] == 0 && j != 0 {
                continue;
            }
            for &i in order.iter() {
                if b_pbc[0] == 0 && i != 0 {
                    continue;
                }
                if j == 0 && k == 0 {
                    continue;
                }
                let mut trial = [0.0f32; 3];
                let mut pos = [0.0f32; 3];
                let mut d2old = 0.0f32;
                let mut d2new = 0.0f32;
                for d in 0..3 {
                    trial[d] = i as f32 * boxm[0][d]
                        + j as f32 * boxm[1][d]
                        + k as f32 * boxm[2][d];
                    pos[d] = if trial[d] < 0.0 {
                        hbox[d].min(-trial[d])
                    } else {
                        (-hbox[d]).max(-trial[d])
                    };
                    d2old += pos[d] * pos[d];
                    let s = pos[d] + trial[d];
                    d2new += s * s;
                }
                if SKEWNESS_MARGIN * d2new < d2old {
                    let mut b_use = true;
                    for (dd, shift) in [i, j, k].iter().enumerate() {
                        if *shift != 0 {
                            let mut d2new_c = 0.0f32;
                            for d in 0..3 {
                                let s = pos[d] + trial[d] - *shift as f32 * boxm[dd][d];
                                d2new_c += s * s;
                            }
                            if d2new_c <= SKEWNESS_MARGIN * d2new {
                                b_use = false;
                            }
                        }
                    }
                    if b_use {
                        if out.len() >= MAX_NTRIC_VEC {
                            continue 'outer;
                        }
                        out.push(trial);
                    }
                }
            }
        }
    }
    out
}

/// `pbc_dx()` from `gromacs/pbcutil/pbc.cpp`.
///
/// For triclinic boxes the plain minimum image is only guaranteed to be the
/// shortest vector when it is within `max_cutoff2`; otherwise the precomputed
/// `tric_vec` shifts are tried as well.
pub fn pbc_dx_gmx(pbc_type: PbcType, boxm: &Matrix, a: &Rvec, b: &Rvec, tric_vec: &[Rvec]) -> Rvec {
    let mut dx = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    if pbc_type == PbcType::No || pbc_type == PbcType::Unset {
        return dx;
    }
    let npbcdim = npbc_dims(pbc_type);
    let hbox = [0.5 * boxm[0][0], 0.5 * boxm[1][1], 0.5 * boxm[2][2]];
    let norm2 = |v: &Rvec| v[0] * v[0] + v[1] * v[1] + v[2] * v[2];

    if crate::frame::is_triclinic(boxm) {
        for i in (0..npbcdim).rev() {
            while dx[i] > hbox[i] {
                for j in (0..=i).rev() {
                    dx[j] -= boxm[i][j];
                }
            }
            while dx[i] <= -hbox[i] {
                for j in (0..=i).rev() {
                    dx[j] += boxm[i][j];
                }
            }
        }
        let mut d2min = norm2(&dx);
        let cutoff = max_cutoff2(pbc_type, boxm);
        if d2min > cutoff {
            let dx_start = dx;
            for t in tric_vec {
                let trial = [
                    dx_start[0] + t[0],
                    dx_start[1] + t[1],
                    dx_start[2] + t[2],
                ];
                let d2trial = norm2(&trial);
                if d2trial < d2min {
                    dx = trial;
                    d2min = d2trial;
                }
            }
        }
    } else {
        for i in 0..npbcdim {
            while dx[i] > hbox[i] {
                dx[i] -= boxm[i][i];
            }
            while dx[i] <= -hbox[i] {
                dx[i] += boxm[i][i];
            }
        }
    }
    dx
}

/// `center_x()`: shifts all atoms so the geometric centre of the selected
/// group ends up at the boxm centre.
/// `calc_pbc_cluster()` from `gmx trjconv`: makes the molecules of the
/// selected cluster whole, then shifts them one by one so that each newly
/// added molecule is the periodic image closest to a molecule that is already
/// part of the cluster.
///
/// `molecules` holds the atom indices of every molecule (one entry per
/// molecule, like `t_topology::mols`).
pub fn calc_pbc_cluster(
    ecenter: ECenter,
    x: &mut [Rvec],
    index_cluster: &[usize],
    molecules: &[Vec<usize>],
    pbc_type: PbcType,
    boxm: &Matrix,
) {
    let nmol = molecules.len();
    if nmol == 0 || index_cluster.is_empty() {
        return;
    }
    let box_center = calc_box_center(ecenter, boxm);
    let tric_vec = triclinic_shift_vectors(pbc_type, boxm);

    // Atom index -> molecule index (the C code binary searches `mols.index`).
    let mut atom_to_mol = vec![usize::MAX; x.len()];
    for (i, mol) in molecules.iter().enumerate() {
        for &a in mol {
            if a < atom_to_mol.len() {
                atom_to_mol[a] = i;
            }
        }
    }
    let mut b_mol = vec![false; nmol];
    let mut b_tmp = vec![false; x.len()];
    for &ai in index_cluster {
        if ai >= x.len() {
            continue;
        }
        b_tmp[ai] = true;
        let m = atom_to_mol[ai];
        if m != usize::MAX {
            b_mol[m] = true;
        }
    }

    let norm2 = |v: &Rvec| v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    let trace = boxm[0][0] + boxm[1][1] + boxm[2][2];
    let mut min_dist2 = 10.0 * trace * trace;
    let mut imol_center: Option<usize> = None;
    let mut cluster: Vec<usize> = Vec::new();
    let mut com = vec![[0.0f32; 3]; nmol];

    for i in 0..nmol {
        let mol = &molecules[i];
        for (k, &j) in mol.iter().enumerate() {
            if b_mol[i] && !b_tmp[j] {
                eprintln!(
                    "Molecule {} marked for clustering but not atom {} in it - check your index!",
                    i + 1,
                    j + 1
                );
                return;
            } else if !b_mol[i] && b_tmp[j] {
                eprintln!(
                    "Atom {} marked for clustering but not molecule {} - this is an internal error...",
                    j + 1,
                    i + 1
                );
                return;
            } else if b_mol[i] {
                if k > 0 {
                    // Make the molecule whole relative to its previous atom.
                    let prev = mol[k - 1];
                    let dx = pbc_dx_gmx(pbc_type, boxm, &x[j], &x[prev], &tric_vec);
                    for d in 0..3 {
                        x[j][d] = x[prev][d] + dx[d];
                    }
                }
                for d in 0..3 {
                    com[i][d] += x[j][d];
                }
            }
        }
        if b_mol[i] {
            let fac = 1.0 / mol.len() as f32;
            for d in 0..3 {
                com[i][d] *= fac;
            }
            // Which marked molecule is closest to the centre of the box?
            let dx = pbc_dx_gmx(pbc_type, boxm, &box_center, &com[i], &tric_vec);
            let r2 = norm2(&dx);
            if r2 < min_dist2 {
                min_dist2 = r2;
                imol_center = Some(i);
            }
            cluster.push(i);
        }
    }

    if cluster.is_empty() {
        eprintln!("No molecules selected in the cluster");
        return;
    }
    let Some(center_mol) = imol_center else {
        eprintln!("No central molecules could be found");
        return;
    };

    let ncluster = cluster.len();
    let mut added: Vec<usize> = vec![center_mol];
    b_mol[center_mol] = false;

    while added.len() < ncluster {
        // Find the closest pair of an already added and a remaining molecule.
        min_dist2 = 10.0 * trace * trace;
        let mut bimin: Option<usize> = None;
        let mut bjmin: Option<usize> = None;
        for &ai in &added {
            for &aj in &cluster {
                if !b_mol[aj] {
                    continue;
                }
                let dx = pbc_dx_gmx(pbc_type, boxm, &com[aj], &com[ai], &tric_vec);
                let r2 = norm2(&dx);
                if r2 < min_dist2 {
                    min_dist2 = r2;
                    bimin = Some(ai);
                    bjmin = Some(aj);
                }
            }
        }
        let (Some(imin), Some(jmin)) = (bimin, bjmin) else {
            break;
        };

        added.push(jmin);
        b_mol[jmin] = false;

        // Shift the new molecule to the image closest to `imin`.
        let dx = pbc_dx_gmx(pbc_type, boxm, &com[jmin], &com[imin], &tric_vec);
        let xtest = [
            com[imin][0] + dx[0],
            com[imin][1] + dx[1],
            com[imin][2] + dx[2],
        ];
        let shift = [
            xtest[0] - com[jmin][0],
            xtest[1] - com[jmin][1],
            xtest[2] - com[jmin][2],
        ];
        for d in 0..3 {
            com[jmin][d] += shift[d];
        }
        for &j in &molecules[jmin] {
            for d in 0..3 {
                x[j][d] += shift[d];
            }
        }
        print!("\rClustering iteration {} of {}...", added.len(), ncluster);
        use std::io::Write;
        let _ = std::io::stdout().flush();
    }
    println!();
}

pub fn center_x(ecenter: ECenter, x: &mut [Rvec], boxm: &Matrix, selected: &[usize]) {
    if selected.is_empty() {
        return;
    }
    let first = x[selected[0]];
    let mut cmin = first;
    let mut cmax = first;
    for &ai in selected {
        for m in 0..3 {
            if x[ai][m] < cmin[m] {
                cmin[m] = x[ai][m];
            } else if x[ai][m] > cmax[m] {
                cmax[m] = x[ai][m];
            }
        }
    }
    let box_center = calc_box_center(ecenter, boxm);
    let mut dx = [0.0f32; 3];
    for m in 0..3 {
        dx[m] = box_center[m] - (cmin[m] + cmax[m]) * 0.5;
    }
    for xi in x.iter_mut() {
        for m in 0..3 {
            xi[m] += dx[m];
        }
    }
}

/// `reset_x_ndim()`: removes the weighted centre of mass of the `ind_cm`
/// atoms from all atoms.  `mass` is indexed by atom and is zero for atoms that
/// should not contribute.
pub fn reset_x_ndim(
    ndim: usize,
    ncm: usize,
    ind_cm: &[usize],
    nreset: usize,
    ind_reset: Option<&[usize]>,
    x: &mut [Rvec],
    mass: &[f64],
) -> Rvec {
    // `reset_x_ndim()` accumulates the centre of mass in `real` (float) and
    // divides by a float total mass; the rounding is observable in the last
    // printed decimal (including the sign of zero) of the output.
    let mut xcm = [0.0f32; 3];
    let mut tm = 0.0f32;
    if ind_cm.is_empty() {
        for i in 0..ncm {
            let mm = mass[i] as f32;
            for m in 0..ndim {
                xcm[m] += mm * x[i][m];
            }
            tm += mm;
        }
    } else {
        for &ai in ind_cm.iter().take(ncm) {
            let mm = mass[ai] as f32;
            for m in 0..ndim {
                xcm[m] += mm * x[ai][m];
            }
            tm += mm;
        }
    }
    if tm != 0.0 {
        for m in 0..ndim {
            xcm[m] /= tm;
        }
    }
    match ind_reset {
        Some(list) => {
            for &i in list.iter().take(nreset) {
                for m in 0..ndim {
                    x[i][m] -= xcm[m];
                }
            }
        }
        None => {
            for i in 0..nreset {
                for m in 0..ndim {
                    x[i][m] -= xcm[m];
                }
            }
        }
    }
    xcm
}

fn do_rotate(a: &mut [Vec<f64>], i: usize, j: usize, k: usize, l: usize, tau: f64, s: f64) {
    let g = a[i][j];
    let h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

/// Jacobi eigen decomposition of a symmetric matrix (`jacobi()` in
/// `gromacs/math/nrjac.cpp`).  Returns eigenvalues and eigenvectors stored in
/// the columns of the second element.
fn jacobi(n: usize, a_in: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut a = a_in.to_vec();
    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    let mut b: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    let mut d: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    let mut z = vec![0.0f64; n];

    for iter in 1..=50 {
        let mut sm = 0.0;
        for ip in 0..n - 1 {
            for iq in ip + 1..n {
                sm += a[ip][iq].abs();
            }
        }
        if sm == 0.0 {
            return (d, v);
        }
        let tresh = if iter < 4 { 0.2 * sm / (n * n) as f64 } else { 0.0 };
        for ip in 0..n - 1 {
            for iq in ip + 1..n {
                let g = 100.0 * a[ip][iq].abs();
                if iter > 4
                    && (d[ip].abs() + g == d[ip].abs())
                    && (d[iq].abs() + g == d[iq].abs())
                {
                    a[ip][iq] = 0.0;
                } else if a[ip][iq].abs() > tresh {
                    let mut h = d[iq] - d[ip];
                    let t = if h.abs() + g == h.abs() {
                        a[ip][iq] / h
                    } else {
                        let theta = 0.5 * h / a[ip][iq];
                        let mut t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                        if theta < 0.0 {
                            t = -t;
                        }
                        t
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    let tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for j in 0..ip {
                        do_rotate(&mut a, j, ip, j, iq, tau, s);
                    }
                    for j in ip + 1..iq {
                        do_rotate(&mut a, ip, j, j, iq, tau, s);
                    }
                    for j in iq + 1..n {
                        do_rotate(&mut a, ip, j, iq, j, tau, s);
                    }
                    for j in 0..n {
                        do_rotate(&mut v, j, ip, j, iq, tau, s);
                    }
                }
            }
        }
        for ip in 0..n {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    (d, v)
}

/// `calc_fit_R()`: quaternion based optimal rotation matrix.
pub fn calc_fit_r(
    ndim: usize,
    natoms: usize,
    w_rls: &[f64],
    xp: &[Rvec],
    x: &[Rvec],
) -> [[f32; 3]; 3] {
    let two_n = 2 * ndim;
    // GROMACS keeps `u`, `vh`, `vk` and `R` in `real` (float) precision and
    // only the eigenvalue problem in double precision.
    let mut u = [[0.0f32; 3]; 3];
    for n in 0..natoms {
        let mn = w_rls[n] as f32;
        if mn != 0.0 {
            for c in 0..ndim {
                let xpc = xp[n][c];
                for r in 0..ndim {
                    let xnr = x[n][r];
                    // The C code keeps `xnr`/`xpc` in double, so the product
                    // is evaluated in double and only the accumulator is
                    // `real` (float).
                    u[c][r] += mn * xnr * xpc;
                }
            }
        }
    }

    let mut omega = vec![vec![0.0f64; two_n]; two_n];
    for r in 0..two_n {
        for c in 0..=r {
            if r >= ndim && c < ndim {
                omega[r][c] = u[r - ndim][c] as f64;
                omega[c][r] = u[r - ndim][c] as f64;
            }
        }
    }

    let (mut d, om) = jacobi(two_n, &omega);
    let mut vh = [[0.0f32; 3]; 3];
    let mut vk = [[0.0f32; 3]; 3];
    let sqrt2 = std::f64::consts::SQRT_2;
    for j in 0..ndim - 1 {
        let mut max_d = -1000.0;
        let mut index = 0usize;
        for i in 0..two_n {
            if d[i] > max_d {
                max_d = d[i];
                index = i;
            }
        }
        d[index] = -10000.0;
        for i in 0..ndim {
            vh[j][i] = (sqrt2 * om[i][index]) as f32;
            vk[j][i] = (sqrt2 * om[i + ndim][index]) as f32;
        }
    }
    if ndim == 3 {
        vh[2] = cprod_f32(vh[0], vh[1]);
        vk[2] = cprod_f32(vk[0], vk[1]);
    } else if ndim == 2 {
        vh[1][0] = -vh[0][1];
        vh[1][1] = vh[0][0];
        vk[1][0] = -vk[0][1];
        vk[1][1] = vk[0][0];
    }

    let mut r = [[0.0f32; 3]; 3];
    for row in 0..ndim {
        for c in 0..ndim {
            for s in 0..ndim {
                r[row][c] += vk[s][row] * vh[s][c];
            }
        }
    }
    for i in ndim..3 {
        r[i][i] = 1.0;
    }
    r
}

fn cprod_f32(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `do_fit_ndim()`.
pub fn do_fit_ndim(ndim: usize, natoms: usize, w_rls: &[f64], xp: &[Rvec], x: &mut [Rvec]) {
    let r = calc_fit_r(ndim, natoms, w_rls, xp, x);
    for j in 0..natoms {
        let old = x[j];
        for row in 0..3 {
            let mut acc = 0.0f32;
            for c in 0..3 {
                acc += r[row][c] * old[c];
            }
            x[j][row] = acc;
        }
    }
}

/// `do_fit()`.
pub fn do_fit(natoms: usize, w_rls: &[f64], xp: &[Rvec], x: &mut [Rvec]) {
    do_fit_ndim(3, natoms, w_rls, xp, x);
}
