// APBS PMGC lapack - Banded linear solver (dpbfa/dpbsl port)
// Port of pmgc/mlinpckd.c Vdpbfa/Vdpbsl (LINPACK dpbfa/dpbsl).
//
// Banded storage follows LINPACK dpbfa conventions: column j of the
// symmetric band matrix occupies abd[(j-1)*lda .. (j-1)*lda+m], where
// abd[(j-1)*lda + m] is the diagonal entry A(j,j) and
// abd[(j-1)*lda + m - l] is A(j-l, j) for l = 1..m.

/// LINPACK ddot with incx == incy == 1, reproducing the unrolled
/// summation order of Vddot (m = n % 5 remainder, then 5 at a time).
fn ddot_unrolled5(n: usize, dx: &[f64], dy: &[f64]) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut dtemp = 0.0;
    let m = n % 5;
    for i in 0..m {
        dtemp += dx[i] * dy[i];
    }
    if n < 5 {
        return dtemp;
    }
    let mut i = m;
    while i < n {
        dtemp += dx[i] * dy[i]
            + dx[i + 1] * dy[i + 1]
            + dx[i + 2] * dy[i + 2]
            + dx[i + 3] * dy[i + 3]
            + dx[i + 4] * dy[i + 4];
        i += 5;
    }
    dtemp
}

/// Factor a symmetric positive definite banded matrix (Cholesky).
/// Returns 0 on success, or the 1-based column index where a
/// non-positive leading minor was detected (matches Vdpbfa's `info`).
pub fn dpbfa(abd: &mut [f64], lda: usize, n: usize, m: usize) -> i32 {
    if lda < m + 1 {
        return -1;
    }
    if n == 0 {
        return 0;
    }

    for j in 1..=n {
        let col_j = (j - 1) * lda;
        let mut s = 0.0f64;
        let ik0 = m + 1;
        let mut jk = j.saturating_sub(m).max(1);
        let mu = (m + 2).saturating_sub(j).max(1);

        if m >= mu {
            let mut ik = ik0;
            for k in mu..=m {
                // t = abd(k,j) - dot(abd(ik.., jk), abd(mu.., j))
                let dx = &abd[(jk - 1) * lda + (ik - 1)..(jk - 1) * lda + (ik - 1) + (k - mu)];
                let dy = &abd[col_j + (mu - 1)..col_j + (mu - 1) + (k - mu)];
                let t = (abd[col_j + (k - 1)] - ddot_unrolled5(k - mu, dx, dy))
                    / abd[(jk - 1) * lda + m];
                abd[col_j + (k - 1)] = t;
                s += t * t;
                ik = ik.saturating_sub(1);
                jk += 1;
            }
        }

        s = abd[col_j + m] - s;
        if s <= 0.0 {
            return j as i32;
        }
        abd[col_j + m] = s.sqrt();
    }
    0
}

/// Solve A * x = b given the dpbfa factorization.
pub fn dpbsl(abd: &[f64], lda: usize, n: usize, m: usize, b: &mut [f64]) {
    if n == 0 {
        return;
    }

    // Solve L * y = b
    for k in 1..=n {
        let lm = (k - 1).min(m);
        let la = m + 1 - lm;
        let lb = k - lm;
        let dx = &abd[(k - 1) * lda + (la - 1)..(k - 1) * lda + (la - 1) + lm];
        let t = ddot_unrolled5(lm, dx, &b[lb - 1..lb - 1 + lm]);
        b[k - 1] = (b[k - 1] - t) / abd[(k - 1) * lda + m];
    }

    // Solve L^T * x = y
    for kb in 1..=n {
        let k = n + 1 - kb;
        let lm = (k - 1).min(m);
        let la = m + 1 - lm;
        let lb = k - lm;
        b[k - 1] /= abd[(k - 1) * lda + m];
        let t = -b[k - 1];
        for i in 0..lm {
            b[lb - 1 + i] += t * abd[(k - 1) * lda + (la - 1) + i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{dpbfa, dpbsl};

    /// A = [[4,1,0],[1,5,2],[0,2,6]] with m = 1, lda = 2.
    #[test]
    fn dpbfa_dpbsl_solve_small_banded_system() {
        let (n, m, lda) = (3usize, 1usize, 2usize);
        // column 0: [_, 4]; column 1: [1, 5]; column 2: [2, 6]
        let mut abd = vec![0.0, 4.0, 1.0, 5.0, 2.0, 6.0];
        let info = dpbfa(&mut abd, lda, n, m);
        assert_eq!(info, 0);

        let mut b = vec![6.0, 17.0, 22.0]; // A * [1, 2, 3]
        dpbsl(&abd, lda, n, m, &mut b);
        assert!((b[0] - 1.0).abs() < 1.0e-12);
        assert!((b[1] - 2.0).abs() < 1.0e-12);
        assert!((b[2] - 3.0).abs() < 1.0e-12);
    }

    #[test]
    fn dpbfa_detects_non_positive_definite() {
        let (n, m, lda) = (2usize, 1usize, 2usize);
        // column 0: [_, -1]  (negative diagonal)
        let mut abd = vec![0.0, -1.0, 0.5, 1.0];
        let info = dpbfa(&mut abd, lda, n, m);
        assert_eq!(info, 1);
    }
}
