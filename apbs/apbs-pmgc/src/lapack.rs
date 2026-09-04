// APBS PMGC lapack - Banded linear solver (dpbfa/dpbsl port)
// Port of LAPACK dpbfa/dpbsl for symmetric positive definite banded matrices

/// Symmetric positive definite banded matrix factorization (Cholesky)
/// Factorizes AB into L * L^T where L is lower triangular banded
/// Returns Ok(()) on success, Err on non-positive-definite
pub fn dpbfa(
    ab: &mut [f64], n: usize, m: usize, lda: usize,
) -> Result<(), usize> {
    let info;

    // Test the input parameters
    if lda < m + 1 {
        return Err(0);
    }

    // Quick return if possible
    if n == 0 {
        return Ok(());
    }

    info = 0;

    for j in 1..=n {
        // Compute the j-th column of L
        let mut t = ab[(j - 1) * lda + m];

        if j > 1 {
            // t = ab(j,m) - sum_{k=1}^{j-1} L(j,k) * L(j,k)^T
            for k in 0..j - 1 {
                let ljk = if k + m + 1 >= j - 1 && k + m + 1 - (j - 1) < lda {
                    ab[(j - 1) * lda + m - (j - 1 - k)]
                } else {
                    0.0
                };
                t -= ljk * ljk;
            }
        }

        if t <= 0.0 {
            return Err(j);
        }

        ab[(j - 1) * lda + m] = t.sqrt();

        // Forward elimination
        if j < n {
            let t_inv = 1.0 / t;
            for i in j..n {
                let mut sum = ab[i * lda + m - (i - j + 1)];
                for k in 0..j - 1 {
                    let ljk = if k + m + 1 >= j - 1 && k + m + 1 - (j - 1) < lda {
                        ab[(j - 1) * lda + m - (j - 1 - k)]
                    } else {
                        0.0
                    };
                    let lik = if k + m + 1 >= i - j + 1 && k + m + 1 - (i - j + 1) < lda {
                        ab[i * lda + m - (i - j) + (j - 1 - k)]
                    } else {
                        0.0
                    };
                    sum -= ljk * lik;
                }
                ab[i * lda + m - (i - j + 1)] = sum * t_inv;
            }
        }
    }

    if info == 0 {
        Ok(())
    } else {
        Err(info)
    }
}

/// Solve A * x = b using banded Cholesky factorization
/// Input: factored AB from dpbfa, right-hand side B
/// Output: solution X
pub fn dpbsl(
    ab: &[f64], n: usize, m: usize, lda: usize,
    b: &mut [f64],
) {
    // Solve L * y = b (forward substitution)
    for j in 1..n {
        let t = ab[(j - 1) * lda + m];
        if t.abs() > 1.0e-20 {
            let mut sum = 0.0;
            for k in 0..j - 1 {
                let ljk = if m >= j - k {
                    ab[(j - 1) * lda + m - (j - 1 - k)]
                } else {
                    0.0
                };
                sum += ljk * b[k];
            }
            b[j - 1] -= sum;
        }
    }

    // Solve L^T * x = y (back substitution)
    for j in (0..n).rev() {
        let t = ab[j * lda + m];
        if t.abs() > 1.0e-20 {
            let mut sum = b[j];
            for k in j + 1..n {
                let ljk = if k < n && m >= k - j {
                    ab[k * lda + m - (k - j)]
                } else {
                    0.0
                };
                sum -= ljk * b[k];
            }
            b[j] = sum / t;
        }
    }
}
