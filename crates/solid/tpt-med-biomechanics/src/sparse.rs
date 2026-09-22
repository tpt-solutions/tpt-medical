//! CSR sparse matrix and conjugate-gradient solver.

/// Compressed sparse row symmetric-pattern matrix (only the upper triangle
/// of the action is needed by CG, but full rows are stored for clarity).
#[derive(Debug, Clone)]
pub struct CsrMatrix {
    /// Number of rows (== number of columns; the system matrix is square).
    pub nrows: usize,
    /// Row offsets, length `nrows + 1`.
    pub indptr: Vec<usize>,
    /// Column indices, length `nnz`.
    pub indices: Vec<usize>,
    /// Values, length `nnz`.
    pub data: Vec<f64>,
}

impl CsrMatrix {
    /// Builds a CSR matrix from (row, col, value) triplets. Duplicate
    /// entries are summed (standard FEM assembly semantics).
    pub fn from_triplets(n: usize, mut triplets: Vec<(usize, usize, f64)>) -> Self {
        triplets.sort_by_key(|&(r, c, _)| (r, c));
        let mut indptr = vec![0usize; n + 1];
        let mut indices = Vec::with_capacity(triplets.len());
        let mut data = Vec::with_capacity(triplets.len());
        let mut last = (usize::MAX, usize::MAX);
        for (r, c, v) in triplets {
            if (r, c) == last {
                *data.last_mut().expect("non-empty") += v;
                continue;
            }
            last = (r, c);
            indices.push(c);
            data.push(v);
            indptr[r + 1] += 1;
        }
        for i in 0..n {
            indptr[i + 1] += indptr[i];
        }
        Self {
            nrows: n,
            indptr,
            indices,
            data,
        }
    }

    /// Matrix–vector product `y = A x`.
    pub fn mul_vec(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..self.nrows {
            let mut sum = 0.0;
            for k in self.indptr[i]..self.indptr[i + 1] {
                sum += self.data[k] * x[self.indices[k]];
            }
            y[i] = sum;
        }
    }

    /// Diagonal (for Jacobi preconditioning). Missing diagonals are 0.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.nrows];
        for i in 0..self.nrows {
            for k in self.indptr[i]..self.indptr[i + 1] {
                if self.indices[k] == i {
                    d[i] += self.data[k];
                }
            }
        }
        d
    }
}

/// Solver outcome statistics.
#[derive(Debug, Clone, Copy)]
pub struct SolveStats {
    /// Iterations performed.
    pub iterations: usize,
    /// Final relative residual norm.
    pub relative_residual: f64,
    /// True if the tolerance was reached within `max_iterations`.
    pub converged: bool,
}

/// Preconditioned conjugate gradient for symmetric positive definite
/// systems, with Jacobi preconditioning.
pub fn conjugate_gradient(
    a: &CsrMatrix,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> (Vec<f64>, SolveStats) {
    let n = a.nrows;
    let mut x = vec![0.0; n];
    let diag = a.diagonal();
    let mut r = b.to_vec();
    let mut z = vec![0.0; n];
    for i in 0..n {
        z[i] = if diag[i].abs() > 1e-30 {
            r[i] / diag[i]
        } else {
            r[i]
        };
    }
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let b_norm = dot(b, b).max(1e-300);

    let mut stats = SolveStats {
        iterations: 0,
        relative_residual: f64::INFINITY,
        converged: false,
    };

    let mut ap = vec![0.0; n];
    for iter in 0..=max_iterations {
        let r_norm = dot(&r, &r);
        stats.iterations = iter;
        stats.relative_residual = (r_norm / b_norm).sqrt();
        if stats.relative_residual.sqrt() < tolerance {
            stats.converged = true;
            break;
        }
        if iter == max_iterations {
            break;
        }
        a.mul_vec(&p, &mut ap);
        let pap = dot(&p, &ap);
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
            z[i] = if diag[i].abs() > 1e-30 {
                r[i] / diag[i]
            } else {
                r[i]
            };
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    (x, stats)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csr_assembles_duplicates_and_multiplies() {
        // [[4,1],[1,3]]
        let a = CsrMatrix::from_triplets(
            2,
            vec![
                (0, 0, 4.0),
                (0, 0, 0.0),
                (0, 1, 1.0),
                (1, 0, 1.0),
                (1, 1, 3.0),
            ],
        );
        let y = {
            let mut y = vec![0.0; 2];
            a.mul_vec(&[1.0, 2.0], &mut y);
            y
        };
        assert_eq!(y, vec![6.0, 7.0]);
    }

    #[test]
    fn cg_solves_spd_system() {
        // K for a 1D chain of springs (stiffness 1): tridiag(-1, 2, -1) with
        // the last diagonal pinned to 1 to remove the rigid mode.
        let n = 50;
        let mut t = Vec::new();
        for i in 0..n {
            if i > 0 {
                t.push((i, i - 1, -1.0));
                t.push((i, i, 2.0));
            } else {
                t.push((i, i, 2.0));
            }
            if i + 1 < n {
                t.push((i, i + 1, -1.0));
            }
        }
        t.push((n - 1, n - 1, -1.0)); // last diagonal: 2 -> 1 (duplicates sum)
        let a = CsrMatrix::from_triplets(n, t);
        let b = vec![1.0; n];
        let (x, stats) = conjugate_gradient(&a, &b, 1e-10, 2 * n);
        assert!(stats.converged, "{stats:?}");
        // Residual check.
        let mut r = vec![0.0; n];
        a.mul_vec(&x, &mut r);
        for i in 0..n {
            r[i] -= b[i];
        }
        let res = dot(&r, &r).sqrt();
        assert!(res < 1e-8, "residual {res}");
    }
}
