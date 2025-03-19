use core::f64;
use extendr_api::prelude::*;
use std::f64::consts::PI;

fn chebychev_polynomial(i: &[f64], y: f64, ya: f64, yb: f64) -> Vec<f64> {

    let a = 2.0 / PI;
    let b = f64::acos((y - ya) / yb);
    let c0 = (y - ya) / yb;
    let c = f64::sqrt(1.0 - c0.powf(2.0));

    let mut result = Vec::with_capacity(i.len());

    for iter in i {
        if (iter - 1.0) == 0.0 {
            result.push(a * f64::cos((iter - 1.0) * b) / (2.0 * c));
        } else {
            result.push(a * f64::cos((iter - 1.0) * b) / c);
        }
    }

    result
}


#[extendr]
fn est_aa_inner_rust(
    ii: usize,
    combn_m_input: Vec<i32>,
    combn_m_nrow: usize,
    g: Vec<f64>,
    precompute_ind: Vec<i32>,
    matrix_ind: Vec<i32>,
    g_dim1: usize,
    g_dim2: usize,
    g_dim3: usize,
) -> RArray<f64, [usize; 2]> {

    use na::DMatrix;
    use nalgebra as na;

    let ii_iter = ii.pow(3);

    // Original dimensions
    // let g_dim1 = 10;
    // let g_dim2 = 16;
    // let g_dim3 = 100000;

    let combn_m = na::DMatrix::from_fn(combn_m_nrow, 3, |row, col| {
        combn_m_input[row + col * combn_m_nrow] as usize
    });

    let precompute_ind = na::DMatrix::from_fn(ii_iter, 3, |row, col| {
        precompute_ind[row + col * ii_iter] as usize
    });

    let matrix_ind = na::DMatrix::from_fn(g_dim1, g_dim2, |row, col| {
        matrix_ind[row + col * g_dim1] as usize
    });

    let g = DMatrix::from_fn(g_dim3, g_dim1 * g_dim2, |i, j| g[j + i * g_dim1 * g_dim2]);
    // Transpose g to take advantage of column-major memory arrangement

    let mut v: Vec<f64> = Vec::with_capacity(combn_m_nrow * ii_iter);

    for j in 0..ii_iter {
        for q in 0..combn_m_nrow {
            let flattened_ind_1 = matrix_ind[(precompute_ind[(j, 0)], combn_m[(q, 0)])];
            let flattened_ind_2 = matrix_ind[(precompute_ind[(j, 1)], combn_m[(q, 1)])];
            let flattened_ind_3 = matrix_ind[(precompute_ind[(j, 2)], combn_m[(q, 2)])];

            let g_1 = g.column(flattened_ind_1);
            let g_2 = g.column(flattened_ind_2);
            let g_3 = g.column(flattened_ind_3);

            let g_sum = g_1.component_mul(&g_2);
            let g_sum2 = g_sum.component_mul(&g_3);

            v.push(g_sum2.sum());
        }
    }

    let v_return = RArray::new_matrix(combn_m_nrow, ii_iter, |row, col| v[row + col * combn_m_nrow]);

    v_return

}


#[extendr]
fn est_cdf_and_mixing_prop_inner_rust(
    y: Vec<f64>,
    II: usize,
    K: usize,
    P: usize,
    L: usize,
    M: usize,
    N_subset: usize,
    g_ind: Vec<i32>,
    precompute_ind: Vec<i32>,
    ind_3rd: Vec<i32>,
    ii: Vec<f64>,
    W: Vec<f64>,
    W_t: Vec<f64>,
    U: Vec<f64>,
    U_t: Vec<f64>,
    cdf_points: Vec<f64>,
    ya: f64,
    yb: f64,
) -> List {

    use na::DMatrix;
    use nalgebra as na;

    let y = na::DMatrix::from_fn(N_subset, M, |row, col| y[row + col * N_subset]);

    let g_ind_rows = g_ind.len() / 2;

    let g_ind = na::DMatrix::from_fn(g_ind_rows, 2, |row, col| {
        g_ind[row + col * g_ind_rows] as usize
    });

    let precompute_ind: Vec<usize> = precompute_ind.iter().map(|&e| e as usize).collect();
    let ind_3rd: Vec<usize> = ind_3rd.iter().map(|&e| e as usize).collect();

    let U = na::DMatrix::from_fn(K, K, |row, col| U[row + col * K]);
    let U_t = na::DMatrix::from_fn(K, K, |row, col| U_t[row + col * K]);

    let W = na::DMatrix::from_fn(K, II, |row, col| W[row + col * K]);
    let W_t = na::DMatrix::from_fn(II, K, |row, col| W_t[row + col * II]);

    let mut a_hat = na::DVector::zeros(II);

    let mut B_hat= Vec::with_capacity(K);
    
    for i in 0..K {
        B_hat.push(na::DMatrix::zeros(II, P));
    }

    let mut CDF = na::DMatrix::from_fn(K, L, |row, col| 0.0);

    for n in 0..N_subset {
        
        let g_precompute = DMatrix::from_vec(
            II,
            M,
            (0..M)
                .flat_map(|m| chebychev_polynomial(&ii, y[(n, m)], ya, yb))
                .collect(),
        );
        // Help from https://stackoverflow.com/questions/56739169/is-there-a-way-to-write-to-a-whole-row-column-of-a-nalgebra-matrix
        // Also needed flat_map() from here: https://stackoverflow.com/questions/49966420/how-do-i-collect-from-a-nested-iterator

        a_hat = a_hat + g_precompute.column_sum();

        let y_comparison = na::DMatrix::from_fn(1, ind_3rd.len(), |_row, col| y[(n, ind_3rd[col])]);
        let y_comparison = y_comparison.iter().enumerate();
        // Needed to break up statements because borrow checker complained

        for k in 0..K {
            let k_U = &U.column(k);
            let k_U_t = &U_t.row(k);

            let mut tau_precompute: Vec<f64> = Vec::with_capacity(g_ind_rows);

            for mm in 0..g_ind_rows {
                let first_outer = &g_precompute.column(g_ind[(mm, 0)]);
                let second_outer = &g_precompute.column(g_ind[(mm, 1)]).transpose();

                let g_center = first_outer * second_outer;
                // Outer product
                // https://github.com/dimforge/nalgebra/blob/main/CHANGELOG.md#removed-3
                // "The free function ::outer has been removed. Use column-vector × row-vector multiplication instead."

                // println!("k_U_t  {:#?}", (k_U_t.nrows(), k_U_t.ncols()));
                // println!("W  {:#?}", (W.nrows(), W.ncols()));
                // println!("g_center  {:#?}", (g_center.nrows(), g_center.ncols()));
                // println!("W_t  {:#?}", (W_t.nrows(), W_t.ncols()));
                // println!("k_U  {:#?}", (k_U.nrows(), k_U.ncols()));

                tau_precompute.push((k_U_t * &W * &g_center * &W_t * k_U)[(0, 0)]);
                // k_U_t: 1 x K
                // W: K x II
                // g_center: II x II
                // W_t: II x K
                // K_U K x 1
            }

            let tauP: Vec<f64> =
                na::DVector::from_fn(P, |iii, _| tau_precompute[precompute_ind[iii]])
                    .data
                    .into();

            let basis_3rd =
                na::DMatrix::from_fn(II, P, |row, col| g_precompute[(row, ind_3rd[col])]);

            let tauP_expanded =
                DMatrix::from_vec(P, II, (0..II).flat_map(|m| tauP.clone()).collect());

            B_hat[k] = B_hat[k].clone() + &tauP_expanded.transpose().component_mul(&basis_3rd);

            for l in 0..L {
                let tauP_filter = y_comparison
                    .clone()
                    .filter_map(|(index, &r)| (r <= cdf_points[l]).then_some(index))
                    .collect::<Vec<_>>();

                if tauP_filter.len() > 0 {
                    let tauP_sum =
                        na::DVector::from_fn(tauP_filter.len(), |row, _col| tauP[tauP_filter[row]]);

                    CDF[(k, l)] = CDF[(k, l)] + tauP_sum.sum();
                }
            }
        }
    }

    let a_hat: Vec<f64> = a_hat.data.into();
    let CDF = RArray::new_matrix(K, L, |row, col| CDF[(row, col)]);

    let B_hat_list = List::from_values(
        (0..K)
            .map(|k| RArray::new_matrix(II, P, |row, col| B_hat[k][(row, col)]))
            .collect::<Vec<_>>(),
    );

    let result = list!(a_hat = a_hat, B_hat_list = B_hat_list, CDF = CDF);

    result

}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod decoyanalysis;
    fn est_aa_inner_rust;
    fn est_cdf_and_mixing_prop_inner_rust;
}
