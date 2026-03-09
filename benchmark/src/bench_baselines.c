/**
 * @file bench_baselines.c
 * @brief BLAS dense and naive loop baseline implementations for mixed state benchmarks
 *
 * Provides:
 * - Gate matrix constants (XMAT, ZMAT, HMAT, CXMAT, SWAPMAT)
 * - Kronecker product and full gate matrix construction
 * - apply_uruh_blas(), apply_uruh_naive()
 * - bench_blas_dense(), bench_naive_loop(), bench_blas_dense_2q()
 */

#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Gate matrices (2x2, column-major)
 * =====================================================================================================================
 */

#define INVSQRT2 0.7071067811865475

const cplx_t XMAT[4] = {0.0, 1.0, 1.0, 0.0};
const cplx_t ZMAT[4] = {1.0, 0.0, 0.0, -1.0};
const cplx_t HMAT[4] = {INVSQRT2, INVSQRT2, INVSQRT2, -INVSQRT2};

/*
 * Two-qubit gate matrices (4x4, column-major)
 *
 * Convention: for gate(q1, q2), matrix index = b_q2 * 2 + b_q1
 * where b_qi is the bit value of qubit qi in the basis state.
 */

/* CX(control, target): flips target when control=1 */
const cplx_t CXMAT[16] = {
    1, 0, 0, 0,   /* col 0: |c=0,t=0> -> |c=0,t=0> */
    0, 0, 0, 1,   /* col 1: |c=1,t=0> -> |c=1,t=1> */
    0, 0, 1, 0,   /* col 2: |c=0,t=1> -> |c=0,t=1> */
    0, 1, 0, 0    /* col 3: |c=1,t=1> -> |c=1,t=0> */
};

/* SWAP(q1, q2): exchanges qubit values */
const cplx_t SWAPMAT[16] = {
    1, 0, 0, 0,   /* col 0: |0,0> -> |0,0> */
    0, 0, 1, 0,   /* col 1: |1,0> -> |0,1> */
    0, 1, 0, 0,   /* col 2: |0,1> -> |1,0> */
    0, 0, 0, 1    /* col 3: |1,1> -> |1,1> */
};

/*
 * Build all ordered qubit pairs (q1 < q2) for an n-qubit system.
 * Declared here for use in bench_blas_dense_2q(); defined in bench_mixed.c.
 */
#define MAX_PAIRS 128  /* sufficient for up to ~16 qubits */
int build_all_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out);

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/** @brief Kronecker product: C = A ⊗ B */
static cplx_t* kron(dim_t ma, dim_t na, const cplx_t *A,
                    dim_t mb, dim_t nb, const cplx_t *B) {
    dim_t mc = ma * mb;
    dim_t nc = na * nb;
    cplx_t *C = malloc(mc * nc * sizeof(cplx_t));
    if (!C) return NULL;

    for (dim_t ja = 0; ja < na; ++ja) {
        for (dim_t ia = 0; ia < ma; ++ia) {
            cplx_t aij = A[ia + ja * ma];
            for (dim_t jb = 0; jb < nb; ++jb) {
                for (dim_t ib = 0; ib < mb; ++ib) {
                    dim_t ic = ib + ia * mb;
                    dim_t jc = jb + ja * nb;
                    C[ic + jc * mc] = aij * B[ib + jb * mb];
                }
            }
        }
    }
    return C;
}

/** @brief Build full gate matrix I⊗...⊗U⊗...⊗I for single-qubit gate on `target` */
static cplx_t* build_full_gate_matrix(qubit_t qubits, const cplx_t gate[4], qubit_t target) {
    qubit_t pos = qubits - 1 - target;
    dim_t dim_left = (dim_t)1 << pos;

    cplx_t *id_left = calloc(dim_left * dim_left, sizeof(cplx_t));
    if (!id_left) return NULL;
    for (dim_t i = 0; i < dim_left; ++i) id_left[i + i * dim_left] = 1.0;

    cplx_t *tmp = kron(dim_left, dim_left, id_left, 2, 2, gate);
    free(id_left);
    if (!tmp) return NULL;
    dim_left *= 2;

    dim_t dim_right = (dim_t)1 << target;
    cplx_t *id_right = calloc(dim_right * dim_right, sizeof(cplx_t));
    if (!id_right) { free(tmp); return NULL; }
    for (dim_t i = 0; i < dim_right; ++i) id_right[i + i * dim_right] = 1.0;

    cplx_t *U = kron(dim_left, dim_left, tmp, dim_right, dim_right, id_right);
    free(tmp);
    free(id_right);
    return U;
}

/** @brief Build full gate matrix for two-qubit gate on qubits q1, q2 */
static cplx_t* build_full_gate_matrix_2q(qubit_t qubits, const cplx_t gate[16],
                                          qubit_t q1, qubit_t q2) {
    dim_t dim = (dim_t)1 << qubits;
    cplx_t *U = calloc(dim * dim, sizeof(cplx_t));
    if (!U) return NULL;

    for (dim_t i = 0; i < dim; ++i) {
        int b1 = (i >> q1) & 1;
        int b2 = (i >> q2) & 1;
        int g_col = b2 * 2 + b1;

        for (int g_row = 0; g_row < 4; ++g_row) {
            cplx_t val = gate[g_row + g_col * 4];
            if (cabs(val) < 1e-15) continue;

            int out_b1 = g_row & 1;
            int out_b2 = (g_row >> 1) & 1;

            dim_t j = i;
            j = (j & ~((dim_t)1 << q1)) | ((dim_t)out_b1 << q1);
            j = (j & ~((dim_t)1 << q2)) | ((dim_t)out_b2 << q2);
            U[j + i * dim] = val;
        }
    }
    return U;
}

/*
 * Compute UρU† using BLAS zgemm.
 *
 * tmp: pre-allocated scratch buffer of dim*dim complex doubles.
 * Passing tmp from outside the hot loop eliminates per-call allocation overhead,
 * keeping the comparison focused on algorithmic cost.
 */
static void apply_uruh_blas(dim_t dim, const cplx_t *U, cplx_t *rho, cplx_t *tmp) {
    const cplx_t alpha = 1.0;
    const cplx_t beta  = 0.0;

    /* tmp = U * rho */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                dim, dim, dim, &alpha, U, dim, rho, dim, &beta, tmp, dim);

    /* rho = tmp * U† */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                dim, dim, dim, &alpha, tmp, dim, U, dim, &beta, rho, dim);
}

/*
 * Compute UρU† using naive triple-nested loops.
 *
 * tmp, out: pre-allocated scratch buffers of dim*dim complex doubles.
 * Same rationale as apply_uruh_blas: allocation lives outside the hot loop.
 */
static void apply_uruh_naive(dim_t dim, const cplx_t *U, cplx_t *rho,
                              cplx_t *tmp, cplx_t *out) {
    /* tmp = U * rho */
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = 0; i < dim; ++i) {
            cplx_t sum = 0.0;
            for (dim_t k = 0; k < dim; ++k)
                sum += U[i + k * dim] * rho[k + j * dim];
            tmp[i + j * dim] = sum;
        }
    }

    /* out = tmp * U† */
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = 0; i < dim; ++i) {
            cplx_t sum = 0.0;
            for (dim_t k = 0; k < dim; ++k)
                sum += tmp[i + k * dim] * conj(U[j + k * dim]);
            out[i + j * dim] = sum;
        }
    }

    memcpy(rho, out, dim * dim * sizeof(cplx_t));
}

/*
 * =====================================================================================================================
 * Benchmark implementations
 * =====================================================================================================================
 */

bench_result_t bench_blas_dense(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "blas_dense";
    result.iterations  = iterations;

    dim_t dim = result.dim;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    cplx_t **U = malloc(qubits * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (qubit_t t = 0; t < qubits; ++t)
        U[t] = build_full_gate_matrix(qubits, gate_mat, t);

    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) {
        for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
        free(U); free(rho);
        return result;
    }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_blas(dim, U[i % qubits], rho, tmp);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) {
        free(tmp); for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
        free(U); free(rho); return result;
    }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            for (qubit_t t = 0; t < qubits; ++t)
                apply_uruh_blas(dim, U[t], rho, tmp);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;
    result.runs           = runs;
    result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * qubits) / (stats.mean / 1000.0) : 0.0;

    free(tmp);
    free(rho);
    for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
    free(U);
    return result;
}

bench_result_t bench_naive_loop(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "naive_loop";
    result.iterations  = iterations;

    dim_t dim = result.dim;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    cplx_t **U = malloc(qubits * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (qubit_t t = 0; t < qubits; ++t)
        U[t] = build_full_gate_matrix(qubits, gate_mat, t);

    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    cplx_t *out = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp || !out) {
        free(tmp); free(out);
        for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
        free(U); free(rho);
        return result;
    }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_naive(dim, U[i % qubits], rho, tmp, out);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) {
        free(tmp); free(out); for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
        free(U); free(rho); return result;
    }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            for (qubit_t t = 0; t < qubits; ++t)
                apply_uruh_naive(dim, U[t], rho, tmp, out);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;
    result.runs           = runs;
    result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * qubits) / (stats.mean / 1000.0) : 0.0;

    free(tmp); free(out);
    free(rho);
    for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
    free(U);
    return result;
}

bench_result_t bench_blas_dense_2q(qubit_t qubits, const char *gate_name,
                                     const cplx_t gate_mat[16],
                                     int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "blas_dense";
    result.iterations  = iterations;

    dim_t dim = result.dim;

    qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
    int num_pairs = build_all_pairs(qubits, q1s, q2s);
    if (num_pairs == 0) return result;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    cplx_t **U = malloc(num_pairs * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (int p = 0; p < num_pairs; ++p)
        U[p] = build_full_gate_matrix_2q(qubits, gate_mat, q1s[p], q2s[p]);

    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) {
        for (int p = 0; p < num_pairs; ++p) free(U[p]);
        free(U); free(rho);
        return result;
    }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_blas(dim, U[i % num_pairs], rho, tmp);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) {
        free(tmp); for (int p = 0; p < num_pairs; ++p) free(U[p]);
        free(U); free(rho); return result;
    }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            for (int p = 0; p < num_pairs; ++p)
                apply_uruh_blas(dim, U[p], rho, tmp);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;
    result.runs           = runs;
    result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * num_pairs) / (stats.mean / 1000.0) : 0.0;

    free(tmp);
    free(rho);
    for (int p = 0; p < num_pairs; ++p) free(U[p]);
    free(U);
    return result;
}
