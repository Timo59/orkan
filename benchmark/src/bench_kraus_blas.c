/**
 * @file bench_kraus_blas.c
 * @brief BLAS Kraus baseline: rho' = sum_k A_k rho A_k^dagger via zgemm.
 *
 * For each position, pre-computes full N x N Kronecker-expanded matrices
 * for all Kraus operators.  Each apply call performs 2K zgemm calls plus
 * accumulation.  Capped at 10 qubits.
 */

#include "bench_kraus.h"
#include <math.h>
#include <complex.h>
#include <string.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

#define KRAUS_BLAS_MAX_QUBITS 10

/* =====================================================================
 * Kronecker product: C = A (x) B  (column-major)
 * ===================================================================== */

static cplx_t *kron(dim_t ma, dim_t na, const cplx_t *A,
                     dim_t mb, dim_t nb, const cplx_t *B) {
    dim_t mc = ma*mb, nc = na*nb;
    cplx_t *C = malloc((size_t)(mc*nc) * sizeof(cplx_t));
    if (!C) return NULL;
    for (dim_t ja = 0; ja < na; ++ja)
        for (dim_t ia = 0; ia < ma; ++ia) {
            cplx_t a = A[ia + ja*ma];
            for (dim_t jb = 0; jb < nb; ++jb)
                for (dim_t ib = 0; ib < mb; ++ib)
                    C[(ib+ia*mb) + (jb+ja*nb)*mc] = a * B[ib+jb*mb];
        }
    return C;
}

/* =====================================================================
 * Build full N x N matrix from a d x d small matrix acting on target qubits
 *
 * For 1-qubit operators: I_{left} (x) A (x) I_{right}
 * ===================================================================== */

static cplx_t *build_full_1q(qubit_t n_qubits, const cplx_t g[4], qubit_t tgt) {
    qubit_t pos = n_qubits - 1 - tgt;
    dim_t dl = (dim_t)1 << pos;
    cplx_t *idl = calloc((size_t)(dl*dl), sizeof(cplx_t));
    if (!idl) return NULL;
    for (dim_t i = 0; i < dl; ++i) idl[i+i*dl] = 1.0;
    cplx_t *tmp = kron(dl, dl, idl, 2, 2, g);
    free(idl);
    if (!tmp) return NULL;
    dl *= 2;
    dim_t dr = (dim_t)1 << tgt;
    cplx_t *idr = calloc((size_t)(dr*dr), sizeof(cplx_t));
    if (!idr) { free(tmp); return NULL; }
    for (dim_t i = 0; i < dr; ++i) idr[i+i*dr] = 1.0;
    cplx_t *U = kron(dl, dl, tmp, dr, dr, idr);
    free(tmp); free(idr);
    return U;
}

/* =====================================================================
 * Context
 * ===================================================================== */

typedef struct {
    cplx_t  *rho;       /* density matrix (dim x dim, column-major) */
    cplx_t  *tmp;       /* scratch for zgemm */
    cplx_t  *acc;       /* accumulator for Kraus sum */
    cplx_t **U_lookup;  /* U_lookup[k * n + tgt] = full expanded Kraus op k at target tgt */
    int      n_ops;
    int      n_qubits;
    dim_t    dim;
} kraus_blas_ctx_t;

/* =====================================================================
 * Apply: rho' = sum_k A_k rho A_k^dagger
 * ===================================================================== */

static void kraus_blas_apply(void *c, const qubit_t *pos) {
    kraus_blas_ctx_t *ctx = (kraus_blas_ctx_t *)c;
    dim_t dim = ctx->dim;
    size_t mat_bytes = (size_t)(dim * dim) * sizeof(cplx_t);
    const cplx_t one = 1.0, zero = 0.0;

    memset(ctx->acc, 0, mat_bytes);

    for (int k = 0; k < ctx->n_ops; ++k) {
        const cplx_t *Ak = ctx->U_lookup[k * ctx->n_qubits + (int)pos[0]];

        /* tmp = A_k * rho */
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    dim, dim, dim, &one, Ak, dim, ctx->rho, dim, &zero, ctx->tmp, dim);
        /* acc += tmp * A_k^dagger */
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                    dim, dim, dim, &one, ctx->tmp, dim, Ak, dim, &one, ctx->acc, dim);
    }

    memcpy(ctx->rho, ctx->acc, mat_bytes);
}

/* =====================================================================
 * Init / cleanup / vtable
 * ===================================================================== */

static void *kraus_blas_init(qubit_t sys_qubits, int n_ops,
                             int tgt_qubits, const cplx_t *ops) {
    if (sys_qubits > KRAUS_BLAS_MAX_QUBITS) return NULL;
    if (tgt_qubits != 1) return NULL;  /* only 1-qubit channels for now */

    int n = (int)sys_qubits;
    dim_t dim = (dim_t)1 << sys_qubits;
    dim_t d = (dim_t)1 << tgt_qubits;

    kraus_blas_ctx_t *ctx = calloc(1, sizeof(kraus_blas_ctx_t));
    if (!ctx) return NULL;
    ctx->dim = dim;
    ctx->n_qubits = n;
    ctx->n_ops = n_ops;

    ctx->rho = malloc((size_t)(dim*dim) * sizeof(cplx_t));
    ctx->tmp = malloc((size_t)(dim*dim) * sizeof(cplx_t));
    ctx->acc = malloc((size_t)(dim*dim) * sizeof(cplx_t));
    if (!ctx->rho || !ctx->tmp || !ctx->acc) goto fail;
    bench_init_hermitian_dense(ctx->rho, dim);

    /* Build Kronecker-expanded matrices: n_ops * n target positions */
    ctx->U_lookup = calloc((size_t)(n_ops * n), sizeof(cplx_t *));
    if (!ctx->U_lookup) goto fail;

    for (int k = 0; k < n_ops; ++k) {
        /* Convert k-th small matrix from row-major to column-major */
        cplx_t g[4];
        g[0] = ops[k*d*d + 0];  g[1] = ops[k*d*d + 2];
        g[2] = ops[k*d*d + 1];  g[3] = ops[k*d*d + 3];

        for (int t = 0; t < n; ++t) {
            ctx->U_lookup[k * n + t] = build_full_1q(sys_qubits, g, (qubit_t)t);
            if (!ctx->U_lookup[k * n + t]) goto fail;
        }
    }
    return ctx;

fail:
    if (ctx->U_lookup) {
        for (int i = 0; i < n_ops * n; ++i) free(ctx->U_lookup[i]);
        free(ctx->U_lookup);
    }
    free(ctx->acc); free(ctx->tmp); free(ctx->rho); free(ctx);
    return NULL;
}

static void kraus_blas_cleanup(void *c) {
    kraus_blas_ctx_t *ctx = (kraus_blas_ctx_t *)c;
    if (ctx->U_lookup) {
        for (int i = 0; i < ctx->n_ops * ctx->n_qubits; ++i)
            free(ctx->U_lookup[i]);
        free(ctx->U_lookup);
    }
    free(ctx->acc); free(ctx->tmp); free(ctx->rho); free(ctx);
}

static size_t kraus_blas_touch_bytes(qubit_t qubits) {
    size_t dim = (size_t)1 << qubits;
    /* Per apply: read rho once, read each A_k once, write acc once,
       plus intermediate reads/writes in zgemm — approximate as 2 * dim^2 */
    return 2 * dim * dim * sizeof(cplx_t);
}

static int kraus_blas_available(qubit_t qubits) {
    return qubits <= KRAUS_BLAS_MAX_QUBITS;
}

const bench_kraus_backend_t kraus_blas_backend = {
    .name        = "BLAS",
    .init        = kraus_blas_init,
    .apply       = kraus_blas_apply,
    .cleanup     = kraus_blas_cleanup,
    .touch_bytes = kraus_blas_touch_bytes,
    .available   = kraus_blas_available
};
