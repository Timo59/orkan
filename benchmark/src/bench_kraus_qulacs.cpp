/**
 * @file bench_kraus_qulacs.cpp
 * @brief Qulacs backend for the Kraus benchmark.
 *
 * Qulacs has no native Kraus channel API at the csim level.
 * Each apply call manually iterates over K Kraus operators:
 * copy rho to scratch, apply dm_single_qubit_dense_matrix_gate
 * (which performs rho -> A rho A^dagger in-place), accumulate
 * via dm_state_add.
 *
 * Each individual conjugation and accumulation is OpenMP-parallelized,
 * but there is no operator-level fusion: K operators means K separate
 * parallel regions for conjugation plus K for accumulation.
 *
 * Capped at 10 qubits.
 */

#ifdef HAVE_QULACS

#include <csim/update_ops_dm.hpp>
#include <csim/stat_ops_dm.hpp>
#include <csim/type.hpp>
#include <cstdlib>
#include <cstring>

#include "bench_kraus.h"

static_assert(sizeof(CTYPE) == sizeof(cplx_t),
    "Qulacs CTYPE size must match cplx_t (double _Complex)");

#define KRAUS_QULACS_MAX_QUBITS 10

/* =====================================================================
 * Context
 * ===================================================================== */

struct kraus_qulacs_ctx_t {
    CTYPE  *rho;        /* density matrix (dim x dim, row-major) */
    CTYPE  *scratch;    /* scratch copy of rho for each operator */
    CTYPE  *acc;        /* accumulator for Kraus sum */
    CTYPE  *ops;        /* n_ops small matrices, each d x d, row-major */
    ITYPE   dim;
    int     n_ops;
    int     tgt_qubits;
};

/* =====================================================================
 * Apply: rho' = sum_k A_k rho A_k^dagger
 * ===================================================================== */

static void kraus_qulacs_apply(void *c, const qubit_t *pos) {
    auto *ctx = static_cast<kraus_qulacs_ctx_t *>(c);
    size_t mat_bytes = static_cast<size_t>(ctx->dim) *
                       static_cast<size_t>(ctx->dim) * sizeof(CTYPE);
    ITYPE d = static_cast<ITYPE>(1) << ctx->tgt_qubits;

    std::memset(ctx->acc, 0, mat_bytes);

    for (int k = 0; k < ctx->n_ops; ++k) {
        /* Copy original rho into scratch */
        std::memcpy(ctx->scratch, ctx->rho, mat_bytes);

        /* Apply A_k rho A_k^dagger in-place on scratch */
        CTYPE *mat = &ctx->ops[k * d * d];
        if (ctx->tgt_qubits == 1) {
            dm_single_qubit_dense_matrix_gate(
                static_cast<UINT>(pos[0]), mat, ctx->scratch, ctx->dim);
        } else {
            UINT tgts[3];
            for (int i = 0; i < ctx->tgt_qubits; ++i)
                tgts[i] = static_cast<UINT>(pos[i]);
            dm_multi_qubit_dense_matrix_gate(
                tgts, static_cast<UINT>(ctx->tgt_qubits),
                mat, ctx->scratch, ctx->dim);
        }

        /* Accumulate: acc += scratch (OpenMP-parallelized) */
        dm_state_add(ctx->scratch, ctx->acc, ctx->dim);
    }

    /* Copy result back */
    std::memcpy(ctx->rho, ctx->acc, mat_bytes);
}

/* =====================================================================
 * Init / cleanup / vtable
 * ===================================================================== */

static void *kraus_qulacs_init(qubit_t sys_qubits, int n_ops,
                               int tgt_qubits, const cplx_t *ops) {
    if (sys_qubits > KRAUS_QULACS_MAX_QUBITS) return nullptr;

    ITYPE dim = static_cast<ITYPE>(1) << sys_qubits;
    ITYPE d = static_cast<ITYPE>(1) << tgt_qubits;
    size_t mat_bytes = static_cast<size_t>(dim) *
                       static_cast<size_t>(dim) * sizeof(CTYPE);

    auto *ctx = static_cast<kraus_qulacs_ctx_t *>(
        std::calloc(1, sizeof(kraus_qulacs_ctx_t)));
    if (!ctx) return nullptr;
    ctx->dim = dim;
    ctx->n_ops = n_ops;
    ctx->tgt_qubits = tgt_qubits;

    ctx->rho     = static_cast<CTYPE *>(std::malloc(mat_bytes));
    ctx->scratch = static_cast<CTYPE *>(std::malloc(mat_bytes));
    ctx->acc     = static_cast<CTYPE *>(std::malloc(mat_bytes));
    if (!ctx->rho || !ctx->scratch || !ctx->acc) goto fail;

    bench_init_hermitian_dense(
        reinterpret_cast<cplx_t *>(ctx->rho), static_cast<dim_t>(dim));

    /* Copy Kraus operators (row-major, matching Qulacs convention) */
    ctx->ops = static_cast<CTYPE *>(
        std::malloc(static_cast<size_t>(n_ops) * d * d * sizeof(CTYPE)));
    if (!ctx->ops) goto fail;
    std::memcpy(ctx->ops, ops,
                static_cast<size_t>(n_ops) * d * d * sizeof(CTYPE));

    return ctx;

fail:
    std::free(ctx->ops);
    std::free(ctx->acc);
    std::free(ctx->scratch);
    std::free(ctx->rho);
    std::free(ctx);
    return nullptr;
}

static void kraus_qulacs_cleanup(void *c) {
    auto *ctx = static_cast<kraus_qulacs_ctx_t *>(c);
    std::free(ctx->ops);
    std::free(ctx->acc);
    std::free(ctx->scratch);
    std::free(ctx->rho);
    std::free(ctx);
}

static size_t kraus_qulacs_touch_bytes(qubit_t qubits) {
    size_t dim = static_cast<size_t>(1) << qubits;
    return 2 * dim * dim * sizeof(CTYPE);
}

static int kraus_qulacs_available(qubit_t qubits) {
    return qubits <= KRAUS_QULACS_MAX_QUBITS;
}

extern "C" {
const bench_kraus_backend_t kraus_qulacs_backend = {
    .name        = "Qulacs",
    .init        = kraus_qulacs_init,
    .apply       = kraus_qulacs_apply,
    .cleanup     = kraus_qulacs_cleanup,
    .touch_bytes = kraus_qulacs_touch_bytes,
    .available   = kraus_qulacs_available
};
}

#endif /* HAVE_QULACS */
