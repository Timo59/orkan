/**
 * @file bench_kraus_aer.cpp
 * @brief Qiskit Aer backend for the Kraus benchmark.
 *
 * Converts Kraus operators to a superoperator via kraus_superop()
 * and vectorize_matrix() at init time.  The timed apply call uses
 * apply_superop_matrix(), matching Aer's production code path.
 */

#ifdef HAVE_AER

#include "simulators/density_matrix/densitymatrix.hpp"
#include "framework/utils.hpp"
#include <complex>
#include <cstdlib>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "bench_kraus.h"

static_assert(sizeof(std::complex<double>) == sizeof(cplx_t),
    "std::complex<double> size must match cplx_t (double _Complex)");

/* =====================================================================
 * Context
 * ===================================================================== */

struct kraus_aer_ctx_t {
    AER::QV::DensityMatrix<double> *dm;
    AER::cvector_t                  superop_vec;
    AER::reg_t                      qubits;
};

/* =====================================================================
 * Apply: superoperator on target qubits
 * ===================================================================== */

static void kraus_aer_apply(void *c, const qubit_t *pos) {
    auto *ctx = static_cast<kraus_aer_ctx_t *>(c);
    /* Update target qubit indices */
    for (size_t i = 0; i < ctx->qubits.size(); ++i)
        ctx->qubits[i] = static_cast<AER::uint_t>(pos[i]);
    ctx->dm->apply_superop_matrix(ctx->qubits, ctx->superop_vec);
}

/* =====================================================================
 * Init / cleanup / vtable
 * ===================================================================== */

static void *kraus_aer_init(qubit_t sys_qubits, int n_ops,
                            int tgt_qubits, const cplx_t *ops) {
    kraus_aer_ctx_t *ctx = nullptr;
    try {
        ctx = new kraus_aer_ctx_t{};

        idx_t d = static_cast<idx_t>(1) << tgt_qubits;

        /* Build Aer matrix objects from our flat row-major operators */
        std::vector<AER::cmatrix_t> kmats(n_ops);
        for (int k = 0; k < n_ops; ++k) {
            kmats[k] = AER::cmatrix_t(d, d);
            for (idx_t i = 0; i < d; ++i)
                for (idx_t j = 0; j < d; ++j) {
                    const auto *src = reinterpret_cast<const std::complex<double> *>(
                        &ops[k * d * d + i * d + j]);
                    kmats[k](i, j) = *src;
                }
        }

        /* Build superoperator and vectorize */
        AER::cmatrix_t superop = AER::Utils::kraus_superop(kmats);
        ctx->superop_vec = AER::Utils::vectorize_matrix(superop);

        /* Target qubit register */
        ctx->qubits.resize(tgt_qubits);
        for (int i = 0; i < tgt_qubits; ++i)
            ctx->qubits[i] = 0;  /* placeholder, updated in apply */

        /* Allocate and init density matrix */
        ctx->dm = new AER::QV::DensityMatrix<double>(sys_qubits);
        ctx->dm->initialize();
#ifdef _OPENMP
        ctx->dm->set_omp_threads(omp_get_max_threads());
#endif
        bench_init_hermitian_dense(
            reinterpret_cast<cplx_t *>(ctx->dm->data()),
            static_cast<idx_t>(1) << sys_qubits);

        return ctx;
    } catch (...) {
        if (ctx) {
            delete ctx->dm;
            delete ctx;
        }
        return nullptr;
    }
}

static void kraus_aer_cleanup(void *c) {
    auto *ctx = static_cast<kraus_aer_ctx_t *>(c);
    delete ctx->dm;
    delete ctx;
}

static size_t kraus_aer_touch_bytes(qubit_t qubits) {
    size_t dim = static_cast<size_t>(1) << qubits;
    return 2 * dim * dim * sizeof(std::complex<double>);
}

static int kraus_aer_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

extern "C" {
const bench_kraus_backend_t kraus_aer_backend = {
    .name        = "Qiskit Aer",
    .init        = kraus_aer_init,
    .apply       = kraus_aer_apply,
    .cleanup     = kraus_aer_cleanup,
    .touch_bytes = kraus_aer_touch_bytes,
    .available   = kraus_aer_available
};
}

#endif /* HAVE_AER */
