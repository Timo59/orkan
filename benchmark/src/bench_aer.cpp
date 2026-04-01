/**
 * @file bench_aer.cpp
 * @brief Qiskit Aer density-matrix backend for the gate benchmark.
 */

#ifdef HAVE_AER

#include "simulators/density_matrix/densitymatrix.hpp"
#include "framework/linalg/matrix_utils/vmatrix_defs.hpp"
#include <complex>
#include <cstdio>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "bench.h"

static_assert(sizeof(std::complex<double>) == sizeof(cplx_t),
    "std::complex<double> size must match cplx_t (double _Complex)");

/* =====================================================================
 * Context: Aer DensityMatrix + resolved apply function
 * ===================================================================== */

typedef void (*aer_apply_fn)(void *ctx, const qubit_t *pos);

struct aer_ctx_t {
    AER::QV::DensityMatrix<double> *dm;
    aer_apply_fn                    apply;
    AER::cvector_t                  vec;     /* pre-built vectorized matrix */
    double                          par;
};

/* =====================================================================
 * Per-gate apply functions — zero branching in hot path
 * ===================================================================== */

static void apply_x(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_x(static_cast<AER::uint_t>(pos[0]));
}
static void apply_y(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_y(static_cast<AER::uint_t>(pos[0]));
}
static void apply_z(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_phase(static_cast<AER::uint_t>(pos[0]), {-1.0, 0.0});
}
static void apply_h(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    AER::reg_t reg = {static_cast<AER::uint_t>(pos[0])};
    a->dm->apply_unitary_matrix(reg, a->vec);
}
static void apply_s(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_phase(static_cast<AER::uint_t>(pos[0]), {0.0, 1.0});
}
static void apply_t(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_phase(static_cast<AER::uint_t>(pos[0]),
                        {M_SQRT1_2, M_SQRT1_2});
}
static void apply_rx(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    AER::reg_t reg = {static_cast<AER::uint_t>(pos[0])};
    a->dm->apply_unitary_matrix(reg, a->vec);
}
static void apply_ry(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    AER::reg_t reg = {static_cast<AER::uint_t>(pos[0])};
    a->dm->apply_unitary_matrix(reg, a->vec);
}
static void apply_rz(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    AER::reg_t reg = {static_cast<AER::uint_t>(pos[0])};
    a->dm->apply_unitary_matrix(reg, a->vec);
}
static void apply_cx(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_cnot(static_cast<AER::uint_t>(pos[0]),
                       static_cast<AER::uint_t>(pos[1]));
}
static void apply_cy(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_cy(static_cast<AER::uint_t>(pos[0]),
                     static_cast<AER::uint_t>(pos[1]));
}
static void apply_cz(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_cphase(static_cast<AER::uint_t>(pos[0]),
                         static_cast<AER::uint_t>(pos[1]), {-1.0, 0.0});
}
static void apply_ccx(void *c, const qubit_t *pos) {
    auto *a = static_cast<aer_ctx_t*>(c);
    a->dm->apply_toffoli(static_cast<AER::uint_t>(pos[0]),
                          static_cast<AER::uint_t>(pos[1]),
                          static_cast<AER::uint_t>(pos[2]));
}

/* =====================================================================
 * Vtable implementation
 * ===================================================================== */

static void *aer_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    aer_ctx_t *ctx = nullptr;
    try {
        ctx = new aer_ctx_t{};  /* value-init: dm=nullptr, apply=nullptr */
        ctx->dm = new AER::QV::DensityMatrix<double>(qubits);
        ctx->dm->initialize();
#ifdef _OPENMP
        ctx->dm->set_omp_threads(omp_get_max_threads());
#endif
        ctx->par = par;

        /* Random Hermitian init */
        bench_init_hermitian_dense(
            reinterpret_cast<cplx_t*>(ctx->dm->data()),
            (idx_t)1 << qubits);

        /* Resolve apply function and pre-build vectorized matrices */
        using VMatrix = AER::Linalg::VMatrix;
        switch (gate) {
        case BG_X:   ctx->apply = apply_x;  break;
        case BG_Y:   ctx->apply = apply_y;  break;
        case BG_Z:   ctx->apply = apply_z;  break;
        case BG_H:   ctx->apply = apply_h;  ctx->vec = VMatrix::H; break;
        case BG_S:   ctx->apply = apply_s;  break;
        case BG_T:   ctx->apply = apply_t;  break;
        case BG_RX:  ctx->apply = apply_rx; ctx->vec = VMatrix::rx(par); break;
        case BG_RY:  ctx->apply = apply_ry; ctx->vec = VMatrix::ry(par); break;
        case BG_RZ:  ctx->apply = apply_rz; ctx->vec = VMatrix::rz(par); break;
        case BG_CX:  ctx->apply = apply_cx;  break;
        case BG_CY:  ctx->apply = apply_cy;  break;
        case BG_CZ:  ctx->apply = apply_cz;  break;
        case BG_CCX: ctx->apply = apply_ccx; break;
        default:
            delete ctx->dm;
            delete ctx;
            return nullptr;
        }
        return ctx;
    } catch (...) {
        if (ctx) {
            delete ctx->dm;  /* safe: dm is nullptr if not yet assigned */
            delete ctx;
        }
        return nullptr;
    }
}

static void aer_apply(void *c, const qubit_t *pos) {
    static_cast<aer_ctx_t*>(c)->apply(c, pos);
}

static void aer_cleanup(void *c) {
    auto *ctx = static_cast<aer_ctx_t*>(c);
    delete ctx->dm;
    delete ctx;
}

static size_t aer_touch_bytes(qubit_t qubits) {
    size_t dim = (size_t)1 << qubits;
    return 2 * dim * dim * sizeof(std::complex<double>);
}

static int aer_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

extern "C" {
const bench_backend_t aer_backend = {
    .name        = "Qiskit Aer",
    .init        = aer_init,
    .apply       = aer_apply,
    .cleanup     = aer_cleanup,
    .touch_bytes = aer_touch_bytes,
    .available   = aer_available
};
}

#endif /* HAVE_AER */
