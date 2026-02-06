// qhipster.c - Quantum gates for pure states (Hilbert space vectors)

#include "gate.h"

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 * Pure states: O(dim) work per gate, need dim >= 4096 (~12 qubits) to offset thread overhead
 */
#define OMP_THRESHOLD 4096

/*
 * =====================================================================================================================
 * Single-qubit gate traversal macro
 * =====================================================================================================================
 *
 * Iterates over all pairs of amplitudes differing only in the target qubit:
 *   - a points to amplitude with target bit = 0
 *   - b points to amplitude with target bit = 1
 *
 * The PAIR_OP(a, b) macro performs the gate-specific operation on each pair.
 *
 * OpenMP parallelization uses collapse(2) to flatten the loop nest, ensuring good
 * parallelism regardless of which qubit is targeted:
 *   - Right qubits (low index): many blocks, small stride
 *   - Left qubits (high index): few blocks, large stride
 * With collapse(2), total iterations = dim/2, evenly distributed across threads.
 */
#define TRAVERSE_PURE_1Q(state, target, PAIR_OP)                               \
do {                                                                           \
    const dim_t dim = POW2((state)->qubits, dim_t);                            \
    const dim_t stride = POW2(target, dim_t);                                  \
    const dim_t step = POW2((target) + 1, dim_t);                              \
    cplx_t *data = (state)->data;                                              \
                                                                               \
    _Pragma("omp parallel for collapse(2) if(dim >= OMP_THRESHOLD)")           \
    for (dim_t i = 0; i < dim; i += step) {                                    \
        for (dim_t j = 0; j < stride; ++j) {                                   \
            cplx_t *a = data + i + j;                                          \
            cplx_t *b = data + i + j + stride;                                 \
            PAIR_OP(a, b);                                                     \
        }                                                                      \
    }                                                                          \
} while(0)

/*
 * =====================================================================================================================
 * Gate operation macros - optimized to minimize floating-point operations
 * =====================================================================================================================
 *
 * Multiplication by special complex values (i, -i, e^{iπ/4}) is implemented using
 * only real arithmetic (swaps, negations, additions) where possible.
 */

// X gate: swap a and b
#define X_OP(a, b) do {         \
    cplx_t tmp = *(a);          \
    *(a) = *(b);                \
    *(b) = tmp;                 \
} while(0)

// Y gate: a' = -i*b, b' = i*a  (zero multiplications)
// (x + yi)*i = -y + xi,  (x + yi)*(-i) = y - xi
#define Y_OP(a, b) do {                             \
    double a_re = creal(*(a)), a_im = cimag(*(a));  \
    double b_re = creal(*(b)), b_im = cimag(*(b));  \
    *(a) = CMPLX(b_im, -b_re);                      \
    *(b) = CMPLX(-a_im, a_re);                      \
} while(0)

// Z gate: b' = -b
#define Z_OP(a, b) do { (void)(a); *(b) = -*(b); } while(0)

// H gate: a' = (a+b)/sqrt(2), b' = (a-b)/sqrt(2)
#define H_OP(a, b) do {                             \
    cplx_t sum = *(a) + *(b);                       \
    cplx_t diff = *(a) - *(b);                      \
    *(a) = sum * M_SQRT1_2;                         \
    *(b) = diff * M_SQRT1_2;                        \
} while(0)

// S gate: b' = i*b  (zero multiplications)
// (x + yi)*i = -y + xi
#define S_OP(a, b) do {                             \
    (void)(a);                                      \
    *(b) = CMPLX(-cimag(*(b)), creal(*(b)));        \
} while(0)

// S-dagger gate: b' = -i*b  (zero multiplications)
// (x + yi)*(-i) = y - xi
#define SDG_OP(a, b) do {                           \
    (void)(a);                                      \
    *(b) = CMPLX(cimag(*(b)), -creal(*(b)));        \
} while(0)

// T gate: b' = e^{iπ/4}*b = (1+i)/sqrt(2) * b
// (x + yi)(1+i)/sqrt(2) = [(x-y) + (x+y)i]/sqrt(2)
#define T_OP(a, b) do {                             \
    (void)(a);                                      \
    double re = creal(*(b)), im = cimag(*(b));      \
    *(b) = CMPLX((re - im) * M_SQRT1_2,             \
                 (re + im) * M_SQRT1_2);            \
} while(0)

// T-dagger gate: b' = e^{-iπ/4}*b = (1-i)/sqrt(2) * b
// (x + yi)(1-i)/sqrt(2) = [(x+y) + (y-x)i]/sqrt(2)
#define TDG_OP(a, b) do {                           \
    (void)(a);                                      \
    double re = creal(*(b)), im = cimag(*(b));      \
    *(b) = CMPLX((re + im) * M_SQRT1_2,             \
                 (im - re) * M_SQRT1_2);            \
} while(0)

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, X_OP);
}

void y_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, Y_OP);
}

void z_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, Z_OP);
}

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, H_OP);
}

void s_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, S_OP);
}

void sdg_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, SDG_OP);
}

void t_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, T_OP);
}

void tdg_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, TDG_OP);
}

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

/**
 * @brief Insert a 0 bit at two positions in value `val`
 *
 * For CNOT, we iterate over indices where both control and target bits are 0.
 * This helper inserts 0 bits at positions lo and hi (lo < hi).
 */
static inline dim_t insertBits2_0(dim_t val, qubit_t lo, qubit_t hi) {
    // Insert 0 at position lo first
    dim_t mask_lo = ((dim_t)1 << lo) - 1;
    val = (val & mask_lo) | ((val & ~mask_lo) << 1);
    // Insert 0 at position hi (already shifted by 1 due to lo insertion)
    dim_t mask_hi = ((dim_t)1 << hi) - 1;
    return (val & mask_hi) | ((val & ~mask_hi) << 1);
}

/*
 * CNOT (Controlled-X) gate: Flips target qubit when control qubit is |1>
 *
 * For each basis state index with control=1, we swap amplitudes where target=0 and target=1.
 * This is a "zero-flop" gate requiring only conditional memory swaps.
 */
void cx_pure(state_t *state, const qubit_t control, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);
    cplx_t *data = state->data;

    // Determine which qubit index is lower/higher for bit insertion
    qubit_t lo = (control < target) ? control : target;
    qubit_t hi = (control < target) ? target : control;

    const dim_t incr_ctrl = POW2(control, dim_t);
    const dim_t incr_tgt = POW2(target, dim_t);
    const dim_t n_base = dim >> 2;  // dim / 4: we iterate over indices with both bits = 0

    #pragma omp parallel for if(dim >= OMP_THRESHOLD)
    for (dim_t k = 0; k < n_base; ++k) {
        // Insert 0 bits at positions lo and hi to get base index
        dim_t base = insertBits2_0(k, lo, hi);

        // Only swap when control=1: swap |ctrl=1,tgt=0> <-> |ctrl=1,tgt=1>
        dim_t idx_c1_t0 = base | incr_ctrl;           // control=1, target=0
        dim_t idx_c1_t1 = base | incr_ctrl | incr_tgt; // control=1, target=1

        cplx_t tmp = data[idx_c1_t0];
        data[idx_c1_t0] = data[idx_c1_t1];
        data[idx_c1_t1] = tmp;
    }
}

/*
 * =====================================================================================================================
 * Rotation gates
 * =====================================================================================================================
 *
 * Rx(θ) = exp(-iθX/2) = cos(θ/2)I - i·sin(θ/2)X
 *       = | cos(θ/2)    -i·sin(θ/2) |
 *         | -i·sin(θ/2)  cos(θ/2)   |
 *
 * Ry(θ) = exp(-iθY/2) = cos(θ/2)I - i·sin(θ/2)Y
 *       = | cos(θ/2)   -sin(θ/2) |
 *         | sin(θ/2)    cos(θ/2) |
 *
 * Rz(θ) = exp(-iθZ/2) = cos(θ/2)I - i·sin(θ/2)Z
 *       = | e^(-iθ/2)     0      |
 *         |    0       e^(iθ/2)  |
 */

void rx_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    cplx_t *data = state->data;

    #pragma omp parallel for collapse(2) if(dim >= OMP_THRESHOLD)
    for (dim_t i = 0; i < dim; i += step) {
        for (dim_t j = 0; j < stride; ++j) {
            cplx_t *a = data + i + j;
            cplx_t *b = data + i + j + stride;

            // a' = c*a - i*s*b,  b' = -i*s*a + c*b
            // -i*(x + yi) = y - xi
            cplx_t a_new = c * (*a) + CMPLX(s * cimag(*b), -s * creal(*b));
            cplx_t b_new = CMPLX(s * cimag(*a), -s * creal(*a)) + c * (*b);
            *a = a_new;
            *b = b_new;
        }
    }
}

void ry_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    cplx_t *data = state->data;

    #pragma omp parallel for collapse(2) if(dim >= OMP_THRESHOLD)
    for (dim_t i = 0; i < dim; i += step) {
        for (dim_t j = 0; j < stride; ++j) {
            cplx_t *a = data + i + j;
            cplx_t *b = data + i + j + stride;

            // a' = c*a - s*b,  b' = s*a + c*b
            cplx_t a_new = c * (*a) - s * (*b);
            cplx_t b_new = s * (*a) + c * (*b);
            *a = a_new;
            *b = b_new;
        }
    }
}

void rz_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    cplx_t *data = state->data;

    #pragma omp parallel for collapse(2) if(dim >= OMP_THRESHOLD)
    for (dim_t i = 0; i < dim; i += step) {
        for (dim_t j = 0; j < stride; ++j) {
            cplx_t *a = data + i + j;
            cplx_t *b = data + i + j + stride;

            // a' = e^(-iθ/2)*a = (c - is)*a
            // b' = e^(iθ/2)*b = (c + is)*b
            // (c - is)*(x + yi) = cx + sy + i(cy - sx)
            // (c + is)*(x + yi) = cx - sy + i(cy + sx)
            double a_re = creal(*a), a_im = cimag(*a);
            double b_re = creal(*b), b_im = cimag(*b);
            *a = CMPLX(c * a_re + s * a_im, c * a_im - s * a_re);
            *b = CMPLX(c * b_re - s * b_im, c * b_im + s * b_re);
        }
    }
}
