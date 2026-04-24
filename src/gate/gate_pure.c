// gate_pure.c - Quantum gates for pure states (Hilbert space vectors)

#include "gate.h"
#include "index.h"

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 * Pure states: O(dim) work per gate, need dim >= 4096 (~12 qubits) to offset thread overhead
 */

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
    const idx_t dim = (idx_t)1 << (state)->qubits;                   \
    const idx_t stride = (idx_t)1 << (target);                       \
    const idx_t step = (idx_t)1 << ((target) + 1);                   \
    cplx_t *data = (state)->data;                                              \
                                                                               \
    _Pragma("omp parallel for collapse(2) if(dim >= OMP_THRESHOLD_PURE)")           \
    for (idx_t i = 0; i < dim; i += step) {                               \
        for (idx_t j = 0; j < stride; ++j) {                              \
            cplx_t *a = data + i + j;                                          \
            cplx_t *b = data + i + j + stride;                                 \
            PAIR_OP(a, b);                                                     \
        }                                                                      \
    }                                                                          \
} while(0)

/*
 * Rotation gate pair-op macros: these capture trig constants (c, s) from the
 * enclosing scope and are used with TRAVERSE_PURE_1Q. Defined adjacent to
 * their use in the rotation gate section below.
 */

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

// Hy gate: a' = (a - b)/√2, b' = (a + b)/√2
#define HY_OP(a, b) do {                            \
    cplx_t diff = *(a) - *(b);                       \
    cplx_t sum  = *(a) + *(b);                       \
    *(a) = diff * M_SQRT1_2;                         \
    *(b) = sum  * M_SQRT1_2;                         \
} while(0)

void hy_pure(state_t *state, const qubit_t target) {
    TRAVERSE_PURE_1Q(state, target, HY_OP);
}
#undef HY_OP

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

// P(θ) gate: a unchanged, b *= e^(iθ)
// b *= e^(iθ) = (c+is)*(x+yi) = (cx-sy) + i(sx+cy)
#define P_OP(a, b) do {                                                      \
    (void)(a);                                                               \
    double re = creal(*(b));                                                 \
    double im = cimag(*(b));                                                 \
    *(b) = CMPLX(c * re - s * im, s * re + c * im);                         \
} while(0)

void p_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta);
    const double s = sin(theta);
    TRAVERSE_PURE_1Q(state, target, P_OP);
}
#undef P_OP


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

/*
 * Rotation pair-op macros: capture c, s from enclosing scope.
 * These are #undef'd after use to avoid polluting the namespace.
 */

// Rx: a' = c*a - i*s*b,  b' = -i*s*a + c*b
// -i*(x + yi) = y - xi
#define RX_OP(a, b) do {                                                     \
    cplx_t a_new = c * (*(a)) + CMPLX(s * cimag(*(b)), -s * creal(*(b)));   \
    cplx_t b_new = CMPLX(s * cimag(*(a)), -s * creal(*(a))) + c * (*(b));   \
    *(a) = a_new;                                                            \
    *(b) = b_new;                                                            \
} while(0)

void rx_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);
    TRAVERSE_PURE_1Q(state, target, RX_OP);
}
#undef RX_OP

// Ry: a' = c*a - s*b,  b' = s*a + c*b
#define RY_OP(a, b) do {                                                     \
    cplx_t a_new = c * (*(a)) - s * (*(b));                                  \
    cplx_t b_new = s * (*(a)) + c * (*(b));                                  \
    *(a) = a_new;                                                            \
    *(b) = b_new;                                                            \
} while(0)

void ry_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);
    TRAVERSE_PURE_1Q(state, target, RY_OP);
}
#undef RY_OP

// Rz: a' = e^(-iθ/2)*a = (c-is)*a,  b' = e^(iθ/2)*b = (c+is)*b
// (c-is)*(x+yi) = cx+sy + i(cy-sx),  (c+is)*(x+yi) = cx-sy + i(cy+sx)
#define RZ_OP(a, b) do {                                                     \
    double a_re = creal(*(a));                                               \
    double a_im = cimag(*(a));                                               \
    double b_re = creal(*(b));                                               \
    double b_im = cimag(*(b));                                               \
    *(a) = CMPLX(c * a_re + s * a_im, c * a_im - s * a_re);                 \
    *(b) = CMPLX(c * b_re - s * b_im, c * b_im + s * b_re);                 \
} while(0)

void rz_pure(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);
    TRAVERSE_PURE_1Q(state, target, RZ_OP);
}
#undef RZ_OP

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

/*
 * CNOT (Controlled-X) gate: Flips target qubit when control qubit is |1>
 *
 * For each basis state index with control=1, we swap amplitudes where target=0 and target=1.
 * This is a "zero-flop" gate requiring only conditional memory swaps.
 */
void cx_pure(state_t *state, const qubit_t control, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    // Determine which qubit index is lower/higher for bit insertion
    qubit_t lo = (control < target) ? control : target;
    qubit_t hi = (control < target) ? target : control;

    const idx_t incr_ctrl = (idx_t)1 << control;
    const idx_t incr_tgt = (idx_t)1 << target;
    const idx_t n_base = dim >> 2;  // dim / 4: we iterate over indices with both bits = 0

    #pragma omp parallel for if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t k = 0; k < n_base; ++k) {
        // Insert 0 bits at positions lo and hi to get base index
        const idx_t base = insertBits2_0(k, lo, hi);

        // Only swap when control=1: swap |ctrl=1,tgt=0> <-> |ctrl=1,tgt=1>
        const idx_t idx_c1_t0 = base | incr_ctrl;           // control=1, target=0
        const idx_t idx_c1_t1 = base | incr_ctrl | incr_tgt; // control=1, target=1

        cplx_t tmp = data[idx_c1_t0];
        data[idx_c1_t0] = data[idx_c1_t1];
        data[idx_c1_t1] = tmp;
    }
}


/*
 * CY (Controlled-Y) gate: Applies Y to target qubit when control qubit is |1>
 *
 * When control=1: a' = -i*b, b' = i*a  (where a=target|0>, b=target|1>)
 * This is a "zero-flop" gate requiring only swaps and sign/component shuffles.
 */
void cy_pure(state_t *state, const qubit_t control, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    qubit_t lo = (control < target) ? control : target;
    qubit_t hi = (control < target) ? target : control;

    const idx_t incr_ctrl = (idx_t)1 << control;
    const idx_t incr_tgt = (idx_t)1 << target;
    const idx_t n_base = dim >> 2;

    #pragma omp parallel for if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t k = 0; k < n_base; ++k) {
        const idx_t base = insertBits2_0(k, lo, hi);

        const idx_t idx_c1_t0 = base | incr_ctrl;
        const idx_t idx_c1_t1 = base | incr_ctrl | incr_tgt;

        // Y gate: a' = -i*b, b' = i*a
        double a_re = creal(data[idx_c1_t0]), a_im = cimag(data[idx_c1_t0]);
        double b_re = creal(data[idx_c1_t1]), b_im = cimag(data[idx_c1_t1]);
        data[idx_c1_t0] = CMPLX(b_im, -b_re);   // -i*b
        data[idx_c1_t1] = CMPLX(-a_im, a_re);    //  i*a
    }
}

/*
 * CZ (Controlled-Z) gate: Applies Z to target qubit when control qubit is |1>
 *
 * When control=1: only the |1> component of the target is negated.
 * This is a "zero-flop" gate requiring only a sign flip.
 */
void cz_pure(state_t *state, const qubit_t control, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    qubit_t lo = (control < target) ? control : target;
    qubit_t hi = (control < target) ? target : control;

    const idx_t incr_ctrl = (idx_t)1 << control;
    const idx_t incr_tgt = (idx_t)1 << target;
    const idx_t n_base = dim >> 2;

    #pragma omp parallel for if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t k = 0; k < n_base; ++k) {
        const idx_t base = insertBits2_0(k, lo, hi);

        // Negate amplitude where both control=1 and target=1
        const idx_t idx_c1_t1 = base | incr_ctrl | incr_tgt;
        data[idx_c1_t1] = -data[idx_c1_t1];
    }
}

/*
 * SWAP gate: Exchanges two qubits by swapping amplitudes where q1 and q2 differ.
 *
 * For each basis state with both bits = 0, we swap |q1=0,q2=1> <-> |q1=1,q2=0>.
 * This is a zero-flop gate requiring only memory swaps, no arithmetic.
 */
void swap_pure(state_t *state, const qubit_t q1, const qubit_t q2) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    const qubit_t lo = (q1 < q2) ? q1 : q2;
    const qubit_t hi = (q1 < q2) ? q2 : q1;
    const idx_t incr_lo = (idx_t)1 << lo;
    const idx_t incr_hi = (idx_t)1 << hi;
    const idx_t n_base = dim >> 2;  // dim / 4: we iterate over indices with both bits = 0

    #pragma omp parallel for if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t k = 0; k < n_base; ++k) {
        const idx_t base = insertBits2_0(k, lo, hi);

        // Swap |q1=0,q2=1> <-> |q1=1,q2=0>
        const idx_t idx01 = base | incr_lo;
        const idx_t idx10 = base | incr_hi;
        cplx_t tmp = data[idx01];
        data[idx01] = data[idx10];
        data[idx10] = tmp;
    }
}

/*
 * CCX (Toffoli) gate: Flips target qubit when both control qubits are |1>
 *
 * For each basis state index with ctrl1=1 and ctrl2=1, we swap amplitudes
 * where target=0 and target=1. This is a zero-flop gate requiring only
 * conditional memory swaps.
 */
void ccx_pure(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    // Sort qubit positions into lo < med < hi for bit insertion
    qubit_t q[3] = {ctrl1, ctrl2, target};
    if (q[0] > q[1]) { qubit_t t = q[0]; q[0] = q[1]; q[1] = t; }
    if (q[1] > q[2]) { qubit_t t = q[1]; q[1] = q[2]; q[2] = t; }
    if (q[0] > q[1]) { qubit_t t = q[0]; q[0] = q[1]; q[1] = t; }
    const qubit_t lo = q[0], med = q[1], hi = q[2];

    const idx_t incr_ctrl1 = (idx_t)1 << ctrl1;
    const idx_t incr_ctrl2 = (idx_t)1 << ctrl2;
    const idx_t incr_tgt   = (idx_t)1 << target;
    const idx_t n_base = dim >> 3;  // dim / 8: iterate over indices with all three bits = 0

    #pragma omp parallel for if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t k = 0; k < n_base; ++k) {
        // Insert 0 bits at positions lo, med, hi to get base index
        const idx_t base = insertBit0(insertBit0(insertBit0(k, lo), med), hi);

        // Only swap when both controls are 1
        const idx_t idx_c11_t0 = base | incr_ctrl1 | incr_ctrl2;
        const idx_t idx_c11_t1 = base | incr_ctrl1 | incr_ctrl2 | incr_tgt;

        cplx_t tmp = data[idx_c11_t0];
        data[idx_c11_t0] = data[idx_c11_t1];
        data[idx_c11_t1] = tmp;
    }
}
