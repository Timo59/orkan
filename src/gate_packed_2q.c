/**
 * @file gate_packed_2q.c
 * @brief Two-qubit gates for mixed states in packed lower-triangular storage
 *
 * Storage: LAPACK packed lower-triangle, column-major, N(N+1)/2 complex elements.
 *
 * All two-qubit gates iterate over quarter-dim blocks indexed by (br, bc) with
 * br >= bc. Each block addresses 4 rows (r00, r01, r10, r11) × 4 columns,
 * producing 6 element-swap pairs. Three pairs always land in the lower triangle;
 * three may cross the diagonal and require conjugation on read/write.
 *
 * Key pitfalls:
 *   - Upper-triangle access: read/write as conj(data[pack_idx(dim, c, r)])
 *   - Diagonal blocks (bc == br): Hermitian pairs share storage; skip one side
 *   - insertBits2_0() requires lo < hi
 */

#include "gate.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems) */
#define OMP_THRESHOLD 64


/**
 * @brief Compute packed array index for element (r, c) where r >= c
 */
static inline gate_idx_t pack_idx(gate_idx_t dim, gate_idx_t r, gate_idx_t c) {
    return c * (2 * dim - c + 1) / 2 + (r - c);
}


/*
 * CX (CNOT): ρ' = CX·ρ·CX, zero-flop permutation gate (CX is self-adjoint).
 *
 * perm(x) = x ^ (1<<target) when bit control is set, else x.
 * ρ'[r,c] = ρ[perm(r), perm(c)].  Involution: each element is fixed or in a 2-cycle.
 *
 * For each block (br, bc) the 6 swap pairs split into:
 *   - Fast path (r00 > c11): all lower-tri, plain swaps
 *   - Slow path: up to 4 pairs may cross the diagonal, requiring conjugation
 *   - Diagonal blocks (bc == br): skip group B (Hermitian pairs share storage)
 */
void cx_packed(state_t * restrict state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = (control > target) ? control : target;
    const qubit_t lo = (control < target) ? control : target;
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

        /* Precompute offsets of columns */
        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* The first two pairs are always in the lower triangle */
            /* |c=1,t=0><c=0,t=0| <--> |c=1,t=1><c=0,t=0| */
            cplx_t tmp = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = data[(r11 - c00) + offset_c00];
            data[(r11 - c00) + offset_c00] = tmp;

            /* |c=1,t=0><c=1,t=0| <--> |c=1,t=1><c=1,t=1| */
            tmp = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = data[(r11 - c11) + offset_c11];
            data[(r11 - c11) + offset_c11] = tmp;

            /* Fast-path: Remaining 4 pairs are in the lower triangle. */
            if (r00 > c11) {
                /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                tmp = data[(r11 - c01) + offset_c01];
                data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                data[(r10 - c01) + offset_c01] = tmp;

                /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                tmp = data[(r11 - c10) + offset_c10];
                data[(r11 - c10) + offset_c10] = data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = tmp;

                /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = data[(r01 - c11) + offset_c11];
                data[(r01 - c11) + offset_c11] = tmp;

                /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                tmp = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = data[(r00 - c11) + offset_c11];
                data[(r00 - c11) + offset_c11] = tmp;
            }
            else {
                /* r10 > c11 implies r10 > c01 */
                if (r10 > c11) {
                    /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = tmp;

                    /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                    tmp = data[(r11 - c01) + offset_c01];
                    data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                    data[(r10 - c01) + offset_c01] = tmp;
                }
                else {
                    /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = conj(data[pack_idx(dim, c11, r10)]);
                    data[pack_idx(dim, c11, r10)] = conj(tmp);

                    /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                    tmp = data[(r11 - c01) + offset_c01];
                    if (r10 > c01) {
                        data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                        data[(r10 - c01) + offset_c01] = tmp;
                    }
                    else {
                        data[(r11 - c01) + offset_c01] = conj(data[pack_idx(dim, c01, r10)]);
                        data[pack_idx(dim, c01, r10)] = conj(tmp);
                    }
                }

                if (bc != br) {
                    /* r00 > c10 implies r01 > c10 */
                    if (r00 > c10) {
                        /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                        /* r00 < c11, otherwise fast path would have been taken */
                        /* (r00, c11) is in the upper triangle => conj(c11, r00) */
                        tmp = data[(r00 - c10) + offset_c10];
                        data[(r00 - c10) + offset_c10] = conj(data[pack_idx(dim, c11, r00)]);
                        data[pack_idx(dim, c11, r00)] = conj(tmp);

                        /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                        /* r00 > c10 also implies r01 > c11 */
                        tmp = data[(r01 - c10) + offset_c10];
                        data[(r01 - c10) + offset_c10] = data[(r01 - c11) + offset_c11];
                        data[(r01 - c11) + offset_c11] = tmp;
                    }
                    else {
                        /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                        /* r00 < c10 implies both elements are in the upper triangle of a lower triangle block */
                        tmp = data[pack_idx(dim, c10, r00)];
                        data[pack_idx(dim, c10, r00)] = data[pack_idx(dim, c11, r00)];
                        data[pack_idx(dim, c11, r00)] = tmp;

                        /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                        /* r00 < c10 implies r01 < c11 */
                        tmp = data[pack_idx(dim, c11, r01)];
                        if (r01 > c10) {
                            data[pack_idx(dim, c11, r01)] = conj(data[(r01 - c10) + offset_c10]);
                            data[(r01 - c10) + offset_c10] = conj(tmp);
                        }
                        else {
                            data[pack_idx(dim, c11, r01)] = data[pack_idx(dim, c10, r01)];
                            data[pack_idx(dim, c10, r01)] = tmp;
                        }
                    }
                }
            }
        }
    }
}

/*
 * CY (controlled-Y): same permutation as CX with per-element phase factors.
 * CY|10⟩ = i|11⟩, CY|11⟩ = -i|10⟩.  ρ'[a,b] = φ(a)·φ*(b)·ρ[perm(a), perm(b)]
 * where φ(00)=φ(01)=1, φ(10)=-i, φ(11)=i.
 *
 * Per-pair phases (same 6 pairs as CX):
 *   A (r10,c00 ↔ r11,c00): ×(-i, i)     E (r10,c10 ↔ r11,c11): plain swap
 *   B (r11,c01 ↔ r10,c01): ×(i, -i)     F (r11,c10 ↔ r10,c11): negate-swap
 *   C,D (r00/r01,c10 ↔ c11): ×(i, -i)
 *
 * Pitfall: when a pair crosses the diagonal, phase and conj() compose
 * differently.  E.g. ×i on lower-tri value v: CMPLX(-Im(v), Re(v)); but
 * ×i on upper-tri stored s (where rho = conj(s)): result stored as
 * CMPLX(Im(s), Re(s)), not CMPLX(-Im(s), Re(s)).
 */
void cy_packed(state_t * restrict state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = (control > target) ? control : target;
    const qubit_t lo = (control < target) ? control : target;
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* Pair A: (r10,c00) ↔ (r11,c00) — always lower tri; phases (-i, i) */
            cplx_t tmp = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = CMPLX(cimag(data[(r11 - c00) + offset_c00]), -creal(data[(r11 - c00) + offset_c00]));
            data[(r11 - c00) + offset_c00] = CMPLX(-cimag(tmp), creal(tmp));

            /* Pair E: (r10,c10) ↔ (r11,c11) — always lower tri; plain swap */
            tmp = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = data[(r11 - c11) + offset_c11];
            data[(r11 - c11) + offset_c11] = tmp;

            /* Fast path: all remaining pairs in lower triangle */
            if (r00 > c11) {
                /* Pair B: (r11,c01) ↔ (r10,c01); phases (i, -i) */
                tmp = data[(r11 - c01) + offset_c01];
                data[(r11 - c01) + offset_c01] = CMPLX(-cimag(data[(r10 - c01) + offset_c01]), creal(data[(r10 - c01) + offset_c01]));
                data[(r10 - c01) + offset_c01] = CMPLX(cimag(tmp), -creal(tmp));

                /* Pair F: (r11,c10) ↔ (r10,c11); negate-swap */
                tmp = data[(r11 - c10) + offset_c10];
                data[(r11 - c10) + offset_c10] = -data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = -tmp;

                /* Pair D: (r01,c10) ↔ (r01,c11); phases (i, -i) */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = CMPLX(-cimag(data[(r01 - c11) + offset_c11]), creal(data[(r01 - c11) + offset_c11]));
                data[(r01 - c11) + offset_c11] = CMPLX(cimag(tmp), -creal(tmp));

                /* Pair C: (r00,c10) ↔ (r00,c11); phases (i, -i) */
                tmp = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = CMPLX(-cimag(data[(r00 - c11) + offset_c11]), creal(data[(r00 - c11) + offset_c11]));
                data[(r00 - c11) + offset_c11] = CMPLX(cimag(tmp), -creal(tmp));
            }
            else {
                if (r10 > c11) {
                    /* Pair F: negate-swap, both lower tri */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = -data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = -tmp;

                    /* Pair B: phases (i, -i), both lower tri */
                    tmp = data[(r11 - c01) + offset_c01];
                    data[(r11 - c01) + offset_c01] = CMPLX(-cimag(data[(r10 - c01) + offset_c01]), creal(data[(r10 - c01) + offset_c01]));
                    data[(r10 - c01) + offset_c01] = CMPLX(cimag(tmp), -creal(tmp));
                }
                else {
                    /* Pair F: (r11,c10) lower, (r10,c11) upper — negate + conj */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = -conj(data[pack_idx(dim, c11, r10)]);
                    data[pack_idx(dim, c11, r10)] = -conj(tmp);

                    /* Pair B: (r11,c01) lower; (r10,c01) may cross diagonal */
                    tmp = data[(r11 - c01) + offset_c01];
                    if (r10 > c01) {
                        data[(r11 - c01) + offset_c01] = CMPLX(-cimag(data[(r10 - c01) + offset_c01]), creal(data[(r10 - c01) + offset_c01]));
                        data[(r10 - c01) + offset_c01] = CMPLX(cimag(tmp), -creal(tmp));
                    }
                    else {
                        data[(r11 - c01) + offset_c01] = CMPLX(cimag(data[pack_idx(dim, c01, r10)]), creal(data[pack_idx(dim, c01, r10)]));
                        data[pack_idx(dim, c01, r10)] = CMPLX(cimag(tmp), creal(tmp));
                    }
                }

                if (bc != br) {
                    if (r00 > c10) {
                        /* Pair C: (r00,c10) lower, (r00,c11) upper */
                        /* r00 < c11, otherwise fast path would have been taken */
                        tmp = data[(r00 - c10) + offset_c10];
                        data[(r00 - c10) + offset_c10] = CMPLX(cimag(data[pack_idx(dim, c11, r00)]), creal(data[pack_idx(dim, c11, r00)]));
                        data[pack_idx(dim, c11, r00)] = CMPLX(cimag(tmp), creal(tmp));

                        /* Pair D: both lower tri (r00 > c10 implies r01 > c11) */
                        tmp = data[(r01 - c10) + offset_c10];
                        data[(r01 - c10) + offset_c10] = CMPLX(-cimag(data[(r01 - c11) + offset_c11]), creal(data[(r01 - c11) + offset_c11]));
                        data[(r01 - c11) + offset_c11] = CMPLX(cimag(tmp), -creal(tmp));
                    }
                    else {
                        /* Pair C: both upper tri */
                        tmp = data[pack_idx(dim, c10, r00)];
                        data[pack_idx(dim, c10, r00)] = CMPLX(cimag(data[pack_idx(dim, c11, r00)]), -creal(data[pack_idx(dim, c11, r00)]));
                        data[pack_idx(dim, c11, r00)] = CMPLX(-cimag(tmp), creal(tmp));

                        /* Pair D: (r01,c11) always upper; (r01,c10) may cross */
                        /* r00 < c10 implies r01 < c11 */
                        tmp = data[pack_idx(dim, c11, r01)];
                        if (r01 > c10) {
                            data[pack_idx(dim, c11, r01)] = CMPLX(cimag(data[(r01 - c10) + offset_c10]), creal(data[(r01 - c10) + offset_c10]));
                            data[(r01 - c10) + offset_c10] = CMPLX(cimag(tmp), creal(tmp));
                        }
                        else {
                            data[pack_idx(dim, c11, r01)] = CMPLX(-cimag(data[pack_idx(dim, c10, r01)]), creal(data[pack_idx(dim, c10, r01)]));
                            data[pack_idx(dim, c10, r01)] = CMPLX(cimag(tmp), -creal(tmp));
                        }
                    }
                }
            }
        }
    }
}

/*
 * CZ: diag(1,1,1,-1), zero-flop diagonal gate — only sign flips.
 *
 * Negate elements where exactly one of (row, col) has both ctrl+tgt bits set:
 *   Group A (row=11): (r11,c00), (r11,c01), (r11,c10) — always lower-tri
 *   Group B (col=11): (r00,c11), (r01,c11), (r10,c11) — may cross diagonal
 *
 * Diagonal blocks (bc == br): skip group B (same storage as A, avoids double-negate).
 */
void cz_packed(state_t *state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = (control > target) ? control : target;
    const qubit_t lo = (control < target) ? control : target;
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* Group A: row = 11 — always in lower triangle */
            data[(r11 - c00) + offset_c00] = -data[(r11 - c00) + offset_c00];
            data[(r11 - c01) + offset_c01] = -data[(r11 - c01) + offset_c01];
            data[(r11 - c10) + offset_c10] = -data[(r11 - c10) + offset_c10];

            /* Group B: col = 11 — skip on diagonal block (Hermitian pairs already negated) */
            if (bc != br) {
                if (r00 > c11) {
                    /* Fast path: all in lower triangle */
                    data[(r00 - c11) + offset_c11] = -data[(r00 - c11) + offset_c11];
                    data[(r01 - c11) + offset_c11] = -data[(r01 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = -data[(r10 - c11) + offset_c11];
                }
                else {
                    /* (r00, c11): always upper tri here — negate stored at (c11, r00) */
                    data[pack_idx(dim, c11, r00)] = -data[pack_idx(dim, c11, r00)];

                    /* (r01, c11) */
                    if (r01 > c11) {
                        data[(r01 - c11) + offset_c11] = -data[(r01 - c11) + offset_c11];
                    } else {
                        data[pack_idx(dim, c11, r01)] = -data[pack_idx(dim, c11, r01)];
                    }

                    /* (r10, c11) */
                    if (r10 > c11) {
                        data[(r10 - c11) + offset_c11] = -data[(r10 - c11) + offset_c11];
                    } else {
                        data[pack_idx(dim, c11, r10)] = -data[pack_idx(dim, c11, r10)];
                    }
                }
            }
        }
    }
}

void cs_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cs: packed not yet implemented");
}

void csdg_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "csdg: packed not yet implemented");
}

void ch_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ch: packed not yet implemented");
}

void chy_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "chy: packed not yet implemented");
}

void ct_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ct: packed not yet implemented");
}

void ctdg_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ctdg: packed not yet implemented");
}

void cp_packed(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cp: packed not yet implemented");
}

void cpdg_packed(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cpdg: packed not yet implemented");
}

/*
 * SWAP: zero-flop permutation gate.  perm(x) = x with bits lo and hi exchanged.
 * ρ'[r,c] = ρ[perm(r), perm(c)].
 *
 * Swap pairs per block:
 *   1: (r01,c00) ↔ (r10,c00)  always lower    4: (r10,c01) ↔ (r01,c10)  may cross
 *   2: (r01,c01) ↔ (r10,c10)  always lower    5: (r00,c01) ↔ (r00,c10)  may cross
 *   3: (r11,c01) ↔ (r11,c10)  always lower    6: (r10,c11) ↔ (r01,c11)  may cross
 *
 * Fast path: r00 > c10 ⟹ all lower-tri.  Slow path: upper-tri partners live
 * in block (bc, br) which is never iterated; both swap directions must be
 * written here with conjugation.
 */
void swap_packed(state_t *state, const qubit_t q1, const qubit_t q2) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = q1 > q2 ? q1 : q2;
    const qubit_t lo = q1 < q2 ? q1 : q2;
    const gate_idx_t incr_hi = (gate_idx_t)1 << hi;
    const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_lo, c10 = c00 | incr_hi, c11 = c10 | incr_lo;

        /* Precompute offsets of columns */
        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;


        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_lo, r10 = r00 | incr_hi, r11 = r10 | incr_lo;

            /* Pair 1: (r01, c00) <-> (r10, c00) — both always in lower triangle */
            cplx_t tmp = data[(r01 - c00) + offset_c00];
            data[(r01 - c00) + offset_c00] = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = tmp;

            /* Pair 2: (r01, c01) <-> (r10, c10) — both always in lower triangle */
            tmp = data[(r01 - c01) + offset_c01];
            data[(r01 - c01) + offset_c01] = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = tmp;

            /* Pair 3: (r11, c01) <-> (r11, c10) — both always in lower triangle */
            tmp = data[(r11 - c01) + offset_c01];
            data[(r11 - c01) + offset_c01] = data[(r11 - c10) + offset_c10];
            data[(r11 - c10) + offset_c10] = tmp;

            /* Fast path: r00 > c10 ⟹ all remaining pairs in lower triangle */
            if (r00 > c10) {
                /* Pair 5: (r00,c01) ↔ (r00,c10) */
                tmp = data[(r00 - c01) + offset_c01];
                data[(r00 - c01) + offset_c01] = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = tmp;

                /* Pair 6: (r01,c11) ↔ (r10,c11) */
                tmp = data[(r01 - c11) + offset_c11];
                data[(r01 - c11) + offset_c11] = data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = tmp;

                /* Pair 4: (r01,c10) ↔ (r10,c01) */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = data[(r10 - c01) + offset_c01];
                data[(r10 - c01) + offset_c01] = tmp;
            }
            else {
                /* Pair 4: (r10, c01) <-> (r01, c10)
                 * (r10, c01) is always lower-tri. (r01, c10) may cross the diagonal:
                 *   r01 > c10: both lower-tri, plain swap.
                 *   r01 <= c10: (r01, c10) is upper-tri; its packed representative is
                 *     conj(data) at (c10, r01). That position belongs to block (bc, br),
                 *     never iterated, so we must write both directions here. */
                tmp = data[(r10 - c01) + offset_c01];
                if (r01 > c10) {
                    data[(r10 - c01) + offset_c01] = data[(r01 - c10) + offset_c10];
                    data[(r01 - c10) + offset_c10] = tmp;
                }
                else {
                    data[(r10 - c01) + offset_c01] = conj(data[pack_idx(dim, c10, r01)]);
                    data[pack_idx(dim, c10, r01)] = conj(tmp);
                }

                /* r00 > c01 implies r10 > c11 */
                if (r00 > c01) {
                    /* Pair 5: (r00,c01) lower ↔ (r00,c10) upper — conj swap */
                    tmp = data[(r00 - c01) + offset_c01];
                    data[(r00 - c01) + offset_c01] = conj(data[pack_idx(dim, c10, r00)]);
                    data[pack_idx(dim, c10, r00)] = conj(tmp);

                    /* Pair 6: (r10,c11) lower ↔ (r01,c11) upper — conj swap */
                    tmp = data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = conj(data[pack_idx(dim, c11, r01)]);
                    data[pack_idx(dim, c11, r01)] = conj(tmp);
                }
                else if (bc != br) {
                    /* Pair 5: (r00,c01) ↔ (r00,c10) — both upper-tri */
                    tmp = data[pack_idx(dim, c01, r00)];
                    data[pack_idx(dim, c01, r00)] = data[pack_idx(dim, c10, r00)];
                    data[pack_idx(dim, c10, r00)] = tmp;

                    /* Pair 6: (r01,c11) ↔ (r10,c11) — both upper-tri */
                    tmp = data[pack_idx(dim, c11, r10)];
                    data[pack_idx(dim, c11, r10)] = data[pack_idx(dim, c11, r01)];
                    data[pack_idx(dim, c11, r01)] = tmp;
                }
            }
        }
    }
}

void ccx_packed(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    GATE_VALIDATE(0, "ccx: packed not yet implemented");
}
