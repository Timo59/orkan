/**
 * @file gate_packed_2q.c
 * @brief Two-qubit gates for mixed states (density matrices in packed lower-triangular storage)
 *
 * Mixed states are stored in LAPACK packed lower-triangular column-major format.
 * For an N×N density matrix, we store N(N+1)/2 complex elements.
 *
 * Two-qubit permutation gates (CX, SWAP) transform ρ → UρU† by permuting rows
 * and columns of the density matrix. For packed Hermitian storage, elements that
 * cross the diagonal require conjugation.
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
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

#define SINGLE_CONTROL_SINGLE_TARGET_PACKED(state, control, target, BLOCK_OP) \
do {  \
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;  \
    const qubit_t hi = (control > target) ? control : target; \
    const qubit_t lo = (control < target) ? control : target; \
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;  \
    const gate_idx_t incr_tgt = (gate_idx_t)1 << target;  \
    cplx_t * restrict data = state->data; \
    const gate_idx_t quarter_dim = dim >> 2;  \
  \
    _Pragma("omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)")       \
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) { \
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi); \
        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt; \
  \
        /* Precompute offsets of columns */ \
        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;  \
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;  \
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;  \
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;  \
  \
        for (gate_idx_t br = bc; br < quarter_dim; ++br) {  \
            const gate_idx_t r00 = insertBits2_0(br, lo, hi); \
            const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt; \
  \
            /* The first two pairs are always in the lower triangle */  \
            /* |c=1,t=0><c=0,t=0| <--> |c=1,t=1><c=0,t=0| */  \
            gate_idx_t idx_a = (r10 - c00) + offset_c00;  \
            gate_idx_t idx_b = (r11 - c00) + offset_c00;  \
            cplx_t tmp = data[idx_a]; \
            data[idx_a] = data[idx_b];  \
            data[idx_b] = tmp;  \
  \
            /* |c=1,t=0><c=1,t=0| <--> |c=1,t=1><c=1,t=1| */  \
            idx_a = (r10 - c10) + offset_c10; \
            idx_b = (r11 - c11) + offset_c11; \
            tmp = data[idx_a];  \
            data[idx_a] = data[idx_b];  \
            data[idx_b] = tmp;  \
} while(0)
/*
 * =====================================================================================================================
 * CNOT gate: ρ → CX · ρ · CX†  (CX = CX† for CNOT)
 * =====================================================================================================================
 *
 * CNOT is a permutation gate that swaps basis states |ctrl=1,tgt=0> <-> |ctrl=1,tgt=1>.
 * For density matrices ρ' = CX·ρ·CX†, the transformation permutes rows and columns:
 *   new_rho[i,j] = rho[perm(i), perm(j)]
 * where perm(x) = x XOR (1<<target) if (x & (1<<control)), else x.
 *
 * This is a "zero-flop" gate in the sense that CNOT involves no matrix multiplication,
 * only permutation of elements. However, for packed Hermitian storage, some elements
 * require conjugation when they cross the diagonal.
 */
/*
 * In-place CNOT for packed storage.
 *
 * CNOT is an involution: perm(perm(x)) = x, where perm(x) = x ^ target_bit
 * when control_bit is set. Since the density matrix transformation is
 * new_rho[r,c] = old_rho[perm(r), perm(c)], and the permutation is its own
 * inverse, every stored element is either a fixed point or part of a 2-cycle.
 *
 * For each stored (r,c) with r >= c, we compute the destination (r',c') =
 * (perm(r), perm(c)) canonicalized to lower-triangle form. If this differs
 * from (r,c), we swap the two packed elements. To avoid double-swapping,
 * we only swap when the source packed index < destination packed index.
 *
 * When the source and destination are on opposite sides of the diagonal
 * (one has r>=c, the other has r<c), the values must be conjugated.
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
            gate_idx_t idx_a = (r10 - c00) + offset_c00;
            gate_idx_t idx_b = (r11 - c00) + offset_c00;
            cplx_t tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* |c=1,t=0><c=1,t=0| <--> |c=1,t=1><c=1,t=1| */
            idx_a = (r10 - c10) + offset_c10;
            idx_b = (r11 - c11) + offset_c11;
            tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Fast-path: Remaining 4 pairs are in the lower triangle. */
            if (r00 > c11) {
                /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                idx_a = (r11 - c01) + offset_c01;
                idx_b = (r10 - c01) + offset_c01;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;

                /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                idx_a = (r11 - c10) + offset_c10;
                idx_b = (r10 - c11) + offset_c11;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;

                /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                idx_a = (r01 - c10) + offset_c10;
                idx_b = (r01 - c11) + offset_c11;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;

                /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                idx_a = (r00 - c10) + offset_c10;
                idx_b = (r00 - c11) + offset_c11;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;
            }
            else {
                /* r10 > c11 implies r10 > c01 */
                if (r10 > c11) {
                    /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                    idx_a = (r11 - c10) + offset_c10;
                    idx_b = (r10 - c11) + offset_c11;
                    tmp = data[idx_a];
                    data[idx_a] = data[idx_b];
                    data[idx_b] = tmp;

                    /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                    idx_a = (r11 - c01) + offset_c01;
                    idx_b = (r10 - c01) + offset_c01;
                    tmp = data[idx_a];
                    data[idx_a] = data[idx_b];
                    data[idx_b] = tmp;
                }
                else {
                    /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                    idx_a = (r11 - c10) + offset_c10;
                    idx_b = pack_idx(dim, c11, r10);
                    tmp = data[idx_a];
                    data[idx_a] = conj(data[idx_b]);
                    data[idx_b] = conj(tmp);

                    /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                    idx_a = (r11 - c01) + offset_c01;
                    tmp = data[idx_a];
                    if (r10 > c01) {
                        idx_b = (r10 - c01) + offset_c01;
                        data[idx_a] = data[idx_b];
                        data[idx_b] = tmp;
                    }
                    else {
                        idx_b = pack_idx(dim, c01, r10);
                        data[idx_a] = conj(data[idx_b]);
                        data[idx_b] = conj(tmp);
                    }
                }

                if (bc != br) {
                    /* r00 > c10 implies r01 > c10 */
                    if (r00 > c10) {
                        /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                        /* r00 < c11, otherwise fast path would have been taken */
                        /* (r00, c11) is in the upper triangle => conj(c11, r00) */
                        idx_a = (r00 - c10) + offset_c10;
                        idx_b = pack_idx(dim, c11, r00);
                        tmp = data[idx_a];
                        data[idx_a] = conj(data[idx_b]);
                        data[idx_b] = conj(tmp);

                        /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                        /* r00 > c10 also implies r01 > c11 */
                        idx_a = (r01 - c10) + offset_c10;
                        idx_b = (r01 - c11) + offset_c11;
                        tmp = data[idx_a];
                        data[idx_a] = data[idx_b];
                        data[idx_b] = tmp;
                    }
                    else {
                        /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                        /* r00 < c10 implies both elements are in the upper triangle of a lower triangle block */
                        idx_a = pack_idx(dim, c10, r00);
                        idx_b = pack_idx(dim, c11, r00);
                        tmp = data[idx_a];
                        data[idx_a] = data[idx_b];
                        data[idx_b] = tmp;

                        /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                        /* r00 < c10 implies r01 < c11 */
                        idx_a = pack_idx(dim, c11, r01);
                        tmp = data[idx_a];
                        if (r01 > c10) {
                            idx_b = (r01 - c10) + offset_c10;
                            data[idx_a] = conj(data[idx_b]);
                            data[idx_b] = conj(tmp);
                        }
                        else {
                            idx_b = pack_idx(dim, c10, r01);
                            data[idx_a] = data[idx_b];
                            data[idx_b] = tmp;
                        }
                    }
                }
            }
        }
    }
}

void cy_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cy: packed not yet implemented");
}

void cz_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cz: packed not yet implemented");
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
 * In-place SWAP for packed storage.
 *
 * SWAP permutes basis states by exchanging bits `lo` and `hi`:
 *   perm(x) = x with bit lo and bit hi swapped
 * For density matrices: new_rho[r,c] = old_rho[perm(r), perm(c)].
 *
 * Like CNOT, SWAP is a zero-flop permutation gate — no arithmetic, only
 * element moves. The permutation maps each 4×4 block (indexed by the two-bit
 * pattern at positions lo,hi) into swap pairs:
 *
 *   Pair 1: (r01, c00) <-> (r10, c00)   always lower triangle
 *   Pair 2: (r01, c01) <-> (r10, c10)   always lower triangle
 *   Pair 3: (r11, c01) <-> (r11, c10)   always lower triangle
 *   Pair 4: (r10, c01) <-> (r01, c10)   may cross diagonal
 *   Pair 5: (r00, c01) <-> (r00, c10)   may cross diagonal
 *   Pair 6: (r10, c11) <-> (r01, c11)   may cross diagonal
 *
 * We iterate only the lower-triangular half of the block grid (br >= bc).
 * When a swap partner falls in the upper triangle (r < c), its packed
 * representative is at the transposed position with a conjugate. These
 * elements belong to block (bc, br) which is never iterated, so both
 * directions of the swap must be written in the current iteration.
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
            gate_idx_t idx_a = (r01 - c00) + offset_c00;
            gate_idx_t idx_b = (r10 - c00) + offset_c00;
            cplx_t tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Pair 2: (r01, c01) <-> (r10, c10) — both always in lower triangle */
            idx_a = (r01 - c01) + offset_c01;
            idx_b = (r10 - c10) + offset_c10;
            tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Pair 3: (r11, c01) <-> (r11, c10) — both always in lower triangle */
            idx_a = (r11 - c01) + offset_c01;
            idx_b = (r11 - c10) + offset_c10;
            tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Fast path: If r00 > c10, all elements are contained in the lower triangle */
            if (r00 > c10) {
                /* |q1=0,q2=0><q1=0,q2=1| <--> |q1=0,q2=0><q1=1,q2=0| */
                idx_a = (r00 - c01) + offset_c01;
                idx_b = (r00 - c10) + offset_c10;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;

                /* |q1=0,q2=1><q1=1,q2=1| <--> |q1=1,q2=0><q1=1,q2=1| */
                idx_a = (r01 - c11) + offset_c11;
                idx_b = (r10 - c11) + offset_c11;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;

                /* |q1=0,q2=1><q1=1,q2=0| <--> |q1=1,q2=0><q1=0,q2=1| */
                idx_a = (r01 - c10) + offset_c10;
                idx_b = (r10 - c01) + offset_c01;
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;
            }
            else {
                /* Pair 4: (r10, c01) <-> (r01, c10)
                 * (r10, c01) is always lower-tri. (r01, c10) may cross the diagonal:
                 *   r01 > c10: both lower-tri, plain swap.
                 *   r01 <= c10: (r01, c10) is upper-tri; its packed representative is
                 *     conj(data) at (c10, r01). That position belongs to block (bc, br),
                 *     never iterated, so we must write both directions here. */
                idx_a = (r10 - c01) + offset_c01;
                tmp = data[idx_a];
                if (r01 > c10) {
                    idx_b = (r01 - c10) + offset_c10;
                    data[idx_a] = data[idx_b];
                    data[idx_b] = tmp;
                }
                else {
                    idx_b = pack_idx(dim, c10, r01);
                    data[idx_a] = conj(data[idx_b]);
                    data[idx_b] = conj(tmp);
                }

                /* r00 > c01 implies r10 > c11 */
                if (r00 > c01) {
                    /* |q1=0,q2=0><q1=0,q2=1| <--> |q1=0,q2=0><q1=1,q2=0| */
                    idx_a = (r00 - c01) + offset_c01;
                    idx_b = pack_idx(dim, c10, r00);
                    tmp = data[idx_a];
                    data[idx_a] = conj(data[idx_b]);
                    data[idx_b] = conj(tmp);

                    /* |q1=0,q2=1><q1=1,q2=1| <--> |q1=1,q2=0><q1=1,q2=1| */
                    idx_a = (r10 - c11) + offset_c11;
                    idx_b = pack_idx(dim, c11, r01);
                    tmp = data[idx_a];
                    data[idx_a] = conj(data[idx_b]);
                    data[idx_b] = conj(tmp);
                }
                else if (bc != br) {
                    /* |q1=0,q2=0><q1=0,q2=1| <--> |q1=0,q2=0><q1=1,q2=0| */
                    idx_a = pack_idx(dim, c01, r00);
                    idx_b = pack_idx(dim, c10, r00);
                    tmp = data[idx_a];
                    data[idx_a] = data[idx_b];
                    data[idx_b] = tmp;

                    /* |q1=0,q2=1><q1=1,q2=1| <--> |q1=1,q2=0><q1=1,q2=1| */
                    idx_a = pack_idx(dim, c11, r10);
                    idx_b = pack_idx(dim, c11, r01);
                    tmp = data[idx_a];
                    data[idx_a] = data[idx_b];
                    data[idx_b] = tmp;
                }
            }
        }
    }
}

void ccx_packed(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    GATE_VALIDATE(0, "ccx: packed not yet implemented");
}
