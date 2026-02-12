// gate.c   - Gate dispatchers with input validation and three-way type dispatch

#include "gate.h"

/*
 * =====================================================================================================================
 * Pure state gate implementations (gate_pure.c)
 * =====================================================================================================================
 */

extern void x_pure(state_t *state, qubit_t target);
extern void y_pure(state_t *state, qubit_t target);
extern void z_pure(state_t *state, qubit_t target);
extern void h_pure(state_t *state, qubit_t target);
extern void s_pure(state_t *state, qubit_t target);
extern void sdg_pure(state_t *state, qubit_t target);
extern void t_pure(state_t *state, qubit_t target);
extern void tdg_pure(state_t *state, qubit_t target);
extern void hy_pure(state_t *state, qubit_t target);

extern void rx_pure(state_t *state, qubit_t target, double theta);
extern void ry_pure(state_t *state, qubit_t target, double theta);
extern void rz_pure(state_t *state, qubit_t target, double theta);
extern void p_pure(state_t *state, qubit_t target, double theta);

extern void cx_pure(state_t *state, qubit_t control, qubit_t target);
extern void cy_pure(state_t *state, qubit_t control, qubit_t target);
extern void cz_pure(state_t *state, qubit_t control, qubit_t target);
extern void cs_pure(state_t *state, qubit_t control, qubit_t target);
extern void csdg_pure(state_t *state, qubit_t control, qubit_t target);
extern void ch_pure(state_t *state, qubit_t control, qubit_t target);
extern void chy_pure(state_t *state, qubit_t control, qubit_t target);
extern void ct_pure(state_t *state, qubit_t control, qubit_t target);
extern void ctdg_pure(state_t *state, qubit_t control, qubit_t target);
extern void cp_pure(state_t *state, qubit_t control, qubit_t target, double theta);
extern void cpdg_pure(state_t *state, qubit_t control, qubit_t target, double theta);
extern void swap_pure(state_t *state, qubit_t q1, qubit_t q2);

extern void ccx_pure(state_t *state, qubit_t ctrl1, qubit_t ctrl2, qubit_t target);

/*
 * =====================================================================================================================
 * Packed mixed state gate implementations (gate_packed.c)
 * =====================================================================================================================
 */

extern void x_packed(state_t *state, qubit_t target);
extern void y_packed(state_t *state, qubit_t target);
extern void z_packed(state_t *state, qubit_t target);
extern void h_packed(state_t *state, qubit_t target);
extern void s_packed(state_t *state, qubit_t target);
extern void sdg_packed(state_t *state, qubit_t target);
extern void t_packed(state_t *state, qubit_t target);
extern void tdg_packed(state_t *state, qubit_t target);
extern void hy_packed(state_t *state, qubit_t target);

extern void rx_packed(state_t *state, qubit_t target, double theta);
extern void ry_packed(state_t *state, qubit_t target, double theta);
extern void rz_packed(state_t *state, qubit_t target, double theta);
extern void p_packed(state_t *state, qubit_t target, double theta);

extern void cx_packed(state_t *state, qubit_t control, qubit_t target);
extern void cy_packed(state_t *state, qubit_t control, qubit_t target);
extern void cz_packed(state_t *state, qubit_t control, qubit_t target);
extern void cs_packed(state_t *state, qubit_t control, qubit_t target);
extern void csdg_packed(state_t *state, qubit_t control, qubit_t target);
extern void ch_packed(state_t *state, qubit_t control, qubit_t target);
extern void chy_packed(state_t *state, qubit_t control, qubit_t target);
extern void ct_packed(state_t *state, qubit_t control, qubit_t target);
extern void ctdg_packed(state_t *state, qubit_t control, qubit_t target);
extern void cp_packed(state_t *state, qubit_t control, qubit_t target, double theta);
extern void cpdg_packed(state_t *state, qubit_t control, qubit_t target, double theta);
extern void swap_packed(state_t *state, qubit_t q1, qubit_t q2);

extern void ccx_packed(state_t *state, qubit_t ctrl1, qubit_t ctrl2, qubit_t target);

/*
 * =====================================================================================================================
 * Tiled mixed state gate implementations (gate_tiled.c, gate_tiled_cx.c, gate_tiled_swap.c)
 * =====================================================================================================================
 */

extern void x_tiled(state_t *state, qubit_t target);
extern void y_tiled(state_t *state, qubit_t target);
extern void z_tiled(state_t *state, qubit_t target);
extern void h_tiled(state_t *state, qubit_t target);
extern void s_tiled(state_t *state, qubit_t target);
extern void sdg_tiled(state_t *state, qubit_t target);
extern void t_tiled(state_t *state, qubit_t target);
extern void tdg_tiled(state_t *state, qubit_t target);
extern void hy_tiled(state_t *state, qubit_t target);

extern void rx_tiled(state_t *state, qubit_t target, double theta);
extern void ry_tiled(state_t *state, qubit_t target, double theta);
extern void rz_tiled(state_t *state, qubit_t target, double theta);
extern void p_tiled(state_t *state, qubit_t target, double theta);

extern void cx_tiled(state_t *state, qubit_t control, qubit_t target);
extern void swap_tiled(state_t *state, qubit_t q1, qubit_t q2);

/*
 * =====================================================================================================================
 * Single-qubit gate dispatchers
 * =====================================================================================================================
 */

#define DISPATCH_1Q(name, state, target)                            \
    GATE_VALIDATE((state) && (state)->data,                         \
                  #name ": null state or data pointer");            \
    GATE_VALIDATE((target) < (state)->qubits,                       \
                  #name ": target qubit out of range");             \
    switch ((state)->type) {                                        \
        case PURE:         name##_pure(state, target);   break;     \
        case MIXED_PACKED: name##_packed(state, target);  break;    \
        case MIXED_TILED:  name##_tiled(state, target);   break;    \
        default: GATE_VALIDATE(0, #name ": unknown state type");    \
    }

void x(state_t *state, const qubit_t target) { DISPATCH_1Q(x, state, target); }
void y(state_t *state, const qubit_t target) { DISPATCH_1Q(y, state, target); }
void z(state_t *state, const qubit_t target) { DISPATCH_1Q(z, state, target); }
void h(state_t *state, const qubit_t target) { DISPATCH_1Q(h, state, target); }
void s(state_t *state, const qubit_t target) { DISPATCH_1Q(s, state, target); }
void sdg(state_t *state, const qubit_t target) { DISPATCH_1Q(sdg, state, target); }
void t(state_t *state, const qubit_t target) { DISPATCH_1Q(t, state, target); }
void tdg(state_t *state, const qubit_t target) { DISPATCH_1Q(tdg, state, target); }
void hy(state_t *state, const qubit_t target) { DISPATCH_1Q(hy, state, target); }

/*
 * =====================================================================================================================
 * Rotation gate dispatchers
 * =====================================================================================================================
 */

#define DISPATCH_ROT(name, state, target, theta)                            \
    GATE_VALIDATE((state) && (state)->data,                                 \
                  #name ": null state or data pointer");                    \
    GATE_VALIDATE((target) < (state)->qubits,                               \
                  #name ": target qubit out of range");                     \
    switch ((state)->type) {                                                \
        case PURE:         name##_pure(state, target, theta);   break;      \
        case MIXED_PACKED: name##_packed(state, target, theta);  break;     \
        case MIXED_TILED:  name##_tiled(state, target, theta);   break;     \
        default: GATE_VALIDATE(0, #name ": unknown state type");            \
    }

void rx(state_t *state, const qubit_t target, const double theta) { DISPATCH_ROT(rx, state, target, theta); }
void ry(state_t *state, const qubit_t target, const double theta) { DISPATCH_ROT(ry, state, target, theta); }
void rz(state_t *state, const qubit_t target, const double theta) { DISPATCH_ROT(rz, state, target, theta); }
void p(state_t *state, const qubit_t target, const double theta) { DISPATCH_ROT(p, state, target, theta); }

/*
 * =====================================================================================================================
 * Two-qubit gate dispatchers
 * =====================================================================================================================
 */

#define DISPATCH_2Q(name, state, control, target)                               \
    GATE_VALIDATE((state) && (state)->data,                                     \
                  #name ": null state or data pointer");                        \
    GATE_VALIDATE((control) < (state)->qubits,                                  \
                  #name ": control qubit out of range");                        \
    GATE_VALIDATE((target) < (state)->qubits,                                   \
                  #name ": target qubit out of range");                         \
    GATE_VALIDATE((control) != (target),                                        \
                  #name ": control and target must differ");                    \
    switch ((state)->type) {                                                    \
        case PURE:         name##_pure(state, control, target);   break;        \
        case MIXED_PACKED: name##_packed(state, control, target);  break;       \
        case MIXED_TILED:  name##_tiled(state, control, target);   break;       \
        default: GATE_VALIDATE(0, #name ": unknown state type");                \
    }

#define DISPATCH_2Q_ROT(name, state, control, target, theta)                        \
    GATE_VALIDATE((state) && (state)->data,                                         \
                  #name ": null state or data pointer");                            \
    GATE_VALIDATE((control) < (state)->qubits,                                      \
                  #name ": control qubit out of range");                            \
    GATE_VALIDATE((target) < (state)->qubits,                                       \
                  #name ": target qubit out of range");                             \
    GATE_VALIDATE((control) != (target),                                            \
                  #name ": control and target must differ");                        \
    switch ((state)->type) {                                                        \
        case PURE:         name##_pure(state, control, target, theta);   break;     \
        case MIXED_PACKED: name##_packed(state, control, target, theta);  break;    \
        case MIXED_TILED:  name##_tiled(state, control, target, theta);   break;    \
        default: GATE_VALIDATE(0, #name ": unknown state type");                    \
    }

/* Variants for 2Q gates without tiled implementation */
#define DISPATCH_2Q_NO_TILED(name, state, control, target)                           \
    GATE_VALIDATE((state) && (state)->data,                                          \
                  #name ": null state or data pointer");                             \
    GATE_VALIDATE((control) < (state)->qubits,                                       \
                  #name ": control qubit out of range");                             \
    GATE_VALIDATE((target) < (state)->qubits,                                        \
                  #name ": target qubit out of range");                              \
    GATE_VALIDATE((control) != (target),                                             \
                  #name ": control and target must differ");                         \
    switch ((state)->type) {                                                         \
        case PURE:         name##_pure(state, control, target);   break;             \
        case MIXED_PACKED: name##_packed(state, control, target);  break;            \
        case MIXED_TILED:  GATE_VALIDATE(0, #name ": not implemented for tiled");    \
        default: GATE_VALIDATE(0, #name ": unknown state type");                     \
    }

#define DISPATCH_2Q_ROT_NO_TILED(name, state, control, target, theta)                    \
    GATE_VALIDATE((state) && (state)->data,                                              \
                  #name ": null state or data pointer");                                 \
    GATE_VALIDATE((control) < (state)->qubits,                                           \
                  #name ": control qubit out of range");                                 \
    GATE_VALIDATE((target) < (state)->qubits,                                            \
                  #name ": target qubit out of range");                                  \
    GATE_VALIDATE((control) != (target),                                                 \
                  #name ": control and target must differ");                             \
    switch ((state)->type) {                                                             \
        case PURE:         name##_pure(state, control, target, theta);   break;          \
        case MIXED_PACKED: name##_packed(state, control, target, theta);  break;         \
        case MIXED_TILED:  GATE_VALIDATE(0, #name ": not implemented for tiled");        \
        default: GATE_VALIDATE(0, #name ": unknown state type");                         \
    }

void cx(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q(cx, state, control, target); }
void cy(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(cy, state, control, target); }
void cz(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(cz, state, control, target); }
void cs(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(cs, state, control, target); }
void csdg(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(csdg, state, control, target); }
void ch(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(ch, state, control, target); }
void chy(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(chy, state, control, target); }
void ct(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(ct, state, control, target); }
void ctdg(state_t *state, const qubit_t control, const qubit_t target) { DISPATCH_2Q_NO_TILED(ctdg, state, control, target); }
void cp(state_t *state, const qubit_t control, const qubit_t target, const double theta) { DISPATCH_2Q_ROT_NO_TILED(cp, state, control, target, theta); }
void cpdg(state_t *state, const qubit_t control, const qubit_t target, const double theta) { DISPATCH_2Q_ROT_NO_TILED(cpdg, state, control, target, theta); }

void swap_gate(state_t *state, const qubit_t q1, const qubit_t q2) {
    GATE_VALIDATE(state && state->data, "swap: null state or data pointer");
    GATE_VALIDATE(q1 < state->qubits, "swap: q1 qubit out of range");
    GATE_VALIDATE(q2 < state->qubits, "swap: q2 qubit out of range");
    GATE_VALIDATE(q1 != q2, "swap: q1 and q2 must differ");
    switch (state->type) {
        case PURE:         swap_pure(state, q1, q2);   break;
        case MIXED_PACKED: swap_packed(state, q1, q2);  break;
        case MIXED_TILED:  swap_tiled(state, q1, q2);   break;
        default: GATE_VALIDATE(0, "swap: unknown state type");
    }
}

/*
 * =====================================================================================================================
 * Three-qubit gate dispatchers
 * =====================================================================================================================
 */

void ccx(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    GATE_VALIDATE(state && state->data, "ccx: null state or data pointer");
    GATE_VALIDATE(ctrl1 < state->qubits, "ccx: ctrl1 qubit out of range");
    GATE_VALIDATE(ctrl2 < state->qubits, "ccx: ctrl2 qubit out of range");
    GATE_VALIDATE(target < state->qubits, "ccx: target qubit out of range");
    GATE_VALIDATE(ctrl1 != ctrl2, "ccx: ctrl1 and ctrl2 must differ");
    GATE_VALIDATE(ctrl1 != target, "ccx: ctrl1 and target must differ");
    GATE_VALIDATE(ctrl2 != target, "ccx: ctrl2 and target must differ");
    switch (state->type) {
        case PURE:         ccx_pure(state, ctrl1, ctrl2, target);   break;
        case MIXED_PACKED: ccx_packed(state, ctrl1, ctrl2, target);  break;
        case MIXED_TILED:  GATE_VALIDATE(0, "ccx: not implemented for tiled");
        default: GATE_VALIDATE(0, "ccx: unknown state type");
    }
}
