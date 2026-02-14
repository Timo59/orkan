/**
 * @file gate_tiled.c
 * @brief Gate registry for tiled density matrix implementations
 *
 * All single-qubit gate implementations have been extracted to standalone files:
 *
 *   Pauli gates:    gate_tiled_x.c, gate_tiled_y.c, gate_tiled_z.c
 *   Clifford gates: gate_tiled_h.c, gate_tiled_hy.c,
 *                   gate_tiled_s.c, gate_tiled_sdg.c,
 *                   gate_tiled_t.c, gate_tiled_tdg.c
 *   Rotation gates: gate_tiled_rx.c, gate_tiled_ry.c,
 *                   gate_tiled_rz.c, gate_tiled_p.c
 *   Two-qubit gates: gate_tiled_cx.c, gate_tiled_swap.c
 *
 * Shared helpers (tile_off, elem_off, dtile_read/write, OMP_THRESHOLD)
 * are in gate_tiled.h.
 */
