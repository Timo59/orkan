/**
 * @file bench_pgfplots.c
 * @brief PGFplots output formatter for the gate benchmark (stub — Phase 3).
 */

#include "bench.h"

void bench_print_pgfplots_agg(const bench_gate_info_t *gate,
                              const bench_options_t *opts,
                              const qubit_t *qubit_list, int n_qubits,
                              const bench_result_t (*results)[BACKEND_COUNT]) {
    (void)gate; (void)opts; (void)qubit_list; (void)n_qubits; (void)results;
}

void bench_print_pgfplots_perq(const bench_gate_info_t *gate,
                               const bench_options_t *opts,
                               qubit_t qubits,
                               const bench_pos_result_t *per_backend[BACKEND_COUNT],
                               int n_pos, int gate_qubits) {
    (void)gate; (void)opts; (void)qubits; (void)per_backend;
    (void)n_pos; (void)gate_qubits;
}
