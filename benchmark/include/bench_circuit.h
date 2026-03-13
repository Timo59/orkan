#ifndef BENCH_CIRCUIT_H
#define BENCH_CIRCUIT_H

/*
 * =====================================================================================================================
 * Circuit-level benchmark infrastructure
 *
 * Defines the circuit abstraction (gate enum, operation struct, circuit struct),
 * circuit constructors (QAOA, VQE-HEA, QV), and runner declarations for each backend.
 * This header is the shared contract for all bench_circuit_*.c/.cpp files.
 * =====================================================================================================================
 */

#include "bench.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Gate type enum
 * =====================================================================================================================
 */

/**
 * @brief Gate type identifiers used in circuit_op_t.
 *
 * Only gates actually used by the three circuit constructors are defined.
 * Dispatch switches must handle all listed cases; add a default branch that
 * prints to stderr in debug builds for unimplemented gate types.
 */
typedef enum {
    GATE_H  = 0,
    GATE_X  = 1,
    GATE_Z  = 2,
    GATE_RX = 3,
    GATE_RY = 4,
    GATE_RZ = 5,
    GATE_CX = 6,
} circuit_gate_t;

/*
 * =====================================================================================================================
 * Circuit operation and circuit structs
 * =====================================================================================================================
 */

/**
 * @brief A single gate operation in a circuit.
 *
 * One-qubit gates use q0; two-qubit gates use q0 (control/first) and q1 (target/second).
 * angle is used by rotation gates; ignored for non-parametric gates.
 */
typedef struct {
    circuit_gate_t gate_type;
    qubit_t        q0;
    qubit_t        q1;      /**< Used only by 2Q gates; 0 otherwise */
    double         angle;   /**< Used only by rotation gates; 0.0 otherwise */
} circuit_op_t;

/**
 * @brief A flat, dynamically-allocated gate sequence.
 */
typedef struct {
    circuit_op_t *ops;       /**< Heap-allocated array of gate operations */
    int           n_ops;     /**< Number of gates currently in the circuit */
    int           cap;       /**< Allocated capacity (internal; do not read) */
    int           n_qubits;  /**< Width of the circuit */
    const char   *name;      /**< Human-readable circuit name (static string) */
} circuit_t;

/*
 * =====================================================================================================================
 * Runner function pointer typedef
 * =====================================================================================================================
 */

/**
 * @brief Function pointer type for circuit runners used by calibrate_iterations.
 *
 * @param c          Circuit to execute.
 * @param ctx        Backend-specific context (cast to the appropriate *_circuit_ctx_t).
 * @param iterations Number of times to execute the full circuit.
 */
typedef void (*runner_fn)(const circuit_t *c, void *ctx, int iterations);

/*
 * =====================================================================================================================
 * Context structs for each runner
 * =====================================================================================================================
 */

/** @brief Context for qlib packed/tiled circuit runner */
typedef struct {
    state_t      *s;
    state_type_t  type;
} qlib_circuit_ctx_t;

#ifdef WITH_QUEST
/** @brief Context for QuEST circuit runner */
typedef struct {
    void *qureg;    /**< Opaque; actually a Qureg, but avoid pulling in QuEST headers */
} quest_circuit_ctx_t;
#endif

#ifdef WITH_QULACS
/** @brief Context for Qulacs circuit runner */
typedef struct {
    void  *rho;     /**< CTYPE* density matrix, opaque here to avoid Qulacs headers */
    size_t dim;     /**< ITYPE dimension = 2^n_qubits */
} qulacs_circuit_ctx_t;
#endif

#ifdef WITH_AER_DM
/** @brief Context for Aer-DM circuit runner */
typedef struct {
    void *dm;       /**< Opaque; AER::QV::DensityMatrix<double>* */
    int   n_qubits;
} aer_dm_circuit_ctx_t;
#endif

/*
 * =====================================================================================================================
 * Circuit lifecycle and append helpers
 * =====================================================================================================================
 */

/**
 * @brief Allocate an empty circuit.
 *
 * @param n_qubits    Width of the circuit.
 * @param initial_cap Initial capacity for the ops array. 0 is valid (ops will be NULL).
 * @param name        Human-readable name (pointer stored, not copied; must outlive circuit).
 * @return            Initialized circuit_t. On allocation failure, ops is NULL and cap is 0.
 */
circuit_t circuit_alloc(int n_qubits, int initial_cap, const char *name);

/**
 * @brief Free a circuit's ops array and zero out the struct.
 */
void circuit_free(circuit_t *c);

/**
 * @brief Append a single-qubit non-parametric gate.
 * @return 0 on success, -1 on allocation failure.
 */
int circuit_append_1q(circuit_t *c, circuit_gate_t g, qubit_t q0);

/**
 * @brief Append a single-qubit rotation gate with an angle parameter.
 * @return 0 on success, -1 on allocation failure.
 */
int circuit_append_1q_angle(circuit_t *c, circuit_gate_t g, qubit_t q0, double angle);

/**
 * @brief Append a two-qubit gate.
 * @return 0 on success, -1 on allocation failure.
 */
int circuit_append_2q(circuit_t *c, circuit_gate_t g, qubit_t q0, qubit_t q1);

/**
 * @brief Append an RZZ composite: CX(q0,q1) . Rz(q1, gamma) . CX(q0,q1).
 * @return 0 on success, -1 on allocation failure.
 */
int circuit_rzz(circuit_t *c, qubit_t q0, qubit_t q1, double gamma);

/*
 * =====================================================================================================================
 * Circuit constructors
 * =====================================================================================================================
 */

/**
 * @brief Construct a QAOA p=1 circuit on a cycle graph.
 *
 * n H gates + n RZZ(gamma) decompositions (3n gates) + n RX(beta) = 5n total gates.
 *
 * @param n     Number of qubits (circuit width).
 * @param gamma Cost Hamiltonian parameter.
 * @param beta  Mixer parameter.
 * @return      Populated circuit. On failure, n_ops == 0.
 */
circuit_t circuit_qaoa_p1(int n, double gamma, double beta);

/**
 * @brief Construct a hardware-efficient VQE ansatz.
 *
 * depth layers of: Ry(pi/4) on all qubits + CX ladder (0->1, 1->2, ..., n-2->n-1).
 * Total gate count: depth * (2n - 1).
 *
 * @param n     Number of qubits.
 * @param depth Number of layers. If < 0, defaults to floor(n/2).
 * @return      Populated circuit. On failure, n_ops == 0.
 */
circuit_t circuit_vqe_hea(int n, int depth);

/**
 * @brief Construct a Quantum Volume circuit.
 *
 * n layers of random permutation + parameterized two-qubit blocks. n must be even.
 * Approximate gate count: 15 * n^2 / 2.
 *
 * @param n    Number of qubits (must be even).
 * @param seed RNG seed for reproducible circuit generation.
 * @return     Populated circuit. On failure or odd n, n_ops == 0.
 */
circuit_t circuit_qv(int n, int seed);

/*
 * =====================================================================================================================
 * Runner functions — circuit-level benchmarks
 * =====================================================================================================================
 */

/**
 * @brief Run circuit benchmark using qlib packed storage.
 */
bench_result_t bench_circuit_qlib_packed(const circuit_t *c, int iterations,
                                         int warmup, int runs);

/**
 * @brief Run circuit benchmark using qlib tiled storage.
 */
bench_result_t bench_circuit_qlib_tiled(const circuit_t *c, int iterations,
                                        int warmup, int runs);

#ifdef WITH_QUEST
/**
 * @brief Run circuit benchmark using QuEST density matrix.
 */
bench_result_t bench_circuit_quest(const circuit_t *c, int iterations,
                                   int warmup, int runs);
#endif

#ifdef WITH_QULACS
/**
 * @brief Run circuit benchmark using Qulacs density matrix.
 */
bench_result_t bench_circuit_qulacs(const circuit_t *c, int iterations,
                                    int warmup, int runs);
#endif

#ifdef WITH_AER_DM
/**
 * @brief Run circuit benchmark using Qiskit-Aer density matrix.
 */
bench_result_t bench_circuit_aer_dm(const circuit_t *c, int iterations,
                                    int warmup, int runs);
#endif

/*
 * =====================================================================================================================
 * Context alloc/free for calibrate_iterations compatibility
 * =====================================================================================================================
 */

/** @brief Allocate a qlib circuit context for the given storage type. */
qlib_circuit_ctx_t *qlib_circuit_ctx_alloc(int n_qubits, state_type_t type);

/** @brief Free a qlib circuit context. */
void qlib_circuit_ctx_free(qlib_circuit_ctx_t *ctx);

#ifdef WITH_QUEST
/** @brief Allocate a QuEST circuit context. Requires bench_quest_init() first. */
quest_circuit_ctx_t *quest_circuit_ctx_alloc(int n_qubits);

/** @brief Free a QuEST circuit context. */
void quest_circuit_ctx_free(quest_circuit_ctx_t *ctx);
#endif

#ifdef WITH_QULACS
/** @brief Allocate a Qulacs circuit context. */
qulacs_circuit_ctx_t *qulacs_circuit_ctx_alloc(int n_qubits);

/** @brief Free a Qulacs circuit context. */
void qulacs_circuit_ctx_free(qulacs_circuit_ctx_t *ctx);
#endif

#ifdef WITH_AER_DM
/** @brief Allocate an Aer-DM circuit context. */
aer_dm_circuit_ctx_t *aer_dm_circuit_ctx_alloc(int n_qubits);

/** @brief Free an Aer-DM circuit context. */
void aer_dm_circuit_ctx_free(aer_dm_circuit_ctx_t *ctx);
#endif

#ifdef __cplusplus
}
#endif

#endif /* BENCH_CIRCUIT_H */
