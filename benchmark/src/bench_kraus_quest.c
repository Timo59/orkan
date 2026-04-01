/**
 * @file bench_kraus_quest.c
 * @brief QuEST backend for the Kraus benchmark.
 *
 * Uses createKrausMap / setKrausMap / mixKrausMap.
 * QuEST internally converts Kraus operators to a superoperator
 * at setKrausMap time; mixKrausMap applies it.
 */

#ifdef HAVE_QUEST

#include "quest/include/quest.h"
#include "bench_kraus.h"

_Static_assert(sizeof(qcomp) == sizeof(cplx_t),
    "QuEST qcomp size must match cplx_t (double _Complex)");

/* =====================================================================
 * Context
 * ===================================================================== */

typedef struct {
    Qureg    qureg;
    KrausMap map;
    int      tgt_qubits;
} kraus_quest_ctx_t;

/* =====================================================================
 * Apply: mixKrausMap at given position
 * ===================================================================== */

static void kraus_quest_apply(void *c, const qubit_t *pos) {
    kraus_quest_ctx_t *ctx = (kraus_quest_ctx_t *)c;
    int tgts[3];
    for (int i = 0; i < ctx->tgt_qubits; ++i)
        tgts[i] = (int)pos[i];
    mixKrausMap(ctx->qureg, tgts, ctx->tgt_qubits, ctx->map);
}

/* =====================================================================
 * Init / cleanup / vtable
 * ===================================================================== */

static int quest_env_initialised = 0;

static void quest_ensure_env(void) {
    if (!quest_env_initialised) {
        initQuESTEnv();
        quest_env_initialised = 1;
    }
}

static void *kraus_quest_init(qubit_t sys_qubits, int n_ops,
                              int tgt_qubits, const cplx_t *ops) {
    quest_ensure_env();

    idx_t d = (idx_t)1 << tgt_qubits;

    kraus_quest_ctx_t *ctx = malloc(sizeof(kraus_quest_ctx_t));
    if (!ctx) return NULL;
    ctx->tgt_qubits = tgt_qubits;

    ctx->qureg = createDensityQureg((int)sys_qubits);
    bench_init_hermitian_dense((cplx_t *)ctx->qureg.cpuAmps,
                                (idx_t)1 << sys_qubits);

    /* Build KrausMap and populate from our flat row-major operators */
    ctx->map = createKrausMap(tgt_qubits, n_ops);

    /* Allocate qcomp*** pointer structure for setKrausMap */
    qcomp ***matrices = malloc((size_t)n_ops * sizeof(qcomp **));
    for (int k = 0; k < n_ops; ++k) {
        matrices[k] = malloc((size_t)d * sizeof(qcomp *));
        for (idx_t i = 0; i < d; ++i) {
            matrices[k][i] = malloc((size_t)d * sizeof(qcomp));
            for (idx_t j = 0; j < d; ++j)
                matrices[k][i][j] = (qcomp)ops[k * d * d + i * d + j];
        }
    }
    setKrausMap(ctx->map, matrices);

    /* Free temporary pointer structure */
    for (int k = 0; k < n_ops; ++k) {
        for (idx_t i = 0; i < d; ++i)
            free(matrices[k][i]);
        free(matrices[k]);
    }
    free(matrices);

    return ctx;
}

static void kraus_quest_cleanup(void *c) {
    kraus_quest_ctx_t *ctx = (kraus_quest_ctx_t *)c;
    destroyKrausMap(ctx->map);
    destroyQureg(ctx->qureg);
    free(ctx);
}

static size_t kraus_quest_touch_bytes(qubit_t qubits) {
    size_t dim = (size_t)1 << qubits;
    return 2 * dim * dim * sizeof(cplx_t);
}

static int kraus_quest_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

const bench_kraus_backend_t kraus_quest_backend = {
    .name        = "QuEST",
    .init        = kraus_quest_init,
    .apply       = kraus_quest_apply,
    .cleanup     = kraus_quest_cleanup,
    .touch_bytes = kraus_quest_touch_bytes,
    .available   = kraus_quest_available
};

#endif /* HAVE_QUEST */
