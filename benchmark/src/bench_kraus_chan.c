/**
 * @file bench_kraus_chan.c
 * @brief Channel definitions: analytic Kraus operators for named channels.
 */

#include "bench_kraus.h"
#include <math.h>
#include <string.h>
#include <strings.h>

/* =====================================================================
 * Depolarizing channel (1 qubit, 4 Kraus operators)
 *
 *   rho -> (1-p) rho + (p/3)(X rho X + Y rho Y + Z rho Z)
 *
 * Kraus operators: { sqrt(1-p) I, sqrt(p/3) X, sqrt(p/3) Y, sqrt(p/3) Z }
 * ===================================================================== */

static void build_depolarize(double p, cplx_t *out) {
    double a = sqrt(1.0 - p);
    double b = sqrt(p / 3.0);

    /* K0 = sqrt(1-p) * I */
    out[0] = a;  out[1] = 0;
    out[2] = 0;  out[3] = a;

    /* K1 = sqrt(p/3) * X */
    out[4] = 0;  out[5] = b;
    out[6] = b;  out[7] = 0;

    /* K2 = sqrt(p/3) * Y */
    out[8]  = 0;          out[9]  = -b * I;
    out[10] = b * I;      out[11] = 0;

    /* K3 = sqrt(p/3) * Z */
    out[12] = b;   out[13] = 0;
    out[14] = 0;   out[15] = -b;
}

/* =====================================================================
 * Channel table
 * ===================================================================== */

const bench_channel_info_t bench_channel_table[BC_COUNT] = {
    [BC_DEPOLARIZE] = {
        .id          = BC_DEPOLARIZE,
        .name        = "Depolarize",
        .n_ops       = 4,
        .tgt_qubits  = 1,
        .parametrised = 1,
        .build       = build_depolarize
    },
};

bench_channel_id_t bench_resolve_channel(const char *name) {
    for (int i = 0; i < BC_COUNT; ++i)
        if (strcasecmp(name, bench_channel_table[i].name) == 0)
            return bench_channel_table[i].id;
    return BC_COUNT;
}
