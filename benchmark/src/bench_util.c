/**
 * @file bench_util.c
 * @brief Non-inline benchmark utilities: PRNG fill, Hermitian init, huge-page allocator.
 */

#include "bench.h"
#include <sys/mman.h>
#if !defined(__APPLE__) && !defined(__FreeBSD__)
#include <sys/random.h>
#endif

/* =====================================================================
 * Random fill — xoshiro256** seeded from OS entropy
 * ===================================================================== */

void bench_fill_random(void *buf, size_t n) {
    uint64_t seed;
#if defined(__APPLE__) || defined(__FreeBSD__)
    arc4random_buf(&seed, sizeof(seed));
#else
    if (getrandom(&seed, sizeof(seed), 0) != (ssize_t)sizeof(seed)) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        seed = (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
    }
#endif

    /* Expand seed via splitmix64 into xoshiro256** state */
    uint64_t z;
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s0 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s1 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s2 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s3 = z ^ (z >> 31);

    uint8_t *p = (uint8_t *)buf;
    size_t remaining = n;
    while (remaining >= 8) {
        uint64_t r = s1 * 5;
        r = ((r << 7) | (r >> 57)) * 9;
        memcpy(p, &r, 8);
        p += 8;
        remaining -= 8;
        uint64_t t = s1 << 17;
        s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t;
        s3 = (s3 << 45) | (s3 >> 19);
    }
    if (remaining > 0) {
        uint64_t r = s1 * 5;
        r = ((r << 7) | (r >> 57)) * 9;
        memcpy(p, &r, remaining);
    }
}

/* =====================================================================
 * Hermitian random initialisation
 * ===================================================================== */

void bench_init_hermitian_dense(cplx_t *rho, dim_t dim) {
    /* Fill entire matrix with random data */
    bench_fill_random(rho, (size_t)dim * (size_t)dim * sizeof(cplx_t));

    /* Enforce Hermiticity: lower triangle stays, upper = conj(lower), diag real */
    for (dim_t i = 0; i < dim; ++i) {
        SETIMAG(rho[i * dim + i], 0.0);
        for (dim_t j = i + 1; j < dim; ++j)
            rho[i * dim + j] = conj(rho[j * dim + i]);
    }
}

void bench_init_hermitian_packed(cplx_t *data, qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    size_t len = (size_t)dim * ((size_t)dim + 1) / 2;

    /* Fill packed lower-triangle storage with random data */
    bench_fill_random(data, len * sizeof(cplx_t));

    /* Zero imaginary part of diagonal elements.
     * In LAPACK column-major packed lower triangle, diagonal element (i,i)
     * is at index: i*(2*dim - i + 1)/2  (the first element of column i). */
    for (dim_t i = 0; i < dim; ++i) {
        size_t idx = (size_t)i * (2 * (size_t)dim - (size_t)i + 1) / 2;
        SETIMAG(data[idx], 0.0);
    }
}

void bench_init_hermitian_tiled(cplx_t *data, qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    dim_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
    size_t n_tile_pairs = (size_t)(n_tiles * (n_tiles + 1) / 2);

    /* Fill all tile storage with random data */
    bench_fill_random(data, n_tile_pairs * TILE_SIZE * sizeof(cplx_t));

    /* Enforce Hermiticity within each diagonal tile.
     * Diagonal tile k (tile-row == tile-col == k) starts at offset
     * k*(k+1)/2 * TILE_SIZE.  Within the tile, elements are row-major:
     * element (lr, lc) is at tile_base + lr * TILE_DIM + lc.
     *
     * The tile represents a TILE_DIM × TILE_DIM submatrix on the main
     * diagonal of the full matrix, so it must itself be Hermitian:
     *   - (lr, lc) with lr > lc: keep random value (lower triangle)
     *   - (lr, lc) with lr < lc: set to conj of (lc, lr)
     *   - (lr, lr): zero imaginary part */
    for (dim_t k = 0; k < n_tiles; ++k) {
        size_t tile_off = (size_t)(k * (k + 1) / 2) * TILE_SIZE;
        cplx_t *tile = data + tile_off;

        for (int lr = 0; lr < TILE_DIM; ++lr) {
            SETIMAG(tile[lr * TILE_DIM + lr], 0.0);
            for (int lc = lr + 1; lc < TILE_DIM; ++lc)
                tile[lr * TILE_DIM + lc] = conj(tile[lc * TILE_DIM + lr]);
        }
    }
}

/* =====================================================================
 * Huge-page memory utilities
 * ===================================================================== */

void *bench_alloc_huge(size_t n) {
#if defined(__linux__)
    if (n >= BENCH_HUGEPAGE_THRESHOLD) {
        void *p = mmap(NULL, n, PROT_READ | PROT_WRITE,
                       MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
        if (p == MAP_FAILED) {
            fprintf(stderr, "bench: mmap failed for %zu bytes\n", n);
            abort();
        }
        madvise(p, n, MADV_HUGEPAGE);
        return p;
    }
    return malloc(n);
#elif defined(__APPLE__)
    if (n >= BENCH_HUGEPAGE_THRESHOLD) {
        void *p = mmap(NULL, n, PROT_READ | PROT_WRITE,
                       MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
        if (p == MAP_FAILED) {
            fprintf(stderr, "bench: mmap failed for %zu bytes\n", n);
            abort();
        }
        return p;
    }
    return malloc(n);
#else
    return malloc(n);
#endif
}

void bench_free_huge(void *p, size_t n) {
#if defined(__linux__) || defined(__APPLE__)
    if (n >= BENCH_HUGEPAGE_THRESHOLD && p != NULL) {
        munmap(p, n);
        return;
    }
#else
    (void)n;
#endif
    free(p);
}
