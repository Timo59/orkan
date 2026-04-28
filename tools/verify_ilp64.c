/**
 * verify_ilp64.c - Verify that BLAS library has 64-bit integer support
 *
 * Compile and run this to verify your BLAS/LAPACK has ILP64 support:
 *
 * On macOS:
 *   gcc -o verify_ilp64 verify_ilp64.c -framework Accelerate -DACCELERATE_NEW_LAPACK -DACCELERATE_LAPACK_ILP64
 *   ./verify_ilp64
 *
 * On Linux (system OpenBLAS):
 *   gcc -o verify_ilp64 verify_ilp64.c -lopenblas -DOPENBLAS_USE64BITINT -I/usr/include
 *   ./verify_ilp64
 *
 * Expected output for ILP64: "✓ ILP64 support detected (sizeof(blasint) = 8)"
 */

#include <stdio.h>
#include <stdint.h>

#ifdef __APPLE__
    #include <Accelerate/Accelerate.h>
    // On Apple with ACCELERATE_LAPACK_ILP64, __LAPACK_int is defined as int64_t
    typedef __LAPACK_int blasint;
#else
    #include <cblas.h>
    // OpenBLAS defines blasint when built with INTERFACE64=1
    #ifndef blasint
        // If blasint is not defined, assume standard int (LP64)
        typedef int blasint;
    #endif
#endif

int main(void) {
    printf("BLAS Integer Type Verification\n");
    printf("================================\n\n");

    // Check the size of the BLAS integer type
    size_t int_size = sizeof(blasint);
    printf("sizeof(blasint) = %zu bytes\n", int_size);

    // Also check standard types for reference
    printf("sizeof(int)     = %zu bytes\n", sizeof(int));
    printf("sizeof(long)    = %zu bytes\n", sizeof(long));
    printf("sizeof(int64_t) = %zu bytes\n\n", sizeof(int64_t));

    // Verify ILP64 support
    if (int_size == 8) {
        printf("✓ ILP64 support detected\n");
        printf("  This BLAS library uses 64-bit integers (ILP64)\n");
        printf("  Compatible with Orkan requirements\n");
        return 0;
    } else if (int_size == 4) {
        printf("✗ LP64 detected (32-bit integers)\n");
        printf("  This BLAS library uses 32-bit integers (LP64)\n");
        printf("  NOT compatible with Orkan requirements\n\n");
        printf("To fix this:\n");
#ifdef __APPLE__
        printf("  - Ensure you compile with -DACCELERATE_LAPACK_ILP64\n");
#else
        printf("  1. Rebuild OpenBLAS with: make INTERFACE64=1\n");
        printf("  2. Or use bundled OpenBLAS: cmake -DUSE_SYSTEM_OPENBLAS=OFF ..\n");
#endif
        return 1;
    } else {
        printf("? Unexpected integer size: %zu bytes\n", int_size);
        printf("  Unable to determine BLAS integer type\n");
        return 2;
    }
}
