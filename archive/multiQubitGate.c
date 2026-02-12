// multiQubitGate.c - Alternative reference implementation for mixed-state gate testing.
//                    Applies a k-local operator on a full hermitian operator via insertBits enumeration.
//                    Archived: not required by the unit test suite.

#include <complex.h>

typedef double complex cplx_t;

/*
 * @brief   Returns the integer obtained by inserting the specified bits at the corresponding positions of the input.
 *
 * @param   i       Input integer
 * @param   n       Number of bits to insert
 * @param   bits    Integer whose rightmost n bits define the bits to insert
 * @param   pos     Ordered (smallest to largest) list of positions at which the bits are inserted (little endian)
 */
static unsigned insertBits(unsigned i, const unsigned n, unsigned bits, const unsigned *pos) {
    // Iterate bits to be inserted
    for (unsigned j = 0; j < n; ++j) {
        // Mask of bits in i right to the current position
        const unsigned right = ((1u << pos[j]) - 1) & i;

        // Copy value of the j-th bit to rightmost bit of a bit set mask
        unsigned set = (bits >> j) & 1;

        // Set all bits right to current position to zero and the bit at position to target bit
        unsigned left = (i >> pos[j]) << 1;
        left |= set;
        left <<= pos[j];

        i = left | right;
    }

    return i;
}

/*
 * @brief   Applies multi-qubit operator on a full hermitian operator in place
 *
 * @param   n   Number of qubits
 * @param   h   Array of size 2**n * 2**n; concatenation of hermitain operator's columns
 * @param   k   Number of qubits the gate acts on non-trivially
 * @param   a   Array of size k; sorted list of indices of qubits acted upon non-trivially (little endian)
 * @param   m   Array of size 2**k * 2**k; of k-local operator in column-major format
 */
void multiQubitGate(const unsigned n, cplx_t *h, const unsigned k, const unsigned *a, const cplx_t *m) {
    const unsigned N = 1 << n;  // Order of the hermitian matrix
    const unsigned K = 1 << k;  // Invariant subspace dimension; number of rows and columns intertwined by the gate
    const unsigned mat_size = K * K;    // Number of entries in each invariant block matrix
    const unsigned M = 1 << (n - k);    // Number of invariant subspaces

    // Iterate the invariant column blocks
    for (unsigned col_block = 0; col_block < M; ++col_block) {

        // Iterate the invariant row blocks
        for (unsigned row_block = 0; row_block < M; ++row_block) {
            cplx_t tmp[mat_size];

            // Fetch elements from invariant matrix blocks to the stack
            for (unsigned j = 0; j < K; ++j) {
                unsigned col = insertBits(col_block, k, j, a);

                for (unsigned i = 0; i < K; ++i) {
                    unsigned row = insertBits(row_block, k, i, a);
                    tmp[i + j * K] = h[row + col * N];
                    h[row + col * N] = 0;
                }
            }

            // Add elements from of conjugation M->M*tmp*M**H of stack matrix to invariant blocks
            for (unsigned j = 0; j < K; ++j) {
                unsigned col = insertBits(col_block, k, j, a);

                for (unsigned i = 0; i < K; ++i) {
                    unsigned row = insertBits(row_block, k, i, a);

                    // H_ij = M_ik tmp_kl (M**H)_lj
                    for (unsigned s = 0; s < K; ++s) {
                        for (unsigned t = 0; t < K; ++t) {
                            h[row + col * N] += m[i + s * K] * tmp[s + t * K] * conj(m[j + t * K]);
                        }
                    }
                }
            }
        }
    }
}
