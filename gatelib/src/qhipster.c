/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "gatelib.h"

/*
 * =================================================================================================
 *                                              Pauli gates
 * =================================================================================================
 */

/*
This function executes the Pauli-X gate on a specified qubit whose state along with the other
qubits' state determine the overall register's pure state.

Input:
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unisgned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    qubit_t qubit: an unsigned integer holding the qubit's position from the right in a bit string
    the Pauli-X gate acts on

Output:
    There is no return value, but the array's elements the pointer in state refers to are exchanged
    depending on the qubit the Pauli-X gate acts on. The function swaps each entry of the array with
    the one whose index in bit representation only differs at the position determined by qubit.
*/

void applyX(state_t* state, qubit_t qubit) {
    dim_t i;                                            // index that determines the starting point
                                                        // of the inner loop changing only those
                                                        // bits in the index's binary representation
                                                        // left to the one specified by qubit

    dim_t j;                                            // index that runs within the inner loop
                                                        // changing only those bits in the index's
                                                        // binary representationindex right to the
                                                        // the one specified by qubit

    dim_t blockDistance = POW2(qubit + 1, dim_t);       // distance in decimal representation of
                                                        // blocks within which the swap of two
                                                        // amplitudes takes place, i.e. encompasses
                                                        // all qubits to the right of the specified
                                                        // qubit and the qubit itself

    dim_t flipDistance = POW2(qubit, dim_t);            // distance in decimal representation of 
                                                        // indices differing only in the value of
                                                        // the specified qubit

    /*
    Starting in the first block, i.e., the block where all bits to the left of qubit in binary
    representation are zero, the integer i for counting the starting index for the block is set to
    zero. After each block, consisting of 2^{qubit - 1} swaps, is processed, it proceeds with
    increasing i by 2^{qubit + 1}. This corresponds to gradually increasing the bits in the index's
    binary representation left to the one specified by qubit while reseting all bits to its right
    to zero.
    */
    for (i = 0; i < state->dimension; i += blockDistance) {
        /*
        The inner loop starts with all bits in the index's binary representation right to the one
        specified by qubit initialized with zero, as handed over by the outer loop. Thereafter j is 
        increased by one in each iteration until all 2^0 + 2^1 + ... + 2^{qubit-1} = 2^qubit - 1
        entries are swapped.
        */
        for (j = i; j < i + flipDistance; ++j) {
            /*
            At each iterational step the value of the j-th entry originating from the pointer stored
            in state is stored to a temporary double complex variable. Subesequently, one sets the
            value the j-th entry's pointer holds to be the value of the entry whose index's binary
            representation only differs in qubit. Finally, the latter pointer is set to the value
            stored under the temporary and j.
            */
            SWAP(state->vector + j, state->vector + (j + flipDistance), cplx_t);
        }
    }
}

void applyY(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            SWAP(state->vector + j, state->vector + (j + flipDistance), cplx_t);
            state->vector[j] *= -I;
            state->vector[j + flipDistance] *= I;
        }
    }
}

void applyZ(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j + flipDistance] *= -1;
        }
    }
}

/*
 * =================================================================================================
 *                                              Clifford gates
 * =================================================================================================
 */

void applyS(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j + flipDistance] *= I;
        }
    }
}

void applySdagger(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j + flipDistance] *= -I;
        }
    }
}

void applyH(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vector[j];
            state->vector[j] = INVSQRT2 * (tmp + state->vector[j + flipDistance]);
            state->vector[j + flipDistance] = INVSQRT2 * (tmp - state->vector[j + flipDistance]);
        }
    }
}

/*
 * =================================================================================================
 *                                              Hadamard-Y gate
 * =================================================================================================
 */

void applyHy(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vector[j];
            state->vector[j] = INVSQRT2 * (tmp - I * state->vector[j + flipDistance]);
            state->vector[j + flipDistance] = INVSQRT2 * (I * tmp - state->vector[j + flipDistance]);
        }
    }
}

/*
 * =================================================================================================
 *                                              T gate
 * =================================================================================================
 */

void applyT(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j + flipDistance] *= (INVSQRT2 + INVSQRT2 * I);
        }
    }
}

void applyTdagger(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j + flipDistance] *= (INVSQRT2 - INVSQRT2 * I);
        }
    }
}

/*
 * =================================================================================================
 *                                              P gate
 * =================================================================================================
 */

void applyP(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j + flipDistance] *= cos(angle) + sin(angle) * I;
        }
    }
}

void applyPdagger(state_t* state, qubit_t qubit, double angle) {
    applyP(state, qubit, -angle);
}

/*
 * =================================================================================================
 *                                              rotation gates
 * =================================================================================================
 */

void applyRX(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    angle /= 2.;
    cplx_t tmp;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vector[j];
            state->vector[j] = cos(angle) * tmp - I * sin(angle) * state->vector[j + flipDistance];
            state->vector[j + flipDistance] = -I * sin(angle) * tmp + cos(angle) \
                                             * state->vector[j + flipDistance];
        }
    }
}

void applyRXdagger(state_t* state, qubit_t qubit, double angle) {
    applyRX(state, qubit, -angle);
}

void applyRY(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    angle /= 2.;
    cplx_t tmp;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vector[j];
            state->vector[j] = cos(angle) * tmp - sin(angle) * state->vector[j + flipDistance];
            state->vector[j + flipDistance] = sin(angle) * tmp + cos(angle) \
                                             * state->vector[j + flipDistance];
        }
    }
}

void applyRYdagger(state_t* state, qubit_t qubit, double angle) {
    applyRY(state, qubit, -angle);
}

void applyRZ(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    angle /= 2.;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vector[j] *= cos(angle) - I * sin(angle);
            state->vector[j + flipDistance] *= cos(angle) + I * sin(angle);
        }
    }
}

void applyRZdagger(state_t* state, qubit_t qubit, double angle) {
    applyRZ(state, qubit, -angle);
}

/*
 * =================================================================================================
 *                                              controlled Pauli gates
 * =================================================================================================
 */

/*
This function executes a Pauli-X gate on the second qubit specified in the input based on whether 
the first qubit specified in the input is in the |1> state and changes the corresponding entries in
the pure state's vector representation.

Input:
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unisgned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    qubit_t control: an unsigned integer holding the qubit's position from the right in a bit string
    that controls the execution of the Pauli-X gate on the target qubit
    qubit_t target: an unsigned integer holding the qubit's position from the right in a bit string
    the Pauli-X gate acts on if the control qubit is in the |1> state
Output:
    There is no output value, but the array's elements the pointer in state refers to are
    exchanged with those whose index in bit representation differs only at the position specified by
    target if the bit at the position specified by control is 1.
*/
void applyCX(state_t* state, qubit_t control, qubit_t target) {
    dim_t i;                                            // index that determines the starting point
                                                        // of the inner loop changing only those
                                                        // bits in the index's binary representation
                                                        // left to the one specified by target

    dim_t j;                                            // index that runs within the inner loop
                                                        // changing only those bits in the index's
                                                        // binary representationindex right to the
                                                        // the one specified by target

    // qubit_t left_qubit = MAX(control, target);           // index of the left qubit

    // qubit_t right_qubit = MIN(control, target);          // index of the right qubit

    dim_t blockDistance = POW2(target + 1, dim_t);
                                                        // distance in decimal representation of
                                                        // blocks within which the swap of two
                                                        // amplitudes takes place, i.e. encompasses
                                                        // all qubits to the right of the specified
                                                        // target and the target itself

    // dim_t minorBlockDistance = POW2(right_qubit + 1, dim_t);

    dim_t flipDistance = POW2(target, dim_t);           // distance in decimal representation of 
                                                        // indices differing only in the value of
                                                        // the specified target

    // dim_t offsetDistance = POW2(control, dim_t);     // distance from the all zero bit to the
                                                        // first bit with a one at the position of
                                                        // the control qubit
    /*
    Starting in the first block, i.e., the block where all bits in the binary representation left to
    the one specified by target are zero, the integer i, which serves as the starting index in the
    block, is set to zero. After each block, consisting of 2^{target - 1} potential swaps, i is
    increased by 2^{target + 1}. This corresponds to gradually increasing the bits left to target
    in the index's binary representation and reseting those right to it to zero until all blocks are
    evaluated when i reaches the integer dimension.
    */
    for (i = 0; i < state->dimension; i += blockDistance) {
        /*
        The inner loop starts with all bits right to target in the index's binary representation
        initialized to zero, as handed over by the outer loop. The integer j then counts through all
        possible integers within block distance, thereby gradually changing the bits right to traget
        in the index's binary representation until no further incrementation in these bits is
        possible.    
        */
        for (j = i; j < i + flipDistance; ++j) {
            /*
            The CNOT operation is executed only if the bit at position control in the current 
            index's binary representation is equal to one. This is done by a bitwise AND operation
            on the current integer and 2^control. Since the latter's bit representation has only
            zeros except for a one at position control from the right, the resulting bit string has
            only zeros except for a one at position control from the right only if j's bit
            representation has a one there.
            */
            if (j & POW2(control ,dim_t)) {
                /*
                Meeting the condition, the value of the j-th entry originating from the pointer
                stored in state is stored to a temporary double complex variable. Subesequently, one
                sets the value the j-th entry's pointer holds to be the value of the entry whose
                index's binary representation only differs in qubit. Finally, the latter pointer is
                set to the value stored under the temporary and j.
                */
                SWAP(state->vector + j, state->vector + (j + flipDistance), cplx_t);
            }
        }
    }
}

void applyCY(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                SWAP(state->vector + j, state->vector + (j + flipDistance), cplx_t);
                state->vector[j] *= -I;
                state->vector[j + flipDistance] *= I;
            }
        }
    }
}

void applyCZ(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vector[j + flipDistance] *= -1;
            }
        }
    }
}

/*
 * =================================================================================================
 *                                              controlled Clifford gates
 * =================================================================================================
 */

void applyCS(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vector[j + flipDistance] *= I;
            }
        }
    }
}

void applyCSdagger(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vector[j + flipDistance] *= -I;
            }
        }
    }
}

void applyCH(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                tmp = state->vector[j];
                state->vector[j] = INVSQRT2 * (tmp + state->vector[j + flipDistance]);
                state->vector[j + flipDistance] = INVSQRT2 * (tmp - state->vector[j + flipDistance]);
            }
        }
    }
}

/*
 * =================================================================================================
 *                                              controlled Hadamard-Y gate
 * =================================================================================================
 */

void applyCHy(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                tmp = state->vector[j];
                state->vector[j] = INVSQRT2 * (tmp - I * state->vector[j + flipDistance]);
                state->vector[j + flipDistance] = INVSQRT2 * (I * tmp - state->vector[j + flipDistance]);
            }
        }
    }
}

/*
 * =================================================================================================
 *                                              controlled T gate
 * =================================================================================================
 */

void applyCT(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vector[j + flipDistance] *= INVSQRT2 + INVSQRT2 * I;
            }
        }
    }
}

void applyCTdagger(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vector[j + flipDistance] *= INVSQRT2 - INVSQRT2 * I;
            }
        }
    }
}

/*
 * =================================================================================================
 *                                              controlled P gate
 * =================================================================================================
 */

void applyCP(state_t* state, qubit_t control, qubit_t target, double angle) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dimension; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vector[j + flipDistance] *= cos(angle) + sin(angle) * I;
            }
        }
    }
}

void applyCPdagger(state_t* state, qubit_t control, qubit_t target, double angle) {
    applyCP(state, control, target, -angle);
}

/*
 * =================================================================================================
 *                                              SWAP gate
 * =================================================================================================
 */

/*
This function swaps two qubits in the system resulting in a rearrangement of the amplitudes in the
state vector, where all entries whose index's binary representation differs at the two specified
positions are exchanged.

Input:
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unisgned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    qubit_t qubit1: unsigned integer holding the position of the qubit from the right that is to be
    swapped with qubit2
    qubit_t qubit2: unsigned integer holding the position of the qubit from the right that is to be
    swapped with qubit1
    ASSUME qubit1 < qubit2
Output:
    There is no output value, but the array's elements the pointer in state refers to are exchanged.
    Each entry whose index in binary representation at position qubit1 differs from that at position
    qubit2 is exchanged with the entry whose index in binary representation differs at position
    qubit1 and qubit2 from the earlier.
    
*/
void applySWAP(state_t* state, qubit_t qubit1, qubit_t qubit2)
{
    dim_t i;                                                        // index that determines the
                                                                    // starting point of the inner
                                                                    // loop, only changing the bits
                                                                    // in the index's binary
                                                                    // representation left to the
                                                                    // left qubit

    dim_t j;                                                        // index for the inner loop,
                                                                    // running over all values that
                                                                    // only differ in the bits right
                                                                    // to the left qubit in the
                                                                    // index's binary representation

    qubit_t left_qubit = MAX(qubit1, qubit2);                       // index of the left qubit that
                                                                    // is to be swapped

    qubit_t right_qubit = MIN(qubit1, qubit2);                      // index of the right qubit that
                                                                    // that is to be swapped

    dim_t majorBlockDistance = POW2(left_qubit + 1, dim_t);         // distance of two indices whose
                                                                    // values differ only by one in 
                                                                    // the bits left to the position
                                                                    // of the left qubit

    dim_t minorBlockDistance = POW2(left_qubit, dim_t);             // distance of two indices whose
                                                                    // values differ only in the bit
                                                                    // at the position of the left
                                                                    // qubit

    dim_t offsetDistance = POW2(right_qubit, dim_t);                // distance of two indices whose
                                                                    // values differ only in the bit
                                                                    // at the position of the right
                                                                    // qubit
                                                                                                      
    dim_t flipDistance = POW2(left_qubit, dim_t) - POW2(right_qubit, dim_t); 
                                                                    // distance of the index whose
                                                                    // bit at the position of the
                                                                    // left and right qubit is 1
                                                                    // and 0, respectively, to the
                                                                    // index whose bit at the
                                                                    // position of the left and
                                                                    // right position is 0 and 1,
                                                                    // respectively
    /*
    Starting with indices whose bits left to the one specified by the left qubit set to zero in the
    index's binary representation, the integer i marks the starting point of the inner loop. After
    each iteration, taking 2^{qubit2} executions of the inner loop, the bits left to the one 
    specified by
    qubit2 are gradually increased by one, which corresponds to an increment of 2^{qubit2 + 1}.
    The loop terminates, when all bits left to the one specified by qubit2 in the index's binary
    representation are set to one.
    */
    for (i = 0; i < state->dimension; i += majorBlockDistance) {
        /*
        The inner loop always starts with all bits right to the one specified by the left qubit in
        the index's binary representation initialised to zero except for the one at the position
        of the right qubit, since it is the first that meets the subsequent condition. After each
        iteration j is increased by one, thereby going through all possible bit strings of the bits
        right to the left qubit. The loop terminates when all bits right to the one specified by the
        left qubit in the index's binary representation are set to one. The bit at the position of
        the left qubit in the index's binary representation remains set to zero.
        */
        for (j = i + offsetDistance; j < i + minorBlockDistance; ++j) {
            /*
            For all indices where the bit at the position of the right qubit in the index's binary
            representation is set to one, a swap is executed. Since the bit at the position of the
            left qubit in the index's binary representation remains set to zero throughout the
            entire loop, the indices that meet the condition have a zero at the left qubit's
            position and a one at the right qubit's position in the index's binary representation.
            Its entry is exchanged with the entry whose index has a one at the left qubit's position
            and a zero at right qubit's position in the index's binary representation. 
            */
            if (j & POW2(right_qubit, dim_t)) {
                /*
                Meeting the condition, the value of the j-th entry originating from the pointer
                stored in state is stored to a temporary double complex variable. Subsequently, one
                sets the value the j-th entry's pointer holds to be the value of the entry whose
                index's binary representation only differs in the left and the right qubit. Finally,
                the latter pointer is set to the value stored under the temporary and j.
                */
                SWAP(state->vector + j, state->vector + (j + flipDistance), cplx_t);
            }
        }
    }
}

/*
This function applies the Toffoli gate to the state.

Input
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unisgned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    qubit_t control1: an unsigned integer holding the first of the two qubits' positions from the
    right in a bit string that control the execution of the Pauli-X gate on the target qubit
    qubit_t control2: an unsigned integer holding the second of the two qubits' positions from the
    right in a bit string that control the execution of the Pauli-X gate on the target qubit
    qubit_t target: an unsigned integer holding the qubit's position from the right in a bit string
    the Pauli-X gate acts on if the control qubits both are in the |1> state
    
*/

/*
 * =================================================================================================
 *                                         Toffoli gate
 * =================================================================================================
 */

void applyToffoli(state_t* state, qubit_t control1, qubit_t control2, qubit_t target)
{
    dim_t i;                                            // index for the outer loop which runs
                                                        // through the blocks, i.e., increseases
                                                        // the bits left to target by one

    dim_t j;                                            // index of the inner loop running through
                                                        // all indices within the block, i.e.,
                                                        // increases bits right to target by one

    dim_t blockDistance = POW2(target + 1, dim_t);      // distance of two integers whose binary
                                                        // representation differs only by one seen
                                                        // from all bits left to target

    dim_t flipDistance = POW2(target, dim_t);           // distance of two integers whose binary
                                                        // representation only differs in the bit 
                                                        // at position target

    /*
    Starting with all bits set to zero the integer i is the starting index for the inner loop and
    is increased by blockDistance after each iteration, thereby incrementing the bits left to target
    by one.
    */
    for (i = 0; i < state->dimension; i += blockDistance)
    {
        /*
        The integer j starts at the same configuration of bits in the index's binary representation
        gradually changing all bits to the right of target until all are set to one
        */
        for (j = i; j < i + flipDistance; ++j)
        {
            /*
            The Pauli-X gate is applied to target only if both control bits in the index's binary 
            representation are set to one.
            */
            if (j & POW2(control1, dim_t) && j & POW2(control2, dim_t)) 
            {
                /*
                Meeting the condition, the value of the j-th entry originating from the pointer
                stored in state is stored to a temporary double complex variable. Subesequently, one
                sets the value the j-th entry's pointer holds to be the value of the entry whose
                index's binary representation only differs in qubit. Finally, the latter pointer is
                set to the value stored under the temporary and j.
                */
                SWAP(state->vector + j, state->vector + (j + flipDistance), cplx_t);
            }
        }
    }
}
