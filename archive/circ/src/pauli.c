/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef PAULI_H
#include "pauli.h"
#endif

/*
 * =====================================================================================================================
 *                                              compound gate constructor
 * =====================================================================================================================
 */

// void initCompgate(compgate_t* compgate, complength_t length, qubit_t qubits) {
//     compgate->components = (pauli_t*) calloc(length * qubits, sizeof(pauli_t));
//     compgate->coefficients = (cplx_t*) calloc(length, sizeof(cplx_t));
//     compgate->length = length;
//     compgate->qubits = qubits;
// }

// void freeCompgate(compgate_t compgate) {
//     free(compgate.components);
//     free(compgate.coefficients);
// }

/* 
 * =====================================================================================================================
 *                                              compound gate fusion
 * =====================================================================================================================
 */

// void fuseCompgates(compgate_t* result, compgate_t compgate1, compgate_t compgate2) {
//     complength_t entry;
//     int8_t imagCounter;
//     for (complength_t i = 0; i < compgate1.length; ++i) {
//         for (complength_t j = 0; j < compgate2.length; ++j) {
//             entry = i * compgate2.length + j;
//             result->coefficients[entry] = compgate1.coefficients[i] * compgate2.coefficients[j];
//             imagCounter = 0;
//             for (qubit_t k = 0; k < compgate1.components->qubits; ++k) {
//                 switch (compgate1.components[i].paulis[k]) {
//                     case ID: {
//                         result->components[entry].paulis[k] = compgate2.components[j].paulis[k];
//                         break;
//                     }
//                     case Z: {
//                         switch (compgate2.components[j].paulis[k]) {
//                             case ID: {
//                                 result->components[entry].paulis[k] = Z;
//                                 break;
//                             }
//                             case Z: {
//                                 result->components[entry].paulis[k] = ID;
//                                 break;
//                             }
//                             case X: {
//                                 result->components[entry].paulis[k] = Y;
//                                 imagCounter += 1;
//                                 break;
//                             }
//                             case Y: {
//                                 result->components[entry].paulis[k] = X;
//                                 imagCounter -= 1;
//                                 break;
//                             }
//                             default: {
//                                 printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                                 exit(1);
//                             }  
//                         }
//                         break;
//                     }
//                     case X: {
//                         switch (compgate2.components[j].paulis[k]) {
//                             case ID: {
//                                 result->components[entry].paulis[k] = X;
//                                 break;
//                             }
//                             case Z: {
//                                 result->components[entry].paulis[k] = Y;
//                                 imagCounter -= 1;
//                                 break;
//                             }
//                             case X: {
//                                 result->components[entry].paulis[k] = ID;
//                                 break;
//                             }
//                             case Y: {
//                                 result->components[entry].paulis[k] = Z;
//                                 imagCounter += 1;
//                                 break;
//                             }
//                             default: {
//                                 printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                                 exit(1);
//                             } 
//                         }
//                         break;
//                     }
//                     case Y: {
//                         switch (compgate2.components[j].paulis[k]) {
//                             case ID: {
//                                 result->components[entry].paulis[k] = Y;
//                                 break;
//                             }
//                             case Z: {
//                                 result->components[entry].paulis[k] = X;
//                                 imagCounter += 1;
//                                 break;
//                             }
//                             case X: {
//                                 result->components[entry].paulis[k] = Z;
//                                 imagCounter -= 1;
//                                 break;
//                             }
//                             case Y: {
//                                 result->components[entry].paulis[k] = ID;
//                                 break;
//                             }
//                             default: {
//                                 printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                                 exit(1);
//                             } 
//                         }
//                         break;                        
//                     }
//                     default: {
//                         printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                         exit(1);
//                     }
//                 }
//             }
//             multiplyByIPower(&(result->coefficients[entry]), imagCounter);
//         }
//     }
// }

/*
 * =====================================================================================================================
 *                      Evolve with Pauli observable with Trotter approximation to first order
 * =====================================================================================================================
 */

/*
 * This function evolves a state with an observable in Trotter approximation to first order
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliObs_t*:            Pauli observable, i.e., a real linear combination of Pauli strings
 *      double angle:           Time of the evolution
 * Output:
 *      This function has no return value, but sequentially evolves the input state with Pauli strings as if the
 *      observable's components commute.
 */

void evolveObsPauliTrotter(state_t* state, const pauliObs_t* obs, double angle) {
    for (complength_t i = 0; i < obs->length; ++i) {
        evolvePauliStr(state, obs->comps + (i * obs->qubits), angle * obs->coeffs[i]);
    }
}

void evolveObsPauliTrotterDagger(state_t* state, const pauliObs_t* obs, double angle) {
    for (complength_t i = 0; i < obs->length; ++i) {
        evolvePauliStr(state,
                       obs->comps + ((obs->length - i - 1) * obs->qubits),
                       -angle * obs->coeffs[obs->length - i - 1]);
    }
}
