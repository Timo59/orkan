// test_gates.c  - Unit tests for quantum gate operations on pure, packed, and tiled states

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef TEST_GATE_H
#include "test_gate.h"
#endif

/*
 * =====================================================================================================================
 * Unity: setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}


/*
 * =====================================================================================================================
 * Test functions: Pure states
 * =====================================================================================================================
 */

void test_X_pure(void) {testSingleQubitGate(x, XMAT);}
void test_Y_pure(void) {testSingleQubitGate(y, YMAT);}
void test_Z_pure(void) {testSingleQubitGate(z, ZMAT);}
void test_H_pure(void) {testSingleQubitGate(h, HMAT);}
void test_Hy_pure(void) {testSingleQubitGate(hy, HYMAT);}
void test_S_pure(void) {testSingleQubitGate(s, SMAT);}
void test_Sdg_pure(void) {testSingleQubitGate(sdg, SDGMAT);}
void test_T_pure(void) {testSingleQubitGate(t, TMAT);}
void test_Tdg_pure(void) {testSingleQubitGate(tdg, TDGMAT);}
void test_Rx_pure(void) {testRotationGate(rx, mat_rx);}
void test_Ry_pure(void) {testRotationGate(ry, mat_ry);}
void test_Rz_pure(void) {testRotationGate(rz, mat_rz);}
void test_P_pure(void) {testRotationGate(p, mat_p);}
void test_CX_pure(void) {testTwoQubitGate(cx, CXMAT);}
void test_CY_pure(void) {testTwoQubitGate(cy, CYMAT);}
void test_CZ_pure(void) {testTwoQubitGate(cz, CZMAT);}
void test_SWAP_pure(void) {testTwoQubitGate(swap_gate, SWAPMAT);}
void test_CCX_pure(void) {testThreeQubitGate(ccx, CCXMAT);}

/*
 * =====================================================================================================================
 * Test functions: Packed mixed states
 * =====================================================================================================================
 */

void test_X_packed(void) {testSingleQubitGateMixed(x, XMAT);}
void test_Y_packed(void) {testSingleQubitGateMixed(y, YMAT);}
void test_Z_packed(void) {testSingleQubitGateMixed(z, ZMAT);}
void test_H_packed(void) {testSingleQubitGateMixed(h, HMAT);}
void test_Hy_packed(void) {testSingleQubitGateMixed(hy, HYMAT);}
void test_S_packed(void) {testSingleQubitGateMixed(s, SMAT);}
void test_Sdg_packed(void) {testSingleQubitGateMixed(sdg, SDGMAT);}
void test_T_packed(void) {testSingleQubitGateMixed(t, TMAT);}
void test_Tdg_packed(void) {testSingleQubitGateMixed(tdg, TDGMAT);}
void test_Rx_packed(void) {testRotationGateMixed(rx, mat_rx);}
void test_Ry_packed(void) {testRotationGateMixed(ry, mat_ry);}
void test_Rz_packed(void) {testRotationGateMixed(rz, mat_rz);}
void test_P_packed(void) {testRotationGateMixed(p, mat_p);}
void test_CX_packed(void) {testTwoQubitGateMixed(cx, CXMAT);}
void test_CY_packed(void) {testTwoQubitGateMixed(cy, CYMAT);}
void test_CZ_packed(void) {testTwoQubitGateMixed(cz, CZMAT);}
void test_SWAP_packed(void) {testTwoQubitGateMixed(swap_gate, SWAPMAT);}
void test_CCX_packed(void) {testThreeQubitGateMixed(ccx, CCXMAT);}

/*
 * =====================================================================================================================
 * Test functions: Tiled mixed states (temporarily disabled)
 * =====================================================================================================================
 */

void test_X_tiled(void) {testSingleQubitGateTiled(x, XMAT);}
void test_Y_tiled(void) {testSingleQubitGateTiled(y, YMAT);}
void test_Z_tiled(void) {testSingleQubitGateTiled(z, ZMAT);}
void test_H_tiled(void) {testSingleQubitGateTiled(h, HMAT);}
void test_Hy_tiled(void) {testSingleQubitGateTiled(hy, HYMAT);}
void test_S_tiled(void) {testSingleQubitGateTiled(s, SMAT);}
void test_Sdg_tiled(void) {testSingleQubitGateTiled(sdg, SDGMAT);}
void test_T_tiled(void) {testSingleQubitGateTiled(t, TMAT);}
void test_Tdg_tiled(void) {testSingleQubitGateTiled(tdg, TDGMAT);}
void test_Rx_tiled(void) {testRotationGateTiled(rx, mat_rx);}
void test_Ry_tiled(void) {testRotationGateTiled(ry, mat_ry);}
void test_Rz_tiled(void) {testRotationGateTiled(rz, mat_rz);}
void test_P_tiled(void) {testRotationGateTiled(p, mat_p);}
void test_CX_tiled(void) {testTwoQubitGateTiled(cx, CXMAT);}
void test_CY_tiled(void) {testTwoQubitGateTiled(cy, CYMAT);}
void test_CZ_tiled(void) {testTwoQubitGateTiled(cz, CZMAT);}
void test_SWAP_tiled(void) {testTwoQubitGateTiled(swap_gate, SWAPMAT);}
void test_CCX_tiled(void) {testThreeQubitGateTiled(ccx, CCXMAT);}

/*
 * =====================================================================================================================
 * Test functions: single_from_mat (tiled mixed states)
 * =====================================================================================================================
 */

static void mat_rx_pi3(cplx_t *mat) {
    const double c = cos(M_PI / 6.0);  /* cos(pi/6) = sqrt(3)/2 */
    const double s = sin(M_PI / 6.0);  /* sin(pi/6) = 1/2       */
    mat[0] =  c + 0.0 * I;
    mat[1] = 0.0 -   s * I;
    mat[2] = 0.0 -   s * I;
    mat[3] =  c + 0.0 * I;
}

void test_single_from_mat_id_tiled(void) { testSingleQubitMatGateTiled(IDMAT); }
void test_single_from_mat_x_tiled(void)  { testSingleQubitMatGateTiled(XMAT);  }
void test_single_from_mat_y_tiled(void)  { testSingleQubitMatGateTiled(YMAT);  }
void test_single_from_mat_h_tiled(void)  { testSingleQubitMatGateTiled(HMAT);  }
void test_single_from_mat_t_tiled(void)  { testSingleQubitMatGateTiled(TMAT);  }

void test_single_from_mat_rx_tiled(void) {
    cplx_t rx_mat[4];
    mat_rx_pi3(rx_mat);
    testSingleQubitMatGateTiled(rx_mat);
}

void test_single_from_mat_xh_tiled(void) { testSingleQubitMatGateTiledDouble(XMAT, HMAT); }
void test_single_from_mat_ht_tiled(void) { testSingleQubitMatGateTiledDouble(HMAT, TMAT); }

/*
 * =====================================================================================================================
 * Unity body
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();

    printf("\n========== PURE Gate Tests ==========\n");
    RUN_TEST(test_X_pure);
    RUN_TEST(test_Y_pure);
    RUN_TEST(test_Z_pure);
    RUN_TEST(test_H_pure);
    RUN_TEST(test_Hy_pure);
    RUN_TEST(test_S_pure);
    RUN_TEST(test_Sdg_pure);
    RUN_TEST(test_T_pure);
    RUN_TEST(test_Tdg_pure);
    RUN_TEST(test_Rx_pure);
    RUN_TEST(test_Ry_pure);
    RUN_TEST(test_Rz_pure);
    RUN_TEST(test_P_pure);
    RUN_TEST(test_CX_pure);
    RUN_TEST(test_CY_pure);
    RUN_TEST(test_CZ_pure);
    RUN_TEST(test_SWAP_pure);
    RUN_TEST(test_CCX_pure);

    printf("\n========== MIXED_PACKED Gate Tests ==========\n");
    RUN_TEST(test_X_packed);
    RUN_TEST(test_Y_packed);
    RUN_TEST(test_Z_packed);
    RUN_TEST(test_H_packed);
    RUN_TEST(test_Hy_packed);
    RUN_TEST(test_S_packed);
    RUN_TEST(test_Sdg_packed);
    RUN_TEST(test_T_packed);
    RUN_TEST(test_Tdg_packed);
    RUN_TEST(test_Rx_packed);
    RUN_TEST(test_Ry_packed);
    RUN_TEST(test_Rz_packed);
    RUN_TEST(test_P_packed);
    RUN_TEST(test_CX_packed);
    RUN_TEST(test_CY_packed);
    RUN_TEST(test_CZ_packed);
    RUN_TEST(test_SWAP_packed);
    RUN_TEST(test_CCX_packed);

    printf("\n========== MIXED_TILED Gate Tests ==========\n");
    RUN_TEST(test_X_tiled);
    RUN_TEST(test_Y_tiled);
    RUN_TEST(test_Z_tiled);
    RUN_TEST(test_H_tiled);
    RUN_TEST(test_Hy_tiled);
    RUN_TEST(test_S_tiled);
    RUN_TEST(test_Sdg_tiled);
    RUN_TEST(test_T_tiled);
    RUN_TEST(test_Tdg_tiled);
    RUN_TEST(test_Rx_tiled);
    RUN_TEST(test_Ry_tiled);
    RUN_TEST(test_Rz_tiled);
    RUN_TEST(test_P_tiled);
    RUN_TEST(test_CX_tiled);
    RUN_TEST(test_CY_tiled);
    RUN_TEST(test_CZ_tiled);
    RUN_TEST(test_SWAP_tiled);
    RUN_TEST(test_CCX_tiled);

    printf("\n========== single_from_mat Tiled Tests ==========\n");
    RUN_TEST(test_single_from_mat_id_tiled);
    RUN_TEST(test_single_from_mat_x_tiled);
    RUN_TEST(test_single_from_mat_y_tiled);
    RUN_TEST(test_single_from_mat_h_tiled);
    RUN_TEST(test_single_from_mat_t_tiled);
    RUN_TEST(test_single_from_mat_rx_tiled);
    RUN_TEST(test_single_from_mat_xh_tiled);
    RUN_TEST(test_single_from_mat_ht_tiled);

    return UNITY_END();
}
