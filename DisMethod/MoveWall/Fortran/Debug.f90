MODULE DEBUG
    USE MESH
    USE MATH
    USE GLOBAL_PARAMETER
    
    IMPLICIT NONE
    
    REAL, ALLOCATABLE, SAVE :: LP(:, :, :)
    
    LOGICAL, SAVE :: ALL_ALLOCATED = .FALSE.
    
    PRIVATE ALL_ALLOCATED
    
    CONTAINS
    
    SUBROUTINE ALLOC_DEBUG
        IMPLICIT NONE
        IF(.NOT. ALL_ALLOCATED) THEN
            ALLOCATE(LP(N1, N2, N3))
            ALL_ALLOCATED = .TRUE.
        END IF
    END SUBROUTINE ALLOC_DEBUG
    
    SUBROUTINE LAPLACE(DP, PB, TT)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: DP(:, :, :), PB(:, :, :), TT
        REAL DPDX2, DPDY2, DPDZ2
        REAL DPDYL(3), DPDX(3), DPDZ(3), DPDY
        REAL DYFAC, DY2FAC
        REAL CROSS, C1, C3
        REAL XC, YC, ZC
        INTEGER I, J, K
        
        CALL ALLOC_DEBUG
        
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            DO I = 1, N1
                XC = (X(I) + X(I - 1)) / 2
                DO J = 2, N2 - 1
                    YC = (Y(J) + Y(J - 1)) / 2
                    
                    DPDX2 = (DP(IP(I), J, K) - 2 * DP(I, J, K) + DP(IM(I), J, K)) / DX / DX
                    DPDZ2 = (DP(I, J, KP(K)) - 2 * DP(I, J, K) + DP(I, J, KM(K))) / DZ / DZ
                    DPDY2 = DP(I, JP(J), K) * DY2H(1, J) + DP(I, J, K) * DY2H(2, J) + DP(I, JM(J), K) * DY2H(3, J)
                    DPDY = DP(I, JP(J), K) * DYH(1, J) + DP(I, J, K) * DYH(2, J) + DP(I, JM(J), K) * DYH(3, J)
                    
                    DYFAC = PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT)  &
                          + PHI2(XC, YC, ZC, TT) * DPHI2DY(XC, YC, ZC, TT)  &
                          + PHI3(XC, YC, ZC, TT) * DPHI3DY(XC, YC, ZC, TT)  &
                          + DPHI1DX(XC, YC, ZC, TT) + DPHI2DY(XC, YC, ZC, TT)&
                          + DPHI3DZ(XC, YC, ZC, TT)
                    
                    DY2FAC = PHI1(XC, YC, ZC, TT) ** 2 + PHI2(XC, YC, ZC, TT) ** 2 + PHI3(XC, YC, ZC, TT) ** 2 &
                           + 2 * PHI2(XC, YC, ZC, TT)
                    
                    DPDX(1) = (DP(IP(I), JP(J), K) - DP(IM(I), JP(J), K)) / 2 / DX
                    DPDX(2) = (DP(IP(I), J, K) - DP(IM(I), J, K)) / 2 / DX
                    DPDX(3) = (DP(IP(I), JM(J), K) - DP(IM(I), JM(J), K)) / 2 / DX
                    C1 = DYH(1, J) * DPDX(1) + DYH(2, J) * DPDX(2) + DYH(3, J) * DPDX(3)
                    
                    DPDZ(1) = (DP(I, JP(J), KP(K)) - DP(I, JP(J), KM(K))) / 2 / DZ
                    DPDZ(2) = (DP(I, J, KP(K)) - DP(I, J, KM(K))) / 2 / DZ
                    DPDZ(3) = (DP(I, JM(J), KP(K)) - DP(I, JM(J), KM(K))) / 2 / DZ
                    C3 = DYH(1, J) * DPDZ(1) + DYH(2, J) * DPDZ(2) + DYH(3, J) * DPDZ(3)
                    
                    CROSS = C1 * 2 * PHI1(XC, YC, ZC, TT) + C3 * 2 * PHI3(XC, YC, ZC, TT)
                    
                    LP(I, J, K) = DPDX2 + DPDY2 + DPDZ2 &
                                + CROSS + DYFAC * DPDY + DY2FAC * DPDY2
                END DO
                J = 1
                YC = (Y(0) + Y(1)) / 2
                DPDX2 = (DP(IP(I), J, K) - 2 * DP(I, J, K) + DP(IM(I), J, K)) / DX / DX
                DPDZ2 = (DP(I, J, KP(K)) - 2 * DP(I, J, K) + DP(I, J, KM(K))) / DZ / DZ
                DPDY2 = ((DP(I, 2, K) - DP(I, 1, K)) / H(2) - PB(I, 1, K)) / DY(1)
                DPDY =  ((DP(I, 2, K) - DP(I, 1, K)) / H(2) + PB(I, 1, K)) / 2
                    
                DYFAC = PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT)  &
                        + PHI2(XC, YC, ZC, TT) * DPHI2DY(XC, YC, ZC, TT)  &
                        + PHI3(XC, YC, ZC, TT) * DPHI3DY(XC, YC, ZC, TT)  &
                        + DPHI1DX(XC, YC, ZC, TT) + DPHI2DY(XC, YC, ZC, TT)&
                        + DPHI3DZ(XC, YC, ZC, TT)
                    
                DY2FAC = PHI1(XC, YC, ZC, TT) ** 2 + PHI2(XC, YC, ZC, TT) ** 2 + PHI3(XC, YC, ZC, TT) ** 2 &
                        + 2 * PHI2(XC, YC, ZC, TT)
                    
                DPDYL(1) = ((DP(IP(I), JP(J), K) - DP(IP(I), J, K)) / H(2) + PB(IP(I), 1, K)) / 2
                DPDYL(3) = ((DP(IM(I), JP(J), K) - DP(IM(I), J, K)) / H(2) + PB(IM(I), 1, K)) / 2
                C1 = (DPDYL(1) - DPDYL(3)) / 2 / DX
                    
                DPDYL(1) = ((DP(I, JP(J), KP(K)) - DP(I, J, KP(K))) / H(2) + PB(I, 1, KP(K))) / 2
                DPDYL(3) = ((DP(I, JP(J), KM(K)) - DP(I, J, KM(K))) / H(2) + PB(I, 1, KM(K))) / 2
                C3 = (DPDYL(1) - DPDYL(3)) / 2 / DZ
                    
                CROSS = C1 * 2 * PHI1(XC, YC, ZC, TT) + C3 * 2 * PHI3(XC, YC, ZC, TT)
                    
                LP(I, J, K) = DPDX2 + DPDY2 + DPDZ2 &
                            + CROSS + DYFAC * DPDY + DY2FAC * DPDY2
                
                J = N2
                YC = (Y(N2) + Y(N2 - 1)) / 2
                DPDX2 = (DP(IP(I), J, K) - 2 * DP(I, J, K) + DP(IM(I), J, K)) / DX / DX
                DPDZ2 = (DP(I, J, KP(K)) - 2 * DP(I, J, K) + DP(I, J, KM(K))) / DZ / DZ
                DPDY2 = (PB(I, 2, K) - (DP(I, N2, K) - DP(I, N2 - 1, K)) / H(N2)) / DY(N2)
                DPDY =  (PB(I, 2, K) + (DP(I, N2, K) - DP(I, N2 - 1, K)) / H(N2)) / 2
                    
                DYFAC = PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT)  &
                        + PHI2(XC, YC, ZC, TT) * DPHI2DY(XC, YC, ZC, TT)  &
                        + PHI3(XC, YC, ZC, TT) * DPHI3DY(XC, YC, ZC, TT)  &
                        + DPHI1DX(XC, YC, ZC, TT) + DPHI2DY(XC, YC, ZC, TT)&
                        + DPHI3DZ(XC, YC, ZC, TT)
                    
                DY2FAC = PHI1(XC, YC, ZC, TT) ** 2 + PHI2(XC, YC, ZC, TT) ** 2 + PHI3(XC, YC, ZC, TT) ** 2 &
                        + 2 * PHI2(XC, YC, ZC, TT)
                    
                DPDYL(1) = (PB(IP(I), 2, K) + (DP(IP(I), N2, K) - DP(IP(I), N2 - 1, K)) / H(N2)) / 2
                DPDYL(3) = (PB(IM(I), 2, K) + (DP(IM(I), N2, K) - DP(IM(I), N2 - 1, K)) / H(N2)) / 2
                C1 = (DPDYL(1) - DPDYL(3)) / 2 / DX
                    
                DPDYL(1) = (PB(I, 2, KP(K)) + (DP(I, N2, KP(K)) - DP(I, N2 - 1, KP(K))) / H(N2)) / 2
                DPDYL(3) = (PB(I, 2, KM(K)) + (DP(I, N2, KM(K)) - DP(I, N2 - 1, KM(K))) / H(N2)) / 2
                C3 = (DPDYL(1) - DPDYL(3)) / 2 / DZ
                    
                CROSS = C1 * 2 * PHI1(XC, YC, ZC, TT) + C3 * 2 * PHI3(XC, YC, ZC, TT)
                    
                LP(I, J, K) = DPDX2 + DPDY2 + DPDZ2 &
                            + CROSS + DYFAC * DPDY + DY2FAC * DPDY2
            END DO
        END DO
    END SUBROUTINE LAPLACE
    
    SUBROUTINE LAPLACES(DPS, PBS, TT)
        IMPLICIT NONE
        
        
    
    SUBROUTINE TEST_PSOLVE(DP, PB, RP, TT)
        IMPLICIT NONE
        REAL, INTENT(IN) :: DP(:, :, :), PB(:, :, :), RP(:, :, :), TT
        
        CALL LAPLACE(DP, PB, TT)
        LP = LP - RP
        PRINT*, MAXVAL(LP), MAXLOC(LP)
    END SUBROUTINE TEST_PSOLVE
    
    END MODULE DEBUG