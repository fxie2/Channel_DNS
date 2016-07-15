MODULE DEBUG
    USE MESH
    USE MATH
    USE GLOBAL_PARAMETER
    
    IMPLICIT NONE
    
    REAL, ALLOCATABLE, SAVE :: LP(:, :, :), LP1(:, :, :), LP2(:, :, :), PHI_FIELD(:, :, :), LPX(:, :, :)
    REAL, ALLOCATABLE, SAVE :: U(:, :, :), V(:, :, :), W(:, :, :), DIV(:, :, :)
    COMPLEX, ALLOCATABLE, SAVE :: LPS(:, :, :)
    
    LOGICAL, SAVE :: ALL_ALLOCATED = .FALSE.
    
    PRIVATE U, V, W, DIV
    PRIVATE GETDIV
    PRIVATE ALL_ALLOCATED
    
    CONTAINS
    
    SUBROUTINE ALLOC_DEBUG
        IMPLICIT NONE
        IF(.NOT. ALL_ALLOCATED) THEN
            ALLOCATE(U(N1, 0:N2+1, N3))
            ALLOCATE(V(N1, 1:N2+1, N3))
            ALLOCATE(W(N1, 0:N2+1, N3))
            ALLOCATE(DIV(N1, N2, N3))
            ALLOCATE(LP(N1, N2, N3))
            ALLOCATE(LP1(N1, N2, N3))
            ALLOCATE(LP2(N1, N2, N3))
            ALLOCATE(LPS(N1, N2, N3))
            ALLOCATE(LPX(N1, N2, N3))
            ALLOCATE(PHI_FIELD(N1, N2, N3))
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
        
        COMPLEX, INTENT(IN) :: DPS(:, :, :), PBS(:, :, :)
        REAL, INTENT(IN) :: TT
        
        REAL ALPHA, BETA
        INTEGER I, J, K, M, N
        REAL XC, YC, ZC
        COMPLEX D2FAC, DFAC, FAC
        COMPLEX PPFAC, PCFAC, PMFAC
        COMPLEX, PARAMETER :: II = (0, 1)
        
        ALPHA = 2 * PI / LX
        BETA = 2 * PI / LZ
        
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            N = K - 1 - N3 / 2
            DO I = 1, N1
                XC = (X(I) + X(I - 1)) / 2
                M = I - 1 - N1 / 2
                DO J = 2, N2 - 1
                    YC = (Y(J) + Y(J - 1)) / 2
                    
                    !D2_DY2 FACTOR
                    D2FAC = 1 + PHI2(XC, YC, ZC, TT) * 2 &
                          + PHI1(XC, YC, ZC, TT) ** 2 &
                          + PHI2(XC, YC, ZC, TT) ** 2 &
                          + PHI3(XC, YC, ZC, TT) ** 2
                    
                    !D_DY FACTOR
                    DFAC = DPHI1DX(XC, YC, ZC, TT) &
                         + DPHI2DY(XC, YC, ZC, TT) &
                         + DPHI3DZ(XC, YC, ZC, TT) &
                         + II * SIN(ALPHA * M * DX) / DX * PHI1(XC, YC, ZC, TT) * 2 &
                         + II * SIN(BETA  * N * DZ) / DZ * PHI3(XC, YC, ZC, TT) * 2 &
                         + PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT) &
                         + PHI2(XC, YC, ZC, TT) * DPHI2DY(XC, YC, ZC, TT) &
                         + PHI3(XC, YC, ZC, TT) * DPHI3DY(XC, YC, ZC, TT)
                    
                    !Y FACTOR
                    FAC = 2 * (COS(ALPHA * M * DX) - 1) / DX / DX + 2 * (COS(BETA * N * DZ) - 1) / DZ / DZ
                    
                    PPFAC = DY2H(1, J) * D2FAC + DYH(1, J) * DFAC
                    PCFAC = DY2H(2, J) * D2FAC + DYH(2, J) * DFAC + FAC
                    PMFAC = DY2H(3, J) * D2FAC + DYH(3, J) * DFAC
                    
                    LPS(I, J, K) = PPFAC * DPS(I, JP(J), K) + PCFAC * DPS(I, J, K) + PMFAC * DPS(I, JM(J), K)
                END DO
                
                !J = 1
                YC = (Y(1) + Y(0)) / 2
                
                D2FAC = 1 + PHI2(XC, YC, ZC, TT) * 2 &
                      + PHI1(XC, YC, ZC, TT) ** 2 &
                      + PHI2(XC, YC, ZC, TT) ** 2 &
                      + PHI3(XC, YC, ZC, TT) ** 2
                    
                DFAC = DPHI1DX(XC, YC, ZC, TT) &
                     + DPHI2DY(XC, YC, ZC, TT) &
                     + DPHI3DZ(XC, YC, ZC, TT) &
                     + II * SIN(ALPHA * M * DX) / DX * PHI1(XC, YC, ZC, TT) * 2 &
                     + II * SIN(BETA  * N * DZ) / DZ * PHI3(XC, YC, ZC, TT) * 2 &
                     + PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT) &
                     + PHI2(XC, YC, ZC, TT) * DPHI2DY(XC, YC, ZC, TT) &
                     + PHI3(XC, YC, ZC, TT) * DPHI3DY(XC, YC, ZC, TT)
                    
                FAC = 2 * (COS(ALPHA * M * DX) - 1) / DX / DX + 2 * (COS(BETA * N * DZ) - 1) / DZ / DZ
                
                PPFAC = 1 / DY(1) / H(2) * D2FAC + 1 / H(2) / 2 * DFAC
                PCFAC =-1 / DY(1) / H(2) * D2FAC - 1 / H(2) / 2 * DFAC + FAC
                PMFAC = 0
                LPS(I, 1, K) = PPFAC * DPS(I, 2, K) + PCFAC * DPS(I, 1, K) - PBS(I, 1, K) / DY(1) * D2FAC + PBS(I, 1, K) / 2 * DFAC
                
                !J = N2
                YC = (Y(N2) + Y(N2 - 1)) / 2
                
                D2FAC = 1 + PHI2(XC, YC, ZC, TT) * 2 &
                      + PHI1(XC, YC, ZC, TT) ** 2 &
                      + PHI2(XC, YC, ZC, TT) ** 2 &
                      + PHI3(XC, YC, ZC, TT) ** 2
                    
                DFAC = DPHI1DX(XC, YC, ZC, TT) &
                     + DPHI2DY(XC, YC, ZC, TT) &
                     + DPHI3DZ(XC, YC, ZC, TT) &
                     + II * SIN(ALPHA * M * DX) / DX * PHI1(XC, YC, ZC, TT) * 2 &
                     + II * SIN(BETA  * N * DZ) / DZ * PHI3(XC, YC, ZC, TT) * 2 &
                     + PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT) &
                     + PHI2(XC, YC, ZC, TT) * DPHI2DY(XC, YC, ZC, TT) &
                     + PHI3(XC, YC, ZC, TT) * DPHI3DY(XC, YC, ZC, TT)
                    
                FAC = 2 * (COS(ALPHA * M * DX) - 1) / DX / DX + 2 * (COS(BETA * N * DZ) - 1) / DZ / DZ
                
                PPFAC = 0
                PCFAC = -1 / H(N2) / DY(N2) * D2FAC + 1 / H(N2) / 2 * DFAC + FAC
                PMFAC =  1 / H(N2) / DY(N2) * D2FAC - 1 / H(N2) / 2 * DFAC
                LPS(I, N2, K) = PCFAC * DPS(I, N2, K) + PMFAC * DPS(I, N2-1, K) + PBS(I, 2, K) / DY(N2) * D2FAC - PBS(I, 2, K) / 2 * DFAC
                
                !IF(M == 0 .AND. N == 0) THEN
                !    LPS(I, 1, K) = 0
                !END IF
                
            END DO
        END DO
    END SUBROUTINE LAPLACES
    
    SUBROUTINE LAPLACE11(DP, PB, TT)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: DP(:, :, :), PB(:, :, :), TT
        REAL DPDX, DPDY, P1
        INTEGER I, J, K
        REAL XC, YC, ZC
        
        CALL ALLOC_DEBUG
        
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            DO I = 1, N1
                XC = (X(I) + X(I - 1)) / 2
                DO J = 2, N2 - 1
                    YC = (Y(J) + Y(J - 1)) / 2
                    
                    DPDX = (DP(IP(I), J, K) - DP(IM(I), J, K)) / 2 / DX / (SIN(DX) / DX)
                    DPDY = DP(I, JP(J), K) * DYH(1, J) + DP(I, J, K) * DYH(2, J) + DP(I, JM(J), K) * DYH(3, J)
                    P1 = PHI1(XC, YC, ZC, TT)
                    LPX(I, J, K) = DPDX + PHI1(XC, YC, ZC, TT) * DPDY
                END DO
                J = 1
                YC = (Y(1) + Y(0)) / 2
                
                DPDX = (DP(IP(I), J, K) - DP(IM(I), J, K)) / 2 / DX / (SIN(DX) / DX)
                DPDY = ((DP(I, 2, K) - DP(I, 1, K)) / H(2) + PB(I, 1, K)) / 2
                P1 = PHI1(XC, YC, ZC, TT)
                LPX(I, 1, K) = DPDX + PHI1(XC, YC, ZC, TT) * DPDY
                
                J = N2
                YC = (Y(N2) + Y(N2-1)) / 2
                DPDX = (DP(IP(I), J, K) - DP(IM(I), J, K)) / 2 / DX / (SIN(DX) / DX)
                DPDY = (PB(I, 2, K) + (DP(I, N2, K) - DP(I, N2 - 1, K)) / H(N2)) / 2
                P1 = PHI1(XC, YC, ZC, TT)
                LPX(I, N2, K) = DPDX + PHI1(XC, YC, ZC, TT) * DPDY
            END DO
        END DO
    END SUBROUTINE LAPLACE11
    
    SUBROUTINE LAPLACE1(DP, PB, TT)
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
                    DPDY2 = DP(I, JP(J), K) * DY2H(1, J) + DP(I, J, K) * DY2H(2, J) + DP(I, JM(J), K) * DY2H(3, J)
                    DPDY = DP(I, JP(J), K) * DYH(1, J) + DP(I, J, K) * DYH(2, J) + DP(I, JM(J), K) * DYH(3, J)
                    
                    DYFAC = PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT)  &
                          + DPHI1DX(XC, YC, ZC, TT)
                    DY2FAC = PHI1(XC, YC, ZC, TT) ** 2
                    
                    DPDX(1) = (DP(IP(I), JP(J), K) - DP(IM(I), JP(J), K)) / 2 / DX
                    DPDX(2) = (DP(IP(I), J, K) - DP(IM(I), J, K)) / 2 / DX
                    DPDX(3) = (DP(IP(I), JM(J), K) - DP(IM(I), JM(J), K)) / 2 / DX
                    C1 = DYH(1, J) * DPDX(1) + DYH(2, J) * DPDX(2) + DYH(3, J) * DPDX(3)
                    
                    CROSS = C1 * 2 * PHI1(XC, YC, ZC, TT)
                    
                    LP1(I, J, K) = DPDX2 + CROSS + DYFAC * DPDY + DY2FAC * DPDY2
                END DO
                J = 1
                YC = (Y(0) + Y(1)) / 2
                DPDX2 = (DP(IP(I), J, K) - 2 * DP(I, J, K) + DP(IM(I), J, K)) / DX / DX
                DPDY2 = ((DP(I, 2, K) - DP(I, 1, K)) / H(2) - PB(I, 1, K)) / DY(1)
                DPDY = ((DP(I, 2, K) - DP(I, 1, K)) / H(2) + PB(I, 1, K)) / 2
                
                DYFAC = PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT) &
                      + DPHI1DX(XC, YC, ZC, TT)
                
                DY2FAC = PHI1(XC, YC, ZC, TT) ** 2
                
                DPDYL(1) = ((DP(IP(I), JP(J), K) - DP(IP(I), J, K)) / H(2) + PB(IP(I), 1, K)) / 2
                DPDYL(3) = ((DP(IM(I), JP(J), K) - DP(IM(I), J, K)) / H(2) + PB(IM(I), 1, K)) / 2
                C1 = (DPDYL(1) - DPDYL(3)) / 2 / DX
                
                CROSS = C1 * 2 * PHI1(XC, YC, ZC, TT)
                
                LP1(I, J, K) = DPDX2 + CROSS + DYFAC * DPDY + DY2FAC * DPDY2
                
                J = N2
                YC = (Y(N2) + Y(N2 - 1)) / 2
                DPDX2 = (DP(IP(I), J, K) - 2 * DP(I, J, K) + DP(IM(I), J, K)) / DX / DX
                DPDY2 = (PB(I, 2, K) - (DP(I, N2, K) - DP(I, N2 - 1, K)) / H(N2)) / DY(N2)
                DPDY =  (PB(I, 2, K) + (DP(I, N2, K) - DP(I, N2 - 1, K)) / H(N2)) / 2
                
                DYFAC = PHI1(XC, YC, ZC, TT) * DPHI1DY(XC, YC, ZC, TT) &
                      + DPHI1DX(XC, YC, ZC, TT)
                
                DY2FAC = PHI1(XC, YC, ZC, TT) ** 2
                
                DPDYL(1) = (PB(IP(I), 2, K) + (DP(IP(I), N2, K) - DP(IP(I), N2 - 1, K)) / H(N2)) / 2
                DPDYL(3) = (PB(IM(I), 2, K) + (DP(IM(I), N2, K) - DP(IM(I), N2 - 1, K)) / H(N2)) / 2
                C1 = (DPDYL(1) - DPDYL(3)) / 2 / DX
                
                CROSS = C1 * 2 * PHI1(XC, YC, ZC, TT)
                
                LP1(I, J, K) = DPDX2 + CROSS + DYFAC * DPDY + DY2FAC * DPDY2
            END DO
        END DO
    END SUBROUTINE LAPLACE1
    
    SUBROUTINE LAPLACE2(DP, PB, TT)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: DP(:, :, :), PB(:, :, :), TT
        REAL DPDY2, DPDYL(3), DPDX(3), DPDZ(3), DPDY
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
                    
                    DPDY2 = DP(I, JP(J), K) * DY2H(1, J) + DP(I, J, K) * DY2H(2, J) + DP(I, JM(J), K) * DY2H(3, J)
                    DY2FAC = (1 + PHI2(XC, YC, ZC, TT)) ** 2
                    
                    LP2(I, J, K) = DY2FAC * DPDY2
                END DO
                J = 1
                YC = (Y(0) + Y(1)) / 2
                !PRINT*, H(2)
                DPDY2 = ((DP(I, 2, K) - DP(I, 1, K)) / H(2) - PB(I, 1, K)) / DY(1)
                DY2FAC = (1 + PHI2(XC, YC, ZC, TT)) ** 2 
                LP2(I, J, K) = DY2FAC * DPDY2
                
                J = N2
                YC = (Y(N2) + Y(N2 - 1)) / 2
                DPDY2 = (PB(I, 2, K) - (DP(I, N2, K) - DP(I, N2 - 1, K)) / H(N2)) / DY(N2)
                DY2FAC = (1 + PHI2(XC, YC, ZC, TT)) ** 2
                LP2(I, J, K) = DY2FAC * DPDY2
            END DO
        END DO
        
    END SUBROUTINE LAPLACE2
    
    SUBROUTINE TEST_PSOLVE(DP, PB, RP, DPS, PBS, DIVS, TT)
        IMPLICIT NONE
        REAL, INTENT(IN) :: DP(:, :, :), PB(:, :, :), TT, RP(:, :, :)
        COMPLEX, INTENT(IN) :: DPS(:, :, :), PBS(:, :, :), DIVS(:, :, :)
        
        COMPLEX, ALLOCATABLE :: LPPS(:, :, :), RPS(:, :, :), PBTS(:, :, :), PTS(:, :, :)
        REAL, ALLOCATABLE :: PT(:, :, :), PBT(:, :, :)
        ALLOCATE(LPPS(N1, N2, N3))
        ALLOCATE(RPS(N1, N2, N3))
        ALLOCATE(PT(N1, N2, N3))
        ALLOCATE(PBT(N1, 2, N3))
        ALLOCATE(PBTS(N1, 2, N3))
        ALLOCATE(PTS(N1, N2, N3))
        
        CALL LAPLACE(DP, PB, TT)
        CALL OUTPUT('LP', LP)
        CALL OUTPUT('RP', RP)
        LP = LP - RP
        CALL OUTPUT('LP_RP', LP)
        PRINT*, MAXVAL(LP), MAXLOC(LP)
        
        CALL FFT(LP, LPPS)
        CALL FFT(RP, RPS)
        CALL LAPLACES(DPS, PBS, TT)
        CALL IFFT(LPS, LP)
        CALL OUTPUT('LPSP', LP)
        LPS = LPS - RPS
        CALL IFFT(LPS, LP)
        CALL OUTPUT('LPS_RPS', LP)
        LP = ABS(LPS)
        CALL OUTPUT('LPS', LP)
        PRINT*, MAXVAL(ABS(LPS)), MAXLOC(ABS(LPS))
        !LPPS = LPS - LPPS
        !CALL IFFT(LPPS, LP)
        !CALL OUTPUT('LPS_LPPS', LP)
        PRINT*, MAXVAL(ABS(LPPS)), MAXLOC(ABS(LPPS))
        RPS = RPS - DIVS
        PRINT*, MAXVAL(ABS(RPS)), MAXLOC(ABS(RPS))
        CALL IFFT(RPS, LP)
        CALL OUTPUT('RPS_DIVS', LP)
        !LPS = LPS - DIVS
        !PRINT*, MAXVAL(ABS(LPS)), MAXLOC(ABS(LPS))
        
        CALL GET_TEST_PCASE(PT, PBT, TT)
        CALL OUTPUT('PT', PT)
        CALL FFT(PT, PTS)
        CALL LAPLACE(PT, PBT, TT)
        !LP(:, 1, :) = 0
        !LP(:, N2, :) = 0
        CALL LAPLACE1(PT, PBT, TT)
        !LP1(:, N2, :) = 0
        CALL OUTPUT('LP1', LP1)
        CALL LAPLACE11(PT, PBT, TT)
        PT = LPX
        CALL LAPLACE11(PT, PBT, TT)
        CALL OUTPUT('LPX', LPX)
        !CALL LAPLACE2(PT, PBT, TT)
        !CALL OUTPUT('LP2', LP2)
        !CALL OUTPUT('LP', LP)
        !CALL IFFT(LPS, LP)
        !CALL OUTPUT('LPSP', LP)
        !CALL EXPORT_PHI
        !CALL FFT(LP, LPPS)
        !CALL FFT(PBT, PBTS)
        !CALL LAPLACES(PTS, PBTS, TT)
        !LPPS = LPS - LPPS
        !PRINT*, MAXVAL(ABS(LPPS)), MAXLOC(ABS(LPPS))
        !DEALLOCATE(LPPS, RPS)
    END SUBROUTINE TEST_PSOLVE
    
    SUBROUTINE TEST_DIV
        IMPLICIT NONE
        
        CALL GET_TEST_VELCASE
        CALL GETDIV(U, V, W, DIV, 1.0)
        CALL OUTPUT('DIV', DIV)
        !CALL OUTPUT('U', U)
        !CALL OUTPUT('V', V)
        !CALL OUTPUT('W', W)
    END SUBROUTINE TEST_DIV
    
    SUBROUTINE GET_TEST_PCASE(PT, PBT, TT)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: TT
        REAL, INTENT(OUT) :: PT(:, :, :), PBT(:, :, :)
        
        INTEGER I, J, K
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    PT(I, J, K) = GETY(I, J, K, TT)
                END DO
            END DO
        END DO
        
        DO K = 1, N3
            DO I = 1, N1
                PBT(I, 1, K) = 1 / (1 + PHI2((X(I) + X(I - 1)) / 2, (Y(0) + Y(1)) / 2, (Z(K) + Z(K - 1)) / 2, TT))
                PBT(I, 2, K) = 1 / (1 + PHI2((X(I) + X(I - 1)) / 2, (Y(N2) + Y(N2 - 1)) / 2, (Z(K) + Z(K - 1)) / 2, TT))
            END DO
        END DO
        
    END SUBROUTINE GET_TEST_PCASE
    
    SUBROUTINE GET_TEST_VELCASE
        IMPLICIT NONE
        
        REAL :: TT = 1.0
        INTEGER I, J, K
        DO K = 1, N3
            DO J = 1, N2+1
                DO I = 1, N1
                    V(I, J, K) = GETY(I, J, K, TT)
                END DO
            END DO
        END DO
        
        U = 0
        W = 0
    END SUBROUTINE GET_TEST_VELCASE
    
    SUBROUTINE OUTPUT(NAME, VAR, TIME)
        IMPLICIT NONE
        INTEGER I, J, K
        CHARACTER(LEN = *), INTENT(IN) :: NAME
        REAL, INTENT(IN) :: VAR(:, :, :)
        INTEGER, INTENT(IN), OPTIONAL :: TIME
        CHARACTER(LEN = 50) PATH
        CHARACTER(LEN = 50), PARAMETER :: DEBUG_PATH = 'E:\DEBUG\'
        
        IF(PRESENT(TIME)) THEN
            WRITE(PATH, *) TIME
            PATH = TRIM(ADJUSTL(DEBUG_PATH))//TRIM(ADJUSTL(NAME))//TRIM(ADJUSTL(PATH))//'.DAT'
        ELSE
            PATH = TRIM(ADJUSTL(DEBUG_PATH))//TRIM(ADJUSTL(NAME))//'.DAT'
        END IF
        
        OPEN(111, FILE = PATH, STATUS = 'REPLACE')
        WRITE(111, *) 'TITLE = '//TRIM(ADJUSTL(NAME))
        WRITE(111, *) 'VARIABLES = "X" "Y" "Z" "'//TRIM(ADJUSTL(NAME))//'"'
        WRITE(111, *) 'ZONE'
        WRITE(111, *) 'I = ', N1, ', J = ', N2, ', K = ', N3
        WRITE(111, *) 'DATAPACKING = BLOCK'
        WRITE(111, *) ((((X(I) + X(I - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    WRITE(111, *) GETY(I, J, K, 1.0)
                END DO
            END DO
        END DO
        WRITE(111, *) ((((Z(K) + Z(K - 1)) / 2, I = 1, N1), J = 1, N2), K = 1, N3)
        IF(NAME .NE. 'U' .AND. NAME .NE. 'W' .AND. NAME .NE. 'DU' .AND. NAME .NE. 'DW') THEN
            WRITE(111, *) (((VAR(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        ELSE
            WRITE(111, *) (((VAR(I, J, K), I = 1, N1), J = 2, N2+1), K = 1, N3)
        END IF
        
        CLOSE(111)
        CALL SYSTEM('preplot '//PATH)
    END SUBROUTINE OUTPUT
    
    SUBROUTINE EXPORT_PHI()
        IMPLICIT NONE
        
        INTEGER I, J, K
        REAL XC, YC, ZC
        REAL TT
        TT = 1.0
        
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            DO J = 1, N2
                YC = (Y(J) + Y(J - 1)) / 2
                DO I = 1, N1
                    XC = (X(I) + X(I - 1)) / 2
                    
                    PHI_FIELD(I, J, K) = PHI1(XC, YC, ZC, TT)
                END DO
            END DO
        END DO
        
        CALL OUTPUT('PHI1', PHI_FIELD)
        
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            DO J = 1, N2
                YC = (Y(J) + Y(J - 1)) / 2
                DO I = 1, N1
                    XC = (X(I) + X(I - 1)) / 2
                    
                    PHI_FIELD(I, J, K) = DPHI1DX(XC, YC, ZC, TT)
                END DO
            END DO
        END DO
        
        CALL OUTPUT('DPHI1DX', PHI_FIELD)
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            DO J = 1, N2
                YC = (Y(J) + Y(J - 1)) / 2
                DO I = 1, N1
                    XC = (X(I) + X(I - 1)) / 2
                    
                    PHI_FIELD(I, J, K) = DPHI1DY(XC, YC, ZC, TT)
                END DO
            END DO
        END DO
        
        CALL OUTPUT('DPHI1DY', PHI_FIELD)
    END SUBROUTINE EXPORT_PHI
    
    SUBROUTINE GETDIV(UU, VV, WW, DDIV, TT)
        IMPLICIT NONE
        REAL, INTENT(IN) :: UU(:, 0:, :), VV(:, :, :), WW(:, 0:, :), TT
        REAL, INTENT(OUT) :: DDIV(:, :, :)
        REAL XC, YC, ZC
        REAL DUDY, DVDY, DWDY
        REAL U1, U2, U3, W1, W2, W3
        INTEGER I, J, K
        
        DO K = 1, N3
            ZC = (Z(K - 1) + Z(K)) / 2
            DO J = 1, N2
                YC = (Y(J - 1) + Y(J)) / 2
                DO I = 1, N1
                    XC = (X(I - 1) + X(I)) / 2
                    !DUDY = ((UU(I, JP(J), K) + UU(IP(I), JP(J), K)) / 2   &
                    !     -  (UU(I, J, K) + UU(IP(I), J, K)) / 2) / H(JP(J))
                    !DVDY = (VV(I, JP(J), K) - VV(I, J, K)) / DY(J)
                    !DWDY = ((WW(I, JP(J), K) + WW(I, JP(J), KP(K))) / 2   &
                    !     -  (WW(I, J, K) + WW(I, J, KP(K))) / 2) / H(JP(J))
                    U1 = (UU(I, JP(J), K) + UU(IP(I), JP(J), K)) / 2
                    U2 = (UU(I, J, K) + UU(IP(I), J, K)) / 2
                    U3 = (UU(I, JM(J), K) + UU(I, JM(J), K)) / 2
                    DUDY = U1 * DYH(1, J) + U2 * DYH(2, J) + U3 * DYH(3, J)
                    W1 = (WW(I, JP(J), K) + WW(I, JP(J), KP(K))) / 2
                    W2 = (WW(I, J, K) + WW(I, J, KP(K))) / 2
                    W3 = (WW(I, JM(J), K) + WW(I, JM(J), KP(K))) / 2
                    DWDY = W1 * DYH(1, J) + W2 * DYH(2, J) + W3 * DYH(3, J)
                    !PRINT*, DY
                    !PRINT*, VV(I, :, K)
                    DVDY = (VV(I, JP(J), K) - VV(I, J, K)) / DY(J)
                    DDIV(I, J, K) = (UU(IP(I), J, K) - UU(I, J, K)) / DX    &
                                  + (VV(I, JP(J), K) - VV(I, J, K)) / DY(J) &
                                  + (WW(I, J, KP(K)) - WW(I, J, K)) / DZ    &
                                  + PHI1(XC, YC, ZC, TT) * DUDY            &
                                  + PHI2(XC, YC, ZC, TT) * DVDY            &
                                  + PHI3(XC, YC, ZC, TT) * DWDY
                END DO
            END DO
        END DO
    END SUBROUTINE GETDIV
    
    END MODULE DEBUG