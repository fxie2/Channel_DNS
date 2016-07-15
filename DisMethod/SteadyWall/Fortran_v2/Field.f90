MODULE FIELD
    USE GLOBAL_PARAMETER
    USE MESH
    USE MATH
    IMPLICIT NONE
    
    REAL, SAVE, ALLOCATABLE :: U(:, :, :), V(:, :, :), W(:, :, :), P(:, :, :)
    REAL, SAVE, ALLOCATABLE :: DU(:, :, :), DV(:, :, :), DW(:, :, :), DP(:, :, :)
    COMPLEX, SAVE, ALLOCATABLE :: DIVS(:, :, :), DPS(:, :, :)
    REAL, SAVE, ALLOCATABLE, TARGET :: R(:, :, :)
    REAL, SAVE, POINTER :: R1(:, :, :), R2(:, :, :), R3(:, :, :), RP(:, :, :)
    REAL, SAVE, ALLOCATABLE :: DIV(:, :, :)
    
    REAL, SAVE :: PGX, PGZ
    REAL, SAVE :: XFLOW, ZFLOW
    REAL, SAVE :: DIVMAX
    REAL, SAVE :: CFLMAX
    
    !SUPPLIMENTAL VECTOR
    !FACTOR VECTOR FOR TD & CTD SOLVER
    REAL, SAVE, ALLOCATABLE :: NMFAC(:), NCFAC(:), NPFAC(:)
    COMPLEX, SAVE, ALLOCATABLE :: PMFAC(:), PCFAC(:), PPFAC(:)
    LOGICAL, SAVE :: FAC_ALLOCATED = .FALSE.
    INTEGER, SAVE :: FAC_LEN
    
    LOGICAL, SAVE :: ALL_ALLOCATED = .FALSE.
    
    PRIVATE ALL_ALLOCATED, DU, DV, DW, DP, R, R1, R2, R3, RP
    PRIVATE NMFAC, NCFAC, NPFAC, FAC_ALLOCATED, FAC_LEN
    PRIVATE DIVS, DPS
    
    CONTAINS
    
    SUBROUTINE ALLOC_FIELD()
        IMPLICIT NONE
        ALLOCATE(U(N1, 0:N2+1, N3))
        ALLOCATE(V(N1, 1:N2+1, N3))
        ALLOCATE(W(N1, 0:N2+1, N3))
        ALLOCATE(P(N1, N2, N3))
        ALLOCATE(DIV(N1, N2, N3))
        ALLOCATE(DU(N1, 0:N2+1, N3))
        ALLOCATE(DV(N1, 1:N2+1, N3))
        ALLOCATE(DW(N1, 0:N2+1, N3))
        ALLOCATE(DP(N1, N2, N3))
        ALLOCATE(DPS(N1, N2, N3))
        ALLOCATE(DIVS(N1, N2, N3))
        ALLOCATE(R(N1, N2, N3))
        ALL_ALLOCATED = .TRUE.
    END SUBROUTINE ALLOC_FIELD
    
    SUBROUTINE DEALLOC_FIELD()
        IMPLICIT NONE
        DEALLOCATE(U, V, W, P, DIV)
        DEALLOCATE(DU, DV, DW, DP)
        DEALLOCATE(R)
        ALL_ALLOCATED = .FALSE.
    END SUBROUTINE DEALLOC_FIELD
    
    SUBROUTINE ALLOC_FAC(N)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: N
        
        IF(FAC_ALLOCATED) THEN
            IF(N > FAC_LEN) THEN
                DEALLOCATE(NMFAC, NCFAC, NPFAC)
                DEALLOCATE(PMFAC, PCFAC, PPFAC)
                ALLOCATE(NMFAC(N))
                ALLOCATE(NCFAC(N))
                ALLOCATE(NPFAC(N))
                ALLOCATE(PMFAC(N))
                ALLOCATE(PCFAC(N))
                ALLOCATE(PPFAC(N))
                FAC_LEN = N
            END IF
        ELSE
            ALLOCATE(NMFAC(N))
            ALLOCATE(NCFAC(N))
            ALLOCATE(NPFAC(N))
            ALLOCATE(PMFAC(N))
            ALLOCATE(PCFAC(N))
            ALLOCATE(PPFAC(N))
            FAC_LEN = N
            FAC_ALLOCATED = .TRUE.
        END IF
        NMFAC = 0
        NCFAC = 0
        NPFAC = 0
        PMFAC = 0
        PCFAC = 0
        PPFAC = 0
    END SUBROUTINE ALLOC_FAC
    
    SUBROUTINE DEALLOC_FAC()
        IMPLICIT NONE
        IF(FAC_ALLOCATED) THEN
            DEALLOCATE(NMFAC, NCFAC, NPFAC)
            DEALLOCATE(PMFAC, PCFAC, PPFAC)
            FAC_LEN = 0
        END IF
    END SUBROUTINE DEALLOC_FAC
    
    SUBROUTINE INIUP()
        IMPLICIT NONE
        
        REAL V1M, V2M, V3M  !MEAN FLOW
        REAL S1, S2, S3     !SLICE AREA
        REAL YH
        REAL RFLOW
        
        INTEGER I, J, K
        
        !CHECK IF ALLOCATED
        IF(.NOT. ALL_ALLOCATED) CALL ALLOC_FIELD()
        
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(U)
        CALL RANDOM_NUMBER(V)
        CALL RANDOM_NUMBER(W)
        
        !UNIFORM DISTRIBUTE IN [-0.5, 0.5)
        U = U - 0.5 
        V = V - 0.5
        W = W - 0.5
        
        !IMPOSE ZERO VELOCITY AT BOUNDARY
        U(:, 0, :) = 0
        V(:, 1, :) = 0
        W(:, 0, :) = 0
        U(:, N2+1, :) = 0
        V(:, N2+1, :) = 0
        W(:, N2+1, :) = 0
        
        !ELIMINATE MEAN QUANTITIES OF RANDOM FLUCTUATIONS
        !U DIRECTION
        DO I = 1, N1
            V1M = 0
            S1 = 0
            DO K = 1, N3
                DO J = 1, N2
                    S1 = S1 + DY(J) * DZ
                    V1M = V1M + U(I, J, K) * DY(J) * DZ
                END DO
            END DO
            V1M = V1M / S1
            U(I, 1:N2, :) = U(I, 1:N2, :) - V1M
        END DO
        
        !V DIRECTION
        DO J = 2, N2
            V2M = 0
            S2 = 0
            DO K = 1, N3
                DO I = 1, N1
                    S2 = S2 + DX * DZ
                    V2M = V2M + V(I, J, K) * DX * DZ
                END DO
            END DO
            V2M = V2M / S2
            V(:, J, :) = V(:, J, :) - V2M
        END DO
        
        !W DIRECTION
        DO K = 1, N3
            V3M = 0
            S3 = 0
            DO J = 1, N2
                DO I = 1, N1
                    S3 = S3 + DY(J) * DX
                    V3M = V3M + W(I, J, K) * DY(J) * DX
                END DO
            END DO
            V3M = V3M / S3
            W(:, 1:N2, K) = W(:, 1:N2, K) - V3M
        END DO
        
        U = U * INIT_TURB_INTENSITY * 2
        V = V * INIT_TURB_INTENSITY * 2
        W = W * INIT_TURB_INTENSITY * 2
        
        !IMPOSE LAMINAR VELOCITY PROFIELS IN U VELOCITIES
        DO J = 1, N2
            YH = (Y(J) + Y(JM(J))) / 2
            U(:, J, :) = U(:, J, :) + 1 - YH * YH
        END DO
        
        RFLOW = 4.0 / 3.0 * LZ !INTEGRAL OF U ALONG Y
        CALL CHECK_FLOW_RATE()
        U = RFLOW / XFLOW * U
        W(:, 1:N2, :) = W(:, 1:N2, :) - ZFLOW / LX / LY
        
        !IMPOSE ZERO-PRESSURE FLUCTUATIONS
        P = 0
        
        !INITIAL MEAN PRESSURE GRADIENT AT LAMINAR FLOW FIELD
        PGX = -2 / RE
        PGZ = 0
    END SUBROUTINE INIUP
    
    SUBROUTINE CHECK_FLOW_RATE()
        IMPLICIT NONE
        
        REAL, ALLOCATABLE :: XMF(:), ZMF(:)
        INTEGER I, J, K
        
        ALLOCATE(XMF(N1))
        ALLOCATE(ZMF(N3))
        
        DO I = 1, N1
            XMF(I) = 0
            DO K = 1, N3
                DO J = 1, N2
                    XMF(I) = XMF(I) + U(I, J, K) * DZ * DY(J)
                END DO
            END DO
        END DO
        
        DO K = 1, N3
            ZMF(K) = 0
            DO J = 1, N2
                DO I = 1, N1
                    ZMF(K) = ZMF(K) + W(I, J, K) * DX * DY(J)
                END DO
            END DO
        END DO
        
        XFLOW = SUM(XMF) / N1
        ZFLOW = SUM(ZMF) / N3
        
    END SUBROUTINE CHECK_FLOW_RATE
        
    SUBROUTINE CHECK_DIV()
        IMPLICIT NONE
        
        CALL GETDIV(U, V, W, DIV, T)
        DIVMAX = MAXVAL(ABS(DIV))
    END SUBROUTINE CHECK_DIV
    
    SUBROUTINE GETDIV(UU, VV, WW, DDIV, TT)
        IMPLICIT NONE
        REAL, INTENT(IN) :: UU(:, 0:, :), VV(:, :, :), WW(:, 0:, :), TT
        REAL, INTENT(OUT) :: DDIV(:, :, :)
        INTEGER I, J, K
        
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    DDIV(I, J, K) = (UU(IP(I), J, K) - UU(I, J, K)) / DX    &
                                  + (VV(I, JP(J), K) - VV(I, J, K)) / DY(J) &
                                  + (WW(I, J, KP(K)) - WW(I, J, K)) / DZ
                END DO
            END DO
        END DO
    END SUBROUTINE GETDIV
    
    SUBROUTINE CHECK_CFL()
        IMPLICIT NONE
        
        REAL CFL
        INTEGER I, J, K
        
        CFLMAX = 0
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    
                    CFL = ABS(U(I, J, K) + U(IP(I), J, K)) * 0.5 / DX   &
                        + ABS(V(I, J, K) + V(I, JP(J), K)) * 0.5 / DY(J)   &
                        + ABS(W(I, J, K) + W(I, J, KP(K))) * 0.5 / DZ
                    
                    CFLMAX = MAX(CFL, CFLMAX)
                END DO
            END DO
        END DO
    END SUBROUTINE CHECK_CFL
    
    SUBROUTINE SOLVEUP()
        IMPLICIT NONE
        CALL GETVEL()
        CALL GETPRE()
        CALL UPDATE_UP()
    END SUBROUTINE SOLVEUP
    
    SUBROUTINE GETVEL()
        IMPLICIT NONE
        
        INTEGER :: ITER = 0
        REAL :: ERR = 0
        
        CALL UPDATE_BC()
        
        CALL GETU()
        CALL GETV()
        CALL GETW()
        CALL FINISH_VEL(ERR)
        
        DU = DU + U
        DV = DV + V
        DW = DW + W
        
    END SUBROUTINE GETVEL
    
    SUBROUTINE UPDATE_BC()
        IMPLICIT NONE
        INTEGER I, K
        
        DU = 0
        DV = 0
        DW = 0
        
    END SUBROUTINE UPDATE_BC
    
    SUBROUTINE GETPRE()
        IMPLICIT NONE
        CALL FORM_RP()
        CALL SOLVE_DP()
    END SUBROUTINE GETPRE
    
    SUBROUTINE GETU()
        IMPLICIT NONE
        CALL FORM_R1()
        CALL SOLVE_UH()
    END SUBROUTINE GETU
    
    SUBROUTINE GETV()
        IMPLICIT NONE
        CALL FORM_R2()
        CALL SOLVE_VH()
    END SUBROUTINE GETV
    
    SUBROUTINE GETW()
        IMPLICIT NONE
        CALL FORM_R3()
        CALL SOLVE_WH()
    END SUBROUTINE GETW
    
    SUBROUTINE FINISH_VEL(ERR)
        IMPLICIT NONE
        REAL, INTENT(OUT) :: ERR
        
        REAL V1, V2, W1, W2, DWVDZ, DWVDY
        REAL U1, U2, W_UP, W_MI, W_DN, DUVDY, DUWDZ, DUWDY
        REAL UERR, VERR, WERR
        INTEGER I, J, K
        
        !FINISH DV
        DO K = 1, N3
            DO J = 2, N2
                DO I = 1, N1
                    
                    !D_W^NP_V^N_DZ
                    V1 = (V(I, J, KM(K)) + V(I, J, K)) / 2
                    V2 = (V(I, J, KP(K)) + V(I, J, K)) / 2
                    W1 = (DW(I, J, K) * DY(JM(J)) + DW(I, JM(J), K) * DY(J)) / 2 / H(J)
                    W2 = (DW(I, J, KP(K)) * DY(JM(J)) + DW(I, JM(J), KP(K)) * DY(J)) / 2 / H(J)
                    DWVDZ = (W2 * V2 - W1 * V1) / DZ
                    
                    DV(I, J, K) = DV(I, J, K) - DT * DWVDZ / 2
                END DO
            END DO
        END DO
        
        !FINISH DU
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    
                    !D_U^N_V^N_DY
                    U1 = (U(I, J, K) * DY(JM(J)) + U(I, JM(J), K) * DY(J)) / 2 / H(J)
                    U2 = (U(I, J, K) * DY(JP(J)) + U(I, JP(J), K) * DY(J)) / 2 / H(JP(J))
                    V1 = (DV(IM(I), J, K) + DV(I, J, K)) / 2
                    V2 = (DV(IM(I), JP(J), K) + DV(I, JP(J), K)) / 2
                    DUVDY = (U2 * V2 - U1 * V1) / DY(J)
                    
                    !D_U^N_W^N_DZ
                    U1 = (U(I, J, K) + U(I, J, KM(K))) / 2
                    U2 = (U(I, J, KP(K)) + U(I, J, K)) / 2
                    W1 = (DW(IM(I), J, K) + DW(I, J, K)) / 2
                    W2 = (DW(I, J, KP(K)) + DW(IM(I), J, KP(K))) / 2
                    DUWDZ = (U2 * W2 - U1 * W1) / DZ
                    
                    DU(I, J, K) = DU(I, J, K)   &
                                - DT * DUVDY / 2   &
                                - DT * DUWDZ / 2
                END DO
            END DO
        END DO

        ERR = 0
        
    END SUBROUTINE FINISH_VEL
    
    SUBROUTINE FORM_R1()
        IMPLICIT NONE
        REAL VISCOS, CROSS, NONLIN, PRESSG
        REAL DUDX(3), DUDZ(3), DUDY
        REAL DUUDX, DUVDY, DUWDZ, DUUDY, DUWDY
        REAL U1, U2, V1, V2, W1, W2, W_UP, W_MI, W_DN
        REAL XC, YC, ZC
        REAL VIS1, VIS2, VIS3, U11, U22, U33, FAC1, FAC2, FAC3
        REAL BC_DN, BC_UP, BCOND
        INTEGER I, J, K
        
        !BIND R1 TO R
        R1 => R
        
        DO K = 1, N3
            ZC = (Z(K) + Z(K - 1)) / 2
            DO J = 1, N2
                YC = (Y(J) + Y(J - 1)) / 2
                DO I = 1, N1
                    XC = X(I - 1)
                    
                    !VISCOS TERM
                    VISCOS = (U(IP(I), J, K) - 2 * U(I, J, K) + U(IM(I), J, K)) / DX / DX   &
                           + (DY2H(1, J) * U(I, JP(J), K) + DY2H(2, J) * U(I, J, K) + DY2H(3, J) * U(I, JM(J), K))  &
                           + (U(I, J, KP(K)) - 2 * U(I, J, K) + U(I, J, KM(K))) / DZ / DZ
                    !FORM Nth-TIME-STEP VISCOS TERM
                    VISCOS = VISCOS / RE * DT
                    
                    !NONLINEAR TERM
                    !D_U^N_U^N_DX
                    U1 = (U(IM(I), J, K) + U(I, J, K)) / 2
                    U2 = (U(IP(I), J, K) + U(I, J, K)) / 2
                    DUUDX = (U2 * U2 - U1 * U1) / DX
                    
                    !D_U^N_V^N_DY & D_U^N_U^N_DY
                    U1 = (U(I, J, K) * DY(JM(J)) + U(I, JM(J), K) * DY(J)) / 2 / H(J)
                    U2 = (U(I, J, K) * DY(JP(J)) + U(I, JP(J), K) * DY(J)) / 2 / H(JP(J))
                    V1 = (V(IM(I), J, K) + V(I, J, K)) / 2
                    V2 = (V(IM(I), JP(J), K) + V(I, JP(J), K)) / 2
                    DUVDY = (U2 * V2 - U1 * V1) / DY(J)
                    
                    !D_U^N_W^N_DZ
                    U1 = (U(I, J, K) + U(I, J, KM(K))) / 2
                    U2 = (U(I, J, KP(K)) + U(I, J, K)) / 2
                    W1 = (W(IM(I), J, K) + W(I, J, K)) / 2
                    W2 = (W(I, J, KP(K)) + W(IM(I), J, KP(K))) / 2
                    DUWDZ = (U2 * W2 - U1 * W1) / DZ
                    
                    NONLIN = DUUDX + DUVDY + DUWDZ
                    
                    !BOUNDARY TERM
                    BC_DN = DY2H(3, 1) * DU(I, 0, K) / RE / 2 &
                          + DU(I, 0, K) * (V(I, 1, K) + V(IM(I), 1, K)) / 2 / DY(1) / 2 &
                          + U(I, 0, K) * (DV(I, 1, K) + DV(IM(I), 1, K)) / 2 / DY(1) / 2
                    
                    BC_UP = DY2H(1, N2) * DU(I, N2+1, K) / RE / 2 &
                          - DU(I, N2+1, K) * (V(I, N2+1, K) + V(IM(I), N2+1, K)) / 2 / DY(N2) / 2 &
                          - U(I, N2+1, K) * (DV(I, N2+1, K) + DV(IM(I), N2+1, K)) / 2 / DY(N2) / 2
                    
                    BCOND = ((N2 - J + 1) - MOD(N2 - J + 1, N2)) / DBLE(N2) * BC_DN &
                          + (J - MOD(J, N2)) / DBLE(N2) * BC_UP
                    
                    !PRESSURE TERM
                    PRESSG = (P(I, J, K) - P(IM(I), J, K)) / DX + PGX
                    
                    !FORM R1 TERM
                    R1(I, J, K) = VISCOS - NONLIN * DT - PRESSG * DT - BCOND * DT
                END DO
            END DO
        END DO
    END SUBROUTINE FORM_R1
    
    SUBROUTINE FORM_R2()
        IMPLICIT NONE
        REAL VISCOS, CROSS, NONLIN, PRESSG
        REAL DVDX(3), DVDZ(3), DVDY
        REAL DUVDX, DVVDY, DWVDZ, DUVDY, DWVDY
        REAL U1, U2, V1, V2, W1, W2
        REAL M21U
        REAL BC_DN, BC_UP, BCOND
        INTEGER I, J, K
        
        !BIND R2 TO R
        R2 => R
        DO K = 1, N3
            DO J = 2, N2
                DO I = 1, N1
                    
                    !VISCOS TERM
                    VISCOS = (V(IP(I), J, K) - 2 * V(I, J, K) + V(IM(I), J, K)) / DX / DX   &
                           + (DY2DY(1, J) * V(I, JP(J), K) + DY2DY(2, J) * V(I, J, K) + DY2DY(3, J) * V(I, JM(J), K))   &
                           + (V(I, J, KP(K)) - 2 * V(I, J, K) + V(I, J, KM(K))) / DZ / DZ
                    
                    VISCOS = VISCOS / RE * DT
                    
                    !NONLINEAR TERM
                    !D_U^N_V^N_DX
                    U1 = (U(I, JM(J), K) * DY(J) + U(I, J, K) * DY(JM(J))) / 2 / H(J)
                    U2 = (U(IP(I), JM(J), K) * DY(J) + U(IP(I), J, K) * DY(JM(J))) / 2 / H(J)
                    V1 = (V(I, J, K) + V(IM(I), J, K)) / 2
                    V2 = (V(I, J, K) + V(IP(I), J, K)) / 2
                    DUVDX = (U2 * V2 - U1 * V1) / DX
                    
                    !D_V^N_V^N_DY & D_U^N_V^N_DY
                    U1 = (U(I, JM(J), K) + U(IP(I), JM(J), K)) / 2
                    U2 = (U(I, J, K) + U(IP(I), J, K)) / 2
                    V1 = (V(I, JM(J), K) + V(I, J, K)) / 2
                    V2 = (V(I, JP(J), K) + V(I, J, K)) / 2
                    DVVDY = (V2 * V2 - V1 * V1) / H(J)
                    
                    !D_W^N_V^N_DZ
                    V1 = (V(I, J, KM(K)) + V(I, J, K)) / 2
                    V2 = (V(I, J, KP(K)) + V(I, J, K)) / 2
                    W1 = (W(I, J, K) * DY(JM(J)) + W(I, JM(J), K) * DY(J)) / 2 / H(J)
                    W2 = (W(I, J, KP(K)) * DY(JM(J)) + W(I, JM(J), KP(K)) * DY(J)) / 2 / H(J)
                    DWVDZ = (W2 * V2 - W1 * V1) / DZ
                    
                    NONLIN = DUVDX + DVVDY + DWVDZ
                    
                    !PRESSURE TERM
                    PRESSG = (P(I, J, K) - P(I, JM(J), K)) / H(J)
                    
                    !BOUNDARY TERM
                    BC_DN = DV(I, 1, K) * DY2DY(3, 2) / RE / 2 &
                          + DV(I, 1, K) * (V(I, 1, K) + V(I, 2, K)) / 2 / H(2) / 2
                    !H = H
                    DY2DY = DY2DY
                    BC_UP = DV(I, N2+1, K) * DY2DY(1, N2) / RE / 2 &
                          - DV(I, N2+1, K) * (V(I, N2+1, K) + V(I, N2, K)) / 2 / H(N2) / 2
                    
                    BCOND = ((N2 - J + 1) - MOD(N2 - J + 1, N2)) / DBLE(N2) * BC_DN &
                          + (J - MOD(J, N2)) / DBLE(N2) * BC_UP
                    
                    !M21 TERM
                    
                    !D_V^N_U^NP_DX
                    U1 = (DY(J) * DU(I, JM(J), K) + DY(JM(J)) * DU(I, J, K)) / 2 / H(J)
                    U2 = (DY(J) * DU(IP(I), JM(J), K) + DY(JM(J)) * DU(IP(I), J, K)) / 2 / H(J)
                    V1 = (V(I, J, K) + V(IM(I), J, K)) / 2
                    V2 = (V(IP(I), J, K) + V(I, J, K)) / 2
                    
                    DUVDX = (U2 * V2 - U1 * V1) / DX
                    
                    M21U = DUVDX / 2
                    
                    !FORM R2 TERM
                    R2(I, J, K) = VISCOS &
                                - NONLIN * DT - PRESSG * DT - BCOND * DT - M21U * DT
                END DO
            END DO
        END DO
    END SUBROUTINE FORM_R2
    
    SUBROUTINE FORM_R3
        IMPLICIT NONE
        REAL VISCOS, CROSS, NONLIN, PRESSG
        REAL DWDX(3), DWDZ(3), DWDY
        REAL DUWDX, DVWDY, DWWDZ, DUWDY, DWWDY
        REAL U1, U2, V1, V2, W1, W2, U_UP, U_MI, U_DN
        REAL M31U, M32V
        REAL BC_DN, BC_UP, BCOND
        INTEGER I, J, K
        
        !BIND R3 TO R
        R3 => R
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    
                    !VISCOS TERM
                    VISCOS = (W(IP(I), J, K) - 2 * W(I, J, K) + W(IM(I), J, K)) / DX / DX   &
                           + (DY2H(1, J) * W(I, JP(J), K) + DY2H(2, J) * W(I, J, K) + DY2H(3, J) * W(I, JM(J), K))  &
                           + (W(I, J, KP(K)) - 2 * W(I, J, K) + W(I, J, KM(K))) / DZ / DZ
                    
                    !FORM Nth-TIME-STEP VISCOS TERM
                    VISCOS = VISCOS / RE * DT
                    
                    !NONLINEAR TERM
                    !D_U^N_W^N_DX
                    U1 = (U(I, J, KP(K)) + U(I, J, K)) / 2
                    U2 = (U(IP(I), J, KP(K)) + U(IP(I), J, K)) / 2
                    W1 = (W(IM(I), J, K) + W(I, J, K)) / 2
                    W2 = (W(IP(I), J, K) + W(I, J, K)) / 2
                    DUWDX = (U2 * W2 - U1 * W1) / DX
                    
                    !D_V^N_W^N_DY & D_W^N_W^N_DY
                    V1 = (V(I, J, K) + V(I, J, KM(K))) / 2
                    V2 = (V(I, JP(J), K) + V(I, JP(J), KM(K))) / 2
                    W1 = (W(I, J, K) * DY(JM(J)) + W(I, JM(J), K) * DY(J)) / 2 / H(J)
                    W2 = (W(I, JP(J), K) * DY(J) + W(I, J, K) * DY(JP(J))) / 2 / H(JP(J))
                    DVWDY = (V2 * W2 - V1 * W1) / DY(J)
                    
                    !D_W^N_W^N_DZ
                    W1 = (W(I, J, KM(K)) + W(I, J, K)) / 2
                    W2 = (W(I, J, KP(K)) + W(I, J, K)) / 2
                    DWWDZ = (W2 * W2 - W1 * W1) / DZ
                    
                    NONLIN = DUWDX + DVWDY + DWWDZ
                    
                    !PRESSURE TERM
                    PRESSG = (P(I, J, K) - P(I, J, KM(K))) / DZ + PGZ
                    
                    !BOUNDARY TERM
                    BC_DN = DW(I, 0, K) * DY2H(3, 1) / RE / 2 &
                          + DW(I, 0, K) * (V(I, 1, K) + V(I, 1, KM(K))) / 2 / DY(1) / 2 &
                          + W(I, 0, K) * (DV(I, 1, K) + DV(I, 1, KM(K))) / 2 / DY(1) / 2
                    
                    BC_UP = DW(I, N2+1, K) * DY2H(1, N2) / RE / 2 &
                          - DW(I, N2+1, K) * (V(I, N2+1, K) + V(I, N2+1, KM(K))) / 2 / DY(N2) / 2 &
                          - W(I, N2+1, K) * (DV(I, N2+1, K) + DV(I, N2+1, KM(K))) / 2 / DY(N2) / 2
                    
                    BCOND = ((N2 - J + 1) - MOD(N2 - J + 1, N2)) / DBLE(N2) * BC_DN &
                          + (J - MOD(J, N2)) / DBLE(N2) * BC_UP
                    
                    !M31 TERM
                    
                    !D_W^N_U^NP_DX
                    U1 = (DU(I, J, KP(K)) + DU(I, J, K)) / 2
                    U2 = (DU(IP(I), J, KP(K)) + DU(IP(I), J, K)) / 2
                    W1 = (W(IM(I), J, K) + W(I, J, K)) / 2
                    W2 = (W(I, J, K) + W(IP(I), J, K)) / 2
                    
                    DUWDX = (U2 * W2 - U1 * W1) / DX
                    
                    M31U = DUWDX / 2
                    
                    !M32 TERM
                    
                    !D_W^N_V^NP_DY
                    V1 = (DV(I, J, K) + DV(I, J, KM(K))) / 2
                    V2 = (DV(I, JP(J), K) + DV(I, JP(J), KM(K))) / 2
                    W1 = (W(I, J, K) * DY(JM(J)) + W(I, JM(J), K) * DY(J)) / 2 / H(J)
                    W2 = (W(I, JP(J), K) * DY(J) + W(I, J, K) * DY(JP(J))) / 2 / H(JP(J))
                    
                    DVWDY = (V2 * W2 - V1 * W1) / H(J)
                    
                    M32V = DVWDY / 2
                    
                    !FORM R3 TERM
                    R3(I, J, K) = VISCOS &
                                - NONLIN * DT - PRESSG * DT - BCOND * DT - M31U * DT - M32V * DT
                END DO
            END DO
        END DO
    END SUBROUTINE FORM_R3
    
    SUBROUTINE SOLVE_UH()
        IMPLICIT NONE
        
        INTEGER I, J, K
        REAL U1, U2, V1, V2, W1, W2, PHI, W_UP, W_MI, W_DN
        REAL VISNM, VISNC, VISNP, VISFAC
        REAL MAXDU
        
        CALL ALLOC_FAC(MAX(N1, N2, N3))
        
        !SOLVE IN Y DIRECTION
        DO K = 1, N3
            DO I = 1, N1
                DO J = 1, N2
                   
                    !PHI * D_V^N_U^NP_DY
                    V1 = (V(IM(I), J, K) + V(I, J, K)) / 2
                    V2 = (V(IM(I), JP(J), K) + V(I, JP(J), K)) / 2
                    NPFAC(J) = V2 / 2 / H(JP(J))
                    NCFAC(J) = (V2 * DY(JP(J)) / 2 / H(JP(J)) - V1 * DY(JM(J)) / 2 / H(J)) / DY(J)
                    NMFAC(J) = -V1 / 2 / H(J)
                    
                    !VISCOS TERM
                    !(1 + PHI_I ** 2) * D2_DY2
                    VISNP = DY2H(1, J)
                    VISNC = DY2H(2, J)
                    VISNM = DY2H(3, J)
                    !DY2H = DY2H
                    NPFAC(J) = NPFAC(J) - VISNP / RE
                    NCFAC(J) = NCFAC(J) - VISNC / RE
                    NMFAC(J) = NMFAC(J) - VISNM / RE
                END DO
                NPFAC = NPFAC / 2 * DT
                NCFAC = NCFAC / 2 * DT + 1
                NMFAC = NMFAC / 2 * DT
                
                CALL TDMA(NMFAC, NCFAC, NPFAC, DU(I, 1 : N2, K), R1(I, 1 : N2, K))
            END DO
        END DO
        
        NPFAC = 0
        NCFAC = 0
        NMFAC = 0
        
        !SOLVE IN X DIRECTION
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    
                    !D_U^N_U^NP_DX
                    U1 = (U(IM(I), J, K) + U(I, J, K)) / 2
                    U2 = (U(IP(I), J, K) + U(I, J, K)) / 2
                    NPFAC(I) = U2 / 2 / DX * DT
                    NCFAC(I) = (U2 - U1) / 2 / DX * DT
                    NMFAC(I) = -U1 / 2 / DX * DT
                    
                    !VISCOS TERM
                    NPFAC(I) = NPFAC(I) - 1 / DX / DX / 2 / RE * DT
                    NCFAC(I) = NCFAC(I) + 1 / DX / DX / RE * DT
                    NMFAC(I) = NMFAC(I) - 1 / DX / DX / 2 / RE * DT
                END DO
                NCFAC(1:N1) = NCFAC(1:N1) + 1
                CALL CTDMA(NMFAC, NCFAC, NPFAC, DU(1:N1, J, K), DU(1:N1, J, K))
            END DO
        END DO
        
        NCFAC = 0
        NPFAC = 0
        NMFAC = 0
        
        !SOLVE IN Z DIRECTION
        DO J = 1, N2
            DO I = 1, N1
                DO K = 1, N3
                    
                    !D_W^N_U^NP_DZ
                    W1 = (W(IM(I), J, K) + W(I, J, K)) / 2
                    W2 = (W(I, J, KP(K)) + W(IM(I), J, KP(K))) / 2
                    NPFAC(K) = W2 / 2 / DZ
                    NCFAC(K) = (W2 - W1) / 2 / DZ
                    NMFAC(K) = -W1 / 2 / DZ
                    
                    !VISCOS TERM
                    NPFAC(K) = NPFAC(K) - 1 / DZ / DZ / RE
                    NCFAC(K) = NCFAC(K) + 2 / DZ / DZ / RE
                    NMFAC(K) = NMFAC(K) - 1 / DZ / DZ / RE
                END DO
                NPFAC = NPFAC / 2 * DT
                NCFAC = NCFAC / 2 * DT + 1
                NMFAC = NMFAC / 2 * DT
                CALL CTDMA(NMFAC, NCFAC, NPFAC, DU(I, J, 1:N3), DU(I, J, 1:N3))
            END DO
        END DO
    END SUBROUTINE SOLVE_UH
    
    SUBROUTINE SOLVE_VH()
        IMPLICIT NONE
        
        INTEGER I, J, K
        REAL U1, U2, V1, V2, W1, W2, PHI
        REAL VISNM, VISNC, VISNP, VISFAC
        
        CALL ALLOC_FAC(MAX(N1, N2, N3))
        
        !SOLVE IN Y DIRECTION
        DO K = 1, N3
            DO I = 1, N1
                DO J = 2, N2
                    
                    !D_V^N_V^NP_DY
                    V1 = (V(I, JM(J), K) + V(I, J, K)) / 2
                    V2 = (V(I, JP(J), K) + V(I, J, K)) / 2
                    NPFAC(J) = V2 / 2 / H(J)
                    NCFAC(J) = (V2 - V1) / 2 / H(J)
                    NMFAC(J) = -V1 / 2 / H(J)
                    
                    !VISCOS TERM
                    !(1 + PHI_I ** 2) * D2_DY2
                    VISFAC = 1. / 2 / RE
                    VISNP = VISFAC * DY2DY(1, J)
                    VISNC = VISFAC * DY2DY(2, J)
                    VISNM = VISFAC * DY2DY(3, J)
                    
                    NPFAC(J) = NPFAC(J) - VISNP
                    NCFAC(J) = NCFAC(J) - VISNC
                    NMFAC(J) = NMFAC(J) - VISNM
                END DO
                NPFAC = NPFAC * DT
                NCFAC = NCFAC * DT + 1
                NMFAC = NMFAC * DT
                
                CALL TDMA(NMFAC(2:N2), NCFAC(2:N2), NPFAC(2:N2), DV(I, 2:N2, K), R2(I, 2:N2, K))
            END DO
        END DO
        
        NPFAC = 0
        NCFAC = 0
        NMFAC = 0
        
        !SOLVE IN X DIRECTION
        DO K = 1, N3
            DO J = 2, N2
                DO I = 1, N1
                    
                    !D_U^N_V^NP_DX
                    U1 = (U(I, JM(J), K) * DY(J) + U(I, J, K) * DY(JM(J))) / 2 / H(J)
                    U2 = (U(IP(I), JM(J), K) * DY(J) + U(IP(I), J, K) * DY(JM(J))) / 2 / H(J)
                    NPFAC(I) = U2 / 2 / DX
                    NCFAC(I) = (U2 - U1) / 2 / DX
                    NMFAC(I) = -U1 / 2 / DX
                    
                    !VISCOS TERM
                    NPFAC(I) = NPFAC(I) - 1 / DX / DX / RE
                    NCFAC(I) = NCFAC(I) + 2 / DX / DX / RE
                    NMFAC(I) = NMFAC(I) - 1 / DX / DX / RE
                END DO
                NPFAC = NPFAC / 2 * DT
                NCFAC = NCFAC / 2 * DT + 1
                NMFAC = NMFAC / 2 * DT
                CALL CTDMA(NMFAC, NCFAC, NPFAC, DV(1:N1, J, K), DV(1:N1, J, K))
            END DO
        END DO
        
        NCFAC = 0
        NPFAC = 0
        NMFAC = 0
        
        !SOLVE IN Z DIRECTION
        DO J = 2, N2
            DO I = 1, N1
                DO K = 1, N3
                    
                    !D_W^N_V^NP_DZ
                    W1 = (W(I, J, K) * DY(JM(J)) + W(I, JM(J), K) * DY(J)) / 2 / H(J)
                    W2 = (W(I, J, KP(K)) * DY(JM(J)) + W(I, JM(J), KP(K)) * DY(J)) / 2 / H(J)
                    NPFAC(K) = W2 / 2 / DZ
                    NCFAC(K) = (W2 - W1) / 2 / DZ
                    NMFAC(K) = -W1 / 2 / DZ
                    
                    !VISCOS TERM
                    NPFAC(K) = NPFAC(K) - 1 / DZ / DZ / RE
                    NCFAC(K) = NCFAC(K) + 2 / DZ / DZ / RE
                    NMFAC(K) = NMFAC(K) - 1 / DZ / DZ / RE
                END DO
                NPFAC = NPFAC / 2 * DT
                NCFAC = NCFAC / 2 * DT + 1
                NMFAC = NMFAC / 2 * DT
                CALL CTDMA(NMFAC, NCFAC, NPFAC, DV(I, J, 1:N3), DV(I, J, 1:N3))
            END DO
        END DO
    END SUBROUTINE SOLVE_VH
    
    SUBROUTINE SOLVE_WH()
        IMPLICIT NONE
        
        INTEGER I, J, K
        REAL U1, U2, V1, V2, W1, W2, PHI, U_UP, U_MI, U_DN
        REAL VISNM, VISNC, VISNP, VISFAC
        
        CALL ALLOC_FAC(MAX(N1, N2, N3))
        
        !SOLVE IN Y DIRECTION
        DO K = 1, N3
            DO I = 1, N1
                DO J = 1, N2
                    
                    !PHI * D_V^N_W^NP_DY
                    V1 = (V(I, J, K) + V(I, J, KM(K))) / 2
                    V2 = (V(I, JP(J), K) + V(I, JP(J), KM(K))) / 2
                    NPFAC(J) = V2 / 2 / H(JP(J))
                    NCFAC(J) = (V2 * DY(JP(J)) / H(JP(J)) - V1 * DY(JM(J)) / H(J)) / 2 / DY(J)
                    NMFAC(J) = -V1 / 2 / H(J)
                    
                    !VISCOS TERM
                    !(1 + PHI_I ** 2) * D2_DY2
                    VISNP = DY2H(1, J)
                    VISNC = DY2H(2, J)
                    VISNM = DY2H(3, J)
                    
                    NPFAC(J) = NPFAC(J) - VISNP / RE
                    NCFAC(J) = NCFAC(J) - VISNC / RE
                    NMFAC(J) = NMFAC(J) - VISNM / RE
                END DO
                NPFAC = NPFAC / 2 * DT
                NCFAC = NCFAC / 2 * DT + 1
                NMFAC = NMFAC / 2 * DT
                
                CALL TDMA(NMFAC, NCFAC, NPFAC, DW(I, 1 : N2, K), R3(I, 1 : N2, K))
            END DO
        END DO

        NPFAC = 0
        NCFAC = 0
        NMFAC = 0
        
        !SOLVE IN X DIRECTION
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    
                    !D_U^N_W^NP_DX
                    U1 = (U(I, J, KP(K)) + U(I, J, K)) / 2
                    U2 = (U(IP(I), J, KP(K)) + U(IP(I), J, K)) / 2
                    NPFAC(I) = U2 / 2 / DX
                    NCFAC(I) = (U2 - U1) / 2 / DX
                    NMFAC(I) = -U1 / 2 / DX
                    
                    !VISCOS TERM
                    NPFAC(I) = NPFAC(I) - 1 / DX / DX / RE
                    NCFAC(I) = NCFAC(I) + 2 / DX / DX / RE
                    NMFAC(I) = NMFAC(I) - 1 / DX / DX / RE
                END DO
                NPFAC = NPFAC / 2 * DT
                NCFAC = NCFAC / 2 * DT + 1
                NMFAC = NMFAC / 2 * DT
                CALL CTDMA(NMFAC, NCFAC, NPFAC, DW(1:N1, J, K), DW(1:N1, J, K))
            END DO
        END DO

        NPFAC = 0
        NCFAC = 0
        NMFAC = 0
        
        !SOLVE IN Z DIRECTION
        DO J = 1, N2
            DO I = 1, N1
                DO K = 1, N3
                    !
                    !!TIME ADVANCE TERM
                    !NCFAC(K) = 1
                    
                    !D_W^N_W^NP_DZ
                    W1 = (W(I, J, KM(K)) + W(I, J, K)) / 2
                    W2 = (W(I, J, KP(K)) + W(I, J, K)) / 2
                    NPFAC(K) = W2 / 2 / DZ * DT
                    NCFAC(K) = (W2 - W1) / 2 / DZ * DT
                    NMFAC(K) = -W1 / 2 / DZ * DT
                    
                    !VISCOS TERM
                    NPFAC(K) = NPFAC(K) - 1 / DZ / DZ / 2 / RE * DT
                    NCFAC(I) = NCFAC(K) + 1 / DZ / DZ / RE * DT
                    NMFAC(I) = NMFAC(K) - 1 / DZ / DZ / 2 / RE * DT
                END DO
                NCFAC(1:N3) = NCFAC(1:N3) + 1
                CALL CTDMA(NMFAC, NCFAC, NPFAC, DW(I, J, 1:N3), DW(I, J, 1:N3))
            END DO
        END DO
    END SUBROUTINE SOLVE_WH
    
    SUBROUTINE FORM_RP()
        IMPLICIT NONE
        
        !BING RP TO R
        RP => R
        CALL GETDIV(DU, DV, DW, RP, T + DT)
        RP = RP / DT
    END SUBROUTINE FORM_RP
    
    SUBROUTINE SOLVE_DP
        IMPLICIT NONE
        REAL :: ALPHA, BETA
        INTEGER I, J, K, M, N
        COMPLEX D2FAC, DFAC, FAC, AVE
        COMPLEX, PARAMETER :: II = (0, 1)
        
        CALL FFT(RP, DIVS)
        ALPHA = 2 * PI / LX
        BETA = 2 * PI / LZ
        DO K = 1, N3
            N = K - 1 - N3 / 2
            DO I = 1, N1
                M = I - 1 - N1 / 2
                DO J = 1, N2
                    
                    !D2_DY2 FACTOR
                    D2FAC = 1
                          
                    !D_DY FACTOR
                    DFAC = 0
                    !Y FACTOR
                    !FAC = -(ALPHA * ALPHA * M * M + BETA * BETA * N * N)
                    FAC = 2 * (COS(ALPHA * M * DX) - 1) / DX / DX + 2 * (COS(BETA * N * DZ) - 1) / DZ / DZ
                        
                    
                    PPFAC(J) = DY2H(1, J) * D2FAC + DYH(1, J) * DFAC
                    PCFAC(J) = DY2H(2, J) * D2FAC + DYH(2, J) * DFAC + FAC
                    PMFAC(J) = DY2H(3, J) * D2FAC + DYH(3, J) * DFAC
                END DO
                
                !PCFAC(1) = -2 / H(2) / (2 * H(1) + H(2)) + FAC
                !PPFAC(1) = 2 / H(2) / (2 * H(1) + H(2))
                !PCFAC(N2) = -2 / H(N2) / (2 * H(N2 + 1) + H(N2)) + FAC
                !PMFAC(N2) = 2 / H(N2) / (H(N2) + 2 * H(N2 + 1))
                PCFAC(N2) = -1 / H(N2) / DY(N2) + FAC
                PMFAC(N2) = 1 / H(N2) / DY(N2)
                PCFAC(1) = -1 / H(2) / DY(1) + FAC
                PPFAC(1) = 1 / H(2) / DY(1)
                
                PMFAC(1) = 0
                PPFAC(N2) = 0
                
                IF(M == 0 .AND. N == 0) THEN
                    PCFAC(1) = 1
                    PPFAC(1) = 0
                    DIVS(I, 1, K) = 0
                END IF
                
                CALL TDMA(PMFAC, PCFAC, PPFAC, DPS(I, :, K), DIVS(I, :, K))
            END DO
        END DO
        
        CALL IFFT(DPS, DP)
    END SUBROUTINE SOLVE_DP
    
    SUBROUTINE UPDATE_UP()
        IMPLICIT NONE
        
        REAL DPGX, DPGZ
        REAL DPDY
        INTEGER I, J, K
        
        !GET N TIME STEP FLOW RATE
        CALL CHECK_FLOW_RATE()
        
        !UPDATE U VELOCITY
        !KEEP CONSTANT MASS FLOW RATE
        DPGX = 0
        DO K = 1, N3
            DO I = 1, N1
                DO J = 1, N2
                    DPGX = DPGX + (DU(I, J, K) * DY(J) * DX * DZ &
                         - DT * (DP(I, J, K) - DP(IM(I), J, K)) * DY(J) * DZ)
                END DO
            END DO
        END DO
        
        DPGX = (DPGX - XFLOW * LX) / LX / LY / LZ / DT

        DO K = 1, N3
            DO I = 1, N1
                DO J = 1, N2
                    U(I, J, K) = DU(I, J, K) - DT * DPGX    &
                               - DT * (DP(I, J, K) - DP(IM(I), J, K)) / DX
                END DO
                !J = 0
                U(I, 0, K) = DU(I, 0, K)
                !J = N2 + 1
                U(I, N2 + 1, K) = DU(I, N2 + 1, K)
            END DO
        END DO
        
        !UPDATE V VELOCITY
        DO K = 1, N3
            DO J = 2, N2
                DO I = 1, N1
                    
                    DPDY = (DP(I, J, K) - DP(I, JM(J), K)) / H(J)
                    V(I, J, K) = DV(I, J, K)    &
                               - DT * DPDY
                END DO
            END DO
        END DO
        
        DO K = 1, N3
            DO I = 1, N1
                V(I, 1, K) = DV(I, 1, K)
                V(I, N2 + 1, K) = DV(I, N2 + 1, K)
            END DO
        END DO
        
        !UPDATE W VELOCITY
        !KEEP 0 FLOW RATE
        DPGZ = 0
        DO K = 1, N3
            DO I = 1, N1
                DO J = 1, N2
                    DPGZ = DPGZ + (DW(I, J, K) * DY(J) * DX * DZ &
                         - DT * (DP(I, J, K) - DP(I, J, KM(K))) * DX * DY(J))
                END DO
            END DO
        END DO
        
        DPGZ = (DPGZ - ZFLOW * LZ) / LX / LY / LZ / DT
        
        DO K = 1, N3
            DO I = 1, N1
                DO J = 1, N2
                    W(I, J, K) = DW(I, J, K) - DT * DPGZ    &
                         - DT * (DP(I, J, K) - DP(I, J, KM(K))) / DZ
                END DO
                !J = 0
                W(I, 0, K) = DW(I, 0, K)
                
                !J = N2 + 1
                W(I, N2 + 1, K) = DW(I, N2 + 1, K)
            END DO
        END DO
        
        !UPDATE PRESSURE
        P = P + DP
        PGX = PGX + DPGX
        PGZ = PGZ + DPGZ
    END SUBROUTINE UPDATE_UP
    
    SUBROUTINE OUTPUT()
        IMPLICIT NONE
        INTEGER I, J, K
        CHARACTER(LEN = 10) NUM
        CHARACTER(LEN = 50) PATH
        
        WRITE(NUM, 110) CURNT_STEP_NUM
110     FORMAT(I10)
        PATH = TRIM(ADJUSTL(SAVE_FILE_PATH))//'U_'//TRIM(ADJUSTL(NUM))//'.DAT'
        OPEN(111, FILE = PATH, FORM = 'BINARY', STATUS = 'REPLACE')
        WRITE(111) (((U(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(111) (((V(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        WRITE(111) (((W(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(111)
        PATH = TRIM(ADJUSTL(SAVE_FILE_PATH))//'P_'//TRIM(ADJUSTL(NUM))//'.DAT'
        OPEN(112, FILE = PATH, FORM = 'BINARY', STATUS = 'REPLACE')
        WRITE(112) (((P(I, J, K), I = 1, N1), J = 1, N2), K = 1, N3)
        CLOSE(112)
    END SUBROUTINE OUTPUT
    
    END MODULE FIELD