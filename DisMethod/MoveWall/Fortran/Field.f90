MODULE FIELD
    USE GLOBAL_PARAMETER
    USE MESH
    IMPLICIT NONE
    
    REAL(8), SAVE, ALLOCATABLE :: U(:, :, :), V(:, :, :), W(:, :, :), P(:, :, :)
    REAL(8), SAVE, ALLOCATABLE :: DU(:, :, :), DV(:, :, :), DW(:, :, :), DP(:, :, :)
    
    REAL(8), SAVE :: PGX, PGY, PGZ
    REAL(8), SAVE :: XFLOW, ZFLOW
    REAL(8), SAVE :: DIVMAX
    REAL(8), SAVE :: CFLMAX
    
    LOGICAL, SAVE :: ALL_ALLOCATED = .FALSE.
    REAL(8), SAVE, ALLOCATABLE :: YSUM(:)
    
    PRIVATE YSUM
    
    CONTAINS
    
    SUBROUTINE ALLOC_FIELD()
        IMPLICIT NONE
        ALLOCATE(U(N1, 0:N2+1, N3))
        ALLOCATE(V(N1, 1:N2+1, N3))
        ALLOCATE(W(N1, 0:N2+1, N3))
        ALLOCATE(P(N1, N2, N3))
        ALLOCATE(DU(N1, 0:N2+1, N3))
        ALLOCATE(DV(N1, 1:N2+1, N3))
        ALLOCATE(DW(N1, 0:N2+1, N3))
        ALLOCATE(DP(N1, N2, N3))
        ALLOCATE(YSUM(N2))
        ALL_ALLOCATED = .TRUE.
    END SUBROUTINE ALLOC_FIELD
    
    SUBROUTINE DEALLOC_FIELD()
        IMPLICIT NONE
        DEALLOCATE(U, V, W, P)
        DEALLOCATE(DU, DV, DW, DP)
        DEALLOCATE(YSUM)
        ALL_ALLOCATED = .FALSE.
    END SUBROUTINE DEALLOC_FIELD
    
    SUBROUTINE INIUP()
        IMPLICIT NONE
        
        REAL V1M, V2M, V3M
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
            DO K = 1, N3
                DO J = 1, N2
                    V1M = V1M + U(I, J, K) * DY(J) * DZ
                END DO
            END DO
            V1M = V1M / LY / LZ
            U(I, :, :) = U(I, :, :) - V1M
        END DO
        
        !V DIRECTION
        DO J = 2, N2
            V2M = 0
            DO K = 1, N3
                DO I = 1, N1
                    V2M = V2M + V(I, J, K) * DX * DZ
                END DO
            END DO
            V2M = V2M / LX / LZ
            V(:, J, :) = V(:, J, :) - V2M
        END DO
        
        !W DIRECTION
        DO K = 1, N3
            V3M = 0
            DO J = 1, N2
                DO I = 1, N1
                    V3M = V3M + W(I, J, K) * DY(J) * DX
                END DO
            END DO
            V3M = V3M / LX / LY
            W(:, :, K) = W(:, :, K) - V3M
        END DO
        
        U = U * INIT_TURB_INTENSITY * 2
        V = V * INIT_TURB_INTENSITY * 2
        W = W * INIT_TURB_INTENSITY * 2
        
        !IMPOSE LAMINAR VELOCITY PROFIELS IN U VELOCITIES
        DO J = 1, N2
            YH = (Y(J) + Y(JP(J))) / 2
            U(:, J, :) = U(:, J, :) + 1 - YH * YH
        END DO
        
        RFLOW = 4 / 3 !INTEGRAL OF U ALONG Y
        CALL CHECK_FLOW_RATE()
        U = RFLOW / XFLOW * U
        W = W - FLOW3
        
        !IMPOSE ZERO-PRESSURE FLUCTUATIONS
        P = 0
        
        !INITIAL MEAN PRESSURE GRADIENT AT LAMINAR FLOW FIELD
        PRX = -2 / RE
        PRZ = 0
    END SUBROUTINE INIUP
    
    SUBROUTINE CHECK_FLOW_RATE()
        IMPLICIT NONE
        
        YSUM = SUM(SUM(U, 1), 2)
        XFLOW = DOT_PRODUCT(YSUM, DY) / N1 / N3
        
        YSUM = SUM(SUM(W, 1), 2)
        ZFLOW = DOT_PRODUCT(YSUM, DY) / N1 / N3
    END SUBROUTINE CHECK_FLOW_RATE
        
    SUBROUTINE CHECK_DIV()
        IMPLICIT NONE
        
        REAL DIV
        REAL XC, YC, ZC
        REAL DUDY, DVDY, DWDY
        INTEGER I, J, K
        
        DIVMAX = 0
        DO K = 1, N3
            ZC = (Z(K) + Z(KP(K))) / 2
            DO J = 1, N2
                YC = (Y(J) + Y(JP(J))) / 2
                DO I = 1, N1
                    XC = (X(I) + X(IP(I))) / 2
                    DUDY = ((U(I, JP(J), K) + U(IP(I), JP(J), K)) / 2   &
                         -  (U(I, J, K) + U(IP(I), J, K)) / 2) / HY(JP(J))
                    DVDY = (V(I, JP(J), K) - V(I, J, K)) / DY(J)
                    DWDY = ((W(I, JP(J), K) + W(I, JP(J), KP(K))) / 2   &
                         -  (W(I, J, K) + W(I, J, KP(K))) / 2) / HY(JP(J))
                    DIV = (U(IP(I), J, K) - U(I, J, K)) / DX    &
                        + (V(I, JP(J), K) - V(I, J, K)) / DY(J) &
                        + (W(I, J, KP(K)) - W(I, J, K)) / DZ    &
                        + PHI1(XC, YC, ZC, T) * DUDY            &
                        + PHI2(XC, YC, ZC, T) * DVDY            &
                        + PHI3(XC, YC, ZC, T) * DWDY
                    DIVMAX = MAX(DIV, DIVMAX)
                END DO
            END DO
        END DO
    END SUBROUTINE CHECK_DIV
    
    SUBROUTINE CHECK_CFL()
        IMPLICIT NONE
        
        REAL CFL
        REAL ETA, ETA0
        INTEGER I, J, K
        
        CFLMAX = 0
        DO K = 1, N3
            DO J = 1, N2
                DO I = 1, N1
                    ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X(I) - UP_WAVE_PSDX * T)  &
                        +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z(K) - UP_WAVE_PSDZ * T)  &
                        -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X(I) - DN_WAVE_NUMX * T)  &
                        -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z(K) - DN_WAVE_NUMZ * T)) / 2
                    
                    ETA0 = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X(I) - UP_WAVE_PSDX * T)  &
                         +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z(K) - UP_WAVE_PSDZ * T)  &
                         -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X(I) - DN_WAVE_NUMX * T)  &
                         -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z(K) - DN_WAVE_NUMZ * T)) / 2
                    
                    CFL = ABS(U(I, J, K) + U(IP(I), J, K)) * 0.5 / DX   &
                        + ABS(V(I, J, K) + V(I, JP(J), K)) * 0.5 / (DY(J) * (1 + ETA) + ETA0)   &
                        + ABS(W(I, J, K) + W(I, J, KP(K))) * 0.5 / DZ
                    
                    CFLMAX = MAX(CFL, CFLMAX)
                END DO
            END DO
        END DO
    END SUBROUTINE CHECK_CFL
    
        
    REAL FUNCTION PHI1(X, Y, Z, T)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADX, DETA0DX
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        DETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        DETA0DX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                +  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        PHI1 = -(Y * DETADX + DETA0DX) / (1 + ETA)
    END FUNCTION PHI1
    
    REAL FUNCTION PHI2(X, Y, Z, T)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: X, Y, Z, T
        
        REAL ETA
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        PHI2 = 1 / (1 + ETA) - 1
    END FUNCTION PHI2
    
    REAL FUNCTION PHI3(X, Y, Z, T)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADZ, DETA0DZ
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        DETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETA0DZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                +  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        PHI3 = -(Y * DETADZ + DETA0DZ) / (1 + ETA)
    END FUNCTION PHI3
    
    REAL FUNCTION PHIT(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL(8), INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADT, DETA0DT
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        DETADT = (UP_WAVE_AMPX * UP_WAVE_PSDX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               +  UP_WAVE_AMPZ * UP_WAVE_PSDZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_PSDX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_PSDZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETA0DT = (UP_WAVE_AMPX * UP_WAVE_PSDX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                +  UP_WAVE_AMPZ * UP_WAVE_PSDZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                +  DN_WAVE_AMPX * DN_WAVE_PSDX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)    &
                +  DN_WAVE_AMPZ * DN_WAVE_PSDZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        PHIT = -(Y * DETADT + DTEA0DT) / (1 + ETA)
    END FUNCTION PHIT
    
    END MODULE FIELD