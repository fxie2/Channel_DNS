MODULE MESH
    USE PROFILE
    IMPLICIT NONE
    
    !VERTICAL POSITION
    REAL, SAVE, ALLOCATABLE :: X(:), Y(:), Z(:)
    
    !CELL LENGTH
    REAL, SAVE, ALLOCATABLE :: H(:), DY2H(:, :), DYH(:, :)
    REAL, SAVE, ALLOCATABLE :: DY(:), DY2DY(:, :), DYDY(:, :)
    REAL, SAVE :: DX, DZ
    
    !WAVE AMPTITUDE
    REAL, SAVE :: UP_WAVE_AMPX, UP_WAVE_AMPZ
    REAL, SAVE :: DN_WAVE_AMPX, DN_WAVE_AMPZ
    REAL, SAVE :: DDT_UP_AMPX, DDT_UP_AMPZ
    REAL, SAVE :: DDT_DN_AMPX, DDT_DN_AMPZ
    
    !RANK
    INTEGER, SAVE, ALLOCATABLE :: IM(:), IP(:)
    INTEGER, SAVE, ALLOCATABLE :: JM(:), JP(:)
    INTEGER, SAVE, ALLOCATABLE :: KM(:), KP(:)
    
    !MESH SCALE IN Y DIRECTION
    REAL, SAVE :: SCALE
    
    LOGICAL, SAVE :: ALL_ALLOCATED = .FALSE.
    
    PRIVATE ALL_ALLOCATED
    
    CONTAINS
    
    SUBROUTINE NEW_MESH()
        IMPLICIT NONE
        IF(.NOT. ALL_ALLOCATED) THEN
            ALLOCATE(X(0:N1))
            ALLOCATE(Y(0:N2))
            ALLOCATE(Z(0:N3))
            ALLOCATE(H(N2+1))
            ALLOCATE(DYH(3, N2))
            ALLOCATE(DY2H(3, N2))
            ALLOCATE(DY(0:N2+1))
            ALLOCATE(DY2DY(3, N2))
            ALLOCATE(DYDY(3, N2))
            ALLOCATE(IM(N1))
            ALLOCATE(IP(N1))
            ALLOCATE(JM(N2))
            ALLOCATE(JP(N2))
            ALLOCATE(KM(N3))
            ALLOCATE(KP(N3))
        ALL_ALLOCATED = .TRUE.
        END IF
    END SUBROUTINE NEW_MESH
    
    SUBROUTINE INIT_MESH()
        IMPLICIT NONE
        
        INTEGER I, J
        
        IF(.NOT. ALL_ALLOCATED) CALL NEW_MESH()
        DX = LX / N1
        DZ = LZ / N3
        X(0) = 0
        Z(0) = 0
        FORALL (I = 1 : N1)
            X(I) = DX * I
        END FORALL
        FORALL (I = 1 : N3)
            Z(I) = DZ * I
        END FORALL
        
        UP_WAVE_AMPX = MAX_UP_AMPX
        UP_WAVE_AMPZ = MAX_UP_AMPZ
        DN_WAVE_AMPX = MAX_DN_AMPX
        DN_WAVE_AMPZ = MAX_DN_AMPZ
        
        !OPEN(100, FILE=MESH_PATH, STATUS = 'OLD')
        !READ(100, *) (Y(I), I = 0, N2)
        !CLOSE(100)

        DO I = 1, N2
            Y(I) = I * 2.0 / N2
        END DO
        Y(0) = 0
        
        !SCALE MESH TO Y IN [-1, 1]
        SCALE = 2 / (Y(N2) - Y(0))
        X = X * SCALE
        Y = Y * SCALE
        Z = Z * SCALE
        LY = 2
        LX = LX * SCALE
        LZ = LZ * SCALE
        Y = Y - 1
        
        DO I = 1, N2
            DY(I) = Y(I) - Y(I-1)
        END DO
        DY(0) = 0
        DY(N2 + 1) = 0
        DO I = 1, N2 + 1
            H(I) = (DY(I-1) + DY(I)) / 2
        END DO
        FORALL (I = 1 : N1)
            IM(I) = I - 1
            IP(I) = I + 1
        END FORALL
        FORALL (I = 1 : N2)
            JM(I) = I - 1
            JP(I) = I + 1
        END FORALL
        FORALL ( I = 1 : N3)
            KM(I) = I - 1
            KP(I) = I + 1
        END FORALL
        IM(1) = N1
        IP(N1) = 1
        KM(1) = N3
        KP(N3) = 1
        DO J = 1, N2
            DY2DY(1, J) = 1 / H(J) / DY(J)
            DY2DY(2, J) = -(1 / DY(J) + 1 / DY(JM(J))) / H(J)
            DY2DY(3, J) = 1 / H(J) / DY(JM(J))
            DYDY(1, J) = DY(J) * DY(JM(J)) / (DY(J) + DY(JM(J))) / DY(J) / DY(J)
            DYDY(2, J) = DY(J) * DY(JM(J)) / (DY(J) + DY(JM(J))) * (1 / DY(JM(J)) / DY(JM(J)) - 1 / DY(J) / DY(J))
            DYDY(3, J) = DY(J) * DY(JM(J)) / (DY(J) + DY(JM(J))) / DY(JM(J)) / DY(JM(J)) * -1.0
            DY2H(1, J) = 2 / (H(J) + H(JP(J))) / H(JP(J))
            DY2H(2, J) = 2 / (H(J) + H(JP(J))) * -(1 / H(JP(J)) + 1 / H(J))
            DY2H(3, J) = 2 / (H(J) + H(JP(J))) / H(J)
            !DY2H(1, J) = 1.0 / H(JP(J)) / DY(J)
            !DY2H(2, J) = -(H(JP(J)) + H(J)) / H(JP(J)) / H(J) / DY(J)
            !DY2H(3, J) = 1.0 / H(J) / DY(J)
            DYH(1, J) = H(J) * H(JP(J)) / (H(J) + H(JP(J))) / H(JP(J)) / H(JP(J))
            DYH(2, J) = H(J) * H(JP(J)) / (H(J) + H(JP(J))) * (1 / H(J) / H(J) - 1 / H(JP(J)) / H(JP(J)))
            DYH(3, J) = H(J) * H(JP(J)) / (H(J) + H(JP(J))) / H(J) / H(J) * -1.0
        END DO
        J = 1
        DY2H(1, J) = 2 / (H(J) + H(JP(J))) / H(JP(J))
        DY2H(2, J) = 2 / (H(J) + H(JP(J))) * -(1 / H(JP(J)) + 1 / H(J))
        DY2H(3, J) = 2 / (H(J) + H(JP(J))) / H(J)
        J = N2
        DY2H(1, J) = 2 / (H(J) + H(JP(J))) / H(JP(J))
        DY2H(2, J) = 2 / (H(J) + H(JP(J))) * -(1 / H(JP(J)) + 1 / H(J))
        DY2H(3, J) = 2 / (H(J) + H(JP(J))) / H(J)        
    END SUBROUTINE INIT_MESH
    
    
    SUBROUTINE DEL_MESH()
        IMPLICIT NONE
        DEALLOCATE(X, Y, Z)
        DEALLOCATE(H, DYH, DY2H)
        DEALLOCATE(DY, DYDY, DY2DY)
        DEALLOCATE(IM, IP, JM, JP, KM, KP)
        ALL_ALLOCATED = .FALSE.
    END SUBROUTINE DEL_MESH
        
    REAL FUNCTION PHI1(X, Y, Z, T)
        IMPLICIT NONE
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADX, DETA0DX
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        DETA0DX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                +  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        PHI1 = -(Y * DETADX + DETA0DX) / (1 + ETA)

    END FUNCTION PHI1
    
    REAL FUNCTION PHI2(X, Y, Z, T)
        IMPLICIT NONE
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        PHI2 = 1 / (1 + ETA) - 1
    END FUNCTION PHI2
    
    REAL FUNCTION PHI3(X, Y, Z, T)
        IMPLICIT NONE
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADZ, DETA0DZ
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETA0DZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                +  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        PHI3 = -(Y * DETADZ + DETA0DZ) / (1 + ETA)
    END FUNCTION PHI3
    
    REAL FUNCTION PHIT(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADT, DETA0DT
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADT = (UP_WAVE_AMPX * UP_WAVE_PSDX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               +  UP_WAVE_AMPZ * UP_WAVE_PSDZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_PSDX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_PSDZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETA0DT = (UP_WAVE_AMPX * UP_WAVE_PSDX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                +  UP_WAVE_AMPZ * UP_WAVE_PSDZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                +  DN_WAVE_AMPX * DN_WAVE_PSDX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)    &
                +  DN_WAVE_AMPZ * DN_WAVE_PSDZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        PHIT = (Y * DETADT + DETA0DT) / (1 + ETA)
    END FUNCTION PHIT
    
    REAL FUNCTION DPHI1DY(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADX
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2

        DPHI1DY = -DETADX / (1 + ETA)
    END FUNCTION DPHI1DY
    
    REAL FUNCTION DPHI2DY(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        
        DPHI2DY = 0
    END FUNCTION DPHI2DY
    
    REAL FUNCTION DPHI3DY(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        REAL ETA, DETADZ
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2

        DPHI3DY = -DETADZ / (1 + ETA)
    END FUNCTION DPHI3DY
    
    REAL FUNCTION DPHI1DX(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        REAL ETA, DETADX, DETA0DX, D2ETADX, D2ETA0DX
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        DETA0DX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                +  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        D2ETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * UP_WAVE_NUMX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                -  DN_WAVE_AMPX * DN_WAVE_NUMX * DN_WAVE_NUMX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / -2
        
        D2ETA0DX = (UP_WAVE_AMPX * UP_WAVE_NUMX * UP_WAVE_NUMX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)   &
                 +  DN_WAVE_AMPX * DN_WAVE_NUMX * DN_WAVE_NUMX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / -2
        
        DPHI1DX = -((Y * D2ETADX + D2ETA0DX) * (1 + ETA) - DETADX * (Y * DETADX + DETA0DX)) / (1 + ETA) / (1 + ETA)
    END FUNCTION DPHI1DX
    
    REAL FUNCTION DPHI3DZ(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        REAL ETA, DETADZ, DETA0DZ, D2ETADZ, D2ETA0DZ
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETA0DZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                +  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2

        D2ETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * UP_WAVE_NUMZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * DN_WAVE_NUMZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / -2
        
        D2ETA0DZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * UP_WAVE_NUMZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)   &
                 +  DN_WAVE_AMPZ * DN_WAVE_NUMZ * DN_WAVE_NUMZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / -2
        
        DPHI3DZ = -((Y * D2ETADZ + D2ETA0DZ) * (1 + ETA) - DETADZ * (Y * DETADZ + DETA0DZ)) / (1 + ETA) / (1 + ETA)
    END FUNCTION DPHI3DZ
    
    REAL FUNCTION GETETA(X, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Z, T
        
        GETETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
               +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
               -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
               -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
    END FUNCTION GETETA
    
    REAL FUNCTION GETETA0(X, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Z, T
        
        GETETA0 = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
                +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
                +  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)  &
                +  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
    END FUNCTION GETETA0
    
    REAL FUNCTION GETX(I, J, K, TT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I, J, K
        REAL, INTENT(IN) :: TT
        
        GETX = I * DX - DX / 2
    END FUNCTION GETX
    
    REAL FUNCTION GETY(I, J, K, TT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I, J, K
        REAL, INTENT(IN) :: TT
        
        REAL ETA, ETA0, XC, YC, ZC
        
        XC = I * DX - DX / 2
        ZC = K * DZ - DZ / 2
        YC = (Y(J) + Y(JM(J))) / 2
        
        ETA = GETETA(XC, ZC, TT)
        ETA0 = GETETA0(XC, ZC, TT)
        
        GETY = YC * (1 + ETA) + ETA0
    END FUNCTION GETY
    
    REAL FUNCTION GETZ(I, J, K, TT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I, J, K
        REAL, INTENT(IN) :: TT
        
        GETZ = K * DZ - DZ / 2
    END FUNCTION GETZ
    
END MODULE MESH