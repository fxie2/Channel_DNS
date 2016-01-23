MODULE MESH
    USE GLOBAL_PARAMETER
    IMPLICIT NONE
    
    !VERTICAL POSITION
    REAL, SAVE, ALLOCATABLE :: X(:), Y(:), Z(:)
    
    !CELL LENGTH
    REAL, SAVE, ALLOCATABLE :: HY(:)
    REAL, SAVE, ALLOCATABLE :: DY(:)
    REAL, SAVE :: DX, DZ
    
    !RANK
    INTEGER, SAVE, ALLOCATABLE :: IM(:), IP(:)
    INTEGER, SAVE, ALLOCATABLE :: JM(:), JP(:)
    INTEGER, SAVE, ALLOCATABLE :: KM(:), KP(:)
    
    !MESH SCALE IN Y DIRECTION
    REAL, SAVE :: SCALE
    
    LOGICAL, SAVE :: ALL_ALLOCATED = .FALSE.
    
    CONTAINS
    
    SUBROUTINE NEW_MESH()
        IMPLICIT NONE
        IF(.NOT. ALL_ALLOCATED) THEN
            ALLOCATE(X(0:N1))
            ALLOCATE(Y(0:N2))
            ALLOCATE(Z(0:N3))
            ALLOCATE(HY(N2+1))
            ALLOCATE(DY(N2))
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
        
        INTEGER I
        
        IF(.NOT. ALL_ALLOCATED) CALL NEW_MESH()
        DX = LX / N1
        DZ = LZ / N3
        X(0) = 0
        Z(0) = 0
        FORALL (I = 1 : N1)
            X(I) = DX * I
            Z(I) = DZ * I
        END FORALL
        OPEN(100, FILE=GRID_FILE_PATH, STATUS = 'OLD')
        READ(100) (Y(I), I = 0, N2)
        CLOSE(100)
        
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
        DO I = 2, N2
            HY(I) = (DY(I-1) + DY(I)) / 2
        END DO
        HY(1) = DY(1) / 2
        HY(N2+1) = DY(N2) / 2
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
        KM(N3) = 1
    END SUBROUTINE INIT_MESH
    
    SUBROUTINE UPDATE_MESH()
        IMPLICIT NONE
    END SUBROUTINE UPDATE_MESH
        
    
    SUBROUTINE DEL_MESH()
        IMPLICIT NONE
        DEALLOCATE(X, Y, Z)
        DEALLOCATE(HY, DY)
        DEALLOCATE(IM, IP, JM, JP, KM, KP)
        ALL_ALLOCATED = .FALSE.
    END SUBROUTINE DEL_MESH
        
    REAL FUNCTION PHI1(X, Y, Z, T)
        IMPLICIT NONE
        REAL, INTENT(IN) :: X, Y, Z, T
        
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
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        PHI2 = 1 / (1 + ETA) - 1
    END FUNCTION PHI2
    
    REAL FUNCTION PHI3(X, Y, Z, T)
        IMPLICIT NONE
        REAL, INTENT(IN) :: X, Y, Z, T
        
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
        
        REAL, INTENT(IN) :: X, Y, Z, T
        
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
        
        PHIT = -(Y * DETADT + DETA0DT) / (1 + ETA)
    END FUNCTION PHIT
    
    REAL FUNCTION DPHI1DY(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        
        REAL ETA, DETADX
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
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
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
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
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        DETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
               -  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        DETA0DX = (UP_WAVE_AMPX * UP_WAVE_NUMX * COS(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                +  DN_WAVE_AMPX * DN_WAVE_NUMX * COS(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / 2
        
        D2ETADX = (UP_WAVE_AMPX * UP_WAVE_NUMX * UP_WAVE_NUMX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)    &
                -  DN_WAVE_AMPX * DN_WAVE_NUMX * DN_WAVE_NUMX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / -2
        
        D2ETA0DX = (UP_WAVE_AMPX * UP_WAVE_NUMX * UP_WAVE_NUMX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)   &
                 -  DN_WAVE_AMPX * DN_WAVE_NUMX * DN_WAVE_NUMX * SIN(DN_WAVE_NUMX * X - DN_WAVE_PSDX * T)) / -2
        
        DPHI1DX = -((Y * D2ETADX + D2ETA0DX) * (1 + ETA) - DETADX * (Y * DETADX + DETA0DX)) / (1 + ETA) / (1 + ETA)
    END FUNCTION DPHI1DX
    
    REAL FUNCTION DPHI3DZ(X, Y, Z, T)
        IMPLICIT NONE
        
        REAL, INTENT(IN) :: X, Y, Z, T
        REAL ETA, DETADZ, DETA0DZ, D2ETADZ, D2ETA0DZ
        
        ETA = (UP_WAVE_AMPX * SIN(UP_WAVE_NUMX * X - UP_WAVE_PSDX * T)  &
            +  UP_WAVE_AMPZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)  &
            -  DN_WAVE_AMPX * SIN(DN_WAVE_NUMX * X - DN_WAVE_NUMX * T)  &
            -  DN_WAVE_AMPZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_NUMZ * T)) / 2
        
        DETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
               -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2
        
        DETA0DZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * COS(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                +  DN_WAVE_AMPZ * DN_WAVE_NUMZ * COS(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / 2

        D2ETADZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * UP_WAVE_NUMZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)    &
                -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * DN_WAVE_NUMZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / -2
        
        D2ETA0DZ = (UP_WAVE_AMPZ * UP_WAVE_NUMZ * UP_WAVE_NUMZ * SIN(UP_WAVE_NUMZ * Z - UP_WAVE_PSDZ * T)   &
                 -  DN_WAVE_AMPZ * DN_WAVE_NUMZ * DN_WAVE_NUMZ * SIN(DN_WAVE_NUMZ * Z - DN_WAVE_PSDZ * T)) / -2
        
        DPHI3DZ = -((Y * D2ETADZ + D2ETA0DZ) * (1 + ETA) - DETADZ * (Y * DETADZ + DETA0DZ)) / (1 + ETA) / (1 + ETA)
    END FUNCTION DPHI3DZ
    
        
END MODULE MESH