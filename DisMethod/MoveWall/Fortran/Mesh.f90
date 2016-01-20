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
        IF(.NOT. ALL_ALLOCATED)
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
    END SUBROUTINE SET_MESH
    
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
    END SUBROUTINE SET_MESH
    
    SUBROUTINE UPDATE_MESH()
        IMPLICIT NONE
        
        
    
    SUBROUTINE DEL_MESH()
        IMPLICIT NONE
        DEALLOCATE(X, Y, Z)
        DEALLOCATE(HY, DY)
        DEALLOCATE(IM, IP, JM, JP, KM, KP)
        ALL_ALLOCATED = .FALSE.
    END SUBROUTINE DEL_MESH