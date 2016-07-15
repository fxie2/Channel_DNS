MODULE MESH
    USE GLOBAL_PARAMETER
    IMPLICIT NONE
    
    !VERTICAL POSITION
    REAL, SAVE, ALLOCATABLE :: X(:), Y(:), Z(:)
    
    !CELL LENGTH
    REAL, SAVE, ALLOCATABLE :: H(:), DY2H(:, :), DYH(:, :)
    REAL, SAVE, ALLOCATABLE :: DY(:), DY2DY(:, :), DYDY(:, :)
    REAL, SAVE :: DX, DZ
    
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
        
        OPEN(100, FILE=GRID_FILE_PATH, STATUS = 'OLD')
        READ(100, *) (Y(I), I = 0, N2)
        CLOSE(100)
        !DO I = 0, N2
        !    Y(I) = I / DBLE(N2) * 2.
        !END DO
        
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
            !DY2H(1, J) = 2 / (H(J) + H(JP(J))) / H(JP(J))
            !DY2H(2, J) = 2 / (H(J) + H(JP(J))) * -(1 / H(JP(J)) + 1 / H(J))
            !DY2H(3, J) = 2 / (H(J) + H(JP(J))) / H(J)
            DY2H(1, J) = 1.0 / H(JP(J)) / DY(J)
            DY2H(2, J) = -(H(JP(J)) + H(J)) / H(JP(J)) / H(J) / DY(J)
            DY2H(3, J) = 1.0 / H(J) / DY(J)
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
        
END MODULE MESH