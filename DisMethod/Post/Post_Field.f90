MODULE FIELD
    USE PROFILE
    USE MESH
    IMPLICIT NONE
    
    !VELOCITY & PRESSURE FIELD
    REAL, SAVE, ALLOCATABLE :: U(:, :, :), V(:, :, :), W(:, :, :), P(:, :, :)
    REAL, SAVE, ALLOCATABLE :: VELABS(:, :, :), TOLPRE(:, :, :)
    
    !VORTEX FIELD
    REAL, SAVE, ALLOCATABLE :: VORX(:, :, :), VORY(:, :, :), VORZ(:, :, :)
    REAL, SAVE, ALLOCATABLE :: VORABS(:, :, :)
    
    !AVERAGE FIELD
    REAL, SAVE, ALLOCATABLE :: VELAVEX(:, :), VELAVEY(:, :), VELAVEZ(:, :)
    REAL, SAVE, ALLOCATABLE :: VORAVEX(:, :), VORAVEY(:, :), VORAVEZ(:, :)
    
    !FLUCTUATE VELOCITY FIELD
    REAL, SAVE, ALLOCATABLE :: DU(:, :, :), DV(:, :, :), DW(:, :, :)
    REAL, SAVE, ALLOCATABLE :: DURMS(:), DVRMS(:), DWRMS(:), DPRMS(:)
    
    !REYNOLDS STRESS
    REAL, SAVE, ALLOCATABLE :: RESX(:, :, :), RESY(:, :, :), RESZ(:, :, :)
    REAL, SAVE, ALLOCATABLE :: RESABS(:, :, :)
    REAL, SAVE, ALLOCATABLE :: RESAVEX(:, :), RESAVEY(:, :), RESAVEZ(:, :)
    
    !TURBULENT KINECTIC
    REAL, SAVE, ALLOCATABLE :: TK(:, :, :)
    
    !DISSAPATION TERM
    REAL, SAVE, ALLOCATABLE :: EPS(:, :, :)
    
    !ALLOC IDENTIFIER
    LOGICAL, SAVE :: ALL_ALLOC = .FALSE.
    
    PRIVATE ALL_ALLOC
    
    CONTAINS
    
    SUBROUTINE ALLOC_FIELD
        IMPLICIT NONE
        IF(.NOT. ALL_ALLOC) THEN
            ALLOCATE(U(N1, 0:N2+1, N3))
            ALLOCATE(V(N1, 1:N2+1, N3))
            ALLOCATE(W(N1, 0:N2+1, N3))
            ALLOCATE(P(N1, N2, N3))
            ALLOCATE(VELABS(N1, N2, N3))
            ALLOCATE(TOLPRE(N1, N2, N3))
            ALLOCATE(VORX(N1, N2, N3))
            ALLOCATE(VORY(N1, N2, N3))
            ALLOCATE(VORZ(N1, N2, N3))
            ALLOCATE(VORABS(N1, N2, N3))
            ALLOCATE(VELAVEX(4, N1))
            ALLOCATE(VELAVEY(4, N2))
            ALLOCATE(VELAVEZ(4, N3))
            ALLOCATE(VORAVEX(4, N1))
            ALLOCATE(VORAVEY(4, N2))
            ALLOCATE(VORAVEZ(4, N3))
            ALLOCATE(DU(N1, N2, N3))
            ALLOCATE(DV(N1, N2, N3))
            ALLOCATE(DW(N1, N2, N3))
            ALLOCATE(DURMS(N2))
            ALLOCATE(DVRMS(N2))
            ALLOCATE(DWRMS(N2))
            ALLOCATE(DPRMS(N2))
            ALLOCATE(RESX(N1, N2, N3))
            ALLOCATE(RESY(N1, N2, N3))
            ALLOCATE(RESZ(N1, N2, N3))
            ALLOCATE(RESABS(N1, N2, N3))
            ALLOCATE(RESAVEX(4, N1))
            ALLOCATE(RESAVEY(4, N2))
            ALLOCATE(RESAVEZ(4, N3))
            ALLOCATE(TK(N1, N2, N3))
            ALLOCATE(EPS(N1, N2, N3))
        
            ALL_ALLOC = .TRUE.
        END IF
    END SUBROUTINE ALLOC_FIELD
    
    SUBROUTINE DEALLOC_FIELD
        IMPLICIT NONE
        IF(ALL_ALLOC) THEN
            DEALLOCATE(U, V, W, P)
            DEALLOCATE(VELABS, TOLPRE)
            DEALLOCATE(VORX, VORY, VORZ, VORABS)
            DEALLOCATE(VELAVEX, VELAVEY, VELAVEZ)
            DEALLOCATE(VORAVEX, VORAVEY, VORAVEZ)
            DEALLOCATE(DU, DV, DW)
            DEALLOCATE(DURMS, DVRMS, DWRMS, DPRMS)
            DEALLOCATE(RESX, RESY, RESZ, RESABS)
            DEALLOCATE(RESAVEX, RESAVEY, RESAVEZ)
            DEALLOCATE(TK, EPS)
            
            ALL_ALLOC = .FALSE.
        END IF
    END SUBROUTINE DEALLOC_FIELD
