PROGRAM MAIN
    USE GLOBAL_PARAMETER
    USE MESH
    USE MATH
    USE DEBUG
    
    REAL, ALLOCATABLE :: P(:, :, :), PB(:, :, :)
    COMPLEX, ALLOCATABLE :: PS(:, :, :), PBS(:, :, :)
    INTEGER I, J, K
    !REAL :: T = 0
    
    CALL INIT_PARAMETERS()
    CALL INIT_MESH
    !ALLOCATE(P(N1, N2, N3))
    !ALLOCATE(PB(N1, 2, N3))
    !!ALLOCATE(LP(N1, N2, N3))
    !ALLOCATE(PS(N1, N2, N3))
    !ALLOCATE(PBS(N1, 2, N3))
    !!ALLOCATE(LPS(N1, N2, N3))
    !
    !DO K = 1, N3
    !    DO J = 1, N2
    !        DO I = 1, N1
    !            P(I, J, K) = SIN(GETX(I, J, K, 0.))
    !        ENDDO
    !    ENDDO
    !ENDDO
    !
    !PB = 0
    !CALL FFT(P, PS)
    !CALL FFT(PB, PBS)
    CALL ALLOC_DEBUG
    !CALL TEST_PSOLVE(P, PB, PS, PBS, 0.)
    CALL TEST_DIV
    END PROGRAM MAIN