INCLUDE "mkl_dfti.f90"
    
MODULE MATH
    USE MKL_DFTI
    IMPLICIT NONE
    
    INTERFACE TDMA
        MODULE PROCEDURE TDMA_REAL
        MODULE PROCEDURE TDMA_COMPLEX
    END INTERFACE
    
    !SUPPLIMENTAL VECTOR
    REAL, SAVE, ALLOCATABLE :: CC(:), BB(:), DD(:)
    REAL, SAVE, ALLOCATABLE :: U(:), V(:)
    INTEGER, SAVE :: LENGTH = -1
    LOGICAL, SAVE :: ALL_ALLOCATED
    
    !FFT/IFFT SUPPLIMENTAL ARRAY
    COMPLEX, SAVE, ALLOCATABLE :: PF_1D(:)
    COMPLEX, SAVE, ALLOCATABLE :: SF_1D(:)
    INTEGER, SAVE :: FFT_DIM(2) = -1, FFT_LEN = -1
    TYPE(DFTI_DESCRIPTOR), SAVE, POINTER :: HANDLER
    INTEGER, SAVE :: STATUS
    LOGICAL, SAVE :: FFT_ALLOCATED = .FALSE.
    
    PRIVATE BB, CC, DD, U, V, LENGTH, ALL_ALLOCATED
    PRIVATE PF_1D, SF_1D, FFT_ALLOCATED
    
    CONTAINS
    
    SUBROUTINE INITVEC(N)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        
        IF(ALL_ALLOCATED) THEN
            IF(N <= LENGTH) THEN
                RETURN
            END IF
            DEALLOCATE(BB, CC, DD, U, V)
            ALLOCATE(BB(N))
            ALLOCATE(CC(N))
            ALLOCATE(DD(N))
            ALLOCATE(U(N))
            ALLOCATE(V(N))
            LENGTH = N
        ELSE
            ALLOCATE(BB(N))
            ALLOCATE(CC(N))
            ALLOCATE(DD(N))
            ALLOCATE(U(N))
            ALLOCATE(V(N))
            LENGTH = N
            ALL_ALLOCATED = .TRUE.
        END IF
    END SUBROUTINE INITVEC
    
    SUBROUTINE DELETEVEC()
        IMPLICIT NONE
        IF(ALL_ALLOCATED) THEN
            DEALLOCATE(BB, CC, DD, U, V)
            LENGTH = 0
        END IF
    END SUBROUTINE DELETEVEC
            
    SUBROUTINE TDMA_REAL(A, B, C, X, R)
        IMPLICIT NONE
        REAL, INTENT(IN) :: A(:), B(:), C(:), R(:)
        REAL, INTENT(OUT) :: X(:)
        
        INTEGER I, N
        
        N = SIZE(X)
        CALL INITVEC(N)
        
        CC(1) = C(1) / B(1)
        DO I = 2, N - 1
            CC(I) = C(I) / (B(I) - A(I) * CC(I - 1))
        END DO
        
        BB = 1
        
        DD(1) = R(1) / B(1)
        DO I = 2, N
            DD(I) = (R(I) - A(I) * DD(I - 1)) / (B(I) - A(I) * CC(I - 1))
        END DO
        
        X(N) = DD(N)
        DO I = N - 1, 1, -1
            X(I) = DD(I) - CC(I) * X(I + 1)
        END DO
    END SUBROUTINE TDMA_REAL
    
    SUBROUTINE TDMA_COMPLEX(A, B, C, X, R)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(:), B(:), C(:), R(:)
        COMPLEX, INTENT(OUT) :: X(:)
        
        INTEGER I, N
        
        N = SIZE(X)
        CALL INITVEC(N)
        
        CC(1) = C(1) / B(1)
        DO I = 2, N - 1
            CC(I) = C(I) / (B(I) - A(I) * CC(I - 1))
        END DO
        
        BB = 1
        
        DD(1) = R(1) / B(1)
        DO I = 2, N
            DD(I) = (R(I) - A(I) * DD(I - 1)) / (B(I) - A(I) * CC(I - 1))
        END DO
        
        X(N) = DD(N)
        DO I = N - 1, 1, -1
            X(I) = DD(I) - CC(I) * X(I + 1)
        END DO
    END SUBROUTINE TDMA_COMPLEX
    
    SUBROUTINE CTDMA(A, B, C, X, R)
        IMPLICIT NONE
        REAL, INTENT(IN) :: A(:), B(:), C(:), R(:)
        REAL, INTENT(OUT) :: X(:)
        
        REAL XN
        INTEGER I, N
        
        N = SIZE(X)
        CALL INITVEC(N)
        
        CC(1) = C(1) / B(1)
        DO I = 2, N - 2
            CC(I) = C(I) / (B(I) - A(I) * CC(I - 1))
        END DO
        
        BB = 1
        
        DD(1) = -A(1) / B(1)
        DO I = 2, N - 2
            DD(I) = -A(I) * DD(I - 1) / (B(I) - A(I) * CC(I - 1))
        END DO
        DD(N - 1) = (-C(N - 1) - A(I) * DD(I - 1)) / (B(I) - A(I) * CC(I - 1))
        
        U(N - 1) = DD(N - 1)
        DO I = N - 2, 1, -1
            U(I) = DD(I) - CC(I) * U(I + 1)
        END DO
        
        DD(1) = R(1) / B(1)
        DO I = 2, N - 1
            DD(I) = (R(I) - A(I) * DD(I - 1)) / (B(I) - A(I) * CC(I - 1))
        END DO
        
        V(N - 1) = DD(N - 1)
        DO I = N - 2, 1, -1
            V(I) = DD(I) - CC(I) * V(I + 1)
        END DO
        
        XN = (R(N) - C(N) * V(1) - A(N) * V(N - 1)) / &
             (C(N) * U(1) + A(N) * U(N - 1) + B(N))
        
        X = U * XN + V
        X(N) = XN
        
    END SUBROUTINE CTDMA
    
    SUBROUTINE FFT_ALLOC(L)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: L(3)
        
        IF(FFT_ALLOCATED) THEN
            IF(SIZE(PF_1D) < L(1) * L(3)) THEN
                DEALLOCATE(PF_1D, SF_1D)
                ALLOCATE(PF_1D(L(1) * L(3)))
                ALLOCATE(SF_1D(L(1) * L(3)))
            END IF
            IF(FFT_DIM(1) .NE. L(1) .OR. FFT_DIM(2) .NE. L(3)) THEN
                STATUS = DFTIFREEDESCRIPTOR(HANDLER)
                STATUS = DFTICREATEDESCRIPTOR(HANDLER, DFTI_DOUBLE, DFTI_COMPLEX, 2, (/L(1), L(3)/))
                STATUS = DFTISETVALUE(HANDLER, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
                STATUS = DFTICOMMITDESCRIPTOR(HANDLER)
            END IF
        ELSE
            ALLOCATE(PF_1D(L(1) * L(3)))
            ALLOCATE(SF_1D(L(1) * L(3)))
            STATUS = DFTICREATEDESCRIPTOR(HANDLER, DFTI_DOUBLE, DFTI_COMPLEX, 2, (/L(1), L(3)/))
            STATUS = DFTISETVALUE(HANDLER, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
            STATUS = DFTICOMMITDESCRIPTOR(HANDLER)
            FFT_ALLOCATED = .TRUE.
        END IF
        
        FFT_DIM(1) = L(1)
        FFT_DIM(2) = L(3)
        FFT_LEN = L(1) * L(3)
    END SUBROUTINE FFT_ALLOC
    
    SUBROUTINE FFT_DEALLOC()
        IMPLICIT NONE
        
        IF(FFT_ALLOCATED) THEN
            DEALLOCATE(PF_1D, SF_1D)
            STATUS = DFTIFREEDESCRIPTOR(HANDLER)
            FFT_ALLOCATED = .FALSE.
        END IF
        
        FFT_DIM = -1
    END SUBROUTINE FFT_DEALLOC
    
    SUBROUTINE FFT(PF, SF)
        IMPLICIT NONE
        REAL, INTENT(IN) :: PF(:, :, :)
        COMPLEX, INTENT(OUT) :: SF(:, :, :)
        INTEGER J
        
        CALL FFT_ALLOC(SHAPE(PF))
        DO J = 1, SIZE(PF, 2)
            PF_1D(1 : FFT_LEN) = RESHAPE(PF(:, J, :), (/FFT_LEN/))
            STATUS = DFTICOMPUTEFORWARD(HANDLER, PF_1D(1 : FFT_LEN), SF_1D(1 : FFT_LEN))
            SF(:, J, :) = RESHAPE(SF_1D, (/SIZE(SF, 1), SIZE(SF, 3)/))
        END DO
        
    END SUBROUTINE FFT
    
    SUBROUTINE IFFT(SF, PF)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: SF(:, :, :)
        REAL, INTENT(OUT) :: PF(:, :, :)
        INTEGER J
        
        CALL FFT_ALLOC(SHAPE(PF))
        DO J = 1, SIZE(PF, 2)
            SF_1D(1 : FFT_LEN) = RESHAPE(SF(:, J, :), (/FFT_LEN/))
            STATUS = DFTICOMPUTEBACKWARD(HANDLER, SF_1D(1 : FFT_LEN), PF_1D(1 : FFT_LEN))
            PF(:, J, :) = REAL(RESHAPE(PF_1D / SIZE(PF, 1) / SIZE(PF, 3), (/SIZE(PF, 1), SIZE(PF, 3)/)))
        END DO
        
    END SUBROUTINE IFFT
    
END MODULE MATH