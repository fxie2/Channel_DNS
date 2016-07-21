include "mkl_dfti.f90"

module Math
    implicit none
    interface FFT
        module procedure FFT_Complex
        module procedure FFT_Real
    end interface
    interface IFFT
        module procedure IFFT_Complex
        module procedure IFFT_Real
    end interface
    interface FCT
        module procedure FCT_Complex
        module procedure FCT_Real
    end interface
    interface IFCT
        module procedure IFCT_Complex
        module procedure IFCT_Real
    end interface

    contains
    
    subroutine FFT_Complex(X, Y)
        use MKL_DFTI
        implicit none
        complex, intent(in)  :: X(:)
        complex, intent(out) :: Y(:)
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, size(X))
        status = DftiCommitDescriptor(handle)
        status = DftiComputeForward(handle, X, Y)
        status = DftiFreeDescriptor(handle)
    end subroutine FFT_Complex
    
    subroutine FFT_Real(X, Y)
        use MKL_DFTI
        implicit none
        real, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_REAL, 1, size(X))
        status = DftiCommitDescriptor(handle)
        status = DftiComputeForward(handle, X, Y)
        status = DftiFreeDescriptor(handle)
    end subroutine FFT_Real
    
    subroutine IFFT_Complex(X, Y)
        use MKL_DFTI
        implicit none
        complex, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, size(X))
        status = DftiCommitDescriptor(handle)
        status = DftiComputeBackward(handle, X, Y)
        status = DftiFreeDescriptor(handle)
    end subroutine IFFT_Complex
    
    subroutine IFFT_Real(X, Y)
        use MKL_DFTI
        implicit none
        real, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, size(X))
        status = DftiCommitDescriptor(handle)
        status = DftiComputeBackward(handle, X*(1,0), Y)
        status = DftiFreeDescriptor(handle)
    end subroutine IFFT_Real

    subroutine FCT_Complex(X, Y)
        implicit none
        complex, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        complex, allocatable :: Xin(:), Xout(:)
        integer N
        N = size(X)
        allocate(Xin(2*N - 2), Xout(2*N - 2))
        Xin(1:N) = X(1:N)
        Xin(N+1:2*N-2) = X(N-1:2:-1)
        call IFFT(Xin, Xout)
        Y = 2 * Xout(1:N)
        deallocate(Xin, Xout)
        Y(1) = Y(1) / 2
        Y(N) = Y(N) / 2
    end subroutine FCT_Complex
    
    subroutine FCT_Real(X, Y)
        implicit none
        real, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        complex, allocatable :: Xin(:), Xout(:)
        integer N
        N = size(X)
        allocate(Xin(2*N - 2), Xout(2*N - 2))
        Xin(1:N) = X(1:N)
        Xin(N+1:2*N-2) = X(N-1:2:-1)
        call IFFT(Xin, Xout)
        Y = 2 * Xout(1:N)
        deallocate(Xin, Xout)
        Y(1) = Y(1) / 2
        Y(N) = Y(N) / 2
    end subroutine FCT_Real

    subroutine IFCT_Complex(X, Y)
        implicit none
        complex, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        complex, allocatable :: Xin(:), Xout(:)
        integer N
        N = size(X)
        allocate(Xin(2*N - 2), Xout(2*N - 2))
        Xin(1) = X(1)
        Xin(2:N-1) = X(2:N-1) / 2
        Xin(N) = X(N)
        Xin(N+1:2*N-2) = X(N-1:2:-1) / 2
        call FFT(Xin, Xout)
        Y = Xout(1:N)
        deallocate(Xin, Xout)
    end subroutine IFCT_Complex
    
    subroutine IFCT_Real(X, Y)
        implicit none
        real, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        complex, allocatable :: Xin(:), Xout(:)
        integer N
        N = size(X)
        allocate(Xin(2*N - 2), Xout(2*N - 2))
        Xin(1) = X(1)
        Xin(2:N-1) = X(2:N-1) / 2
        Xin(N) = X(N)
        Xin(N+1:2*N-2) = X(N-1:2:-1) / 2
        call FFT(Xin, Xout)
        Y = Xout(1:N)
        deallocate(Xin, Xout)
    end subroutine IFCT_Real
end module Math