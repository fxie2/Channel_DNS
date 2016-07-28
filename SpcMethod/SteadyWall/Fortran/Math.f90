include "mkl_dfti.f90"

module Math
    use MKL_DFTI
    implicit none
    interface FFT
        module procedure FFT_Complex
        module procedure FFT_Real
        module procedure FFT2_Complex
        module procedure FFT2_Real
    end interface
    interface IFFT
        module procedure IFFT_Complex
        module procedure IFFT_Real
        module procedure IFFT2_Complex
        module procedure IFFT2_Real
    end interface
    interface FFTX
        module procedure FFT2X_Complex
        module procedure FFT2X_Real
    end interface
    interface IFFTX
        module procedure IFFT2X_Complex
        module procedure IFFT2X_Real
    end interface
    interface FCT
        module procedure FCT_Complex
        module procedure FCT_Real
    end interface
    interface IFCT
        module procedure IFCT_Complex
        module procedure IFCT_Real
    end interface
    interface XZTrans
        module procedure XZTrans_C
        module procedure XZTrans_R
    end interface
    interface IXZTrans
        module procedure IXZTrans_C
        module procedure IXZTrans_R
    end interface
    interface YTrans
        module procedure YTrans_C
        module procedure YTrans_R
    end interface
    interface IYTrans
        module procedure IYTrans_C
        module procedure IYTrans_R
    end interface

    !Variables
    !FFT handlers
    type(DFTI_DESCRIPTOR), pointer :: handle_1D
    type(DFTI_DESCRIPTOR), pointer :: handle_2D
    type(DFTI_DESCRIPTOR), pointer :: handle_2DX
    !FFT Arrays
    complex, allocatable :: X_1D(:), X_1DX(:)
    
    !FCT Arrays
    complex, allocatable :: Xin(:), Xout(:)
    
    !FCL Arrays
    complex, allocatable :: U_trans_y(:, :, :), V_trans_y(:, :, :), W_trans_y(:, :, :)
    complex, allocatable :: U_extend(:, :, :), V_extend(:, :, :), W_extend(:, :, :)
    real, allocatable :: U_real(:, :, :), V_real(:, :, :), W_real(:, :, :)
    
    !Trans Arrays
    complex, allocatable :: US_temp(:, :, :)
    
    !Ode solver Arrays
    real, allocatable :: infac(:, :), bdfac(:, :)
    complex, allocatable :: b(:)
    
    contains
    
    subroutine MathKernel_init(NX, NY, NZ)
        implicit none
        integer, intent(in) :: NX, NY, NZ
        integer :: status
        if(sizeof(1.0) == 8) then
            !Create handler for 1D FFT/IFFT
            status = DftiCreateDescriptor(handle_1D, DFTI_DOUBLE, DFTI_COMPLEX, 1, (NY - 1) * 2)
            status = DftiCommitDescriptor(handle_1D)
            !Create handler for 2D FFT/IFFT
            status = DftiCreateDescriptor(handle_2D, DFTI_DOUBLE, DFTI_COMPLEX, 2, (/NX, NZ/))
            status = DftiCommitDescriptor(handle_2D)
            !Create handler for 2D Extend FFT/IFFT
            status = DftiCreateDescriptor(handle_2DX, DFTI_DOUBLE, DFTI_COMPLEX, 2, (/NX * 3 / 2, NZ * 3 / 2/))
            status = DftiCommitDescriptor(handle_2DX)
        else
            !Create handler for 1D FFT/IFFT
            status = DftiCreateDescriptor(handle_1D, DFTI_SINGLE, DFTI_COMPLEX, 1, (NY - 1) * 2)
            status = DftiCommitDescriptor(handle_1D)
            !Create handler for 2D FFT/IFFT
            status = DftiCreateDescriptor(handle_2D, DFTI_SINGLE, DFTI_COMPLEX, 2, (/NX, NZ/))
            status = DftiCommitDescriptor(handle_2D)
            !Create handler for 2D Extend FFT/IFFT
            status = DftiCreateDescriptor(handle_2DX, DFTI_SINGLE, DFTI_COMPLEX, 2, (/NX * 3 / 2, NZ * 3 / 2/))
            status = DftiCommitDescriptor(handle_2DX)
        end if
        
        allocate(X_1D(NX * NZ))
        allocate(X_1DX(NX * NZ * 9 / 4))
        allocate(Xin(2*NY - 2))
        allocate(Xout(2*NY - 2))
        allocate(U_trans_y(NY, NX, NZ))
        allocate(V_trans_y(NY, NX, NZ))
        allocate(W_trans_y(NY, NX, NZ))
        allocate(U_extend(NY, NX / 2 * 3, NZ / 2 * 3))
        allocate(V_extend(NY, NX / 2 * 3, NZ / 2 * 3))
        allocate(W_extend(NY, NX / 2 * 3, NZ / 2 * 3))
        allocate(U_real(NY, NX / 2 * 3, NZ / 2 * 3))
        allocate(V_real(NY, NX / 2 * 3, NZ / 2 * 3))
        allocate(W_real(NY, NX / 2 * 3, NZ / 2 * 3))
        allocate(US_temp(NY, NX, NZ))
        allocate(infac(3, NY - 2))
        allocate(bdfac(NY, 2))
        allocate(b(NY))
    end subroutine MathKernel_init
    
    subroutine MathKernel_finalize
        implicit none
        integer :: status
        !Delete handlers for FFT/IFFT
        status = DftiFreeDescriptor(handle_1D)
        status = DftiFreeDescriptor(handle_2D)
        status = DftiFreeDescriptor(handle_2DX)
        
        deallocate(Xin, Xout)
        deallocate(X_1D, X_1DX)
        deallocate(U_trans_y, V_trans_y, W_trans_y)
        deallocate(U_extend, V_extend, W_extend)
        deallocate(U_real, V_real, W_real)
        deallocate(US_temp)
        deallocate(infac, bdfac, b)
    end subroutine MathKernel_finalize
    
    subroutine FFT_Complex(X, Y)
        implicit none
        complex, intent(in)  :: X(:)
        complex :: Y(:)
        integer :: status
        Y = X
        status = DftiComputeForward(handle_1D, Y)
    end subroutine FFT_Complex
    
    subroutine FFT_Real(X, Y)
        implicit none
        real, intent(in) :: X(:)
        complex :: Y(:)
        integer :: status
        Y = X
        status = DftiComputeForward(handle_1D, Y)
    end subroutine FFT_Real
    
    subroutine FFT2_Complex(X, Y)
        implicit none
        complex, intent(in)     :: X(:, :)
        complex, intent(out)    :: Y(:, :)
        integer :: status
        X_1D = reshape(X, shape(X_1D))
        status = DftiComputeForward(handle_2D, X_1D)
        Y = reshape(X_1D, shape(Y))
    end subroutine FFT2_Complex
    
    subroutine FFT2_Real(X, Y)
        implicit none
        real, intent(in)        :: X(:, :)
        complex, intent(out)    :: Y(:, :)
        integer :: status
        X_1D = reshape(X, shape(X_1D))
        status = DftiComputeForward(handle_2D, X_1D)
        Y = reshape(X_1D, shape(Y))
    end subroutine FFT2_Real
    
    subroutine FFT2X_Complex(X, Y)
        implicit none
        complex, intent(in)     :: X(:, :)
        complex, intent(out)    :: Y(:, :)
        integer :: status
        X_1DX = reshape(X, shape(X_1DX))
        status = DftiComputeForward(handle_2DX, X_1DX)
        Y = reshape(X_1DX, shape(Y))
    end subroutine FFT2X_Complex
    
    subroutine FFT2X_Real(X, Y)
        implicit none
        real, intent(in)     :: X(:, :)
        complex, intent(out)    :: Y(:, :)
        integer :: status
        X_1DX = reshape(X, shape(X_1DX))
        status = DftiComputeForward(handle_2DX, X_1DX)
        Y = reshape(X_1DX, shape(Y))
    end subroutine FFT2X_Real
    
    subroutine IFFT2X_Complex(X, Y)
        implicit none
        complex, intent(in)     :: X(:, :)
        complex, intent(out)    :: Y(:, :)
        integer :: status
        X_1DX = reshape(X, shape(X_1DX))
        status = DftiComputeBackward(handle_2DX, X_1DX)
        Y = reshape(X_1DX, shape(Y)) / size(X, 1) / size(X, 2)
    end subroutine IFFT2X_Complex
    
    subroutine IFFT2X_Real(X, Y)
        implicit none
        complex, intent(in)     :: X(:, :)
        real, intent(out)    :: Y(:, :)
        integer :: status
        X_1DX = reshape(X, shape(X_1DX))
        status = DftiComputeBackward(handle_2DX, X_1DX)
        Y = reshape(real(X_1DX), shape(Y)) / size(X, 1) / size(X, 2)
    end subroutine IFFT2X_Real

    subroutine IFFT_Complex(X, Y)
        implicit none
        complex, intent(in) :: X(:)
        complex :: Y(:)
        integer :: status
        Y = X
        status = DftiComputeBackward(handle_1D, Y)
        Y = Y / size(Y, 1)
    end subroutine IFFT_Complex
    
    subroutine IFFT_Real(X, Y)
        implicit none
        real, intent(in) :: X(:)
        complex :: Y(:)
        integer :: status
        Y = X
        status = DftiComputeBackward(handle_1D, Y)
        Y = Y / size(Y, 1)
    end subroutine IFFT_Real
    
    subroutine IFFT2_Complex(X, Y)
        implicit none
        complex, intent(in) :: X(:, :)
        complex, intent(out) :: Y(:, :)
        integer :: status
        X_1D = reshape(X, shape(X_1D))
        status = DftiComputeBackward(handle_2D, X_1D)
        Y = reshape(X_1D, shape(Y)) / size(Y, 1) / size(Y, 2)
    end subroutine IFFT2_Complex
    
    subroutine IFFT2_Real(X, Y)
        implicit none
        complex, intent(in) :: X(:, :)
        real, intent(out) :: Y(:, :)
        integer :: status
        X_1D = reshape(X, shape(X_1D))
        status = DftiComputeBackward(handle_2D, X_1D)
        Y = reshape(real(X_1D), shape(Y)) / size(Y, 1) / size(Y, 2)
    end subroutine IFFT2_Real
    
    subroutine FFTshift(X)
        implicit none
        complex, intent(inout) :: X(:, :)
        integer :: N1, N2
        N1 = size(X, 1)
        N2 = size(X, 2)
        X = cshift(X, -N1/2, 1)
        X = cshift(X, -N2/2, 2)
    end subroutine FFTshift
    
    subroutine IFFTshift(X)
        implicit none
        complex, intent(inout) :: X(:, :)
        integer :: N1, N2
        N1 = size(X, 1)
        N2 = size(X, 2)
        X = cshift(X, N2/2, 2)
        X = cshift(X, N1/2, 1)
    end subroutine IFFTshift
    
    subroutine FCT_Complex(X, Y)
        implicit none
        complex, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        integer N
        N = size(X)
        Xin(1:N) = X(1:N)
        Xin(N+1:2*N-2) = X(N-1:2:-1)
        call IFFT(Xin, Xout)
        Y = 2 * Xout(1:N)
        Y(1) = Y(1) / 2
        Y(N) = Y(N) / 2
    end subroutine FCT_Complex
    
    subroutine FCT_Real(X, Y)
        implicit none
        real, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        integer N
        N = size(X)
        Xin(1:N) = X(1:N)
        Xin(N+1:2*N-2) = X(N-1:2:-1)
        call IFFT(Xin, Xout)
        Y = 2 * Xout(1:N)
        Y(1) = Y(1) / 2
        Y(N) = Y(N) / 2
    end subroutine FCT_Real

    subroutine IFCT_Complex(X, Y)
        implicit none
        complex, intent(in) :: X(:)
        complex, intent(out) :: Y(:)
        integer N
        N = size(X)
        Xin(1) = X(1)
        Xin(2:N-1) = X(2:N-1) / 2
        Xin(N) = X(N)
        Xin(N+1:2*N-2) = X(N-1:2:-1) / 2
        call FFT(Xin, Xout)
        Y = Xout(1:N)
    end subroutine IFCT_Complex
    
    subroutine IFCT_Real(X, Y)
        implicit none
        complex, intent(in) :: X(:)
        real, intent(out) :: Y(:)
        integer N
        N = size(X)
        Xin(1) = X(1)
        Xin(2:N-1) = X(2:N-1) / 2
        Xin(N) = X(N)
        Xin(N+1:2*N-2) = X(N-1:2:-1) / 2
        call FFT(Xin, Xout)
        Y = real(Xout(1:N))
    end subroutine IFCT_Real
    
    subroutine FCL(U, V, W)
        implicit none
        complex, intent(in) :: U(:, :, :), V(:, :, :)
        complex, intent(out) :: W(:, :, :)
        integer nx, ny, nz
        integer i, j
        
        nx = size(U, 2)
        ny = size(U, 1)
        nz = size(U, 3)
        
        do j = 1, nz
            do i = 1, nx
                call IFCT(U(:, i, j), U_trans_y(:, i, j))
            end do
        end do
        do j = 1, nz
            do i = 1, nx
                call IFCT(V(:, i, j), V_trans_y(:, i, j))
            end do
        end do
        
        U_extend = 0
        V_extend = 0
        U_extend(:, nx / 4 + 1 : nx / 4 + nx, nz / 4 + 1 : nz / 4 + nz) = U_trans_y
        V_extend(:, nx / 4 + 1 : nx / 4 + nx, nz / 4 + 1 : nz / 4 + nz) = V_trans_y
        U_real = 0
        V_real = 0
        do i = 1, ny
            call IFFTshift(U_extend(i, :, :))
            call IFFTX(U_extend(i, :, :), U_real(i, :, :))
        end do
        do i = 1, ny
            call IFFTshift(V_extend(i, :, :))
            call IFFTX(V_extend(i, :, :), V_real(i, :, :))
        end do
        W_real = U_real * V_real
        W_extend = 0
        do i = 1, ny
            call FFTX(W_real(i, :, :), W_extend(i, :, :))
            call FFTshift(W_extend(i, :, :))
        end do
        W_trans_y = W_extend(:, nx / 4 + 1 : nx / 4 + nx, nz / 4 + 1 : nz / 4 + nz)
        do j = 1, nz
            do i = 1, nx
                call FCT(W_trans_y(:, i, j), W(:, i, j))
            end do
        end do
        W = W * 9 / 4
    end subroutine FCL
    
    subroutine XZTrans_C(UP, US)
        implicit none
        complex, intent(in) :: UP(:, :, :)
        complex, intent(out) :: US(:, :, :)
        integer ny, i
        
        ny = size(UP, 1)
        
        do i = 1, ny
            call FFT(UP(i, :, :), US(i, :, :))
            call FFTshift(US(i, :, :))
        end do
    end subroutine XZTrans_C
    
    subroutine XZTrans_R(UP, US)
        implicit none
        real, intent(in) :: UP(:, :, :)
        complex, intent(out) :: US(:, :, :)
        integer ny, i
        
        ny = size(UP, 1)
        
        do i = 1, ny
            call FFT(UP(i, :, :), US(i, :, :))
            call FFTshift(US(i, :, :))
        end do
    end subroutine XZTrans_R
    
    subroutine IXZTrans_C(US, UP)
        implicit none
        complex, intent(in) :: US(:, :, :)
        complex, intent(out) :: UP(:, :, :)
        integer ny, i
        
        ny = size(UP, 1)
        US_temp = US
        
        do i = 1, ny
            call IFFTshift(US_temp(i, :, :))
            call IFFT(US_temp(i, :, :), UP(i, :, :))
        end do
    end subroutine IXZTrans_C
    
    subroutine IXZTrans_R(US, UP)
        implicit none
        complex, intent(in) :: US(:, :, :)
        real, intent(out) :: UP(:, :, :)
        integer ny, i
        
        ny = size(UP, 1)
        US_temp = US
        
        do i = 1, ny
            call IFFTshift(US_temp(i, :, :))
            call IFFT(US_temp(i, :, :), UP(i, :, :))
        end do
    end subroutine IXZTrans_R
    
    subroutine YTrans_C(UP, US)
        implicit none
        complex, intent(in) :: UP(:, :, :)
        complex, intent(out) :: US(:, :, :)
        integer nx, nz, i, j
        
        nx = size(UP, 2)
        nz = size(UP, 3)
        
        do j = 1, nz
            do i = 1, nx
                call FCT(UP(:, i, j), US(:, i, j))
            end do
        end do
    end subroutine YTrans_C
    
    subroutine YTrans_R(UP, US)
        implicit none
        real, intent(in) :: UP(:, :, :)
        complex, intent(out) :: US(:, :, :)
        integer nx, nz, i, j
        
        nx = size(UP, 2)
        nz = size(UP, 3)
        
        do j = 1, nz
            do i = 1, nx
                call FCT(UP(:, i, j), US(:, i, j))
            end do
        end do
    end subroutine YTrans_R
    
    subroutine IYTrans_C(US, UP)
        implicit none
        complex, intent(in) :: US(:, :, :)
        complex, intent(out) :: UP(:, :, :)
        integer nx, nz, i, j
        
        nx = size(UP, 2)
        nz = size(UP, 3)
        
        do j = 1, nz
            do i = 1, nx
                call IFCT(US(:, i, j), UP(:, i, j))
            end do
        end do
    end subroutine IYTrans_C
    
    subroutine IYTrans_R(US, UP)
        implicit none
        complex, intent(in) :: US(:, :, :)
        real, intent(out) :: UP(:, :, :)
        integer nx, nz, i, j
        
        nx = size(UP, 2)
        nz = size(UP, 3)
        
        do j = 1, nz
            do i = 1, nx
                call IFCT(US(:, i, j), UP(:, i, j))
            end do
        end do
    end subroutine IYTrans_R
    
    subroutine GTrans(UP, US)
        implicit none
        real, intent(in) :: UP(:, :, :)
        complex, intent(out) :: US(:, :, :)
        
        call XZTrans(UP, US)
        call YTrans(US, US)
    end subroutine GTrans
    
    subroutine IGTrans(US, UP)
        implicit none
        complex, intent(in) :: US(:, :, :)
        real, intent(out) :: UP(:, :, :)
        
        call IYTrans(US, US_temp)
        call IXZTrans(US_temp, UP)
    end subroutine IGTrans
    
    subroutine D_DX(U, DU, alpha)
        implicit none
        complex, intent(in) :: U(:, :, :)
        complex, intent(out) :: DU(:, :, :)
        real :: alpha
        integer nx, i
        
        nx = size(U, 2)
        forall(i = 1:nx)
            DU(:, i, :) = (0, 1) * alpha * U(:, i, :) * (i - nx / 2 - 1)
        end forall
    end subroutine D_DX
    
    subroutine D_DZ(U, DU, beta)
        implicit none
        complex, intent(in) :: U(:, :, :)
        complex, intent(out) :: DU(:, :, :)
        real :: beta
        integer nz, i
        
        nz = size(U, 3)
        forall(i = 1:nz)
            DU(:, :, i) = (0, 1) * beta * U(:, :, i) * (i - nz / 2 - 1)
        end forall
    end subroutine D_DZ
    
    subroutine D_DY(U, DU)
        implicit none
        complex, intent(in) :: U(:, :, :)
        complex, intent(out) :: DU(:, :, :)
        integer nx, nz, i, j, k
        
        nx = size(U, 2)
        nz = size(U, 3)
        
        do j = 1, nz
            do i = 1, nx
                call D_DY1D(U(:, i, j), DU(:, i, j))
            end do
        end do
        
    end subroutine D_DY
    
    pure subroutine D_DY1D(F, DF)
        implicit none
        complex, intent(in) :: F(:)
        complex, intent(out) :: DF(:)
        integer p, i, j
        
        p = size(F, 1) - 1
        DF = 0
        do i = 1, p, 2
            DF(1) = DF(1) + i * F(i + 1)
        end do
        do j = 2, p
            do i = j, p, 2
                DF(j) = DF(j) + 2 * i * F(i + 1)
            end do
        end do
    end subroutine D_DY1D
    
    subroutine ode_solve(k, x, s, boundary_type, bc1, bc2)
        implicit none
        real, intent(in) :: k
        complex, intent(in) :: s(:), bc1, bc2
        complex, intent(out) :: x(:)
        integer, intent(in) :: boundary_type
        integer p, i, fac
        
        p = size(s) - 1
        do i = 1, p - 3
            infac(1, i) = k / 4 / i / (i + 1)
            infac(2, i) = -(1 + k / 2 / ((i + 1) ** 2 - 1))
            infac(3, i) = k / 4 / (i + 1) / (i + 2)
        end do
        infac(1, 1) = k / 4
        infac(3, p - 4 : p - 3) = 0
        infac(1, p - 2) = k / 4 / (p - 1) / (p - 2)
        infac(2, p - 2) = -1
        infac(3, p - 2) = 0
        infac(1, p - 1) = k / 4 / p / (p - 1)
        infac(2, p - 1) = -1
        infac(3, p - 1) = 0
        if(boundary_type == 1) then
            do i = 0, p
                bdfac(i + 1, 1) = 1 - mod(i, 2) * 2
            end do
            bdfac(:, 2) = 1
        else
            do i = 0, p
                bdfac(i + 1, 1) = (1 - mod(i + 1, 2) * 2) * i * i
                bdfac(i + 1, 1) = i * i
            end do
        end if
        b(1) = -(s(1) / 4 - s(3) / 6 + s(5) / 24)
        do i = 3, p - 4
            b(i - 1) = -(s(i - 1) / (4 * i * (i - 1)) - s(i + 1) / (2 * (i * i - 1))    &
                     + s(i + 3) / (4 * i * (i + 1)))
        end do
        b(p - 4) = -(s(p - 4) / (4 * (p - 3) * (p - 4)) - s(p - 2) / (2 * (p - 3) ** 2 - 1))
        b(p - 3) = -(s(p - 3) / (4 * (p - 2) * (p - 3)) - s(p - 1) / (2 * (p - 2) ** 2 - 1))
        b(p - 2) = -s(p - 2) / (4 * (p - 1) * (p - 2))
        b(p - 1) = -s(p - 1) / (4 * (p - 1) * p)
        b(p) = bc1
        b(p + 1) = bc2
        if(abs(k) > epsilon(1.0)) then
            do i = 1, p - 5
                fac = bdfac(i, 1) / infac(1, i)
                bdfac(i, 1) = bdfac(i, 1) - fac * infac(1, i)
                bdfac(i + 2, 1) = bdfac(i + 2, 1) - fac * infac(2, i)
                bdfac(i + 4, 1) = bdfac(i + 4, 1) - fac * infac(3, i)
                b(p) = b(p) - fac * b(i)
                fac = bdfac(i, 2) / infac(1, i)
                bdfac(i, 2) = bdfac(i, 2) - fac * infac(1, i)
                bdfac(i + 2, 2) = bdfac(i + 2, 2) - fac * infac(2, i)
                bdfac(i + 4, 2) = bdfac(i + 4, 2) - fac * infac(3, i)
                b(p + 1) = b(p + 1) - fac * b(i)
            end do
            do i = p - 4, p - 1
                fac = bdfac(i, 1) / infac(1, i)
                bdfac(i, 1) = bdfac(i, 1) - fac * infac(1, i)
                bdfac(i + 2, 1) = bdfac(i + 2, 1) - fac * infac(2, i)
                b(p) = b(p) - fac * b(i)
                fac = bdfac(i, 2) / infac(1, i)
                bdfac(i, 2) = bdfac(i, 2) - fac * infac(1, i)
                bdfac(i + 2, 2) = bdfac(i + 2, 2) - fac * infac(2, i)
                b(p + 1) = b(p + 1) - fac * b(i)
            end do
            x(p) = (bdfac(p + 1, 2) * b(p) - bdfac(p + 1, 1) * b(p + 1))    &
                / (bdfac(p, 1) * bdfac(p + 1, 2) - bdfac(p + 1, 1) * bdfac(p, 2))
            x(p + 1) = (bdfac(p, 1) * b(p + 1) - bdfac(p, 2) * b(p))    &
                / (bdfac(p, 1) * bdfac(p + 1, 2) - bdfac(p + 1, 1) * bdfac(p, 2))
            do i = p - 1, p - 4, -1
                x(i) = (b(i) - infac(2, i) * x(i + 2)) / infac(1, i)
            end do
            do i = p - 5, 1, -1
                x(i) = (b(i) - infac(3, i) * x(i + 4) - infac(2, i) * x(i + 2)) / infac(1, i)
            end do
        else if(boundary_type == 1) then
            x(3:p + 1) = -b(1:p - 1)
            b(p) = b(p) - sum(x(3:p + 1) * bdfac(3:p + 1, 1))
            b(p + 1) = b(p + 1) - sum(x(3:p + 1) * bdfac(3:p + 1, 2))
            x(1) = (b(p) * bdfac(2, 2) - b(p + 1) * bdfac(2, 1)) &
                 / (bdfac(1, 1) * bdfac(2, 2) - bdfac(2, 1) * bdfac(1, 2))
            x(2) = (bdfac(1, 1) * b(p + 1) - bdfac(1, 2) * b(p)) &
                 / (bdfac(1, 1) * bdfac(2, 2) - bdfac(2, 1) * bdfac(1, 2))
        else
            x(3:p + 1) = -b(1:p - 1)
            x(2) = (b(p) - sum(x(3:p + 1) * bdfac(3:p + 1, 1))) / bdfac(2, 1)
            x(1) = 0
        end if
    end subroutine ode_solve
    
end module Math