include "mkl_dfti.f90"
module FieldMatrix
    use Global_parameter
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
    
    !Real field
    real(8), save, allocatable :: U(:, :, :), V(:, :, :), W(:, :, :)
    real(8), save, allocatable :: P(:, :, :)
    !Spectral field
    !_S for s step, _1 for s-1 step, _2 for s-2 step, _NEW for s+1 step
    
    !Velocity field
    complex, save, allocatable :: US(:, :, :), VS(:, :, :), WS(:, :, :)
    complex, save, allocatable :: U1(:, :, :), V1(:, :, :), W1(:, :, :)
    complex, save, allocatable :: U2(:, :, :), V2(:, :, :), W2(:, :, :)
    complex, save, allocatable :: UNEW(:, :, :), VNEW(:, :, :), WNEW(:, :, :)
    
    !Rotation field
    complex, save, allocatable :: rot_x(:, :, :), rot_y(:, :, :), rot_z(:, :, :)
    
    !Temp space
    complex, save, allocatable :: temp(:, :, :)
    
    !Non-linear term
    complex, save, allocatable :: FX(:, :, :), FY(:, :, :), FZ(:, :, :)
    complex, save, allocatable :: FX1(:, :, :), FY1(:, :, :), FZ1(:, :, :)
    complex, save, allocatable :: FX2(:, :, :), FY2(:, :, :), FZ2(:, :, :)
    
    !G term
	complex, save, allocatable :: G(:, :, :), G1(:, :, :), G2(:, :, :)
    contains
    subroutine alloc_field()
        implicit none
        
        allocate(U(NY, NX, NZ))
        allocate(V(NY, NX, NZ))
        allocate(W(NY, NX, NZ))
        allocate(P(NY, NX, NZ))
        
        allocate(US(NY, NX, NZ))
        allocate(VS(NY, NX, NZ))
        allocate(WS(NY, NX, NZ))
        allocate(U1(NY, NX, NZ))
        allocate(V1(NY, NX, NZ))
        allocate(W1(NY, NX, NZ))
        allocate(U2(NY, NX, NZ))
        allocate(V2(NY, NX, NZ))
        allocate(W2(NY, NX, NZ))
        allocate(UNEW(NY, NX, NZ))
        allocate(VNEW(NY, NX, NZ))
        allocate(WNEW(NY, NX, NZ))
        
        allocate(FX(NY, NX, NZ))
        allocate(FY(NY, NX, NZ))
        allocate(FZ(NY, NX, NZ))
        allocate(FX1(NY, NX, NZ))
        allocate(FY1(NY, NX, NZ))
        allocate(FZ1(NY, NX, NZ))
        allocate(FX2(NY, NX, NZ))
        allocate(FY2(NY, NX, NZ))
        allocate(FZ2(NY, NX, NZ))
        
        allocate(G(NY, NX, NZ))
        allocate(G1(NY, NX, NZ))
        allocate(G2(NY, NX, NZ))
        allocate(rot_x(NY, NX, NZ))
        allocate(rot_y(NY, NX, NZ))
        allocate(rot_z(NY, NX, NZ))
        allocate(temp(NY, NX, NZ))
    end subroutine alloc_field
    
    subroutine dealloc_field()
        implicit none
        deallocate(U, V, W, P)
        deallocate(US, VS, WS)
        deallocate(U1, V1, W1, U2, V2, W2, UNEW, VNEW, WNEW)
        deallocate(FX, FY, FZ, FX1, FY1, FZ1, FX2, FY2, FZ2)
        deallocate(rot_x, rot_y, rot_z, temp)
        deallocate(G, G1, G2)
    end subroutine dealloc_field
        
    subroutine init_field()
        implicit none
        integer iter_x, iter_y
        if(allocated(U) == .false.) call alloc_field()
        U = 0
        V = 0
        W = 0
        P = 0
        !Velocity initialize
        forall(iter_y = 1:NY)
            U(iter_y, :, :) = 1 - cos((iter_y - 1) * PI / (NY - 1)) ** 2
        end forall
        U1 = U
        U2 = U
        V1 = V
        V2 = V
        W1 = W
        W2 = W
        call GTrans(U, US)
        call GTrans(V, VS)
        call GTrans(W, WS)
        call Lamb
        FX1 = FX
        FX2 = FX
        FY1 = FY
        FY2 = FY
        FZ1 = FZ
        FZ2 = FZ
        call Gterm
        G1 = G
        G2 = G
    end subroutine init_field
    
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
    
    subroutine GTrans(UP, US)
        use MKL_DFTI
        implicit none
        real, intent(in) :: UP(NY, NX, NZ)
        complex, intent(out) :: US(NY, NX, NZ)
        complex temp_XZ(NX, NZ), temp_Y(NY)
        integer iter_x, iter_y, iter_z
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, [NX, NZ])
        status = DftiCommitDescriptor(handle)
        do iter_y = 1, NY
            temp_XZ(:, :) = UP(iter_y, :, :)
            status = DftiComputeForward(handle, temp_XZ(:,1))
            temp_XZ = cshift(temp_XZ, -NX/2, 2)
            temp_XZ = cshift(temp_XZ, -NZ/2, 1)
            US(iter_y, :, :) = temp_XZ
        end do
        status = DftiFreeDescriptor(handle)
        do iter_z = 1, NZ
            do iter_x = 1, NX
                temp_Y = US(:, iter_x, iter_z)
                call FCT(temp_Y, temp_Y)
                US(:, iter_x, iter_z) = temp_Y
            end do
        end do
    end subroutine GTrans
    
    subroutine IGTrans(US, UP)
        use MKL_DFTI
        implicit none
        complex, intent(in) :: US(NY, NX, NZ)
        real, intent(out) :: UP(NY, NX, NZ)
        complex temp_XZ(NX, NZ), temp_Y_C(NY)
        real    temp_Y(NY)
        integer iter_x, iter_y, iter_z
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, [NX, NZ])
        status = DftiCommitDescriptor(handle)
        do iter_y = 1, NY
            temp_XZ = US(iter_y, :, :)
            temp_XZ = cshift(temp_XZ, NZ/2, 1)
            temp_XZ = cshift(temp_XZ, NX/2, 2)
            status = DftiComputeBackward(handle, temp_XZ(:,1))
            UP(iter_y, :, :) = real(temp_XZ)
        end do
        status = DftiFreeDescriptor(handle)
        do iter_z = 1, NZ
            do iter_x = 1, NX
                temp_Y = UP(:, iter_x, iter_z)
                call IFCT(temp_Y, temp_Y_C)
                UP(:, iter_x, iter_z) = real(temp_Y_C)
            end do
        end do
    end subroutine IGTrans

    subroutine YTrans(UP, US)
        implicit none
        complex, intent(in) :: UP(NY, NX, NZ)
        complex, intent(out) :: US(NY, NX, NZ)
        complex temp_Y(NY)
        integer iter_x, iter_z
        do iter_z = 1, NZ
            do iter_x = 1, NX
                temp_Y = US(:, iter_x, iter_z)
                call FCT(temp_Y, temp_Y)
                US(:, iter_x, iter_z) = temp_Y
            end do
        end do
    end subroutine YTrans
    
    subroutine IYTrans(US, UP)
        implicit none
        complex, intent(in) :: US(NY, NX, NZ)
        complex, intent(out) :: UP(NY, NX, NZ)
        complex temp_Y(NY)
        integer iter_x, iter_z
        do iter_z = 1, NZ
            do iter_x = 1, NX
                temp_Y = UP(:, iter_x, iter_z)
                call IFCT(temp_Y, temp_Y)
                UP(:, iter_x, iter_z) = temp_Y
            end do
        end do
    end subroutine IYTrans
    
    subroutine XZTrans(UP, US)
        use MKL_DFTI
        implicit none
        complex, intent(in) :: UP(:, :, :)
        complex, intent(out) :: US(:, :, :)
        complex, allocatable :: temp_XZ(:, :)
        integer iter_y, N_X, N_Y, N_Z
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        N_X = size(UP, 2)
        N_Y = size(UP, 1)
        N_Z = size(UP, 3)
        allocate(temp_XZ(N_X, N_Z))
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, [N_X, N_Z])
        status = DftiCommitDescriptor(handle)
        do iter_y = 1, N_Y
            temp_XZ(:, :) = UP(iter_y, :, :)
            status = DftiComputeForward(handle, temp_XZ(:,1))
            temp_XZ = cshift(temp_XZ, -N_X/2, 2)
            temp_XZ = cshift(temp_XZ, -N_Z/2, 1)
            US(iter_y, :, :) = temp_XZ
        end do
        status = DftiFreeDescriptor(handle)
        deallocate(temp_XZ)
    end subroutine XZTrans
    
    subroutine IXZTrans(US, UP)
        use MKL_DFTI
        implicit none
        complex, intent(in) :: US(:, :, :)
        complex, intent(out) :: UP(:, :, :)
        complex, allocatable :: temp_XZ(:, :)
        integer iter_y, N_X, N_Y, N_Z
        type(DFTI_DESCRIPTOR), pointer :: handle
        integer :: status
        N_X = size(UP, 2)
        N_Y = size(UP, 1)
        N_Z = size(UP, 3)
        allocate(temp_XZ(N_X, N_Z))
        status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, [N_X, N_Z])
        status = DftiCommitDescriptor(handle)
        do iter_y = 1, N_Y
            temp_XZ = US(iter_y, :, :)
            temp_XZ = cshift(temp_XZ, N_Z/2, 1)
            temp_XZ = cshift(temp_XZ, N_X/2, 2)
            status = DftiComputeBackward(handle, temp_XZ(:,1))
            UP(iter_y, :, :) = real(temp_XZ)
        end do
        status = DftiFreeDescriptor(handle)
        deallocate(temp_XZ)
    end subroutine IXZTrans    
    
    subroutine FCL(U, V, W)
        implicit none
        complex, intent(in) :: U(NY, NX, NZ), V(NY, NX, NZ)
        complex, intent(out):: W(NY, NX, NZ)
        complex     U_trans_y(NY, NX, NZ), V_trans_y(NY, NX, NZ), W_trans_y(NY, NX, NZ)
        complex     U_extend(NY, NX32, NZ32), V_extend(NY, NX32, NZ32), W_extend(NY, NX32, NZ32)
        complex     U_real(NY, NX32, NZ32),   V_real(NY, NX32, NZ32),   W_real(NY, NX32, NZ32)
        complex     temp_Y_in(NY), temp_Y_out(NY)
        integer     iter_x, iter_y, iter_z
        
        call IYTrans(U, U_trans_y)
        call IYTrans(V, V_trans_y)
        
        U_extend = 0
        V_extend = 0
        U_extend(:, NX/4+1 : NX/4+NX, NZ/4+1 : NZ/4+NZ) = U_trans_y
        V_extend(:, NX/4+1 : NX/4+NX, NZ/4+1 : NZ/4+NZ) = V_trans_y
        
        call IXZTrans(U_extend, U_real)
        call IXZTrans(V_extend, V_real)
        W_real = U_real * V_real
        
        call XZTrans(W_real, W_extend)
        W_trans_y = W_extend(:, NX/4+1 : NX/4+NX, NZ/4+1 : NZ/4+NZ)
        
        call YTrans(W_trans_y, W)
        W = 9/4 * W
    end subroutine FCL
    
    subroutine Lamb
        implicit none
        
        call D_DY(WS, rot_x)
        call D_DZ(VS, temp, beta)
        rot_x = rot_x - temp
        call D_DZ(US, rot_y, beta)
        call D_DX(WS, temp, alpha)
        rot_y = rot_y - temp
        call D_DX(VS, rot_z, alpha)
        call D_DY(US, temp)
        rot_z = rot_z - temp
        
        call FCL(VS, rot_z, FX)
        call FCL(WS, rot_y, temp)
        FX = FX - temp
        call FCL(WS, rot_x, FY)
        call FCL(US, rot_z, temp)
        FY = FY - temp
        call FCL(US, rot_y, FZ)
        call FCL(VS, rot_x, temp)
        FZ = FZ - temp
    end subroutine Lamb
    
    subroutine Gterm
        implicit none
        
        call D_DX(rot_z, G, alpha)
        call D_DZ(rot_x, temp, beta)
        G = (G - temp) / re + FY
    end subroutine Gterm
    
end module FieldMatrix