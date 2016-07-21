module FieldMatrix
    use Global_parameter
    implicit none
    
    !Real field
    real(8), save, allocatable :: U(:, :, :), V(:, :, :), W(:, :, :)
    real(8), save, allocatable :: P(:, :, :)
    !Spectral field
    !_S for s step, _1 for s-1 step, _2 for s-2 step, _NEW for s+1 step
    
    !Velocity field
    complex(8), save, allocatable :: US(:, :, :), VS(:, :, :), WS(:, :, :)
    complex(8), save, allocatable :: U1(:, :, :), V1(:, :, :), W1(:, :, :)
    complex(8), save, allocatable :: U2(:, :, :), V2(:, :, :), W2(:, :, :)
    complex(8), save, allocatable :: UNEW(:, :, :), VNEW(:, :, :), WNEW(:, :, :)
    
    !Non-linear term
    complex(8), save, allocatable :: FX(:, :, :), FY(:, :, :), FZ(:, :, :)
    complex(8), save, allocatable :: FX1(:, :, :), FY1(:, :, :), FZ1(:, :, :)
    complex(8), save, allocatable :: FX2(:, :, :), FY2(:, :, :), FZ2(:, :, :)
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
    end subroutine alloc_field
    
    subroutine dealloc_field()
        implicit none
        deallocate(U, V, W, P)
        deallocate(US, VS, WS)
        deallocate(U1, V1, W1, U2, V2, W2, UNEW, VNEW, WNEW)
        deallocate(FX, FY, FZ, FX1, FY1, FZ1, FX2, FY2, FZ2)
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
    end subroutine init_field
    
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
end module FieldMatrix