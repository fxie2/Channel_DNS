module FieldMatrix
    use Global_parameter
    use Math
    implicit none
    
    !Real field
    real, save, allocatable :: U(:, :, :), V(:, :, :), W(:, :, :)
    real, save, allocatable :: P(:, :, :)
    !Spectral field
    !_S for s step, _1 for s-1 step, _2 for s-2 step, _NEW for s+1 step
    
    !Velocity field
    complex, save, pointer :: US(:, :, :), VS(:, :, :), WS(:, :, :)
    complex, save, pointer :: U1(:, :, :), V1(:, :, :), W1(:, :, :)
    complex, save, pointer :: U2(:, :, :), V2(:, :, :), W2(:, :, :)
    complex, save, pointer :: UNEW(:, :, :), VNEW(:, :, :), WNEW(:, :, :)
    complex, save, allocatable, target :: US_R(:, :, :), VS_R(:, :, :), WS_R(:, :, :)
    complex, save, allocatable, target :: U1_R(:, :, :), V1_R(:, :, :), W1_R(:, :, :)
    complex, save, allocatable, target :: U2_R(:, :, :), V2_R(:, :, :), W2_R(:, :, :)
    complex, save, allocatable, target :: UNEW_R(:, :, :), VNEW_R(:, :, :), WNEW_R(:, :, :)
    
    !Pressure field
    complex, save, allocatable :: PS(:, :, :), s(:, :, :) !source term for solve
    
    !Rotation field
    complex, save, allocatable :: rot_x(:, :, :), rot_y(:, :, :), rot_z(:, :, :)
    
    !Temp space
    complex, save, allocatable :: temp(:, :, :)
    
    !Non-linear term
    complex, save, pointer :: FX(:, :, :), FY(:, :, :), FZ(:, :, :)
    complex, save, pointer :: FX1(:, :, :), FY1(:, :, :), FZ1(:, :, :)
    complex, save, pointer :: FX2(:, :, :), FY2(:, :, :), FZ2(:, :, :)
    complex, save, allocatable, target :: FX_R(:, :, :), FY_R(:, :, :), FZ_R(:, :, :)
    complex, save, allocatable, target :: FX1_R(:, :, :), FY1_R(:, :, :), FZ1_R(:, :, :)
    complex, save, allocatable, target :: FX2_R(:, :, :), FY2_R(:, :, :), FZ2_R(:, :, :)
    
    !G term
	complex, save, pointer :: G(:, :, :), G1(:, :, :), G2(:, :, :)
	complex, save, allocatable, target :: G_R(:, :, :), G1_R(:, :, :), G2_R(:, :, :)
    
    !Boundary term
    complex, save, allocatable :: boundary(:, :, :), upper_bd(:, :), lower_bd(:, :)
    contains
    subroutine alloc_field()
        implicit none
        
        allocate(U(NY, NX, NZ))
        allocate(V(NY, NX, NZ))
        allocate(W(NY, NX, NZ))
        allocate(P(NY, NX, NZ))
        
        allocate(US_R(NY, NX, NZ))
        allocate(VS_R(NY, NX, NZ))
        allocate(WS_R(NY, NX, NZ))
        allocate(PS(NY, NX, NZ))
        allocate(U1_R(NY, NX, NZ))
        allocate(V1_R(NY, NX, NZ))
        allocate(W1_R(NY, NX, NZ))
        allocate(U2_R(NY, NX, NZ))
        allocate(V2_R(NY, NX, NZ))
        allocate(W2_R(NY, NX, NZ))
        allocate(UNEW_R(NY, NX, NZ))
        allocate(VNEW_R(NY, NX, NZ))
        allocate(WNEW_R(NY, NX, NZ))
        
        allocate(FX_R(NY, NX, NZ))
        allocate(FY_R(NY, NX, NZ))
        allocate(FZ_R(NY, NX, NZ))
        allocate(FX1_R(NY, NX, NZ))
        allocate(FY1_R(NY, NX, NZ))
        allocate(FZ1_R(NY, NX, NZ))
        allocate(FX2_R(NY, NX, NZ))
        allocate(FY2_R(NY, NX, NZ))
        allocate(FZ2_R(NY, NX, NZ))
        
        allocate(G_R(NY, NX, NZ))
        allocate(G1_R(NY, NX, NZ))
        allocate(G2_R(NY, NX, NZ))
        allocate(rot_x(NY, NX, NZ))
        allocate(rot_y(NY, NX, NZ))
        allocate(rot_z(NY, NX, NZ))
        allocate(boundary(NY, NX, NZ))
        allocate(upper_bd(NX, NZ))
        allocate(lower_bd(NX, NZ))
        allocate(temp(NY, NX, NZ))
        allocate(s(NY, NX, NZ))
    end subroutine alloc_field
    
    subroutine dealloc_field()
        implicit none
        deallocate(U, V, W, P)
        deallocate(US_R, VS_R, WS_R, PS)
        deallocate(U1_R, V1_R, W1_R, U2_R, V2_R, W2_R, UNEW_R, VNEW_R, WNEW_R)
        deallocate(FX_R, FY_R, FZ_R, FX1_R, FY1_R, FZ1_R, FX2_R, FY2_R, FZ2_R)
        deallocate(rot_x, rot_y, rot_z)
        deallocate(boundary, upper_bd, lower_bd)
        deallocate(G_R, G1_R, G2_R)
        deallocate(temp, s)
        call MathKernel_finalize
    end subroutine dealloc_field
        
    subroutine init_field()
        implicit none
        integer iter_x, iter_y
        real turb
        if(allocated(U) == .false.) call alloc_field()
        call MathKernel_init(NX, NY, NZ)
        US => US_R
        VS => VS_R
        WS => WS_R
        U1 => U1_R
        V1 => V1_R
        W1 => W1_R
        U2 => U2_R
        V2 => V2_R
        W2 => W2_R
        UNEW => UNEW_R
        VNEW => VNEW_R
        WNEW => WNEW_R
        FX => FX_R
        FY => FY_R
        FZ => FZ_R
        FX1 => FX1_R
        FY1 => FY1_R
        FZ1 => FZ1_R
        FX2 => FX2_R
        FY2 => FY2_R
        FZ2 => FZ2_R
        G => G_R
        G1 => G1_R
        G2 => G2_R
        
        U = 0
        V = 0
        W = 0
        P = 0
        !Velocity initialize
        call RANDOM_SEED()
        do iter_y = 1, NY
            call RANDOM_NUMBER(turb)
            turb = turb - 0.5
            print*, turb
            U(iter_y, :, :) = 1 - cos((iter_y - 1) * PI / (NY - 1)) ** 2 + 1e-5 * turb
        end do
        call GTrans(U, US)
        call GTrans(V, VS)
        call GTrans(W, WS)
        U1 = US
        U2 = US
        V1 = VS
        V2 = VS
        W1 = WS
        W2 = WS
        call Lamb
        FX(1, NX / 2 + 1, NZ / 2 + 1) = FX(1, NX / 2 + 1, NZ / 2 + 1) + dpdx * NX * NZ
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
    
    subroutine Calculate
        implicit none
        
        call Mass_correct
        call NonLiner_step
        call Pressure_step
        call Viscos_step
        call Div_correct
        call Update_field
    end subroutine Calculate
    
    subroutine Mass_correct
        implicit none
    end subroutine Mass_correct
    
    subroutine NonLiner_step
        implicit none
        
        call Lamb
        FX(1, NX / 2 + 1, NZ / 2 + 1) = FX(1, NX / 2 + 1, NZ / 2 + 1) + dpdx * NX * NZ
        
        UNEW = 3 * US - 1.5 * U1 + 1. / 3 * U2 + dt * (3 * FX - 3 * FX1 + FX2)
        VNEW = 3 * VS - 1.5 * V1 + 1. / 3 * V2 + dt * (3 * FY - 3 * FY1 + FY2)
        WNEW = 3 * WS - 1.5 * W1 + 1. / 3 * W2 + dt * (3 * FZ - 3 * FZ1 + FZ2)
    end subroutine NonLiner_step
    
    subroutine Pressure_step
        implicit none
        integer i, k
        
        call Gterm
        boundary = 3 * G - 3 * G1 + G2
        call IYTrans(boundary, boundary)
        upper_bd = boundary(1, :, :)
        lower_bd = boundary(NY, :, :)
        call D_DX(UNEW, s, alpha)
        call D_DY(VNEW, temp)
        s = s + temp
        call D_DZ(WNEW, temp, beta)
        s = s + temp
        s = s / dt
        do k = 1, NZ
            do i = 1, NX
                call ode_solve(alpha ** 2 * (i - 1 - NXH) ** 2 + beta ** 2 * (k - 1 - NZH) ** 2,&
                               PS(:, i, k), s(:, i, k), 2, lower_bd(i, k), upper_bd(i, k))
            end do
        end do
        call D_DX(PS, temp, alpha)
        UNEW = UNEW - dt * temp
        call D_DY(PS, temp)
        VNEW = VNEW - dt * temp
        call D_DZ(PS, temp, beta)
        WNEW = WNEW - dt * temp
    end subroutine Pressure_step
    
    subroutine Viscos_step
        implicit none
        integer i, k
        real fac
        
        do k = 1, NZ
            do i = 1, NX
                fac = alpha ** 2 * (i - 1 - NXH) ** 2 + beta ** 2 * (k - 1 - NZH) ** 2 + 11. / 6 * re / dt
                call ode_solve(fac, UNEW(:, i, k), -re / dt * UNEW(:, i, k), 1, (0, 0), (0, 0))
                call ode_solve(fac, VNEW(:, i, k), -re / dt * VNEW(:, i, k), 1, (0, 0), (0, 0))
                call ode_solve(fac, WNEW(:, i, k), -re / dt * WNEW(:, i, k), 1, (0, 0), (0, 0))
            end do
        end do
    end subroutine Viscos_step
    
    subroutine Div_correct
        implicit none
    end subroutine Div_correct
    
    subroutine Update_field
        implicit none
        complex, pointer :: pt(:, :, :)
        
        pt => U2
        U2 => U1
        U1 => US
        US => UNEW
        UNEW => pt
        pt => V2
        V2 => V1
        V1 => VS
        VS => VNEW
        VNEW => pt
        pt => W2
        W2 => W1
        W1 => WS
        WS => WNEW
        WNEW => pt
        pt => FX2
        FX2 => FX1
        FX1 => FX
        FX => pt
        pt => FY2
        FY2 => FY1
        FY1 => FY
        FY => pt
        pt => FZ2
        FZ2 => FZ1
        FZ1 => FZ
        FZ => pt
        pt => G2
        G2 => G1
        G1 => G
        G => pt
        
        call IGTrans(PS, P)
        call IGTrans(US, U)
        call IGTrans(VS, V)
        call IGTrans(WS, W)
    end subroutine Update_field
    
    subroutine Output
        implicit none
    end subroutine Output
    
end module FieldMatrix