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