module Global_parameter
    implicit none
    integer :: NX, NY, NZ
    integer :: NXH,NYH,NZH
    integer :: NX32, NY32, NZ32
    real,save    :: re
    real,save    :: alpha, beta
    real,save    :: t, dt, start_time, end_time
    integer,save :: step_num
    real,save    :: dpdx, mass_flux
    real,parameter  ::  PI = ATAN(1.0) * 4
    !parameter(NX = 128, NY = 129, NZ = 128)
    contains
    subroutine init_parm
        implicit none
        call init_parm_default
        NXH = NX / 2
        NYH = NY / 2
        NZH = NZ / 2
        NX32 = NX * 3 / 2
        NY32 = NY * 3 / 2
        NZ32 = NZ * 3 / 2
        t = start_time
        step_num = int((end_time - start_time) / dt)
    end subroutine init_parm
    
    subroutine init_parm_user()
        implicit none
        print*, "Please input the mesh size : [NX, NY, NZ] "
        read(*,*) NX, NY, NZ
        print*, "Please input the reynolds number : "
        read(*,*) re
        print*, "Please input the wave number in X and Z direction : "
        read(*,*) alpha, beta
        print*, "Please input the time space, start time and end time : "
        read(*,*) dt, start_time, end_time
        print*, "Please input the pressure gradient and mass flux (0 for disabled) : "
        read(*,*) dpdx, mass_flux
    end subroutine init_parm_user
    
    subroutine init_parm_file()
        implicit none
    end subroutine init_parm_file
    
    subroutine init_parm_default
        implicit none
        NX = 8
        NY = 17
        NZ = 8
        re = 1000
        alpha = 1.0
        beta = 1.0
        dt = 0.01
        start_time = 0
        end_time = dt * 100
        dpdx = 2. / re
        mass_flux = 0.
    end subroutine init_parm_default
    
    end module Global_parameter