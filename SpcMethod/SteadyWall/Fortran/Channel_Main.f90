program Channel_DNS
    use Global_parameter
    use FieldMatrix
    implicit none
    call init_parm()
    call init_field()
    do while(t < end_time)
        call Calculate
        call Output
        t = t + dt
    end do
    print*, maxval(U)
    call dealloc_field()
    end program Channel_DNS