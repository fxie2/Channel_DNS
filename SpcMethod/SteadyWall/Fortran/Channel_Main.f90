program Channel_DNS
    use Global_parameter
    use FieldMatrix
    implicit none
    call init_parm()
    call alloc_field()
    call init_field()
    read*
    call dealloc_field()
    end program Channel_DNS