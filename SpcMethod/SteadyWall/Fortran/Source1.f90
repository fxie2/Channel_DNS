    program basic_sp_complex_dft_2d

    use MKL_DFTI, forget => DFTI_SINGLE, DFTI_SINGLE => DFTI_SINGLE_R

    ! Sizes of 2D transform
    integer, parameter :: N1 = 7
    integer, parameter :: N2 = 13

    ! Arbitrary harmonic to test the FFT
    integer, parameter :: H1 = 1
    integer, parameter :: H2 = 1

    ! need single precision
    integer, parameter :: WP = selected_real_kind(6,37)

    ! Execution status
    integer :: status = 0, ignored_status

    complex(WP), allocatable :: x (:,:)

    type(DFTI_DESCRIPTOR), POINTER :: hand

    hand => null()

    print *,"Example basic_sp_complex_dft_2d"
    print *,"Forward and backward single-precision complex-to-complex ",       &
        " in-place 2D transform"
    print *,"Configuration parameters:"
    print *,"DFTI_PRECISION      = DFTI_SINGLE"
    print *,"DFTI_FORWARD_DOMAIN = DFTI_COMPLEX"
    print *,"DFTI_DIMENSION      = 2"
    print '(" DFTI_LENGTHS        = /"I0","I0"/" )', N1, N2

    print *,"Create DFTI descriptor"
    status = DftiCreateDescriptor(hand, DFTI_SINGLE, DFTI_COMPLEX, 2,[N1,N2])
    if (0 /= status) goto 999

    print *,"Commit DFTI descriptor"
    status = DftiCommitDescriptor(hand)
    if (0 /= status) goto 999

    print *,"Allocate array for input data"
    allocate (x(N1, N2))

    print *,"Initialize input for forward transform"
    call init(x, N1, N2, H1, H2)

    print *,"Compute forward transform"
    status = DftiComputeForward(hand, x(:,1))
    if (0 /= status) goto 999

100 continue

    print *,"Release the DFTI descriptor"
    ignored_status = DftiFreeDescriptor(hand)

    print*, x(1,:)
    print*
    print*, x(2,:)
    
    if (allocated(x))  then
        print *,"Deallocate data array"
        deallocate(x)
    endif

999 print '("  Error, status = ",I0)', status
    goto 100

    contains

    ! Initialize array with harmonic /H1, H2/
    subroutine init(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    complex(WP), parameter :: I_TWOPI = (0,6.2831853071795864769_WP)

    forall (k1=1:N1, k2=1:N2)
        x(k1, k2) = k1 + k2
    end forall
    
    end subroutine init
    end program basic_sp_complex_dft_2d
