program test_fft

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    integer :: plan_fft, plan_ifft
    integer :: nextpow2, nsnap, nomega, ntimes, ntimes1
    integer, parameter :: fftw_mode = FFTW_ESTIMATE

    real(8), dimension(:), allocatable       :: datat
    complex(8), dimension(:), allocatable    :: dataf
    real(8), dimension(:), allocatable       :: datat1
    complex(8), dimension(:), allocatable    :: dataf1

    nsnap = 2
    nextpow2 = 2
    do while (nextpow2 < nsnap) 
        nextpow2 = nextpow2 * 2
    end do

    nomega = nextpow2 + 1
    ntimes = nextpow2 * 2
    ntimes1 = nextpow2 * 2

    write(6,*) 'ntimes :', ntimes
    write(6,*) 'ntimes1:', ntimes1
    write(6,*) 'nomega :', nomega

    allocate(datat(1:ntimes))
    allocate(dataf(1:nomega))
    allocate(datat1(1:ntimes1))
    allocate(dataf1(1:nomega))

    call dfftw_plan_dft_r2c_1d(plan_fft, ntimes, datat, dataf, fftw_mode)
    write(6,*) 'ntimes :', ntimes
    write(6,*) 'ntimes1:', ntimes1

    !call dfftw_plan_dft_c2r_1d(plan_ifft, ntimes, dataf1, datat1, fftw_mode)

    datat = 0
    datat(1) = 1
    write(6,'(100(f5.1))') datat(:)
    
    call dfftw_execute_dft_r2c(plan_fft, datat, dataf)
    write(6,'(100("(", f5.1, f5.1, ")"))') dataf(:)

    !call dfftw_execute_dft_c2r(plan_ifft, dataf, datat)
    !write(6,'(100(f5.1))') datat(:)

    !call dfftw_destroy_plan(plan_fft)
    !call dfftw_destroy_plan(plan_ifft)

end program
