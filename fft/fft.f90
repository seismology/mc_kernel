program test_fft

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    integer(8)  :: plan_fft, plan_ifft
    integer     :: rank, istride, ostride, nextpow2, nsnap, nomega, ntimes
    integer     :: ntracesperstep
    integer     :: i

    integer, parameter :: fftw_mode = FFTW_ESTIMATE

    real(kind=8), dimension(:,:), allocatable       :: datat
    complex(kind=8), dimension(:,:), allocatable    :: dataf
    real(kind=8), dimension(:,:), allocatable       :: datat1
    complex(kind=8), dimension(:,:), allocatable    :: dataf1

    nsnap = 2
    ntracesperstep = 2

    nextpow2 = 2
    do while (nextpow2 < nsnap) 
        nextpow2 = nextpow2 * 2
    end do

    nomega = nextpow2 + 1
    ntimes = nextpow2 * 2

    write(6,*) 'ntimes:', ntimes
    write(6,*) 'nomega:', nomega

    rank    = 1
    istride = 1
    ostride = 1

    call dfftw_plan_many_dft_r2c(plan_fft, rank, ntimes, ntracesperstep, datat, &
                                 ntracesperstep, istride, ntimes, dataf, &
                                 ntracesperstep, ostride, nomega, fftw_mode)
    call dfftw_plan_many_dft_c2r(plan_ifft, rank, ntimes, ntracesperstep, dataf, &
                                 ntracesperstep, istride, nomega, datat, &
                                 ntracesperstep, ostride, ntimes, fftw_mode)

    allocate(datat(1:ntimes, 1:ntracesperstep))
    allocate(dataf(1:nomega, 1:ntracesperstep))
    allocate(datat1(1:ntimes, 1:ntracesperstep))
    allocate(dataf1(1:nomega, 1:ntracesperstep))


    datat = 0
    !datat(2,:) = 1
    datat1 = 0
    datat1(2,:) = 1

    do i=1, ntracesperstep
       write(6,'(100(f5.1))') datat(:,i)
    enddo

    do i=1, ntracesperstep
       write(6,'(100(f5.1))') datat1(:,i)
    enddo

    ! use specific interfaces including the buffer arrays to make sure the
    ! compiler does not change the order of execution (stupid, but known issue)
    call dfftw_execute_dft_r2c(plan_fft, datat, dataf)
    call dfftw_execute_dft_r2c(plan_fft, datat1, dataf1)


    do i=1, ntracesperstep
       write(6,'(100("(", f5.1, f5.1, ")"))') dataf(:,i)
    enddo

    do i=1, ntracesperstep
       write(6,'(100("(", f5.1, f5.1, ")"))') dataf1(:,i)
    enddo

    datat = 0

    ! use specific interfaces including the buffer arrays to make sure the
    ! compiler does not change the order of execution (stupid, but known issue)
    call dfftw_execute_dft_c2r(plan_ifft, dataf, datat)
    call dfftw_execute_dft_c2r(plan_ifft, dataf1, datat1)

    ! normalization, see
    ! http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
    datat = datat / ntimes
    datat1 = datat1 / ntimes

    do i=1, ntracesperstep
       write(6,'(100(f5.1))') datat(:,i)
    enddo

    do i=1, ntracesperstep
       write(6,'(100(f5.1))') datat1(:,i)
    enddo

    call dfftw_destroy_plan(plan_fft)
    call dfftw_destroy_plan(plan_ifft)

end program
