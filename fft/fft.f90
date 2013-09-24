program test_fft

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    integer(8)  :: plan_fft, plan_ifft
    integer     :: rank, istride, ostride, nextpow2, nsnap, nomega, ntimes
    integer     :: ntracesperstep
    integer     :: i

    integer, parameter :: fftw_mode = FFTW_ESTIMATE

    ! IMPORTANT: 
    ! these buffer arrays may be used as work arrays in the fft and
    ! hence the original data can be destroyed. Needs verfication for the
    ! specific application!
    real(kind=8), dimension(:,:), allocatable       :: datat
    complex(kind=8), dimension(:,:), allocatable    :: dataf

    rank    = 1
    istride = 1
    ostride = 1

    nsnap = 5
    ntracesperstep = 4

    nextpow2 = 2
    do while (nextpow2 < nsnap) 
        nextpow2 = nextpow2 * 2
    end do

    nomega = nextpow2 + 1
    ntimes = nextpow2 * 2

    write(6,*) 'ntimes:', ntimes
    write(6,*) 'nomega:', nomega

    allocate(datat(1:ntimes, 1:ntracesperstep))
    allocate(dataf(1:nomega, 1:ntracesperstep))

    call dfftw_plan_many_dft_r2c(plan_fft, rank, ntimes, ntracesperstep, datat, &
                                 ntracesperstep, istride, ntimes, dataf, &
                                 ntracesperstep, ostride, nomega, fftw_mode)
    call dfftw_plan_many_dft_c2r(plan_ifft, rank, ntimes, ntracesperstep, dataf, &
                                 ntracesperstep, istride, nomega, datat, &
                                 ntracesperstep, ostride, ntimes, fftw_mode)

    datat = 0
    datat(1,:) = 1

    do i=1, ntracesperstep
       write(6,'(100(f5.1))') datat(:,i)
    enddo

    call dfftw_execute_dft_r2c(plan_fft, datat, dataf)

    do i=1, ntracesperstep
       write(6,'(100("(", f5.1, f5.1, ")"))') dataf(:,i)
    enddo

    datat = 0

    call dfftw_execute_dft_c2r(plan_ifft, dataf, datat)

    datat = datat / ntimes

    do i=1, ntracesperstep
       write(6,'(100(f5.1))') datat(:,i)
    enddo

    call dfftw_destroy_plan(plan_fft)
    call dfftw_destroy_plan(plan_ifft)

end program
