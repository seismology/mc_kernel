program test_fft

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    integer :: plan_fftf
    integer :: rank, istride, ostride, nextpow2, nsnap, nomega, ntimes
    integer :: ntracesperstep

    real(kind=8), dimension(:,:), allocatable       :: datat
    complex(kind=8), dimension(:,:), allocatable    :: dataf

    rank    = 1
    istride = 1
    ostride = 1

    nsnap = 1000
    ntracesperstep = 100

    nextpow2 = 2
    do while (nextpow2 < nsnap) 
        nextpow2 = nextpow2 * 2
    end do

    nomega = nextpow2 + 1
    ntimes = nextpow2 * 2

    allocate(datat(1:ntimes, 1:ntracesperstep))
    allocate(dataf(1:nomega, 1:ntracesperstep))

    call dfftw_plan_many_dft_r2c(plan_fftf, rank, ntimes, ntracesperstep, datat, &
                                 ntracesperstep, istride, ntimes, dataf, &
                                 ntracesperstep, ostride, nomega, FFTW_ESTIMATE)

    call dfftw_execute(plan_fftf)

    call dfftw_destroy_plan(plan_fftf)

end program
