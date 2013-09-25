module test_fft

  use fft
  use ftnunit
  implicit none
  public
contains

subroutine test_fft_dirac
    integer     :: nomega, ntimes, ntraces
    integer     :: i

    type(rfft_type) :: fftt

    real(kind=8), dimension(:,:), allocatable       :: datat
    complex(kind=8), dimension(:,:), allocatable    :: dataf, dataf_ref

    ntimes = 4
    ntraces = 2

    call fftt%init(ntimes, ntraces)

    ntimes = fftt%get_ntimes()
    nomega = fftt%get_nomega()

    allocate(datat(1:ntimes, 1:ntraces))
    allocate(dataf(1:nomega, 1:ntraces))

    datat = 0
    datat(1,:) = 1

    call fftt%rfft(datat, dataf)

    allocate(dataf_ref(1:nomega, 1:ntraces))
    dataf_ref = 1
    
    do i=1, ntraces
        call assert_comparable_real1d(real(dataf(:,i), 4), real(dataf_ref(:,i), 4), &
                1e-5, 'fft of dirac = 1 + 0j')
        call assert_comparable_real1d(real(imagpart(dataf(:,i)), 4), &
                real(imagpart(dataf_ref(:,i)), 4), 1e-5, 'fft of dirac = 1 + 0j')
    enddo

    call fftt%freeme()
end subroutine

end module
