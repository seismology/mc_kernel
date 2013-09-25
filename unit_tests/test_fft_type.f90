module test_fft

  use fft
  use ftnunit
  implicit none
  public
contains

subroutine test_fft_dirac
    integer     :: nomega, ntimes, ntraces
    integer     :: i, j

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

    ! dirac
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

    ! phase shifted dirac
    datat = 0
    datat(2,:) = 1

    call fftt%rfft(datat, dataf)
    
    do i=1, ntraces
        call assert_comparable_real1d(real(abs(dataf(:,i)), 4), &
                    real(dataf_ref(:,i), 4), 1e-5, &
                    'fft of dirac = 1 + 0j * phase shift')
    enddo


    call fftt%freeme()
end subroutine

subroutine test_fft_inverse
    integer     :: nomega, ntimes, ntraces
    integer     :: i, j

    type(rfft_type) :: fftt

    real(kind=8), dimension(:,:), allocatable       :: datat, datat_ref
    complex(kind=8), dimension(:,:), allocatable    :: dataf

    ntimes = 4
    ntraces = 1

    call fftt%init(ntimes, ntraces)

    ntimes = fftt%get_ntimes()
    nomega = fftt%get_nomega()

    allocate(datat(1:ntimes, 1:ntraces))
    allocate(datat_ref(1:ntimes, 1:ntraces))
    allocate(dataf(1:nomega, 1:ntraces))

    ! some random data
    datat = 0
    datat(1,:) = 1
    datat(2,:) = 2
    datat(3,:) = 0
    datat(4,:) = 1

    datat_ref(:,:) = datat(:,:)

    call fftt%rfft(datat, dataf)
    call fftt%irfft(dataf, datat_ref)

    do i=1, ntraces
        call assert_comparable_real1d(real(datat(:,i), 4) + 1, real(datat_ref(:,i), 4) + 1, &
                1e-5, 'ifft(fft(datat)) = datat')
    enddo

    call fftt%freeme()
end subroutine

end module
