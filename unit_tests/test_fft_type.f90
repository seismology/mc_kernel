!=========================================================================================
module test_fft

  use global_parameters
  use fft, only: rfft_type, taperandzeropad
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_fft_dirac
    integer     :: nomega, ntimes, ntraces
    integer     :: i, j

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat
    complex(kind=dp), dimension(:,:), allocatable    :: dataf, dataf_ref

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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_inverse
    integer     :: nomega, ntimes, ntraces
    integer     :: i, j

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat, datat_ref
    complex(kind=dp), dimension(:,:), allocatable    :: dataf

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
        call assert_comparable_real1d(real(datat(:,i), 4) + 1, &
                                      real(datat_ref(:,i), 4) + 1, &
                                      1e-5, 'ifft(fft(datat)) = datat')
    enddo

    call fftt%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

subroutine test_fft_convolve
    integer     :: nomega, ntimes, ntraces, ntimes_fft
    integer     :: i, j

    type(rfft_type) :: fftt

    real(kind=sp), dimension(:,:), allocatable       :: datat1, datat2
    real(kind=dp), dimension(:,:), allocatable       :: dataconv, dataconv_ref
    complex(kind=dp), dimension(:,:), allocatable    :: dataf1, dataf2

    ntimes = 6
    ntraces = 1

    call fftt%init(ntimes, ntraces)
    
    ntimes_fft = fftt%get_ntimes()
    nomega = fftt%get_nomega()

    allocate(datat1(1:ntimes, 1:ntraces))
    allocate(datat2(1:ntimes, 1:ntraces))
    allocate(dataconv(1:ntimes_fft, 1:ntraces))
    allocate(dataconv_ref(1:ntimes_fft, 1:ntraces))

    allocate(dataf1(1:nomega, 1:ntraces))
    allocate(dataf2(1:nomega, 1:ntraces))
   
    datat1(1,:) = 0
    datat1(2,:) = 1
    datat1(3,:) = 1
    datat1(4,:) = 1
    datat1(5,:) = 1
    datat1(6,:) = 0

    datat2 = datat1

    dataconv_ref(:,1) = [0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0]

    call fftt%rfft(taperandzeropad(datat1, ntimes_fft), dataf1)
    call fftt%rfft(taperandzeropad(datat2, ntimes_fft), dataf2)
    call fftt%irfft(dataf1 * dataf2, dataconv)

    do i=1, ntraces
        call assert_comparable_real1d(real(dataconv(:,i), sp) + 1, &
                                      real(dataconv_ref(:,i), sp) + 1, &
                                      1e-5, 'convolution of rectangular pulses')
    enddo

end subroutine


subroutine test_fft_taperandzeropad
    integer, parameter :: len_orig = 32, len_padded = 64
    real(kind=sp) :: data_orig(len_orig,1), data_padded(len_padded,1)
    integer       :: i

    data_orig = 1

    data_padded = taperandzeropad(data_orig, len_padded)

    call assert_equal(len_orig/2, count(abs(data_padded-1)<1e-8), &
                     'length of flat part in tapering window is 50%')

    call assert_comparable_real(real(sum(data_padded(len_orig+1:,1))), 0.0, &
                                  1e-8, 'zeropadded part is really zero')

end subroutine

end module
!=========================================================================================
