!=========================================================================================
module test_fft_type

  use global_parameters
  use fft, only: rfft_type, taperandzeropad
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_fft_dirac
    integer     :: nomega, ntimes, ntraces
    integer     :: i

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat
    complex(kind=dp), dimension(:,:), allocatable    :: dataf, dataf_ref

    ntimes = 4
    ntraces = 2

    call fftt%init(ntimes, 1, ntraces)

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
        call assert_comparable_real1d(real(dataf(:,i), sp), real(dataf_ref(:,i), sp), &
                1e-5, 'fft of dirac = 1 + 0j')
        call assert_comparable_real1d(real(aimag(dataf(:,i)), sp), &
                real(aimag(dataf_ref(:,i)), sp), 1e-5, 'fft of dirac = 1 + 0j')
    enddo

    ! phase shifted dirac
    datat = 0
    datat(2,:) = 1

    call fftt%rfft(datat, dataf)
    
    do i=1, ntraces
        call assert_comparable_real1d(real(abs(dataf(:,i)), sp), &
                    real(dataf_ref(:,i), 4), 1e-5, &
                    'fft of dirac = 1 + 0j * phase shift')
    enddo


    call fftt%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_sine
    integer     :: nomega, ntimes, ntraces
    integer     :: i

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat
    complex(kind=dp), dimension(:,:), allocatable    :: dataf, dataf_ref
    real(kind=dp), dimension(:,:), allocatable       :: T
    real(kind=dp), dimension(:), allocatable         :: F
    real(kind=dp)                                    :: dt, df
    real(kind=sp)                                    :: f_max, f_ref

    ntimes = 100
    ntraces = 1

    dt = 0.1
    call fftt%init(ntimes, 1, ntraces, dt)

    ntimes = fftt%get_ntimes()
    nomega = fftt%get_nomega()

    allocate(datat(1:ntimes, 1:ntraces))
    allocate(dataf(1:nomega, 1:ntraces))
    allocate(T(1:ntimes, 1:ntraces))
    allocate(F(1:nomega))


    do i = 1, ntimes
       T(i,1) = dt * i
    end do

    !df = 1. / (dt*fftt%get_ntimes())

    !do i = 1, nomega
    !   F(i) = df * (i-1)
    !end do
    df = fftt%get_df()
    F = fftt%get_f()

    ! Sine with period 1s
    f_ref = 1.0
    datat = sin(2*pi*T * f_ref)

    call fftt%rfft(taperandzeropad(datat, fftt%get_ntimes()), dataf)

    f_max = F( maxloc( abs(dataf(:,1)), 1 ) )

    call assert_comparable_real(f_max, f_ref, &
                                real(dF/2), 'Maximum amplitude at 1 Hz (Sine)')

    ! Cosine with period 1s
    f_ref = 1.0
    datat = cos(2*pi*T * f_ref)

    call fftt%rfft(taperandzeropad(datat, fftt%get_ntimes()), dataf)

    f_max = F( maxloc( abs(dataf(:,1)), 1 ) )

    call assert_comparable_real(f_max, f_ref, &
                                real(dF/2), 'Maximum amplitude at 1 Hz (Cosine)')

    ! Sine with period 2s
    f_ref = 0.5
    datat = sin(2*pi*T * f_ref)

    call fftt%rfft(taperandzeropad(datat, fftt%get_ntimes()), dataf)

    f_max = F( maxloc( abs(dataf(:,1)), 1 ) )

    call assert_comparable_real(f_max, f_ref, &
                                real(dF/2), 'Maximum amplitude at 0.5 Hz')

    ! Sine with period 0.5s
    f_ref = 2
    datat = sin(2*pi*T * f_ref)

    call fftt%rfft(taperandzeropad(datat, fftt%get_ntimes()), dataf)

    f_max = F( maxloc( abs(dataf(:,1)), 1 ) )

    call assert_comparable_real(f_max, f_ref, &
                                real(dF/2), 'Maximum amplitude at 2 Hz')

    call fftt%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_inverse
    integer     :: nomega, ntimes, ntraces
    integer     :: i

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat, datat_ref
    complex(kind=dp), dimension(:,:), allocatable    :: dataf

    ntimes = 4
    ntraces = 1

    call fftt%init(ntimes, 1, ntraces)

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
    call fftt%irfft(dataf, datat)

    do i=1, ntraces
        call assert_comparable_real1d(real(datat(:,i), 4) + 1, &
                                      real(datat_ref(:,i), 4) + 1, &
                                      1e-5, 'ifft(fft(datat)) = datat')
    enddo

    call fftt%freeme()
end subroutine test_fft_inverse
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_inverse_md
    integer     :: nomega, ntimes, ntraces
    integer     :: i

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat, datat_ref
    complex(kind=dp), dimension(:,:), allocatable    :: dataf

    ntimes = 16
    ntraces = 10

    call fftt%init(ntimes, 1, ntraces)

    ntimes = fftt%get_ntimes()
    nomega = fftt%get_nomega()

    allocate(datat(1:ntimes, 1:ntraces))
    allocate(datat_ref(1:ntimes, 1:ntraces))
    allocate(dataf(1:nomega, 1:ntraces))

    ! some random data
    datat = 0
    datat(5,:)  = 1
    datat(6,:)  = 2
    datat(7,:)  = 0
    datat(8,:)  = 1
    datat(9,:)  = 1
    datat(10,:) = 2
    datat(11,:) = 0
    datat(12,:) = 1

    datat_ref(:,:) = datat(:,:)

    call fftt%rfft(datat_ref, dataf)
    call fftt%irfft(dataf, datat)

    do i=1, ntraces
        call assert_comparable(datat(:,i) + 1, &
                               datat_ref(:,i) + 1, &
                               1d-5, 'ifft(fft(datat)) = datat')
    enddo

    call fftt%freeme()
end subroutine test_fft_inverse_md
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_inverse_taz
    integer             :: nomega, ntimes
    integer, parameter  :: ntimes_orig = 16, ntraces = 10
    integer             :: i

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat, datat_ref
    complex(kind=dp), dimension(:,:), allocatable    :: dataf


    call fftt%init(ntimes_orig, 1, ntraces)

    ntimes = fftt%get_ntimes()
    nomega = fftt%get_nomega()

    allocate(datat_ref(1:ntimes_orig, 1:ntraces))
    allocate(datat(1:ntimes, 1:ntraces))
    allocate(dataf(1:nomega, 1:ntraces))

    ! some random data
    datat_ref = 0
    datat_ref(5,:)  = 1
    datat_ref(6,:)  = 1
    datat_ref(7,:)  = 1
    datat_ref(8,:)  = 1
    datat_ref(9,:)  = 1
    datat_ref(10,:) = 1
    datat_ref(11,:) = 1
    datat_ref(12,:) = 1

    datat(1:ntimes_orig,:) = datat_ref(:,:)

    call fftt%rfft(taperandzeropad(datat_ref, ntimes), dataf)
    call fftt%irfft(dataf, datat)

    do i=1, ntraces
        call assert_comparable(datat(1:ntimes_orig,i) + 1, &
                               datat_ref(:,i) + 1, &
                               1d-5, 'ifft(fft(datat)) = datat')
    enddo

    call fftt%freeme()
end subroutine test_fft_inverse_taz
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_convolve
    integer     :: nomega, ntimes, ntraces, ntimes_fft
    integer     :: i

    type(rfft_type) :: fftt

    real(kind=dp), dimension(:,:), allocatable       :: datat1, datat2
    real(kind=dp), dimension(:,:), allocatable       :: dataconv, dataconv_ref
    complex(kind=dp), dimension(:,:), allocatable    :: dataf1, dataf2

    ntimes = 6
    ntraces = 1

    call fftt%init(ntimes, ndim = 1, ntraces=ntraces)
    
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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_fft_taperandzeropad
    integer, parameter :: len_orig = 32, len_padded = 64
    real(kind=dp) :: data_orig(len_orig,1), data_padded(len_padded,1)

    data_orig = 1

    data_padded = taperandzeropad(data_orig, len_padded)

    call assert_comparable(data_orig(9:24,1), data_padded(9:24,1), 1d-8, &
                           'Central part of tapered signal remains the same')

    call assert_equal(len_orig/2, count(abs(data_padded-1)<1e-8), &
                     'length of flat part in tapering window is 50%')

    call assert_comparable_real(real(sum(data_padded(len_orig+1:,1))), 0.0, &
                                  1e-8, 'zeropadded part is really zero')

end subroutine
!-----------------------------------------------------------------------------------------

end module test_fft_type
!=========================================================================================
