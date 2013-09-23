!=====================
  module filtering
!=====================
  use global_parameters
  use data_fft
  implicit none
  public :: define_filter,filter_field,filter_vector,time_shift_field,time_shift_vector
  private
  contains

!-------------------------------------------------------------------------------------
subroutine define_filter(filter_type2,per)
  use data_mesh, only : mynum

  character(len=3),intent(in) :: filter_type2
  real(kind=8)  :: decay,shift_fact,power
  real(kind=realkind), intent(in), optional :: per(2)
  integer :: iomega
  real(kind=realkind) :: filter_period2,omegamid,omegamin2,omegamax2

  if (present(per)) then 
    omegamin=2.*pi/per(2)
    omegamax=2.*pi/per(1)
    filter_period2=per(1)
  else
    filter_period2=filter_period
    omegamin2=omegamin
    omegamax2=omegamax
endif

  filt_power=8
  decay=3.5
  shift_fact=1.
  power = filt_power
  omegamid=(omegamax2-omegamin2)/2.+omegamin2

   if (filter_type2=='low') then ! Low pass
      do iomega=0, nomega
         specfilt(iomega)= 1./(1.+cmplx(0.,1.)*(omega(iomega)/omegamax2))**power
      end do

  elseif (filter_type2=='hig') then  !High pass
      do iomega=0, nomega
         specfilt(iomega)= 1./(1.+cmplx(0.,1.)*(omegamin2/(omega(iomega))  / &
                              (1.+cmplx(0.,1.)*omegamin2/(omega(iomega)) ))) **power
      end do

  elseif (filter_type2=='buh') then  ! Butterworth high pass
      do iomega=0, nomega
         specfilt(iomega)= 1./(1.+ cmplx(0.,1.)* ( omegamin2/(omega(iomega))) **(power/2) )
      end do
     power=1.

  elseif (filter_type2=='bul') then  ! Butterworth low pass
      do iomega=0, nomega
         specfilt(iomega)= 1./(1.+ cmplx(0.,1.)*( omega(iomega)/omegamax2) **(power/2) )
      end do
     power=1.

  elseif (filter_type2=='bub') then  ! Butterworth band reject
      do iomega=0, nomega
         specfilt(iomega)= 1./(1.+ cmplx(0.,1.)*( (omegamax2-omegamin2)* (omega(iomega)-omegamid)  / &
                                                   ( (omega(iomega)-omegamid)**2 - omegamax2**2 )  ) **(power/2) )
      end do
     power=1.

    elseif (filter_type2=='ban') then ! Bandpass
      do iomega=0, nomega
         specfilt(iomega) = 1. / ( 1. + &
              cmplx(0.,1.)*(omegamin2/omega(iomega)+omega(iomega)/omegamax2)**power)
     enddo

     elseif (filter_type2=='ba2') then ! Exponential band pass
        do iomega=0, nomega
           specfilt(iomega) = 1./ ( 1. + exp( (omega(iomega) - omegamax2)/(0.01d0 ) ) ) - &
                                         1. / ( 1. + exp( (omega(iomega)-omegamin2)/(0.01d0 ) ))
      enddo
     power=1.

     elseif (filter_type2=='gal') then ! Gaussian low pass
      do iomega=0, nomega
         specfilt(iomega)  = cmplx(0.,1.)*exp(-( 0.5*omega(iomega)/omegamax2 )**2 )
     enddo 
     power=1.

     elseif (filter_type2=='gah') then  ! Gaussian high pass
      do iomega=0, nomega
         specfilt(iomega)  = 1.d0-cmplx(0.,1.)*exp(-( 0.5*omega(iomega)/omegamin2)**2 )
     enddo 
     power=1.

     elseif (filter_type2=='gab') then  ! Gaussian band reject
      do iomega=0, nomega
         specfilt(iomega)  = 1.d0-cmplx(0.,1.)*exp(-(    ( (omega(iomega)-omegamid)**2 - omegamax2**2  ) /  &
                                                     (omega(iomega)-omegamid)*(omegamax2-omegamin2)        )**2 )
     enddo 
     power=1.

    else
       write(6,*)'Filter',filter_type2,'undefined!'
       stop
  endif

  write(6,*)'                                   '
  write(6,*)'      DEFINE FILTER  FOR ',trim(filter_type2)
  write(6,*)'size omega, filter period:',size(omega),filter_period2
  write(6,*)'omegamin/max:',omegamin2,omegamax2
  write(6,*)'period min/max:',2.*pi/omegamin2,2.*pi/omegamax2

if (mynum==0) then
 open(unit=7000,file="Info/specfilt_"//trim(filter_type2)//".dat")
 open(unit=9000,file="Info/powerspectrum_filter_period"//trim(filter_type2)//".dat")
 open(unit=9000,file="Info/powerspectrum_filter_frequency"//trim(filter_type2)//".dat")
  do iomega=1,nomega
     write(7000,*)2.*pi/omega(iomega),real(specfilt(iomega))
     write(9000,*)2.*pi/omega(iomega),(abs(specfilt(iomega)))**2
     write(8000,*)1./(2.*pi/omega(iomega)),(abs(specfilt(iomega)))**2
  enddo
 close(8000);close(7000); close(9000)
endif

end subroutine define_filter
!-------------------------------------------------------------------------------------


!dk filter_field----------------------------------------------------------
  subroutine filter_field(uspec)
  use data_mesh, only : npts
  complex(kind=8), dimension(0:nomega,1:npts),intent(inout) :: uspec
  integer :: ipt

  do ipt = 1, npts
     uspec(0:nomega,ipt) = specfilt(0:nomega)*uspec(0:nomega,ipt)
  enddo

  end subroutine filter_field
!--------------------------------------------------------------------------


!dk filter_timeseries----------------------------------------------------------
  subroutine filter_vector(uspec)
  complex(kind=8), dimension(0:nomega),intent(inout) :: uspec

  uspec(0:nomega) = specfilt(0:nomega)*uspec(0:nomega)

  end subroutine filter_vector
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine time_shift_field(uspec,ndim,dtshift)
  use data_fft, only : nomega, omega
  use data_mesh, only : npts,mynum

  integer, intent(in) :: ndim
  real(kind=realkind), intent(in) :: dtshift
  complex(kind=8), dimension(0:nomega,1:npts,1:ndim),intent(inout) :: uspec
  complex(kind=8), dimension(0:nomega) :: timeshift_fourier
  integer :: idim,i

   timeshift_fourier(0:nomega) = exp(cmplx(0.,1.)*omega(0:nomega)*dble(dtshift))
   do idim=1,ndim
      do i=1,npts
         uspec(0:nomega,i,idim) = uspec(0:nomega,i,idim) * timeshift_fourier(0:nomega)
      enddo
   enddo
   if (mynum==0 ) write(6,*)'Shifted field by [s]',dtshift

end subroutine time_shift_field
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine time_shift_vector(uspec,dtshift)
  use data_fft, only : nomega, omega
  use data_mesh, only : mynum

  real(kind=realkind), intent(in) :: dtshift
  complex(kind=8), dimension(0:nomega),intent(inout) :: uspec
  complex(kind=8), dimension(0:nomega) :: timeshift_fourier

   timeshift_fourier(0:nomega) = exp(cmplx(0.,1.)*omega(0:nomega)*dble(dtshift))
   uspec(0:nomega) = uspec(0:nomega) * timeshift_fourier(0:nomega)
   if (mynum==0 ) write(6,*)'Shifted vector by [s]',dtshift

end subroutine time_shift_vector
!--------------------------------------------------------------------------


!========================
  end module filtering
!=======================
