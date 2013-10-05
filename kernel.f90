module kernel

use global_parameters, only : sp, dp
contains 

function calc_misfit_kernel(timeseries)
   real(kind=sp), intent(in)               :: timeseries(:,:)
   real(kind=sp)                           :: calc_misfit_kernel(size(timeseries,2))

   calc_misfit_kernel = sum(timeseries,1)

end function

function cut_timewindow(t, f, timewindow)
   real(kind=sp), intent(in)  :: t(:), f(:,:)
   real(kind=sp), intent(in)  :: timewindow(2)
   real(kind=sp), allocatable :: cut_timewindow(:,:)
   real(kind=sp), allocatable :: cut_timewindow_temp(:,:)
   integer                    :: ntimes 
   integer                    :: ntraces
   
   ntimes = size(t)
   ntraces = size(f,2)
   dt = t(2) - t(1)

   if (nintimewindow.le.0) then
      write(*,*) 'length of time window is negative'
      write(*,*) 'Beginning: ', timewindow(1), '; end: ', timewindow(2)
      write(*,*) 'nintimewindow:', nintimewindow
      stop
   end if

   allocate(cut_timewindow_temp(ntimes, ntraces))
   iintimewindow = 0
   do i = 1, ntimes
      if (t(i).lt.timewindow(2).and.t(i).gt.timewindow(1)) then
         iintimewindow = iintimewindow + 1
         cut_timewindow_temp(iintimewindow,:) = f(i,:)
      end if
      if (t(i).gt.timewindow(2)) exit
   end do

   if (allocated(cut_timewindow)) deallocate(cut_timewindow)
   allocate(cut_timewindow(iintimewindow, ntraces))

   cut_timewindow(:,:) = cut_timewindow_temp(1:iintimewindow, :)
   

end function

end module
