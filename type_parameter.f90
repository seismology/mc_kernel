module type_parameter
    use global_parameters,               only : sp, dp, pi, deg2rad, verbose, lu_out
    use source_class,                    only : src_param_type
    use kernel,                          only : kernelspec_type
    use receiver_class,                  only : rec_param_type
    use filtering,                       only : filter_type
    use commpi,                          only : pbroadcast_char, pbroadcast_int, pbroadcast_dble, &
                                                pbroadcast_dble_arr, master
    implicit none    

    type parameter_type
        type(src_param_type)                 :: source
        type(rec_param_type),   allocatable  :: receiver(:)
        type(kernelspec_type),  allocatable  :: kernel(:)
        type(filter_type),      allocatable  :: filter(:)
        integer                              :: nrec

        real(kind=dp)                        :: allowed_error
        character(len=512)                   :: fwd_dir
        character(len=512)                   :: bwd_dir
        character(len=512)                   :: source_file
        character(len=512)                   :: receiver_file
        character(len=512)                   :: filter_file
        character(len=512)                   :: mesh_file
        character(len=1)                     :: component
        character(len=4)                     :: model_parameter
        character(len=32)                    :: whattodo
        integer                              :: nsim_fwd, nsim_bwd
        integer                              :: nkernel
        integer                              :: nelems_per_task
        integer                              :: npoints_per_step
        integer                              :: max_iter
        contains
           procedure, pass                   :: read_parameters
           procedure, pass                   :: read_receiver
           procedure, pass                   :: read_source
           procedure, pass                   :: read_kernel
           procedure, pass                   :: read_filter
    end type

!    type kernelspec_type
!        real, dimension(2)                   :: time_window
!        real, dimension(4)                   :: corner_frequencies
!        integer                              :: filter_type
!        character                            :: misfit_type
!        character                            :: model_parameter
!        integer                              :: receiver_index
!        integer                              :: src_index
!        type(rec_param_type), pointer        :: receiver
!        !pointer                              :: filter
!    end type

contains

!------------------------------------------------------------------------------
subroutine read_parameters(this, input_file)
   class(parameter_type)           :: this
   character(len=*), intent(in)    :: input_file
   integer                         :: iinparam_basic=500, ioerr
   character(len=256)              :: line
   character(len=256)              :: keyword, keyvalue


   if (master) then
     if (verbose > 0) write(lu_out,'(A)', advance='no') '    Reading inparam_basic...'
     open(unit=iinparam_basic, file=trim(input_file), status='old', action='read',  iostat=ioerr)
     if (ioerr /= 0) then
        print *, 'ERROR: Check input file ''', trim(input_file), '''! Is it still there?' 
        stop
     end if

     do
        read(iinparam_basic, fmt='(a256)', iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
     
        read(line,*) keyword, keyvalue 
      
        parameter_to_read : select case(trim(keyword))
        case('ALLOWED_ERROR')
           read(keyvalue, *) this%allowed_error

        case('FWD_DIR')
           this%fwd_dir = keyvalue

        case('BWD_DIR')
           this%bwd_dir = keyvalue

        case('SOURCE_FILE')
           this%source_file = keyvalue

        case('RECEIVER_FILE')
           this%receiver_file = keyvalue

        case('FILTER_FILE')
           this%filter_file = keyvalue

        case('MESH_FILE')
           this%mesh_file = keyvalue

        case('MAXIMUM_ITERATIONS')
           read(keyvalue, *) this%max_iter

        case('ELEMENTS_PER_TASK')
           read(keyvalue, *) this%nelems_per_task

        case('POINTS_PER_MC_STEP')
           read(keyvalue, *) this%npoints_per_step

        case('WHAT_TO_DO')
           this%whattodo = keyvalue

        end select parameter_to_read

     end do
  end if

  
  ! broadcast values to other processors
  call pbroadcast_char(this%fwd_dir, 0) 
  call pbroadcast_char(this%bwd_dir, 0) 
  call pbroadcast_dble(this%allowed_error, 0)
  call pbroadcast_char(this%source_file, 0) 
  call pbroadcast_char(this%receiver_file, 0)
  call pbroadcast_char(this%filter_file, 0)
  call pbroadcast_char(this%mesh_file, 0)
  call pbroadcast_int(this%max_iter, 0)
  call pbroadcast_int(this%nelems_per_task, 0)
  call pbroadcast_int(this%npoints_per_step, 0)
  call pbroadcast_char(this%whattodo, 0)

  call flush(6)

end subroutine read_parameters
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine read_source(this)
   class(parameter_type)          :: this
   real(kind=dp)                  :: Mij_dyncm(6), latd, lond, depth
   character(len=16)              :: junk
   character(len=16)              :: event_name
   integer, parameter             :: lu_source=1000

   if (master) then
       write(lu_out,*)'  reading source from file ', trim(this%source_file)
       open(unit=lu_source, file=trim(this%source_file), status='old')
       read(lu_source,*) junk
       read(lu_source,*) junk, event_name
       read(lu_source,*) junk
       read(lu_source,*) junk
       read(lu_source,*) junk, latd
       read(lu_source,*) junk, lond  
       read(lu_source,*) junk, depth
       read(lu_source,*) junk, Mij_dyncm(1) !Mrr
       read(lu_source,*) junk, Mij_dyncm(2) !Mtt
       read(lu_source,*) junk, Mij_dyncm(3) !Mpp
       read(lu_source,*) junk, Mij_dyncm(4) !Mrt
       read(lu_source,*) junk, Mij_dyncm(5) !Mrp
       read(lu_source,*) junk, Mij_dyncm(6) !Mtp
       close(lu_source)
   end if

   call pbroadcast_dble(latd, 0)
   call pbroadcast_dble(lond, 0)
   call pbroadcast_dble(depth, 0)
   call pbroadcast_dble_arr(Mij_dyncm, 0)
   

   call this%source%init(lat = latd, lon = lond, mij = Mij_dyncm*1.E-7)
   call flush(6)

end subroutine read_source
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine read_receiver(this)
   class(parameter_type)         :: this
   integer, parameter            :: lu_receiver = 1001
   integer                       :: irec, firstkernel, lastkernel
   integer                       :: ikernel, recnkernel
   real(kind=dp)                 :: timewindow(2), reclatd, reclond
   character(len=16)             :: recname, kernelname, filtername
   character(len=80)             :: fmtstring

   if (master) then
       write(lu_out,*)'  reading receivers from file ', trim(this%receiver_file)
       open(unit=lu_receiver, file=trim(this%receiver_file), status='old')
       read(lu_receiver,*) this%nrec
       read(lu_receiver,*) this%model_parameter, this%component
   end if

   call pbroadcast_int(this%nrec, 0)
   call pbroadcast_char(this%model_parameter, 0)
   call pbroadcast_char(this%component, 0)

   if (master) then
       fmtstring = '("  Using ", I5, " receivers")'
       write(lu_out, fmtstring) this%nrec
   
       fmtstring = '("  Kernel for parameter ", A, " on component ", A)'
       write(lu_out, fmtstring) this%model_parameter, this%component
   end if

   allocate(this%receiver(this%nrec))

   firstkernel = 1
   lastkernel  = 0
   do irec = 1, this%nrec
      if (master) then
          read(lu_receiver, *) recname, reclatd, reclond, recnkernel
      end if
      call pbroadcast_char(recname, 0)
      call pbroadcast_dble(reclatd, 0)
      call pbroadcast_dble(reclond, 0)
      call pbroadcast_int(recnkernel, 0)

      lastkernel = lastkernel + recnkernel

      call this%receiver(irec)%init(name        = recname        , &
                                    lat         = reclatd        , &
                                    lon         = reclond        , &
                                    component   = this%component , &
                                    nkernel     = recnkernel     , &  
                                    firstkernel = firstkernel    , &
                                    lastkernel  = lastkernel       )
      firstkernel = firstkernel + recnkernel

      fmtstring = '("  Receiver ", A, ", lat: ", F8.3, ", lon: ", F8.3)'
      write(lu_out, fmtstring) trim(recname), reclatd, reclond
      if (master) then
          do ikernel = 1, recnkernel
             read(lu_receiver, *) !kernelname, filtername, timewindow
          end do
      end if

      call this%receiver(irec)%rotate_receiver( this%source%trans_rot_mat )
   end do
   if (master) then 
       close(lu_receiver)
   end if
   write(lu_out,*) ' In total ', this%nkernel, ' Kernels to calculate'
   call flush(lu_out)

   this%nkernel = lastkernel

end subroutine read_receiver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine read_kernel(this, sem_data, filter)
   use readfields, only            : semdata_type
   use filtering,  only            : filter_type
   class(parameter_type)          :: this
   type(semdata_type), intent(in), optional :: sem_data
   type(filter_type), target, intent(in), optional  :: filter(:)
   integer                        :: lu_kernel
   integer                        :: irec, nfilter
   integer                        :: ikernel, ifilter
   integer                        :: lu_receiver
   real(kind=dp)                  :: timewindow(2), junk
   character(len=1)               :: component
   character(len=4)               :: misfit_type
   character(len=32)              :: recname, kernelname, filtername, kernel_shortname
   character(len=80)              :: fmtstring

   if (master) then
       write(6,*)'  reading kernels from file ', trim(this%receiver_file)
       open(newunit=lu_receiver, file=trim(this%receiver_file), status='old')
       read(lu_receiver,*) 
       read(lu_receiver,*)
   end if

   allocate(this%kernel(this%nkernel))

   if (.not. master) nfilter     = size(filter)

   do irec = 1, this%nrec
      if (master) then
          read(lu_receiver, *) ! recname, reclatd, reclond, recnkernel

          fmtstring = '("  Receiver ", A, ", has ", I3, " kernel, first and last:", 2(I5))'
          write(lu_out,fmtstring) trim(this%receiver(irec)%name),  this%receiver(irec)%nkernel,&
                           this%receiver(irec)%firstkernel, this%receiver(irec)%lastkernel
      end if

      do ikernel = this%receiver(irec)%firstkernel, this%receiver(irec)%lastkernel

         if (master) read(lu_receiver, *) kernel_shortname, filtername, misfit_type, timewindow
         call pbroadcast_char(kernel_shortname, 0)
         call pbroadcast_char(filtername, 0)
         call pbroadcast_char(misfit_type, 0)
         call pbroadcast_dble_arr(timewindow, 0)

         kernelname = trim(this%receiver(irec)%name)//'_'//trim(kernel_shortname)

         if (.not.master) then
             do ifilter = 1, nfilter
                if (trim(filtername).eq.trim(filter(ifilter)%name)) exit
             end do
             if (ifilter == nfilter + 1) then
                print *, 'Could not find filter ', trim(filtername), ', which was requested'
                print *, 'by kernel ', trim(kernelname)
                print *, 'Available filters: ', [(filter(ifilter)%name, ifilter = 1, nfilter)]
                stop
             end if

             call this%kernel(ikernel)%init(name            = kernelname                , &
                                            time_window     = timewindow                , &
                                            filter          = filter(ifilter)           , &
                                            misfit_type     = misfit_type               , &  
                                            model_parameter = this%model_parameter      , &
                                            seis            = sem_data%veloseis(:,irec) , &
                                            dt              = sem_data%dt               , &
                                            timeshift_fwd   = sem_data%timeshift_fwd    )

         else
             this%kernel(ikernel)%name = kernelname
         end if


      end do
   end do
   if (master) close(lu_receiver)

   this%nkernel = this%receiver(this%nrec)%lastkernel

   call flush(6)
end subroutine read_kernel
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine read_filter(this, nomega, df)
   use filtering, only        : filter_type
   class(parameter_type)     :: this
   integer, intent(in), optional        :: nomega
   real(kind=dp), intent(in), optional  :: df
   integer                              :: lu_filter
   integer                              :: ifilter, nfilter
   character(len=32)                    :: filtername, filtertype
   character(len=80)                    :: fmtstring
   real(kind=dp)                        :: freqs(4)

   write(lu_out,*)'  reading filters from file ', trim(this%filter_file)
   if (master) then
       open(newunit=lu_filter, file=trim(this%filter_file), status='old')
       read(lu_filter, *) nfilter
   end if

   call pbroadcast_int(nfilter, 0)

   allocate(this%filter(nfilter))
  
   fmtstring = '("  Creating ", I5, " filters")'
   write(lu_out,fmtstring) nfilter

   fmtstring = '("  Creating filter ", A, " of type ", A/, "   freqs: ", 4(F8.3))'
   do ifilter = 1, nfilter
      if (master) read(lu_filter, *) filtername, filtertype, freqs
      call pbroadcast_char(filtername, 0)
      call pbroadcast_char(filtertype, 0)
      call pbroadcast_dble_arr(freqs, 0)

      write(lu_out,fmtstring) trim(filtername), trim(filtertype), freqs
      
      if(.not.master) then
          call this%filter(ifilter)%create(filtername, df, nomega, filtertype, freqs)
      end if

   end do

   if (master) close(lu_filter)
   call flush(lu_out)
end subroutine
!------------------------------------------------------------------------------

end module
