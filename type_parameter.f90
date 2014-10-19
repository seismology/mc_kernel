!=========================================================================================
module type_parameter
    use global_parameters,  only : sp, dp, pi, deg2rad, verbose, lu_out, master, testing
    use source_class,       only : src_param_type
    use kernel,             only : kernelspec_type
    use receiver_class,     only : rec_param_type
    use filtering,          only : filter_type
    use commpi,             only : pbroadcast_char, pbroadcast_int, pbroadcast_dble, &
                                   pbroadcast_dble_arr, pbroadcast_log, pbarrier, &
                                   pbroadcast_char_arr, pabort 
    use simple_routines,    only : lowtrim
    implicit none    

    type parameter_type
        type(src_param_type)                 :: source
        type(rec_param_type),   allocatable  :: receiver(:)
        type(kernelspec_type),  allocatable  :: kernel(:)
        type(filter_type),      allocatable  :: filter(:)
        integer                              :: nrec

        real(kind=dp)                        :: allowed_error
        real(kind=dp)                        :: allowed_relative_error = 1d-10

        character(len=512)                   :: fwd_dir
        character(len=512)                   :: bwd_dir
        character(len=512)                   :: source_file
        character(len=512)                   :: receiver_file
        character(len=512)                   :: filter_file
        character(len=512)                   :: mesh_file
        character(len=512)                   :: mesh_file_type
        character(len=512)                   :: mesh_file_vert
        character(len=512)                   :: mesh_file_face
        character(len=512)                   :: output_file = 'kerner'
        character(len=1)                     :: component
        character(len=32)                    :: whattodo
        character(len=32)                    :: int_type
        character(len=32)                    :: dump_type
        character(len=32)                    :: fftw_plan = 'MEASURE'
        character(len=32)                    :: strain_type_fwd
        integer                              :: nsim_fwd, nsim_bwd
        integer                              :: nkernel
        integer                              :: nfilter
        integer                              :: nelems_per_task
        integer                              :: npoints_per_step
        integer                              :: max_iter
        integer                              :: strain_buffer_size
        integer                              :: displ_buffer_size
        logical                              :: parameters_read      = .false.
        logical                              :: receiver_read        = .false.
        logical                              :: source_read          = .false.
        logical                              :: filter_read          = .false.
        logical                              :: kernel_read          = .false.
        logical                              :: detailed_convergence = .false.
        logical                              :: deconv_stf           = .false.
        logical                              :: write_smgr           = .true.
        logical                              :: quasirandom          = .true.
        contains
           procedure, pass                   :: read_parameters
           procedure, pass                   :: read_receiver
           procedure, pass                   :: read_source
           procedure, pass                   :: read_kernel
           procedure, pass                   :: read_filter
    end type

contains

!-----------------------------------------------------------------------------------------
subroutine read_parameters(this, input_file_in)
   class(parameter_type)           :: this
   character(len=*), intent(in), optional :: input_file_in
   character(len=256)              :: input_file
   integer                         :: lu_inparam_basic, ioerr, narg
   character(len=256)              :: line
   character(len=256)              :: keyword, keyvalue

   call pbarrier

   if (master) then

     if (present(input_file_in)) then
         input_file = input_file_in
     else
         narg = command_argument_count()
         if (narg<1) then
             print *, 'ERROR: Argument "parameter file" is missing'
             call pabort
         end if
         call get_command_argument(1, input_file)
     end if

     write(lu_out,'(" Reading ", A, "...")') trim(input_file)
     open(newunit=lu_inparam_basic, file=trim(input_file), status='old', &
          action='read', iostat=ioerr)
     if (ioerr /= 0) then
        print *, 'ERROR: Check input file ''', trim(input_file), '''! Is it still there?' 
        call pabort
     end if

     do
        read(lu_inparam_basic, fmt='(a256)', iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
     
        read(line,*) keyword, keyvalue 
      
        parameter_to_read : select case(trim(keyword))
        case('ALLOWED_ERROR')
           read(keyvalue, *) this%allowed_error

        case('ALLOWED_RELATIVE_ERROR')
           read(keyvalue, *) this%allowed_relative_error

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

        case('MESH_FILE_TYPE')
           this%mesh_file_type = keyvalue

        case('MESH_FILE_ABAQUS')
           this%mesh_file = keyvalue

        case('MESH_FILE_VERTICES')
           this%mesh_file_vert = keyvalue

        case('MESH_FILE_FACETS')
           this%mesh_file_face = keyvalue

        case('OUTPUT_FILE')
           this%output_file = keyvalue

        case('DUMP_TYPE')
           this%dump_type = keyvalue

        case('STRAIN_BUFFER_SIZE')
           read(keyvalue, *) this%strain_buffer_size

        case('DISPL_BUFFER_SIZE')
           read(keyvalue, *) this%displ_buffer_size

        case('MAXIMUM_ITERATIONS')
           read(keyvalue, *) this%max_iter

        case('ELEMENTS_PER_TASK')
           read(keyvalue, *) this%nelems_per_task

        case('POINTS_PER_MC_STEP')
           read(keyvalue, *) this%npoints_per_step

        case('WRITE_DETAILED_CONVERGENCE')
           read(keyvalue, *) this%detailed_convergence

        case('USE_QUASIRANDOM_NUMBERS')
           read(keyvalue, *) this%quasirandom

        case('DECONVOLVE_STF')
           read(keyvalue, *) this%deconv_stf

        case('FFTW_PLAN')
           read(keyvalue, *) this%fftw_plan

        case('WHAT_TO_DO')
           this%whattodo = keyvalue

        case('INT_TYPE')
           this%int_type = keyvalue

        case('WRITE_SEISMOGRAMS')
           read(keyvalue, *) this%write_smgr

        end select parameter_to_read
        if (verbose>0) write(lu_out, "('  ', A32,' ',A)")  keyword, trim(keyvalue)
     end do
     close(lu_inparam_basic)

  else
     write(lu_out,*) 'Waiting for master to distribute inparam_basic...'

  end if

  
  ! broadcast values to other processors
  call pbroadcast_char(this%fwd_dir, 0) 
  call pbroadcast_char(this%bwd_dir, 0) 
  call pbroadcast_dble(this%allowed_error, 0)
  call pbroadcast_dble(this%allowed_relative_error, 0)
  call pbroadcast_char(this%source_file, 0) 
  call pbroadcast_char(this%receiver_file, 0)
  call pbroadcast_char(this%filter_file, 0)
  call pbroadcast_char(this%mesh_file_type, 0)
  call pbroadcast_char(this%output_file, 0)
  call pbroadcast_char(this%dump_type, 0)
  call pbroadcast_int(this%strain_buffer_size, 0)
  call pbroadcast_int(this%displ_buffer_size, 0)
  call pbroadcast_int(this%max_iter, 0)
  call pbroadcast_int(this%nelems_per_task, 0)
  call pbroadcast_int(this%npoints_per_step, 0)
  call pbroadcast_log(this%detailed_convergence, 0)
  call pbroadcast_log(this%deconv_stf, 0)
  call pbroadcast_log(this%write_smgr, 0)
  call pbroadcast_log(this%quasirandom, 0)
  call pbroadcast_char(this%fftw_plan, 0)
  call pbroadcast_char(this%whattodo, 0)
  call pbroadcast_char(this%int_type, 0)

  if (lowtrim(this%mesh_file_type)=='abaqus') then
     call pbroadcast_char(this%mesh_file, 0)
  else
     call pbroadcast_char(this%mesh_file_vert, 0)
     call pbroadcast_char(this%mesh_file_face, 0)
  end if

  write(lu_out,*)
  call flush(lu_out)

  this%parameters_read = .true.

end subroutine read_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_source(this)
   class(parameter_type)          :: this
   real(kind=dp)                  :: Mij_dyncm(6), latd, lond, depth
   character(len=16)              :: junk
   character(len=16)              :: event_name
   integer, parameter             :: lu_source=1000

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_source!'
       stop
   end if

   call pbarrier

   if (master) then
       write(lu_out,*) 'Reading source from file ', trim(this%source_file)
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

   else
      write(lu_out,*) 'Waiting for master to distribute source...'

   end if

   call pbroadcast_dble(latd, 0)
   call pbroadcast_dble(lond, 0)
   call pbroadcast_dble(depth, 0)
   call pbroadcast_dble_arr(Mij_dyncm, 0)
   

   call this%source%init(lat=latd, lon=lond, mij=Mij_dyncm*1.E-7, depth=depth)

   write(lu_out,*)
   call flush(lu_out)

   this%source_read = .true.

end subroutine read_source
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_receiver(this)
   use simple_routines, only      : to_lower 
   use kernel,          only      : tabulate_kernels
   class(parameter_type)         :: this
   integer, parameter            :: lu_receiver = 1001
   integer                       :: irec, firstkernel, lastkernel
   integer                       :: ikernel, recnkernel
   real(kind=dp)                 :: reclatd, reclond
   character(len=16)             :: recname, trash(5), model_parameter
   character(len=80)             :: fmtstring

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_receiver!'
       call pabort
   end if

   if (.not.this%source_read) then
       write(6,*) 'Source parameters not read. Call read_source before read_receiver!'
       call pabort
   end if

   this%strain_type_fwd = 'straintensor_trace'

   call pbarrier

   if (master) then
       write(lu_out,*)'Reading receivers from file ', trim(this%receiver_file)
       open(unit=lu_receiver, file=trim(this%receiver_file), status='old')
       read(lu_receiver,*) this%nrec
       read(lu_receiver,*) this%component

   else
      write(lu_out,*) 'Waiting for master to distribute receivers...'

   end if

   call pbroadcast_int(this%nrec, 0)
   call pbroadcast_char(this%component, 0)


   fmtstring = '("  Using ", I5, " receivers")'
   write(lu_out, fmtstring) this%nrec

   fmtstring = '("  Kernel on component ", A)'
   write(lu_out, fmtstring) this%component

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
             ! Just get the model parameter
             read(lu_receiver, *) trash(1:5), model_parameter
             
             if ( (model_parameter.ne.'lambda') .and.      &
                  (model_parameter.ne.'vp')          ) then
               this%strain_type_fwd = 'straintensor_full'
             end if
                  
          end do
       end if


       call this%receiver(irec)%rotate_receiver( this%source%trans_rot_mat )
   end do

   call pbroadcast_char(this%strain_type_fwd, 0)

   if (master) then 
       close(lu_receiver)
   end if

   this%nkernel = lastkernel

   write(lu_out,*) ' In total ', this%nkernel, ' Kernels to calculate'

   this%receiver_read = .true.

   write(lu_out,*)
   call flush(lu_out)

end subroutine read_receiver
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_kernel(this, sem_data, filter)
   use readfields, only            : semdata_type
   use filtering,  only            : filter_type
   class(parameter_type)          :: this
   type(semdata_type), intent(in), optional :: sem_data
   type(filter_type), target, intent(in), optional  :: filter(:)
   integer                        :: irec, nfilter
   integer                        :: ikernel, ifilter
   integer                        :: lu_receiver
   real(kind=dp)                  :: timewindow(2)
   character(len=4)               :: misfit_type
   character(len=16)              :: model_parameter
   character(len=32)              :: kernelname, filtername, kernel_shortname
   character(len=80)              :: fmtstring

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_kernel!'
       call pabort
   end if

   if (.not.this%source_read) then
       write(6,*) 'Source parameters not read. Call read_source before read_kernel!'
       call pabort
   end if

   if (.not.this%receiver_read) then
       write(6,*) 'Receiver parameters not read. Call read_receiver before read_kernel!'
       call pabort
   end if

   if (.not.this%filter_read) then
       write(6,*) 'Filter parameters not read. Call read_filter before read_kernel!'
       call pabort
   end if

   call pbarrier

   if (master) then
       write(6,*)'Reading kernels from file ', trim(this%receiver_file)
       open(newunit=lu_receiver, file=trim(this%receiver_file), status='old')
       read(lu_receiver,*) 
       read(lu_receiver,*) this%component

   else
      write(lu_out,*) 'Waiting for master to distribute kernels...'

   end if

   allocate(this%kernel(this%nkernel))

   if ((.not.master).or.testing) then
       nfilter = size(filter)

       if (this%deconv_stf) then
          do ifilter = 1, nfilter
              call filter(ifilter)%add_stfs(sem_data%stf_fwd, sem_data%stf_bwd)
          end do
       end if
   else
       nfilter = 0
   end if

   do irec = 1, this%nrec
      if (master) then
          read(lu_receiver, *)
      end if
      
      ! Start with the assumption that this receiver needs no base kernels at all
      this%receiver(irec)%needs_basekernel = .false.

      fmtstring = '("  Receiver ", A, ", has ", I3, " kernel, first and last:", 2(I5))'
      write(lu_out,fmtstring) trim(this%receiver(irec)%name),  this%receiver(irec)%nkernel,   &
                              this%receiver(irec)%firstkernel, this%receiver(irec)%lastkernel

      do ikernel = this%receiver(irec)%firstkernel, this%receiver(irec)%lastkernel

         if (master) read(lu_receiver, *) kernel_shortname, filtername, & 
                                          misfit_type, timewindow, model_parameter

         call pbroadcast_char(model_parameter, 0)
         call pbroadcast_char(kernel_shortname, 0)
         call pbroadcast_char(filtername, 0)
         call pbroadcast_char(misfit_type, 0)
         call pbroadcast_dble_arr(timewindow, 0)
         call pbroadcast_char(model_parameter, 0)

         kernelname = trim(this%receiver(irec)%name)//'_'//trim(kernel_shortname)

         if ((.not.master).or.(testing)) then
             if (.not.(present(sem_data).and.present(filter))) then
                write(*,*) 'ERROR: Only master may call read_kernel without sem_data and filter'
                call pabort()
             end if
             do ifilter = 1, nfilter
                if (trim(filtername).eq.trim(filter(ifilter)%name)) exit
             end do
             if (ifilter == nfilter + 1) then
                print *, 'Could not find filter ', trim(filtername), ', which was requested'
                print *, 'by kernel ', trim(kernelname)
                print *, 'Available filters: ', [(filter(ifilter)%name, ifilter = 1, nfilter)]
                call pabort
             end if

             call this%kernel(ikernel)%init(name            = kernelname                , &
                                            time_window     = timewindow                , &
                                            filter          = filter(ifilter)           , &
                                            misfit_type     = misfit_type               , &
                                            model_parameter = model_parameter           , &
                                            seis            = sem_data%veloseis(:,irec) , &
                                            dt              = sem_data%dt               , &
                                            timeshift_fwd   = sem_data%timeshift_fwd    , &
                                            write_smgr      = this%write_smgr)

             ! Update which base kernels are needed for this receiver
             this%receiver(irec)%needs_basekernel =   &
                                  this%receiver(irec)%needs_basekernel &
                             .or. this%kernel(ikernel)%needs_basekernel
         else
             this%kernel(ikernel)%name = kernelname
         end if


      end do ! ikernel
 
      ! TODO: This might allow a finer control over what has to be read. The final idea:
      ! If any kernel of this receiver needs more than a lambda base kernel,
      ! the whole strain tensor has to be read in. If only lambda, the trace is enough.
      ! However, this would require some involved changes in readfields.f90. 
      if (any(this%receiver(irec)%needs_basekernel(2:6))) then
         this%receiver(irec)%strain_type = 'straintensor_full'
         !this%strain_type_fwd = 'straintensor_full'
      else
         this%receiver(irec)%strain_type = 'straintensor_trace'
      end if

   end do ! irec
   if (master) close(lu_receiver)

   this%nkernel = this%receiver(this%nrec)%lastkernel

   write(lu_out,*)
   call flush(lu_out)

   this%kernel_read = .true.

end subroutine read_kernel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_filter(this, nomega, df)
   use filtering, only        : filter_type
   class(parameter_type)     :: this
   integer, intent(in), optional        :: nomega
   real(kind=dp), intent(in), optional  :: df
   integer                              :: lu_filter
   integer                              :: ifilter
   character(len=32)                    :: filtername, filtertype
   character(len=80)                    :: fmtstring
   real(kind=dp)                        :: freqs(4)

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_filter!'
       call pabort
   end if

   write(lu_out,*)'Reading filters from file ', trim(this%filter_file)
   if (master) then
       open(newunit=lu_filter, file=trim(this%filter_file), status='old')
       read(lu_filter, *) this%nfilter

   else
      write(lu_out,*) 'Waiting for master to distribute filters...'

   end if

   call pbroadcast_int(this%nfilter, 0)

   allocate(this%filter(this%nfilter))
  
   fmtstring = '("  Creating ", I5, " filters")'
   write(lu_out,fmtstring) this%nfilter

   fmtstring = '("  Creating filter ", A, " of type ", A/, "   freqs: ", 4(F8.3))'
   do ifilter = 1, this%nfilter
      if (master) read(lu_filter, *) filtername, filtertype, freqs
      call pbroadcast_char(filtername, 0)
      call pbroadcast_char(filtertype, 0)
      call pbroadcast_dble_arr(freqs, 0)

      write(lu_out,fmtstring) trim(filtername), trim(filtertype), freqs
      
      if((.not.master).or.(testing)) then
          if (.not.(present(df).and.present(nomega))) then
             write(*,*) 'ERROR: Only master may call read_filter without nomega and df'
             call pabort()
          end if
          call this%filter(ifilter)%create(filtername, df, nomega, filtertype, freqs)
      end if

   end do

   if (master) close(lu_filter)

   write(lu_out,*)
   call flush(lu_out)

   this%filter_read = .true.
end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
