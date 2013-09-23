!========================
module parameters
!========================
 
    use global_parameters
    use commpi
    
    implicit none
  
  
    public :: read_kernel_input, define_sem_kernel_meshes
    public :: read_parameters
    !public :: save_simulation_info
!    public :: ioparamtype, kernelparamtype, ncparamtype, paramtype
    !private

  contains 


!-------------------------------------------------------------------------------------------------
subroutine read_parameters(params)
    use parameter_types
    type(paramtype), intent(inout) :: params

    call check_fwd_and_bwd(params%io%fwd(1), params%io%bwd(1))

    call read_source(params%src_fwd, params%src_bwd)

    call read_params_from_solver(params%src_fwd, params%src_bwd, params%io%fwd(1), params%io%bwd(1))

    call read_kernel_input(params)

end subroutine read_parameters


!-------------------------------------------------------------------------------------------------
subroutine read_kernel_input(params)
!<  Read input_kernel.dat and time windows 
 
    use parameter_types
    type(paramtype), intent(inout) :: params
    integer                        :: i, ierr, idim
    character(len=20)              :: stf_type
 
    
    open(unit=12,file='input_kernel.dat',status='old',iostat=ierr)
    read(12,*)  !calc_or_load_kernels
    read(12,'(a8)') !kernmesh_type
    read(12,*) params%k%thr
    read(12,*) params%k%phr
    read(12,*) params%k%phi0
    read(12,*) params%k%char_wfkern
    read(12,*) params%k%twf1
    read(12,*) params%k%twf2
    read(12,*) !td_fd
    read(12,*) params%filter%what
    read(12,*) params%filter%type
    read(12,*) params%filter%period_low, params%filter%period_hi
    params%filter%fmin = 1. / params%filter%period_hi
    params%filter%fmax = 1. / params%filter%period_low
    params%filter%omegamin = 2.*pi / params%filter%period_hi
    params%filter%omegamax = 2.*pi / params%filter%period_low
    !   endif
    ! Debunked according to TNM
    !if (filter_type=='gau' .or. filter_type=='low') filter_period=filter_period_low
    !if (filter_type=='hig') filter_period=filter_period_hi
 
    
    read(12,*) params%k%compute_src_kernels
    read(12,*) params%k%do_rho
    read(12,*) params%k%do_lam
    read(12,*) params%k%do_mu
    read(12,*) params%k%do_vp
    read(12,*) params%k%do_vs
    read(12,*) params%k%do_imped
    read(12,*) params%k%save_snaps
    read(12,*) params%k%times_medium
    read(12,*) params%k%norm_kernels
    read(12,*) !dump_vtk
    read(12,*) !dump_avs
    read(12,*) !dump_bin
    read(12,*) !dump_ascii
    read(12,*) !params%src%nsim_req !Not necessary, derived from dir
    read(12,*) !(params%src%Mij(i),i=1,6)
    close(12)

    !if (params%src%nsim.ne.params%io%nsim) then
    !    if (mynum.eq.0) write(6,*) 'Forward was not run with multiple simulations, but', &
    !    ' this is requested in input_kernel.dat'
    !    stop
    !end if
    !if (mynum==0) write(6,*)'Number of simulations:', params%src%nsim
 
    call read_time_windows(params%k%twf1, params%k%twf2, params%k%ntw, &
                           params%k%tf1, params%k%tf2, params%k%phase_name)
 
    if (mynum==0) then 
       write(6,*)
       write(6,*)' receiver theta,phi:', params%k%thr, params%k%phr
       write(6,*)' waveform time window t1,t2:',params%k%twf1, params%k%twf2
    endif
    call flush(6)
 
    params%k%thr=params%k%thr*pi/180.d0
    params%k%phr=params%k%phr*pi/180.d0
 
 
    !params%k%dim_fwd = 3
    !if (params%src_fwd%type1=='monopole') params%k%dim_fwd = 2
 
    !if (params%src_fwd%nsim.eq.4) then 
    !   params%src_fwd%srctype = ['mzz','mxx','mxz','mxy']
    !elseif (params%src_fwd%nsim.eq.2) then
    !   params%src_fwd%srctype = ['vertforce','xforce']
    !else
    !if (params%src_fwd%nsim.eq.1) then
    !   params%src%srctype = trim(params%src%type2fwd)
    !   params%src%exctype = trim(params%src%type1fwd)
    ! !  params%src%Mij = params%src%magnitude
    !endif
 
    !!@TODO: IS THIS RIGHT?
    !params%src%Mij = params%src%Mij / params%src%magnitude

 end subroutine read_kernel_input
!-------------------------------------------------

subroutine read_source(fwd, bwd)
    use parameter_types
    type(sourceparamtype), intent(inout) :: fwd, bwd
    integer                              :: ifileparam, nval, ioerr
    character(len=256)                   :: line
    character(len=256), dimension(10)    :: vals

    open(unit=ifileparam, file='input_source.dat', status='old', action='read',  iostat=ioerr)
    if (ioerr.ne.0) stop 'Check input file ''input_source.dat''! Is it still there?' 
    do
        read(ifileparam,fmt='(a256)',iostat=ioerr) line
        if (ioerr.lt.0) exit
        if (len(trim(line)).lt.1.or.line(1:1).eq.'#') cycle
        call read_line(line,vals,nval)
    
        parameter_to_read : select case(trim(vals(1)))
        case('SOURCE_TYPE')
            select case(trim(vals(2))) 
            case('moment')
                if (fwd%nsim.ne.4) then
                    write(6,*) 'Forward simulation needs to be run with -s moment'
                    stop
                else
                    allocate(fwd%srctype(4))
                    allocate(fwd%exctype(4))
                    fwd%srctype = ['mzz','mxx','mxz','mxy']
                    fwd%exctype = ['monopole','monopole','dipole  ','quadpole']
                end if

            case('forces')
                if (fwd%nsim.ne.2) then 
                    write(6,*) 'Forward simulation needs to be run with -s forces'
                    stop
                else
                    allocate(fwd%srctype(2))
                    allocate(fwd%exctype(2))
                    fwd%srctype = ['vertforce', 'xforce   ']
                    fwd%exctype = ['monopole', 'dipole  ']
                end if

            case('single')
                if (fwd%nsim.ne.1) then
                    write(6,*) 'Forward simulation needs to be run with -s forces'
                    stop
                else
                    allocate(fwd%srctype(1))
                    allocate(fwd%exctype(1))
                    ! Is determined by Solver then
                    !read(ifileparam,*) fwd%exc_type 
                    !read(ifileparam,*) fwd%srctype 
                end if
            case default 
                write(6,*) 'ERROR: unrecognized forward source type!'
                stop
            end select

        case('CMT_SRC')
            if (fwd%nsim.eq.4) then
                write(vals(2:7) ,*) fwd%mij
            end if

        case('FORCES_SRC')
            if (fwd%nsim.eq.2) then
                write(vals(2:4),*) fwd%forces
            end if

        case('SCALAR_MOMENT_SRC')
            if (fwd%nsim.eq.1) then
                write(vals(2), *) fwd%magnitude
            end if
    
        case('RECEIVER_TYPE')
            select case(trim(vals(2))) 
            case('moment')
                if (fwd%nsim.ne.4) then
                    write(6,*) 'Forward simulation needs to be run with -s moment'
                    stop
                else
                    allocate(bwd%srctype(4))
                    allocate(bwd%exctype(4))
                    bwd%srctype = ['mzz','mxx','mxz','mxy']
                    bwd%exctype = ['monopole','monopole','dipole  ','quadpole']
                end if

            case('forces')
                if (bwd%nsim.ne.2)then 
                    write(6,*) 'Forward simulation needs to be run with -s forces'
                    stop
                else
                    allocate(bwd%srctype(2))
                    allocate(bwd%exctype(2))
                    bwd%srctype = ['vertforce', 'xforce   ']
                    bwd%exctype = ['monopole', 'dipole  ']
                end if

            case('single')
                if (bwd%nsim.ne.1)then 
                    write(6,*) 'Forward simulation needs to be run with -s forces'
                    stop
                else
                    allocate(bwd%srctype(1))
                    allocate(bwd%exctype(1))
                    ! Is determined by Solver then
                    !read(ifileparam,*) bwd%exc_type 
                    !read(ifileparam,*) bwd%srctype 
                end if
            case default 
                write(6,*) 'ERROR: unrecognized forward source type!'
                stop
            end select
        
        case('CMT_REC')
            if (bwd%nsim.eq.4) then
                write(vals(2:7) ,*) bwd%mij
            end if

        case('FORCES_REC')
            if (bwd%nsim.eq.2) then
                write(vals(2:4) ,*) bwd%forces
            end if

        case('SCALAR_MOMENT_REC')
            if (bwd%nsim.eq.1) then
                write(vals(2) ,*) bwd%magnitude
            end if
    
        end select parameter_to_read
    end do

    end subroutine

!-------------------------------------------------
    subroutine read_line(line,vals,nval)
    !< Reads string and gives back array of its separate words. Taken from K.Sigloch
    implicit none

    integer, parameter :: NMAX = 30

    character(len=256), intent(in) :: line
    integer, intent(out)  :: nval
    character(len=*), dimension(NMAX), intent(out) :: vals
    
    integer :: i,wcnt,lcnt
    character :: symbol
    character(len=256) :: word
    logical :: in_word
    
    ! init    
    wcnt=0 
    in_word=.FALSE.
    word=''

    do i=1,256
      symbol=line(i:i)
      if (symbol.ne.' ') then
         if(.not.in_word) then ! beginning of new word
            in_word=.TRUE.
            wcnt=wcnt+1
            if(wcnt.gt.NMAX) go to 100
            lcnt=1
            word(lcnt:lcnt)=symbol
         else                  ! in the middle of a word
            lcnt=lcnt+1      
            word(lcnt:lcnt)=symbol
         endif
      else  ! blank found
         if(in_word) then      ! previous symbol was end of word
           in_word=.FALSE.     
           vals(wcnt)=word
           word=''  ! wipe blank before filling anew  
         endif 
      endif
    end do

    nval=wcnt
    return

100 print *, 'ERROR in read_line: cannot read more than ',NMAX, ' words'
    stop

    end subroutine


!-------------------------------------------------
subroutine read_params_from_solver(src_fwd, src_bwd, fwd, bwd)
    use parameter_types
    use nc_routines, only: nc_read_att_char, nc_read_att_real, nc_read_att_int    
    type(sourceparamtype), intent(inout) :: src_fwd, src_bwd
    type(ncparamtype), intent(in)        :: fwd, bwd
    
    call nc_read_att_real(src_fwd%magnitude, &
                          'scalar source magnitude', fwd)
    call nc_read_att_real(src_fwd%shift_fact, &
                          'source shift factor in sec', fwd)
    call nc_read_att_real(src_fwd%shift_fact_deltat, &
                          'source shift factor for deltat', fwd)
    call nc_read_att_real(src_fwd%shift_fact_deltat_coarse, &
                          'source shift factor for deltat_coarse', fwd)
    call nc_read_att_real(src_fwd%shift_fact_seis_dt, &
                          'source shift factor for seis_dt', fwd)
 
    call nc_read_att_real(src_bwd%magnitude, &
                          'scalar source magnitude', bwd)
    call nc_read_att_real(src_bwd%shift_fact, &
                          'source shift factor in sec', bwd)
    call nc_read_att_real(src_bwd%shift_fact_deltat, &
                          'source shift factor for deltat', bwd)
    call nc_read_att_real(src_bwd%shift_fact_deltat_coarse, &
                          'source shift factor for deltat_coarse', bwd)
    call nc_read_att_real(src_bwd%shift_fact_seis_dt, &
                          'source shift factor for seis_dt', bwd)
 
 
    call nc_read_att_char(src_fwd%type1, 'excitation type', fwd)
    call nc_read_att_char(src_fwd%type2, 'source type', fwd)
    call nc_read_att_real(src_fwd%depth, 'source depth in km', fwd)
    call nc_read_att_char(src_fwd%stf_type, 'source time function', fwd)


    call nc_read_att_char(src_bwd%type1, 'excitation type', bwd)
    call nc_read_att_char(src_bwd%type2, 'source type', bwd)
    call nc_read_att_real(src_bwd%depth, 'source depth in km', bwd)
    call nc_read_att_char(src_bwd%stf_type, 'source time function', bwd)
    
    if (src_bwd%nsim.eq.1) then
       src_bwd%srctype = trim(src_bwd%type2)
       src_bwd%exctype = trim(src_bwd%type1)
    endif

    if (src_fwd%nsim.eq.1) then
       src_fwd%srctype = trim(src_fwd%type2)
       src_fwd%exctype = trim(src_fwd%type1)
    endif
end subroutine read_params_from_solver
!-------------------------------------------------


subroutine check_fwd_and_bwd(fwd, bwd)
!< Checks whether forward and backward wavefield were calculated with a reasonable 
!! combination of parameters (i.e. same background model, same point number etc.) 
!! Does not actually set any parameters
    use parameter_types
    use nc_routines, only: nc_read_att_char, nc_read_att_real, nc_read_att_int    
    type(ncparamtype), intent(in)  :: fwd, bwd

    logical                        :: testfail=.false.
    character(len=200)             :: bgmodel_fwd, bgmodel_bwd
    real                           :: period_fwd, period_bwd
    integer                        :: ibeg_fwd, ibeg_bwd, iend_fwd, iend_bwd
    real                           :: dt_fwd, dt_bwd
    integer                        :: ndumps_fwd, ndumps_bwd

    if (mynum.eq.0) write(6,*) 'Doing some tests on compatibility', &
        ' of forward and backward field'

    call nc_read_att_char(bgmodel_fwd, 'background model', fwd)
    call nc_read_att_char(bgmodel_bwd, 'background model', bwd)
    if (bgmodel_fwd.eq.bgmodel_bwd) then
        if (mynum.eq.0) write(6,*) ' PASSED: Same background models'
    else
        if (mynum.eq.0) write(6,*) ' FAILED: Background models differ'
        testfail = .true.
    end if


    call nc_read_att_real(period_fwd, 'dominant source period', fwd)
    call nc_read_att_real(period_bwd, 'dominant source period', bwd)
    
    if (period_fwd.eq.period_bwd) then
        if (mynum.eq.0) write(6,*) ' PASSED: Same source period'
    else
        if (mynum.eq.0) write(6,*) ' FAILED: Source periods differ'
        testfail = .true.
    end if

    call nc_read_att_int(ibeg_fwd, 'ibeg', fwd)
    call nc_read_att_int(iend_fwd, 'iend', fwd)
    call nc_read_att_int(ibeg_bwd, 'ibeg', bwd)
    call nc_read_att_int(iend_bwd, 'iend', bwd)
    if ((iend_fwd.eq.iend_bwd).and.(ibeg_fwd.eq.ibeg_bwd)) then
        if (mynum.eq.0) write(6,*) ' PASSED: Same ibeg and iend'
    else
        if (mynum.eq.0) write(6,*) ' FAILED: ibeg and iend differ'
        testfail = .true.
    end if

    call nc_read_att_real(dt_fwd, 'time step in sec', fwd)
    call nc_read_att_real(dt_bwd, 'time step in sec', bwd)

    if (dt_fwd.eq.dt_bwd) then
        if (mynum.eq.0) write(6,*) ' PASSED: Same time step'
    else
        if (mynum.eq.0) write(6,*) ' FAILED: time steps differ'
        testfail = .true.
    end if

    call nc_read_att_int(ndumps_fwd, 'number of strain dumps', fwd)
    call nc_read_att_int(ndumps_bwd, 'number of strain dumps', bwd)

    if (ndumps_fwd.eq.ndumps_bwd) then
        if (mynum.eq.0) write(6,*) ' PASSED: Same number of strain dumps'
    else
        if (mynum.eq.0) write(6,*) ' FAILED: number of strain dumps differs'
        testfail = .true.
    end if


    if (testfail) then
        if (mynum.eq.0) write(6,*) ' ERROR: One or more tests failed. Please check '
        stop
    end if
end subroutine

!-------------------------------------------------
subroutine check_input_parameters(k, src_fwd, src_bwd)
!< Put some tests into this separate subroutine. Am not actually aware, whether they make
!! sense. Only TNM knows...

    use parameter_types
    type(kernelparamtype), intent(inout) :: k
    type(sourceparamtype), intent(in)    :: src_fwd, src_bwd


! define component of the receiver-seismogram

!! TNM Nov 2010......... kinda stupid: xforce and yforce, being the 
! same simulations... should be handled jointly, NOT like below.... plus rotations are necessary!
    if (k%do_vp .or. k%do_vs .or. k%do_imped  & 
        .and. .not. k%times_medium) then 
       if (mynum==0) write(6,*)'Need medium parameters for wavespeed/impedance!'
       k%times_medium=.true.
    endif
    
    if (k%do_vp .and. .not. k%do_lam) then
       if (mynum==0) write(6,*)'Need bulk kernel for compressional vel. kernel!'
       k%do_lam=.true.
    endif
 
    if (k%do_vs .and. .not. ( k%do_lam .or. k%do_mu ) ) then 
       if (mynum==0) write(6,*)'need bulk & shear kernels for shear vel. kernel!'
       k%do_lam=.true.
       k%do_mu=.true.
    endif
    
    if (k%do_imped .and. .not. &
       (k%do_lam .or. k%do_mu .or. k%do_rho) ) then 
       if (mynum==0) write(6,*)'need all primitive kernels for impedance kernel'
       k%do_lam=.true.
       k%do_mu=.true.
       k%do_rho=.true.
    endif

    !> @TODO: The following block serves a rather opaque purpose...
    if (src_bwd%type1=='monopole') then
        if (src_fwd%type1=='monopole') then 
            k%reccomp=2
        else
            k%reccomp=3
        endif
    elseif (src_bwd%type2=='xforce') then 
        k%reccomp=1
    elseif (src_bwd%type2=='yforce') then 
        if (src_fwd%type1=='monopole') then 
            write(6,*) 'This combination does not make sense:'
            write(6,*)'monopole fwd source will not have any signal on the transverse bwd field...'
            stop
        else
            k%reccomp=2
        endif
    else 
        k%reccomp=2
    endif
end subroutine
!-------------------------------------------------

!-------------------------------------------------
subroutine read_time_windows(twf1, twf2, ntw, t1, t2, phase_name)
    !> Beginning and end of wevefield kernel window
    real, intent(inout)            :: twf1, twf2 
    !> Number of time windows
    integer, intent(out)           :: ntw
    !> Beginning and end of all time windows
    real, intent(out), allocatable :: t1(:), t2(:)

    character(len=10), intent(out), allocatable :: phase_name(:)
    character(len=4) :: appmywin
    logical taup_exists
    integer itw


    inquire(file='input_timewindows_taup.dat',exist=taup_exists)
    if (taup_exists) then
       if (mynum==0) write(6,*) 'reading input_timewindows_taup.dat'
       open(unit=13,file='input_timewindows_taup.dat') 
    else
       if (mynum==0) write(6,*) 'reading input_timewindows.dat'
       open(unit=13,file='input_timewindows.dat')
    endif
    read(13,*) ntw ! number of time windows
    if (mynum==0) write(6,*)'Number of time windows:',ntw
    if (ntw==0) then
       if (mynum==0) write(6,*)'No time window? doing waveform kernels only.'
       ntw=1
       allocate(t1(ntw),t2(ntw))
       t1(ntw)=twf1 
       t2(ntw)=twf2
       !char_wfkern='allt'
    else
       allocate(t1(ntw),t2(ntw))
       allocate(phase_name(ntw))
       do itw=1,ntw
          if (taup_exists) then 
             read(13,*)t1(itw),t2(itw),phase_name(itw)
          else
             read(13,*)t1(itw),t2(itw)
             call define_io_appendix(appmywin,itw)
             phase_name(itw)=appmywin
          endif
          if (mynum==0) write(6,*)'time window #,from/to/phase',itw , &
              t1(itw),t2(itw),trim(phase_name(itw))

          if (t1(itw)<twf1) then 
             if (mynum==0) write(6,*)'time window starts before wf kernel!'
             twf1=t1(itw)
          elseif (t2(itw)>twf2) then 
             if (mynum==0) write(6,*)'time window extends beyond wf kernel!'
             twf2=t2(itw)
          endif
       enddo
    endif
    close(13)
end subroutine

!-------------------------------------------------
subroutine define_sem_kernel_meshes(io, fwd_nsim, bwd_nsim)
    use parameter_types
    type(ioparamtype), intent(inout)  :: io
    integer, intent(out)              :: fwd_nsim, bwd_nsim
    character(len=200)                :: filename
    logical                           :: moment, force, fwdexist, bwdexist
    integer                           :: ifilesem = 98
    character(len=512), dimension(10) :: vals
    character(len=512)                :: line
    integer                           :: nval, ioerr

    open(unit=ifilesem,file="input_sem.dat")
    do 
        read(ifilesem,fmt='(a512)',iostat=ioerr) line
        if (ioerr.lt.0) exit
        if ((len_trim(line).lt.1).or.(line(1:1).eq.'#')) cycle
        call read_line(line,vals,nval)
        select case(trim(vals(1)))
        case('FORWARD_DIR')
            io%dir_fwdmesh = trim(vals(2))
        case('BACKWARD_DIR')
            io%dir_bwdmesh = trim(vals(2))
        case('KERNEL_MESH')
            io%ext_mesh_name = trim(vals(2))
        case('DATA_DIR')
            io%data_dir = trim(vals(2))
        end select
    end do
    !read(98,*) io%dir_fwdmesh
    !read(98,*) io%dir_bwdmesh
    !read(98,*) io%ext_mesh_name
    !read(98,*) io%data_dir
    close(ifilesem)
    io%lffwd = index(io%dir_fwdmesh,' ')-1
    io%lfbwd = index(io%dir_bwdmesh,' ')-1
    io%lfext = index(io%ext_mesh_name,' ')-1
 
    if (mynum==0) then 
        write(6,*)'forward calc :',trim(io%dir_fwdmesh)
        write(6,*)'backward calc:',trim(io%dir_bwdmesh)
        write(6,*)'external mesh:',trim(io%ext_mesh_name)
        write(6,*)
    endif

    filename = trim(io%dir_fwdmesh)
    inquire(file=filename, exist=fwdexist)
    if (.not.fwdexist) then 
        if (mynum.eq.0) write(6,*) 'ERROR: Forward directory ', trim(io%dir_fwdmesh), &
                                   ' does not exist.'
        stop
    end if

    filename = trim(io%dir_bwdmesh)
    inquire(file=filename, exist=bwdexist)
    if (.not.bwdexist) then 
        if (mynum.eq.0) write(6,*) 'ERROR: Backward directory ', trim(io%dir_bwdmesh), &
                                   ' does not exist.'
        stop
    end if

    filename = trim(io%dir_fwdmesh)//'MZZ'
    inquire(file=filename, exist=moment)
    filename = trim(io%dir_fwdmesh)//'PZ'
    inquire(file=filename, exist=force)
    
    if (moment) then
        if (mynum.eq.0) write(6,*) 'Forward wavefield was calculated with multiple ', &
                                   'simulations. Summing up according to CMT.'
        fwd_nsim = 4
    elseif (force) then
        fwd_nsim = 2
    else
        fwd_nsim = 1
    end if

    filename = trim(io%dir_bwdmesh)//'MZZ'
    inquire(file=filename, exist=moment)
    filename = trim(io%dir_bwdmesh)//'PZ'
    inquire(file=filename, exist=force)
    
    if (moment) then
        if (mynum.eq.0) write(6,*) 'Backward wavefield was calculated with multiple ', &
                                   'simulations. Summing up according to CMT.'
        bwd_nsim = 4
    elseif (force) then
        bwd_nsim = 2
    else
        bwd_nsim = 1
    end if
    !filename = trim(io%dir_bwdmesh)//'MZZ'
    !inquire(file=filename, exist=multisim)
    !if (multisim) then 
    !    if (mynum.eq.0) write(6,*) 'ERROR: Backward wavefield was calculated with multiple simulations. This makes no sense.'
    !    stop
    !end if

    io%nsim_fwd= fwd_nsim
    io%nsim_bwd= bwd_nsim
    allocate(io%fwd(fwd_nsim))
    allocate(io%bwd(bwd_nsim))

    if (fwd_nsim.eq.4) then 
        io%fwd(1)%meshdir = trim(io%dir_fwdmesh)//'MZZ/'
        io%fwd(2)%meshdir = trim(io%dir_fwdmesh)//'MXX_P_MYY/'
        io%fwd(3)%meshdir = trim(io%dir_fwdmesh)//'MXZ_MYZ/'
        io%fwd(4)%meshdir = trim(io%dir_fwdmesh)//'MXY_MXX_M_MYY/'
    elseif (fwd_nsim.eq.2) then
        io%fwd(1)%meshdir = trim(io%dir_fwdmesh)//'PZ/'
        io%fwd(2)%meshdir = trim(io%dir_fwdmesh)//'PX/'
    else
        io%fwd(1)%meshdir = trim(io%dir_fwdmesh)
    endif

    if (bwd_nsim.eq.4) then 
        io%bwd(1)%meshdir = trim(io%dir_bwdmesh)//'MZZ/'
        io%bwd(2)%meshdir = trim(io%dir_bwdmesh)//'MXX_P_MYY/'
        io%bwd(3)%meshdir = trim(io%dir_bwdmesh)//'MXZ_MYZ/'
        io%bwd(4)%meshdir = trim(io%dir_bwdmesh)//'MXY_MXX_M_MYY/'
    elseif (bwd_nsim.eq.2) then
        io%bwd(1)%meshdir = trim(io%dir_bwdmesh)//'PZ/'
        io%bwd(2)%meshdir = trim(io%dir_bwdmesh)//'PX/'
    else
        io%bwd(1)%meshdir = trim(io%dir_bwdmesh)
    endif

end subroutine define_sem_kernel_meshes
!------------------------------------------------

!!-------------------------------------------------
!subroutine save_simulation_info(thr,phr)
!
!!use data_arrays
!!use data_mesh
!implicit none
!
!real(kind=realkind), intent(in) :: thr,phr
!character(len=6) :: nnproc
!integer :: i,iwin,nnnproc
!character(len=20) :: stf_type
!real(kind=realkind) :: src_depthbwd
!
!  if (mynum==0) write(6,*)'saving Data/simulation.info ........' 
!!af changed all read into write for simulation info
!!tnm changed all af-changed write into read for simulation info.
!
!   open(unit=45,file='Data/simulation.info',status='unknown')
!!   read(45,*)nnproc,nnnproc
!   write(45,15)'nproc=',nproc
!   write(45,15)'npts=',npts
!
!12 format(a20,a15,a15)
!13 format(a20,1pe11.3)
!14 format(a20,a15)
!15 format(a20,i10)
!16 format(a20,3(a5))
!17 format(a20,2(1pe11.3))
!18 format(a20,i10,2(1pe11.3))
!19 format(a20,i10,1pe11.3)
!20 format(a20,l10)
!21 format(a20,3(l10))
!28 format(a20,i10,2(i10),2(1pe11.3))
!
!   write(45,12)'fwdsrctype= ',src_type1fwd,src_type2fwd
!   write(45,13)'fwdsrcdepth=',src_depth
!   write(45,14)'fwdstf_type= ',stf_type
!   write(45,15)'dim_fwd= ',dim_fwd
!
!   open(unit=20000,file='sourceparams_bwd.dat')
!   read(20000,*)realjunk
!   read(20000,*)src_type1bwd
!   read(20000,*)src_type2bwd
!   read(20000,*)
!   read(20000,*)
!   read(20000,*)src_depthbwd
!   read(20000,*)realjunk
!   read(20000,*)realjunk
!   read(20000,*)stf_type
!   close(20000)
!
!   dim_bwd = 3; if (src_type1bwd=='monopole') dim_bwd = 2; 
!
!   write(45,12)'bwdsrctype= ',src_type1bwd,src_type2bwd
!   write(45,13)'bwdsrcdepth=',src_depthbwd
!   write(45,14)'bwdstf_type= ',stf_type
!   write(45,15)'dim_bwd= ',dim_bwd
!
!   write(45,13)'thetar=',thr*180./pi
!   write(45,13)'phir=',phr*180./pi
!   write(45,20)'save_wfkern=',save_wfkern
!   write(45,21)'rho_lam_mu=',do_rho,do_lam,do_mu
!   write(45,21)'vp_vs_imped=',do_vp,do_vs,do_imped
!   write(45,13)'dt=',dt
!   write(45,17)'twf1_twf2=',twf1,twf2
!   write(45,15)'nwindows=',nwin
!   do i=1,nwin
!      write(45,28)'iwin_t1_t2=',i,begdumps(i),enddumps(i),time(begdumps(i)),time(enddumps(i))
!   enddo
!   write(45,15)'ntw=',ntw
!   do i=1,ntw
!      write(45,18)'ttwin_t1_t2=',i,t1(i),t2(i)
!   enddo
!   write(45,15)'nwfkernels=',num_it
!   do iwin=1,nwin
!      do i=begdumps(iwin),enddumps(iwin)
!         write(45,19)'kernid_time=',i,time(i)
!      enddo
!   enddo
!   close(45)
!   write(6,*)'saved simulation info into Data/simulation.info'
!
!
!
!end subroutine save_simulation_info
!-------------------------------------------------

!========================
end module parameters
!========================
