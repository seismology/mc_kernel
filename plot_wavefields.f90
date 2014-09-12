!=========================================================================================
module plot_wavefields_mod

    implicit none

contains

!-----------------------------------------------------------------------------------------
subroutine plot_wavefields()

    use global_parameters,          only : sp, dp, pi, deg2rad, init_random_seed
    use inversion_mesh,             only : inversion_mesh_data_type
    use readfields,                 only : semdata_type
    use type_parameter,             only : parameter_type
    use fft,                        only : rfft_type, taperandzeropad
    use filtering,                  only : timeshift_type

    type(inversion_mesh_data_type)      :: inv_mesh
    type(parameter_type)                :: parameters
    type(semdata_type)                  :: sem_data
    type(rfft_type)                     :: fft_data
    type(timeshift_type)                :: timeshift_fwd, timeshift_bwd

    integer                             :: nelems, ntimes, nomega, nrec
    integer                             :: idump, ndumps, irec, icomp
    integer                             :: nvertices, ndim
    real(kind=dp),    allocatable       :: co_points(:,:)
    real(kind=dp),    allocatable       :: fw_field(:,:,:)
    real(kind=dp),    allocatable       :: bw_field(:,:,:)
    complex(kind=dp), allocatable       :: fw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: bw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: conv_field_fd(:,:)
    real(kind=dp),    allocatable       :: conv_field(:,:)
    real(kind=dp),    allocatable       :: modelcoeffs(:,:)
    real(kind=dp)                       :: df
    character(len=64)                   :: fmtstring
    character(len=64)                   :: cname
    character(len=64)                   :: rname

    write(*,*) '***************************************************************'
    write(*,*) ' Read input files for parameters, source and receivers'
    write(*,*) '***************************************************************'
    call parameters%read_parameters()
    call parameters%read_source()
    call parameters%read_receiver()

    nrec = size(parameters%receiver)


    write(*,*) '***************************************************************'
    write(*,*) ' Initialize and open AxiSEM wavefield files'
    write(*,*) '***************************************************************'
    call sem_data%set_params(parameters%fwd_dir,     &
                             parameters%bwd_dir,     &
                             parameters%strain_buffer_size, & 
                             parameters%displ_buffer_size, & 
                             parameters%strain_type)
    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    ! Read seismogram
    call sem_data%load_seismogram_rdbm(parameters%receiver, parameters%source)
!    call sem_data%load_seismogram(parameters%receiver, parameters%source)

    ndumps = sem_data%ndumps
    ndim   = sem_data%get_ndim()

    write(*,*) '***************************************************************'
    write(*,*) ' Read inversion mesh'
    write(*,*) '***************************************************************'
    call inv_mesh%read_abaqus_mesh(parameters%mesh_file,parameters%int_type)

    nvertices = inv_mesh%get_nvertices()
    nelems    = inv_mesh%get_nelements()
    fmtstring = '(A, I8, A, I8)'
    print fmtstring, '  nvertices: ',  nvertices, ', nelems: ', nelems


    write(*,*) '***************************************************************'
    write(*,*) ' Initialize FFT'
    write(*,*) '***************************************************************'
    call fft_data%init(ndumps, ndim=ndim, ntraces=nvertices, dt=sem_data%dt)
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    df     = fft_data%get_df()
    fmtstring = '(A, I8, A, I8)'
    print fmtstring, '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
    fmtstring = '(A, F8.3, A, F8.3, A)'
    print fmtstring, '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

    allocate(co_points(3, nvertices))
    co_points = inv_mesh%get_vertices()

    write(*,*) ' Read in forward field'
    allocate(fw_field(ndumps, ndim, nvertices))
    allocate(modelcoeffs(3,size(co_points)))
    fw_field(:,:,:) = sem_data%load_fw_points(dble(co_points), parameters%source, coeffs=modelcoeffs)

    write(*,*) ' FFT forward field'
    allocate(fw_field_fd(nomega, ndim, nvertices))
    call fft_data%rfft(taperandzeropad(fw_field, ntimes, 2), fw_field_fd)
    deallocate(fw_field)

    write(*,*) ' Timeshift forward field'
    call timeshift_fwd%init(fft_data%get_f(), sem_data%timeshift_fwd)
    call timeshift_fwd%apply(fw_field_fd)
    allocate(fw_field(ntimes, ndim, nvertices))

    call fft_data%irfft(fw_field_fd, fw_field)
    call timeshift_fwd%freeme()

    print *, ' Initialize XDMF file'
    call inv_mesh%init_node_data(ndumps*ndim)

    write(*,*) ' Dump forward field to XDMF file'
    do idump = 1, ndumps
        if (mod(idump, 100)==0) &
            write(*,*) '  Passing dump ', idump, ' to inversion mesh datatype'
        !Test of planar wave , works
        !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
        !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
        do icomp=1,ndim
           write(cname,'("comp_",I0.2)') icomp
           call inv_mesh%set_node_data_snap(real(fw_field(idump,icomp,:), kind=sp), &
                                            idump + ndumps*(icomp-1), &
                                            'fwd_'//trim(cname))
        end do
    end do
    deallocate(fw_field)

    print *, ' Save XDMF file'
    call inv_mesh%dump_node_data_xdmf(trim(parameters%output_file)//'_wavefield_fwd')
    call inv_mesh%free_node_and_cell_data()




    do irec = 1, nrec

        write(rname,'("_",I0.2)') irec

        write(*,*) ' Read in backward field of receiver', irec
        allocate(bw_field(ndumps, ndim, nvertices))
        bw_field(:,:,:) = sem_data%load_bw_points(dble(co_points), &
                                                  parameters%receiver(irec))

        allocate(bw_field_fd  (nomega, ndim, nvertices))
        call fft_data%rfft(taperandzeropad(bw_field, ntimes, 2), bw_field_fd)
        deallocate(bw_field)
        
        write(*,*) ' Timeshift backward field'
        call timeshift_bwd%init(fft_data%get_f(), sem_data%timeshift_bwd)
        call timeshift_bwd%apply(bw_field_fd)
        allocate(bw_field(ntimes, ndim, nvertices))
        
        call fft_data%irfft(bw_field_fd, bw_field)
        call timeshift_bwd%freeme()

        call inv_mesh%init_node_data(ndumps*ndim)

        write(*,*) ' Dump backward field to XDMF file'
        do idump = 1, ndumps
            if (mod(idump, 100)==0) &
                write(*,*) '  Passing dump ', idump, ' to inversion mesh datatype'
            !Test of planar wave , works
            !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
            !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)

            do icomp = 1,ndim
               write(cname,'("comp_",I0.2)') icomp
               call inv_mesh%set_node_data_snap(real(bw_field(idump,icomp,:), kind=sp), &
                                                idump + ndumps*(icomp-1), &
                                                'bwd_'//trim(parameters%receiver(irec)%name)//'_'//trim(cname))
            end do
        end do
        deallocate(bw_field)

        print *, ' Save XDMF file'
        call inv_mesh%dump_node_data_xdmf(trim(parameters%output_file)//'_wavefield_'//trim(rname)//'_bwd')
        call inv_mesh%free_node_and_cell_data()

        write(*,*) ' Convolve wavefields'
        allocate(conv_field_fd(nomega, nvertices))

        ! sum over dimension 2 necessary for vs kernels
        conv_field_fd = sum(fw_field_fd * bw_field_fd, 2)
        deallocate(bw_field_fd)
        
        allocate(conv_field(ntimes, nvertices))
        call fft_data%irfft(conv_field_fd, conv_field)
        deallocate(conv_field_fd)

        call inv_mesh%init_node_data(ndumps)

        write(*,*) ' Dump convolved fields to XDMF file'
        do idump = 1, ndumps
           if (mod(idump, 100)==0) write(*,*) ' Passing dump ', idump, ' of convolved wavefield'
           call inv_mesh%set_node_data_snap(real(conv_field(idump,:), kind=sp), &
                                            idump, &
                                            'convolved_'//trim(parameters%receiver(irec)%name))
        end do 
        deallocate(conv_field)

        write(*,*) ' Save XDMF file'
        call inv_mesh%dump_node_data_xdmf(trim(parameters%output_file)//'_wavefield_'//trim(rname)//'_conv')
        call inv_mesh%free_node_and_cell_data()


    end do ! irec

    deallocate(fw_field_fd)

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Free memory of inversion mesh datatype'
    write(*,*) '***************************************************************'
    call inv_mesh%freeme

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Close AxiSEM wavefield files'
    write(*,*) '***************************************************************'
    call sem_data%close_files()

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Free memory of FFT datatype'
    write(*,*) '***************************************************************'
    call fft_data%freeme()

    write(*,*)
    write(*,*) ' Finished!'

end subroutine plot_wavefields
!-----------------------------------------------------------------------------------------

end module plot_wavefields_mod
!=========================================================================================

