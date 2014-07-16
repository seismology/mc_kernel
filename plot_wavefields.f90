!=========================================================================================
module plot_wavefields_mod

contains

!-----------------------------------------------------------------------------------------
subroutine plot_wavefields()

    use global_parameters,  only : sp, dp, pi, deg2rad, verbose, init_random_seed
    use inversion_mesh,     only : inversion_mesh_data_type
    use readfields,         only : semdata_type
    use type_parameter,     only : parameter_type
    use fft,                only : rfft_type, taperandzeropad

    implicit none
    type(inversion_mesh_data_type)      :: inv_mesh
    type(parameter_type)                :: parameters
    type(semdata_type)                  :: sem_data
    type(rfft_type)                     :: fft_data

    integer                             :: nelems, ntimes, nomega, nrec
    integer                             :: idump, ndumps, irec
    integer                             :: nvertices
    real(kind=dp),    allocatable       :: co_points(:,:)
    real(kind=dp),    allocatable       :: fw_field(:,:,:)
    real(kind=dp),    allocatable       :: bw_field(:,:,:)
    complex(kind=dp), allocatable       :: fw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: bw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: conv_field_fd(:,:,:)
    real(kind=dp),    allocatable       :: conv_field(:,:,:)
    real(kind=dp)                       :: df
    character(len=64)                   :: fmtstring

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
                             parameters%buffer_size, & 
                             parameters%model_param)
    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    call sem_data%load_seismogram(parameters%receiver, parameters%source)

    ndumps = sem_data%ndumps


    write(*,*) '***************************************************************'
    write(*,*) ' Read inversion mesh'
    write(*,*) '***************************************************************'
    !call inv_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')
    !call inv_mesh%read_abaqus_mesh('unit_tests/tetrahedron.inp')
    !call inv_mesh%read_abaqus_mesh('unit_tests/flat_triangles.inp')
    call inv_mesh%read_abaqus_mesh(parameters%mesh_file,parameters%inttype)

    nvertices = inv_mesh%get_nvertices()
    nelems    = inv_mesh%get_nelements()
    fmtstring = '(A, I8, A, I8)'
    print fmtstring, '  nvertices: ',  nvertices, ', nelems: ', nelems


    write(*,*) '***************************************************************'
    write(*,*) ' Initialize FFT'
    write(*,*) '***************************************************************'
    call fft_data%init(ndumps, ndim=1, ntraces=nvertices, dt=sem_data%dt)
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    df     = fft_data%get_df()
    fmtstring = '(A, I8, A, I8)'
    print fmtstring, '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
    fmtstring = '(A, F8.3, A, F8.3, A)'
    print fmtstring, '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

    print *, 'Initialize XDMF file'
    allocate(co_points(3, nvertices))
    co_points = inv_mesh%get_vertices()
    call inv_mesh%init_node_data(ndumps + 2*ndumps*nrec)

    write(*,*) ' Read in forward field'
    allocate(fw_field(ndumps, 1, nvertices))
    fw_field(:,:,:) = sem_data%load_fw_points(dble(co_points), parameters%source)

    
    ! Dump forward field to XDMF file
    do idump = 1, ndumps
        if (mod(idump, 100)==0) &
            write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
        !Test of planar wave , works
        !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
        !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
        call inv_mesh%set_node_data_snap(real(fw_field(idump,1,:), kind=sp), &
                                         idump, 'fwd_wavefield')
    end do

    write(*,*) ' FFT forward field'
    allocate(fw_field_fd(nomega, 1, nvertices))
    call fft_data%rfft(taperandzeropad(fw_field, ntimes), fw_field_fd)
    deallocate(fw_field)

    do irec = 1, nrec
        write(*,*) ' Read in backward field of receiver', irec
        allocate(bw_field(ndumps, 1, nvertices))
        bw_field(:,:,:) = sem_data%load_bw_points(dble(co_points), &
                                                  parameters%receiver(irec))
        ! Dump backward field to XDMF file
        do idump = 1, ndumps
            if (mod(idump, 100)==0) &
                write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
            !Test of planar wave , works
            !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
            !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
            call inv_mesh%set_node_data_snap(real(bw_field(idump,1,:), kind=sp), &
                                             idump+(ndumps*irec), &
                                             'bwd_'//trim(parameters%receiver(irec)%name))
        end do
        write(*,*) ' FFT backward field'
        allocate(bw_field_fd  (nomega, 1, nvertices))
        call fft_data%rfft(taperandzeropad(bw_field, ntimes), bw_field_fd)
        deallocate(bw_field)

        write(*,*) ' Convolve wavefields'
        allocate(conv_field(ntimes, 1, nvertices))
        allocate(conv_field_fd(nomega, 1, nvertices))
        conv_field_fd = fw_field_fd * bw_field_fd
        call fft_data%irfft(conv_field_fd, conv_field)
        deallocate(bw_field_fd)
        deallocate(conv_field_fd)

    !   write(*,*) ' Apply filter'
    !   fw_field_fd = gabor20%apply_2d(fw_field_fd)

    !   write(*,*) ' iFFT product of fields'
    !   allocate(conv_field   (ntimes, nptperstep))
    !   call fft_data%irfft(fw_field_fd*bw_field_fd, conv_field)
    !
    !   !write(*,*) ' Write pickpoint to disc'
    !   !write(41,*) fw_field(:,10000)
    !   !write(42,*) bw_field(:,10000)
    !   !write(43,*) conv_field(:,10000)
    !   !read(*,*)


        do idump = 1, ndumps
           if (mod(idump, 100)==0) write(*,*) ' Passing dump ', idump, ' of convolved wavefield'
           call inv_mesh%set_node_data_snap(real(conv_field(idump,1,:), kind=sp), &
                                            idump + ndumps*2 + (irec-1)*ndumps, &
                                            'convolved_'//trim(parameters%receiver(irec)%name))
        end do 
        deallocate(conv_field)

    end do ! irec

    write(*,*)
    write(*,*) ' Writing data to disk'
    call inv_mesh%dump_node_data_xdmf('wavefield')
    call inv_mesh%freeme()

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

