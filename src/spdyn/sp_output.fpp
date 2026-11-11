!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_output_md
!> @brief   definition of output
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Yasuhiro Matsunaga (YM),
!!          Ryuhei Harada (RH), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_output_mod

  use sp_dynvars_mod
  use sp_output_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_constraints_mod
  use molecules_str_mod
  use fileio_control_mod
  use fileio_rst_mod
  use fileio_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_out_info
    character(MaxFilename) :: dcdfile    = ''
    character(MaxFilename) :: dcdvelfile = ''
    character(MaxFilename) :: rstfile    = ''
    character(MaxFilename) :: pdbfile    = ''
  end type s_out_info

  ! valiables
  real(wip),        private, allocatable :: tmp_coord1(:,:)
  real(wip),        private, allocatable :: tmp_coord2(:,:)

  ! subroutines
  public  :: show_ctrl_output
  public  :: read_ctrl_output
  public  :: setup_output_md
  public  :: open_output
  public  :: close_output
  public  :: output_md

  private :: output_restart_md
  private :: write_trajectory_dcd
  private :: write_trajectory_dcdvel
  private :: reduce_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_output
  !> @brief        show OUTPUT section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_output(show_all)

    ! formal arguments
    logical,                 intent(in)    :: show_all


    if (show_all) then

      write(MsgOut,'(A)') '[OUTPUT]'
      write(MsgOut,'(A)') '# dcdfile    = sample.dcd    # DCD trajectory file'
      write(MsgOut,'(A)') '# dcdvelfile = sample.dcd    # DCD velocity file'
      write(MsgOut,'(A)') '# rstfile    = sample.rst    # restart file'
      write(MsgOut,'(A)') '# pdbfile    = sample.pdb    # PDB file'
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine show_ctrl_output
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_output
  !> @brief        read OUTPUT section in the control file
  !! @authors      TM, JJ
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   out_info : OUTPUT section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_output(handle, out_info)

    ! parameters
    character(*),            parameter     :: Section = 'Output'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_out_info),        intent(inout) :: out_info


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string(handle, Section, 'dcdfile',   out_info%dcdfile)
    call read_ctrlfile_string(handle, Section, 'dcdvelfile',out_info%dcdvelfile)
    call read_ctrlfile_string(handle, Section, 'rstfile',   out_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'pdbfile',   out_info%pdbfile)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Output> Output Files'
      if (out_info%dcdfile .ne. '') &
        write(MsgOut,*) ' dcdfile    = ', trim(out_info%dcdfile)
      if (out_info%dcdvelfile .ne. '') &
        write(MsgOut,*) ' dcdvelfile = ', trim(out_info%dcdvelfile)
      if (out_info%rstfile .ne. '') &
        write(MsgOut,*) ' rstfile    = ', trim(out_info%rstfile)
      if (out_info%pdbfile .ne. '') &
        write(MsgOut,*) ' pdbfile    = ', trim(out_info%pdbfile)
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_md
  !> @brief        setup output information for MD
  !! @authors      TM
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_md(out_info, dynamics, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_output),          intent(inout) :: output


    output%verbose = dynamics%verbose

    if (dynamics%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') &
        call error_msg('Setup_Output_Md> Error: dcdfile name is not'//&
                  ' specified in [OUTPUT] (crdout_period > 0 in [DYNAMICS])')
      output%dcdout  = .true.
      output%dcdfile = out_info%dcdfile
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') &
        call error_msg('Setup_Output_Md> Error: dcdvelfile name is not'//&
                  ' specified in [OUTPUT] (velout_period > 0 in [DYNAMICS])')
      output%dcdvelout  = .true.
      output%dcdvelfile = out_info%dcdvelfile
    end if

    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Md> Error: rstfile name is not'//&
                       'specified in [OUTPUT] (rstout_period > 0 in [DYNAMICS])')
      output%rstout  = .true.
      output%rstfile = out_info%rstfile
    end if

    if (out_info%pdbfile .ne. '') then
      output%pdbout  = .true.
      output%pdbfile = out_info%pdbfile
    end if

    return

  end subroutine setup_output_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_output
  !> @brief        open output file
  !! @authors      TM
  !! @param[inout] outout : information of output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_output(output)

    ! formal arguments
    type(s_output),          intent(inout) :: output

    ! local variables
    integer                  :: file


    ! open dcdfile
    !
    if (output%dcdout) then

      if (main_rank) then

        call open_binary_file(file, output%dcdfile, IOFileOutputNew)
        output%dcdunit = file

      end if

    end if

    ! open dcdvelfile
    !
    if (output%dcdvelout) then

       if (main_rank) then

        call open_binary_file(file, output%dcdvelfile, IOFileOutputNew)
        output%dcdvelunit = file

      end if

    end if

    ! open rstfile
    !   Just check whether rstfile is already existing.
    !
    if (output%rstout) then

       if (main_rank) then

        call open_file(file, output%rstfile, IOFileOutputNew)
        call close_file(file)

      end if

    end if

    ! open pdbfile
    !   Just check whether pdbfile is already existing.
    !
    if (output%pdbout) then
      if (main_rank) then
        call open_file(file, output%pdbfile, IOFileOutputNew)
        call close_file(file)
      end if
    end if

    return

  end subroutine open_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_output
  !> @brief        close output files
  !! @authors      TM
  !! @param[in]    output : information of output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_output(output)

    ! formal arguments
    type(s_output),          intent(in)    :: output


    ! close dcdfile
    !
    if (output%dcdout) &
      call close_file(output%dcdunit)

    ! close dcdvelfile
    !
    if (output%dcdvelout) &
      call close_file(output%dcdvelunit)


    return

  end subroutine close_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_md
  !> @brief        output trajectory and restart data for MD
  !! @authors      JJ, TM, NT
  !! @param[inout] output      : output information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_md(output, dynamics, boundary, pairlist, ensemble, &
                       constraints, dynvars, domain, enefunc)

    ! formakl arguments
    type(s_output),          intent(inout) :: output
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(in)    :: constraints
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    ! local variables
    real(wip)                 :: scale_v, gamma_t, dt
    real(wip), allocatable    :: tmp_coord(:,:,:)
    integer                   :: i, ix, istep

    integer,          pointer :: ncell, natom(:)


    ncell        => domain%num_cell_local
    natom        => domain%num_atom

    gamma_t      =  ensemble%gamma_t*AKMA_PS
    dt           =  dynamics%timestep/AKMA_PS
    istep        =  dynvars%step

    ! Output coordinates (Wrap_all is not done)
    !
    if (dynamics%crdout_period > 0 .and. istep > 0) then

      if (mod(istep,dynamics%crdout_period) == 0) then

        if (.not. allocated(tmp_coord)) &
          allocate(tmp_coord(3, MaxAtom, ncell))

        tmp_coord(1:3,1:MaxAtom,1:ncell) = 0.0_wip

        if (dynamics%integrator == IntegratorLEAP) then
          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = domain%coord_ref(1:3,ix,i)
            end do
          end do
        else
          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = domain%coord(1:3,ix,i)
            end do
          end do
        end if

          call write_trajectory_dcd(output,            &
                                    boundary,          &
                                    domain,            &
                                    istep,             &
                                    dynamics%nsteps,   &
                                    dynamics%crdout_period, &
                                    dynamics%timestep, &
                                    tmp_coord)

      end if
    end if


    ! Output velocities (Wrap_all is not done)
    !
    if (dynamics%velout_period > 0 .and. istep > 0) then
      if (mod(istep,dynamics%velout_period) == 0) then
        if (.not. allocated(tmp_coord)) &
          allocate(tmp_coord(3, MaxAtom, ncell))

        tmp_coord(1:3,1:MaxAtom,1:ncell) = 0.0_wip

        if (dynamics%integrator == IntegratorLEAP) then
          if (ensemble%tpcontrol == TpcontrolLangevin) then
            scale_v = sqrt(1.0_wip+0.5_wip*gamma_t*dt)
          else
            scale_v = 1.0_wip
          end if

          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = 0.5_wip * &
                                    (domain%velocity_ref(1:3,ix,i) + &
                                     domain%velocity    (1:3,ix,i))*scale_v
            end do
          end do
        else
          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = domain%velocity(1:3,ix,i)
            end do
          end do
        end if

          call write_trajectory_dcdvel(output,            &
                                       boundary,          &
                                       domain,            &
                                       istep,             &
                                       dynamics%nsteps,   &
                                       dynamics%velout_period, &
                                       dynamics%timestep, &
                                       tmp_coord)

      end if
    end if


    ! Output restart data
    !
    if (dynamics%rstout_period > 0 .and. istep > 0) then
      if (mod(istep,dynamics%rstout_period) == 0) then

          call output_restart_md(output, domain, dynamics, boundary, dynvars)

      end if
    end if

    return

  end subroutine output_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_md
  !> @brief        output restart data for MD
  !! @authors      NT
  !! @param[in]    output   : output information
  !! @param[in]    domain   : domain information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_md(output, domain, dynamics, boundary, dynvars)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_domain),  target, intent(in)    :: domain
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars

    ! local variables
    type(s_rst)              :: rst
    integer                  :: i, ix, istep
    
    real(wip),       pointer :: coord(:,:,:), vel(:,:,:)
    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    coord     => domain%coord
    vel       => domain%velocity
    istep     =  dynvars%step


    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! setup restart information
    !
    call alloc_rst(rst, RestartAtom, int(natom_all))

    rst%rstfile_type = RstfileTypeMd
    rst%iseed                  = dynamics%iseed_init_velocity
    rst%num_atoms              = natom_all
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars%thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars%barostat_momentum(1:3)

    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    rst%coord(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! reduce velocities
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = vel(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    rst%velocity(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! output restart information
    !
    if (main_rank) &
      call output_rst(output%rstfile, rst)

    if (istep == dynamics%nsteps) &
      call dealloc_rst_all(rst)

    return

  end subroutine output_restart_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trajectory_dcd
  !> @brief        write coordinates in DCD format with domain decomposition
  !! @authors      RH, TM, JJ
  !! @param[inout] output    : output information
  !! @param[in]    boundary  : boundary information
  !! @param[in]    domain    : dmain information
  !! @param[in]    istep     : dynamics step
  !! @param[in]    nstep     : number of dynamics steps
  !! @param[in]    outperiod : coordinate output period
  !! @param[in]    timestep  : time step in PS
  !! @param[in]    coord     : coordinates per domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trajectory_dcd(output, boundary, domain, istep, nstep, &
                                  outperiod, timestep, coord)

    ! parameters
    character(4),            parameter     :: HEADER  = 'CORD'
    integer(4),              parameter     :: NTITLE  = 4
    real(wip),               parameter     :: TIMEFAC = 48.88821

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),  target, intent(in)    :: domain
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    integer,                 intent(in)    :: outperiod
    real(wip),               intent(in)    :: timestep
    real(wip),               intent(in)    :: coord(:,:,:)

    ! local variable
    real(dp)                 :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: file
    integer                  :: i, ix
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date

    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    file      =  output%dcdunit

    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do
    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    if (.not. main_rank) return
  

    ! write header
    !
    if (output%out_dcd_header) then

      ts_namd = real(timestep * 1000_wip / TIMEFAC, sp)

      icntrl(:)  = 0
      icntrl(1)  = nstep / outperiod            ! total # of frames   (NSET)
      icntrl(2)  = istep                        ! first time step     (NPRIV)
      icntrl(3)  = outperiod                    ! output period       (NSAVC)
      icntrl(4)  = nstep                        ! number of time step (NSTEP)
      icntrl(10) = transfer(ts_namd,icntrl(10)) ! length of time step (DELTA)
      icntrl(11) = 1                            ! flag for with unit-cell
      icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

      call fdate (date)
      call getlog(name)

      title(1) = 'REMARKS CREATED BY GENESIS                                                      '
      title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
      title(3) = '                                                                                '
      title(4) = '                                                                                '

      write(file) HEADER, icntrl(1:20)
      write(file) NTITLE,(title(i),i=1,NTITLE)
      write(file) natom_all

      output%out_dcd_header = .false.

    end if

    ! write coordinates
    !
    dcd_cell(1:6) = 0.0_dp
    dcd_cell(1)   = boundary%box_size_x_ref
    dcd_cell(3)   = boundary%box_size_y_ref
    dcd_cell(6)   = boundary%box_size_z_ref
    write(file) dcd_cell(:)

    write(file) real(tmp_coord1(1,:),sp)
    write(file) real(tmp_coord1(2,:),sp)
    write(file) real(tmp_coord1(3,:),sp)

    return

  end subroutine write_trajectory_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trajectory_dcdvel
  !> @brief        write velocities in DCD format with domain decomposition
  !! @authors      RH, TM, YM
  !! @param[inout] output    : output information
  !! @param[in]    boundary  : boundary information
  !! @param[in]    domain    : dmain information
  !! @param[in]    istep     : dynamics step
  !! @param[in]    nstep     : number of dynamics steps
  !! @param[in]    outperiod : coordinate output period
  !! @param[in]    timestep  : time step in PS
  !! @param[in]    velocity  : velocities per domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trajectory_dcdvel(output, boundary, domain, istep, nstep, &
                                     outperiod, timestep, velocity)

    ! parameters
    character(4),            parameter     :: HEADER  = 'VELD'
    integer(4),              parameter     :: NTITLE  = 4
    real(wip),               parameter     :: TIMEFAC = 48.88821

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    integer,                 intent(in)    :: outperiod
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(in)    :: timestep
    real(wip),               intent(in)    :: velocity(:,:,:)

    ! local variable
    real(wip)                :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: file
    integer                  :: i, ix
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date

    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g


    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = velocity(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    if (.not. main_rank) return


    file    = output%dcdvelunit

    ! write header
    !
    if (output%out_dcdvel_header) then

      ts_namd = real(timestep * 1000_wip / TIMEFAC,sp)

      icntrl(:)  = 0
      icntrl(1)  = nstep / outperiod            ! total # of frames   (NSET)
      icntrl(2)  = istep                        ! first time step     (NPRIV)
      icntrl(3)  = outperiod                    ! output period       (NSAVC)
      icntrl(4)  = nstep                        ! number of time step (NSTEP)
      icntrl(10) = transfer(ts_namd,icntrl(10)) ! length of time step (DELTA)
      icntrl(11) = 1                            ! flag for with unit-cell
      icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

      call fdate (date)
      call getlog(name)

      title(1) = 'REMARKS CREATED BY GENESIS                                                      '
      title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
      title(3) = '                                                                                '
      title(4) = '                                                                                '

      write(file) HEADER, icntrl(1:20)
      write(file) NTITLE,(title(i),i=1,NTITLE)
      write(file) natom_all

      output%out_dcdvel_header = .false.

    end if

    ! write coordinates
    !
    dcd_cell(1:6) = 0.0_wip
    dcd_cell(1)   = boundary%box_size_x_ref
    dcd_cell(3)   = boundary%box_size_y_ref
    dcd_cell(6)   = boundary%box_size_z_ref
    write(file) dcd_cell(:)

    write(file) real(tmp_coord1(1,:),sp)
    write(file) real(tmp_coord1(2,:),sp)
    write(file) real(tmp_coord1(3,:),sp)

    return

  end subroutine write_trajectory_dcdvel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_coord
  !> @brief        reduce coord
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_coord(coord, temporary, natom)

    ! formal arguments
    real(wip),               intent(inout) :: coord(:,:)
    real(wip),               intent(inout) :: temporary(:,:)
    integer,                 intent(in)    :: natom

#ifdef HAVE_MPI_GENESIS

    ! local variables
    integer                  :: j
    integer                  :: ncycle, icycle, nlen, ixx


    ! Reduce coordinate
    !
    temporary(1:3,1:natom) = 0.0_wip
    do j = 1, natom
      temporary(1,j) = coord(1,j)
      temporary(2,j) = coord(2,j)
      temporary(3,j) = coord(3,j)
    end do

    ncycle = (natom - 1) / mpi_drain + 1
    nlen = mpi_drain
    ixx  = 1

    do icycle = 1, ncycle
      if (icycle == ncycle) nlen = natom - (ncycle-1) * mpi_drain
      call mpi_reduce(temporary(1,ixx), coord(1,ixx), 3*nlen,    &
                      mpi_wip_real, mpi_sum, 0, mpi_comm_country, &
                      ierror)
      ixx = ixx + nlen
    end do

#endif

    return

  end subroutine reduce_coord

end module sp_output_mod
