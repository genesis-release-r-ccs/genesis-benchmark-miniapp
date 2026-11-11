!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  SPDYN
!! @brief   Molecular dynamics Simulation of BioMolecules using
!!          Spacial Decomposition Scheme
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Yuji Sugita (YS), 
!!          Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

program spdyn

  use sp_setup_mpi_mod
  use sp_md_vverlet_mod
  use sp_md_leapfrog_mod
  use sp_dynamics_mod
  use sp_setup_spdyn_mod
  use sp_control_mod
  use sp_energy_mod
  use sp_energy_pme_mod
  use sp_communicate_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_output_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use hardwareinfo_mod
  use timers_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

  ! local variables
  character(MaxFilename)      :: ctrl_filename
  type(s_ctrl_data)           :: ctrl_data
  type(s_molecule)            :: molecule
  type(s_enefunc)             :: enefunc
  type(s_dynvars)             :: dynvars
  type(s_pairlist)            :: pairlist
  type(s_boundary)            :: boundary
  type(s_constraints)         :: constraints
  type(s_ensemble)            :: ensemble
  type(s_dynamics)            :: dynamics
  type(s_output)              :: output
  type(s_domain)              :: domain
  type(s_comm)                :: comm
  integer                     :: omp_get_max_threads

#ifdef HAVE_MPI_GENESIS
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world, my_world_rank, ierror)
  call mpi_comm_size(mpi_comm_world, nproc_world,   ierror)
  main_rank = (my_world_rank == 0)
#else
  my_world_rank = 0
  nproc_world   = 1
  main_rank     = .true.
#endif

#ifdef OMP
  nthread = omp_get_max_threads()
#else
  nthread = 1
#endif
 
  ! show usage
  !
  call usage(ctrl_filename)


  ! check run mode in control file
  !
  call check_genesis_mode(ctrl_filename)


  ! run genesis
  !
  call domain_decomposition_genesis(ctrl_filename)

#ifdef HAVE_MPI_GENESIS
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_finalize(ierror)
#endif

  stop

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_genesis_mode
  !> @brief        check that the control file is for MD
  !! @authors      TM
  !! @param[in]    ctrl_filename    : control file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_genesis_mode(ctrl_filename)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename


    if (find_ctrlfile_section(ctrl_filename, 'ALCHEMY')  .or. &
        find_ctrlfile_section(ctrl_filename, 'REMD')     .or. &
        find_ctrlfile_section(ctrl_filename, 'RPATH')    .or. &
        find_ctrlfile_section(ctrl_filename, 'MINIMIZE')) then
      call error_msg('Check_Genesis_Mode> Only MD is available in mini APP.')
    end if

    if (.not. find_ctrlfile_section(ctrl_filename, 'DYNAMICS')) &
      call error_msg('Check_Genesis_Mode> ERROR: DYNAMICS section is not found.')

    return

  end subroutine check_genesis_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_decomposition_genesis
  !> @brief        run genesis using domain decomposition scheme
  !! @authors      TM
  !! @param[in]    ctrl_filename   : control file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_decomposition_genesis(ctrl_filename)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename

#ifdef USE_GPU
    integer :: my_device_id
#endif


    ! set timer
    !
    call timer(TimerTotal, TimerOn)

    ! [Step0] Architecture & Compiler information
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP0] Architecture and Compiler Information'
      write(MsgOut,'(A)') ' '

      call hw_information

#ifdef USE_GPU
    else
      ! assign GPU
      call assign_gpu(my_device_id)
#endif /* USE_GPU */
    end if

    ! [Step1] Read control file
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP1] Read Control Parameters'
      write(MsgOut,'(A)') ' '
    end if

      call control_md  (ctrl_filename, ctrl_data)


    ! [Step2] Setup MPI
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2] Setup MPI'
      write(MsgOut,'(A)') ' '
    end if

      call setup_mpi_md  (ctrl_data%ene_info)

    ! [Step3] Set relevant variables and structures 
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP3] Set Relevant Variables and Structures'
      write(MsgOut,'(A)') ' '
    end if

      call setup_spdyn_md(ctrl_data, output, molecule, enefunc, pairlist,    &
                         dynvars, dynamics, constraints, ensemble, boundary, &
                         domain, comm)


    ! [Step4] Compute single point energy for molecules
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP4] Compute Single Point Energy for Molecules'
      write(MsgOut,'(A)') ' '
    end if

    call compute_energy(domain, enefunc, pairlist, boundary, domain%coord, &
                        .false., .true., .true., .true., & 
                        enefunc%nonb_limiter,            &
                        dynvars%energy,                  &
                        domain%atmcls_pbc,               &
                        domain%translated,               &
                        domain%force,                    &
                        domain%force_long,               &
                        domain%force_omp,                &
                        domain%force_pbc,                &
                        domain%virial_cellpair,          &
                        dynvars%virial,                  &
                        dynvars%virial_long,             &
                        dynvars%virial_extern)

    call output_energy(dynvars%step, enefunc, dynvars%energy)


#ifdef HAVE_MPI_GENESIS
    call mpi_barrier(mpi_comm_world, ierror)
#endif
    call timer(TimerDynamics, TimerOn)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Molecular Dynamics Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_md(output, domain, enefunc, dynvars, dynamics, pairlist, &
                  boundary, constraints, ensemble, comm)
    
    call timer(TimerDynamics, TimerOff)
#ifdef GENESIS_GPU_TIMER_VERBOSE
    call gpu_print_kernel_timer()
#endif


    ! [Step6] Deallocate arrays
    !
    if (main_rank) then
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '[STEP6] Deallocate Arrays'
      write(MsgOut,'(A)') ' '
    end if

    call dealloc_pme(enefunc)
    call dealloc_constraints_all(constraints)
    call dealloc_boundary_all   (boundary)
    call dealloc_pairlist_all   (pairlist)
    call dealloc_enefunc_all    (enefunc)
    call dealloc_domain_all     (domain)

    call timer(TimerTotal, TimerOff)


    ! output process time
    !
    call output_time
    call output_time_gpu

    return

  end subroutine domain_decomposition_genesis

end program spdyn
