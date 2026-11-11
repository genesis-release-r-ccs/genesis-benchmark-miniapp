!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_setup_spdyn
!> @brief   setup variables and structures in MD (DD) simulaton
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_setup_spdyn_mod

  use sp_control_mod
  use sp_restart_mod
  use sp_output_mod
  use sp_input_mod
  use sp_dynamics_mod
  use sp_dynvars_mod
  use sp_domain_mod
  use sp_ensemble_mod
  use sp_enefunc_table_mod
  use sp_restraints_mod
  use sp_constraints_mod
  use sp_boundary_mod
  use sp_pairlist_mod
  use sp_enefunc_mod
  use sp_energy_mod
  use sp_energy_pme_mod
  use sp_update_domain_mod
  use sp_communicate_mod
  use sp_output_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_enefunc_fit_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use molecules_mod
  use molecules_str_mod
  use fitting_mod
  use fitting_str_mod
  use fileio_localres_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use fileio_crd_mod
  use fileio_rst_mod
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_spdyn_md

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_md
  !> @brief        setup variables and structures in MD simulation
  !! @authors      JJ, HO
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_md(ctrl_data, output, molecule, enefunc, pairlist,    &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_localres)         :: localres

    ! read input files
    !
    call input_md(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                  pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                  localres)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, psf, ref, &
                          prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)

    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if

    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        rst, boundary)
  
      call setup_domain(ctrl_data%ene_info,  &
                        ctrl_data%cons_info, &
                        boundary, molecule, enefunc, constraints, domain)

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info,  &
                        par, prmtop, grotop, &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc, comm)

    call setup_fitting_spdyn(ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    ! setup experiments
    !

    call dealloc_restraints_all(restraints)
    call dealloc_localres(localres, LocalRestraint)

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if


    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info,   &
                        molecule, dynamics)

    ! mass repartitioning
    !
    if (dynamics%hydrogen_mr) &
      call setup_mass_repartitioning(dynamics, constraints, domain, &
                                     enefunc)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)


    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, domain, &
                           enefunc, constraints)

    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)


    ! set output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)


    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Md> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    ! pressure degree of freedom for group pressure
    !
    if (ensemble%group_tp) then
      if (constraints%rigid_bond) then
        call compute_group_deg_freedom(domain, constraints)
      else
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Setup_Spdyn_Md> group_tp = No is assigned with rigid_bond = No'
        ensemble%group_tp = .false.
      end if
    end if

    return

  end subroutine setup_spdyn_md

end module sp_setup_spdyn_mod
