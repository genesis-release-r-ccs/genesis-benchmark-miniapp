!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_gpu_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_gpu_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use timers_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond_tbl_lnr_gpu
  public  :: compute_energy_nonbond_tbl_ljpme_gpu
  public  :: compute_force_nonbond_tbl_lnr_gpu
  public  :: compute_force_nonbond_tbl_ljpme_gpu
  public  :: compute_energy_nonbond_notbl_gpu
  public  :: compute_force_nonbond_notbl_gpu

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_gpu
  !> @brief        calculate nonbonded energy with gpu
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_gpu( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,ixx,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_linear_univ( &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_tbl_lnr_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_ljpme_gpu
  !> @brief        calculate nonbonded energy with gpu (LJ-PME)
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_ljpme_gpu( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_ljpme_univ(  &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_tbl_ljpme_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_gpu
  !> @brief        calculate nonbonded energy with gpu
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_gpu( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,ixx,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_notable_univ(      &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_notbl_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_gpu( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_table_linear_univ( &
         coord_pbc, force(:,:,:,1), ene_virial,              &            
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check, domain%num_atom_domain, MaxAtom,      &
         MaxAtomCls, num_atom_cls,                           &            
         ncell_local, ncell_bound, ncell_max, cutoff_int,    &           
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_tbl_lnr_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_ljpme_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_ljpme_gpu( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_table_ljpme_univ(  &
         coord_pbc, force(:,:,:,1), ene_virial,              &  
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check, domain%num_atom_domain, MaxAtom,      &
         MaxAtomCls, num_atom_cls,                           &  
         ncell_local, ncell_bound, ncell_max, cutoff_int,    & 
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_tbl_ljpme_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_gpu( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_notable_univ(      &
         coord_pbc, force(:,:,:,1), ene_virial,              &  
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check, domain%num_atom_domain, MaxAtom,      &
         MaxAtomCls, num_atom_cls,                           &  
         ncell_local, ncell_bound, ncell_max, cutoff_int,    & 
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_notbl_gpu

end module sp_energy_nonbond_gpu_mod
