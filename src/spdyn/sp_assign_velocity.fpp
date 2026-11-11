!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_assign_velocity_mod
!> @brief   generate initial velocities and remove trans- and rotational motions
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_assign_velocity_mod

  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_dynvars_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! subroutines
  public  :: initial_velocity
  public  :: stop_translation
  private :: reduce_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_velocity
  !> @brief        generate initial velocities
  !! @authors      JJ
  !! @param[in]    firstt     : initial temperature
  !! @param[in]    natom      : number of atom
  !! @param[in]    num_degree : degree of freedom
  !! @param[in]    ncell      : total number of cells in each domain
  !! @param[in]    id_g2l     : map of atom index (global -> local)
  !! @param[in]    inv_mass   : inverse mass in each domain
  !! @param[inout] iseed      : random number seed
  !! @param[inout] velocity   : velocities generated
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initial_velocity(firstt, natom, num_degree, ncell, id_g2l, &
                              inv_mass, iseed, velocity)

    ! formal arguments
    real(wip),               intent(in)    :: firstt
    integer,                 intent(in)    :: natom(:)
    integer(iintegers),      intent(in)    :: num_degree
    integer,                 intent(in)    :: ncell
    integer(int2),           intent(in)    :: id_g2l(:,:)
    real(wip),               intent(in)    :: inv_mass(:,:)
    integer,                 intent(inout) :: iseed
    real(wip),               intent(inout) :: velocity(:,:,:)

    ! local variables
    real(wip)                :: alpha1, beta1, alpha2, beta2, alpha3, beta3
    real(wip)                :: sigma, myu, kBT
    integer                  :: ix, icel


    if (main_rank) then
      write(DynvarsOut,'(A)') 'Initial_Velocity> Generate initial velocities'
      write(DynvarsOut,'(A20,I10)')   '  iseed           = ', iseed
      write(DynvarsOut,'(A20,F10.3)') '  temperature     = ', firstt
      write(DynvarsOut,'(A)') ' '
    end if

    ! generate initial velocities 
    !
    kBT = KBOLTZ * firstt
    myu = 0.0_wip

    do icel = 1, ncell

      do ix = 1, natom(icel)

        alpha1 = sqrt( -2.0_wip * log(real(random_get_legacy(iseed),wip)))
        beta1  = cos ( 2.0_wip * PI * real(random_get_legacy(iseed),wip) )

        alpha2 = sqrt( -2.0_wip * log(real(random_get_legacy(iseed),wip)))
        beta2  = cos ( 2.0_wip * PI * real(random_get_legacy(iseed),wip) )

        alpha3 = sqrt( -2.0_wip * log(real(random_get_legacy(iseed),wip)))
        beta3  = cos ( 2.0_wip * PI * real(random_get_legacy(iseed),wip) )

        sigma = sqrt(kBT*inv_mass(ix,icel))
        velocity(1,ix,icel) = sigma * alpha1 * beta1 + myu
        velocity(2,ix,icel) = sigma * alpha2 * beta2 + myu
        velocity(3,ix,icel) = sigma * alpha3 * beta3 + myu

      end do

    end do
   
    return 

  end subroutine initial_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    stop_translation
  !> @brief        remove translational motion about the center of mass
  !! @authors      JJ
  !! @param[in]    ncell      : number of cell
  !! @param[in]    natom      : number of atom in each cell
  !! @param[in]    stop_trans : flag for stop translational motion
  !! @param[in]    mass       : mass in each domain
  !! @param[inout] velocity   : velocities in each domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine stop_translation(ncell, natom, stop_trans, mass, velocity)

    ! formal arguments
    integer,                 intent(in)    :: ncell
    integer,                 intent(in)    :: natom(:)
    logical,                 intent(in)    :: stop_trans
    real(wip),               intent(in)    :: mass(:,:)
    real(wip),               intent(inout) :: velocity(:,:,:)

    ! local variables
    real(dp )                :: total_mass, vcom(3)
    integer                  :: i, ix


    if (.not. stop_trans) return

    ! calculate the center of mass of velocities
    !
    vcom(1:3)  = 0.0_dp 
    total_mass = 0.0_dp 
    do i = 1, ncell
      do ix = 1, natom(i)
        vcom(1:3)  = vcom(1:3) + velocity(1:3,ix,i)*mass(ix,i)
        total_mass = total_mass + mass(ix,i)
      end do
    end do

    call reduce_com(vcom, total_mass)
    
    vcom(1:3) = vcom(1:3) / total_mass

    ! remove translational motion
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        velocity(1:3,ix,i) = velocity(1:3,ix,i) - real(vcom(1:3),wip)
      end do
    end do

    return

  end subroutine stop_translation

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_com
  !> @brief
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_com(val1, val2)

    ! formal arguments
    real(dp),               intent(inout) :: val1(:), val2

    ! local variables
    real(dp)                :: before_reduce(4), after_reduce(4)


    before_reduce(1:3) = val1(1:3)
    before_reduce(4)   = val2

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(before_reduce, after_reduce, 4, mpi_real8, &
                       mpi_sum, mpi_comm_city, ierror)
#else
    after_reduce(1:4) = before_reduce(1:4)
#endif

    val1(1:3)    = after_reduce(1:3)
    val2         = after_reduce(4)

    return

  end subroutine reduce_com

end module sp_assign_velocity_mod
