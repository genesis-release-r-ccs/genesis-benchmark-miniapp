!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dihedral_libs_mod
!> @brief   libralies for dihedral angles
!! @authors Yasuhiro Matsunaga (YM), Chigusa Kobayashi (CK)
!! @date    2014/4/16-2014/4/20 (YM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module dihedral_libs_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  public  :: calculate_dihedral

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_dihedral
  !> @brief        calculate dihedral angle
  !! @authors      CK
  !! @param[in]    aindex  : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] cos_dih : cosin of dihedral angles
  !! @param[inout] sin_dih : sin of dihedral angles
  !! @param[inout] grad    : gradient of dihedral angles
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

    ! formal arguments
    integer,   intent(in)    :: aindex(:)
    real(wp),  intent(in)    :: coord(:,:)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(1:9)
    real(wp),  intent(inout) :: v(1:3,1:3)

    ! local variable
    integer                  :: i, j
    real(wp)                 :: dij(1:3), djk(1:3), dlk(1:3)
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2, vtmp
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl
      

    dij(1:3) = coord(1:3,aindex(1)) - coord(1:3,aindex(2))
    djk(1:3) = coord(1:3,aindex(2)) - coord(1:3,aindex(3))
    dlk(1:3) = coord(1:3,aindex(4)) - coord(1:3,aindex(3))

    aijk(1) = dij(2)*djk(3) - dij(3)*djk(2)
    aijk(2) = dij(3)*djk(1) - dij(1)*djk(3)
    aijk(3) = dij(1)*djk(2) - dij(2)*djk(1)

    ajkl(1) = dlk(2)*djk(3) - dlk(3)*djk(2)
    ajkl(2) = dlk(3)*djk(1) - dlk(1)*djk(3)
    ajkl(3) = dlk(1)*djk(2) - dlk(2)*djk(1)

    raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
    rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)

    inv_raijk2  = 1.0_wp / raijk2
    inv_rajkl2  = 1.0_wp / rajkl2

    inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)

    cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl

    rjk     = sqrt(djk(1)*djk(1) + djk(2)*djk(2) + djk(3)*djk(3))
    inv_rjk = 1.0_wp/rjk

    tmp(1)  = aijk(1)*dlk(1) + aijk(2)*dlk(2) + aijk(3)*dlk(3)
    sin_dih = tmp(1) * rjk * inv_raijkl

    dotpro_ijk = dij(1)*djk(1) + dij(2)*djk(2) + dij(3)*djk(3)
    dotpro_jkl = djk(1)*dlk(1) + djk(2)*dlk(2) + djk(3)*dlk(3)

    tmp(1) = rjk*inv_raijk2
    tmp(2) = rjk*inv_rajkl2

    tmp(3) =  dotpro_ijk*inv_raijk2*inv_rjk
    tmp(4) =  dotpro_jkl*inv_rajkl2*inv_rjk

    grad(1:3) =  tmp(1)*aijk(1:3)
    grad(4:6) = -tmp(3)*aijk(1:3) + tmp(4)*ajkl(1:3)
    grad(7:9) =                   - tmp(2)*ajkl(1:3)

    do i = 1, 3
      do j = i + 1, 3
      vtmp = - (grad(j)*dij(i) + grad(j+3)*djk(i) + grad(j+6)*dlk(i))
      v(j, i) = vtmp
      v(i, j) = vtmp
      end do
      vtmp = - (grad(i)*dij(i) + grad(i+3)*djk(i) + grad(i+6)*dlk(i))
      v(i, i) = vtmp
    end do

    return
  end subroutine calculate_dihedral

end module dihedral_libs_mod
