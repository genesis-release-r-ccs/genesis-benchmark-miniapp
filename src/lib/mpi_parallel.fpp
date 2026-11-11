!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   mpi_parallel
!> @brief   utilities for mpi related programs                 
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module mpi_parallel_mod

#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! variables for MPI in MD
  integer, public :: mpi_comm_country  = 0
  integer, public :: mpi_comm_city     = 0
  integer, public :: nproc_world       = 1
  integer, public :: nproc_country     = 1
  integer, public :: nproc_city        = 1
  integer, public :: nprocx            = 1
  integer, public :: nprocy            = 1
  integer, public :: nprocz            = 1
  integer, public :: my_world_rank     = 0
  integer, public :: my_country_rank   = 0
  integer, public :: my_city_rank      = 0
  integer, public :: my_x_rank         = 0
  integer, public :: my_y_rank         = 0
  integer, public :: my_z_rank         = 0
  integer, public :: mpi_drain         = 50000
  integer, public :: grid_commx        = 0               
  integer, public :: grid_commy        = 0               
  integer, public :: grid_commz        = 0               
  integer, public :: my_xy_rank        = 0
  integer, public :: nprocxy           = 1
  integer, public :: grid_commxy       = 0               
  integer, public :: nthread           = 0
  integer, public :: ierror            = 0

  logical, public :: main_rank         = .false.
  logical, public :: real_calc         = .true.
  logical, public :: reciprocal_calc   = .true.

#ifdef USE_GPU
  integer, public :: my_node_local_rank = 0
#endif /* USE_GPU */

  ! subroutines

contains

end module mpi_parallel_mod
