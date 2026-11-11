!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_output_str_md
!> @brief   structure of output information
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_output_str_mod

  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_output
    integer                :: dcdunit
    integer                :: dcdvelunit

    character(MaxFilename) :: dcdfile
    character(MaxFilename) :: dcdvelfile
    character(MaxFilename) :: rstfile
    character(MaxFilename) :: pdbfile

    logical                :: out_dcd_header    = .true.
    logical                :: out_dcdvel_header = .true.

    logical                :: dcdout     = .false.
    logical                :: dcdvelout  = .false.
    logical                :: rstout     = .false.
    logical                :: pdbout     = .false.
                    


    logical                :: verbose

  end type s_output

end module sp_output_str_mod
