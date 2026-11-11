!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fitting_mod
!> @brief   module for fitting
!! @authors Yuji Sugita (YS), Norio Takase (NT), Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fitting_mod

  use fitting_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use messages_mod
  use fileio_control_mod
  use mpi_parallel_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_fit_info
    integer    :: fitting_method    = FittingMethodNO
    integer    :: fitting_atom      = 1
    logical    :: mass_weight       = .false.
    logical    :: force_no_fitting  = .false.
  end type s_fit_info

  public  :: read_ctrl_fitting_md

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_fitting_md
  !> @brief        read FITTING section in the control file
  !! @authors      NT, CK
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   fit_info : FITTING section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_fitting_md(handle, fit_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Fitting'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_fit_info),        intent(inout) :: fit_info


    ! read parameters 
    !
    fit_info%fitting_method = FittingMethodTR_ROT

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'fitting_method', &
                               fit_info%fitting_method, FittingMethodTypes)

    call read_ctrlfile_integer(handle, Section, 'fitting_atom', &
                               fit_info%fitting_atom)

    call read_ctrlfile_logical(handle, Section, 'mass_weight',    &
                               fit_info%mass_weight)

    call read_ctrlfile_logical(handle, Section, 'force_no_fitting',    &
                               fit_info%force_no_fitting)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Fitting_MD> Parameters of Fitting'
      write(MsgOut,'(A20,A10)') &
                      '  fitting method  = ', &
                      FittingMethodTypes(fit_info%fitting_method)
      write(MsgOut,'(A20,A5,I0)') &
                      '  fitting atom    = ', 'group', fit_info%fitting_atom
      if (fit_info%mass_weight) then
        write(MsgOut,'(A30)')  &
                      '  mass_weight     =        yes'
      else
        write(MsgOut,'(A30)')  &
                      '  mass_weight     =         no'
      end if
      if (fit_info%force_no_fitting) then
        write(MsgOut,'(A30)')  &
                      ' force_no_fitting =        yes'
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_fitting_md

end module fitting_mod
