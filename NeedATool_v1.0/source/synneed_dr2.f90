!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday,
!                          Matthias Bartelmann, Hans K. Eriksen,
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------

! OCTOBER 2008 - D.Pietrobon
!
! When using this module please refere to
!
! "Spherical needlets for cosmic microwave background data analysis"
! Marinucci et al.
! 2008MNRAS.383..539M - arxiv:0707.0844
!
!  Thanks to P.Cabella for suggestions on the implementation.
      program syneed

!  This program computes the needlets of a given temperature map, storing bl2 coefficients in file.
!
!  If required, direct space needlet is built.

! Healpix modules 
        USE healpix_types
! ---
        use extension, only: nArguments, getArgument 
!        use paramfile_io
!!$        use misc_utils,only: assert, fatal_error!, strlowcase

! --- my modules
        USE needlets_mod
!        USE need_ioput

        IMPLICIT NONE

! insert parsing module and necessary declaration variables

        INTEGER(i4b)                :: n_args, i
        CHARACTER(LEN=FILENAMELEN)  :: arg
!        logical(lgt)                :: do_double
        CHARACTER(len=FILENAMELEN)  :: paramfile
        CHARACTER(LEN=10)           :: lcode

        LOGICAL :: interactive = .TRUE.

!        TYPE(input_keys) :: input_v
!        TYPE(io_vars) :: out_v

! count arguments, should be 0, 1 or 2
!        n_args = nArguments()

        WRITE(*,*) " "
        WRITE(*,*) "  *****************************"
        WRITE(*,*) "   Welcome to "//TRIM(SYNCODE)//" "//TRIM(MYVERSION)
        WRITE(*,*) "  *****************************"
        WRITE(*,*) " "

        lcode = SYNCODE

!!$        n_args = COMMAND_ARGUMENT_COUNT()
!!$        call assert(n_args <= 2,' Usage: '//TRIM(SYNCODE)//' [-s|--single|-d|--double] [parameter_file_name]')

        paramfile = ''
        if (nArguments() == 1) call getArgument(1, paramfile) 
        call parsing_hpx(paramfile, lcode)
!!$        SELECT CASE(n_args)
!!$           CASE (0)
!!$              interactive = .TRUE.
!!$              CALL get_input_pars_interactive(out_v,lcode)
!!$           CASE (1)
!!$              interactive = .FALSE.
!!$              CALL getArgument(1, arg)
!!$              arg = TRIM(ADJUSTL(arg))
!!$              paramfile = arg
!!$              CALL parsing(paramfile)
!!$              CALL keys2values(out_v,lcode)
!!$           CASE DEFAULT
!!$              CALL fatal_error(SYNCODE//': Invalid argument: '//trim(arg))
!!$        END SELECT


!        lcode = trim(strlowcase(CODE))
!        call assert(n_args <= 2,' Usage: '//trim(lcode)//' [-s|--single|-d|--double] [parameter_file_name]')

!        parafile = ''
!        do_double = .false.

! parse arguments
!        do i=1, n_args
!           call getArgument(i, arg)
!           arg = trim(adjustl(arg))
!           if (arg(1:1) == '-') then
!              select case (trim(arg))
!              case ('-d', '--double')
!                 do_double = .true.
!              case ('-s', '--single')
!                 do_double = .false.
!              case default
!                 call fatal_error(CODE//': Invalid argument: '//trim(arg))
!              end select
!           else
!              parafile = arg
!           endif
!        enddo

! start calculations
!        if (do_double) then
!           call syn_sub_d(parafile)
!        else
!           call syn_sub_s(parafile)
!        endif

!        print*, " input finished!"

!        STOP

        print*, ""
        print*, " Starting computation..."
        print*, ""

        CALL synneed_sub

 !!$        CALL write_params_used(lcode)

        STOP

      END PROGRAM syneed

! ************************************************************************************************
