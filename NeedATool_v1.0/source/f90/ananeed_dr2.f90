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
      PROGRAM Ananeed

!  This program computes the needlets of a given temperature map, storing bl2 coefficients in file.
!
!  If required, direct space needlet is built.

! Healpix modules 
        USE healpix_types
! ---
        use extension, only: nArguments, getArgument
!        use misc_utils,only: assert, fatal_error!, strlowcase

! --- my modules
        USE needlets_mod
!        USE need_ioput

        IMPLICIT NONE

        CHARACTER(len=FILENAMELEN)  :: paramfile
        CHARACTER(LEN=10)           :: lcode

        WRITE(*,*) " "
        WRITE(*,*) "  ***********************************"
        WRITE(*,*) "   Welcome to "//TRIM(ANACODE)//" "//TRIM(MYVERSION)
        WRITE(*,*) "  ***********************************"
        WRITE(*,*) " "

        lcode = ANACODE

         paramfile = ''
        if (nArguments() == 1) call getArgument(1, paramfile) 
        call parsing_hpx(paramfile, lcode)

        print*, ""
        print*, " Starting computation..."
        print*, ""

        CALL ananeed_sub

        print*, 'Done. Good Bye.'

      END PROGRAM ananeed

! ************************************************************************************************
