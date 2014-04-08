! June 2010 - making the input more user friendly, using getsize_fits healpix routine to get arguments of the input map.
!
! Modification to take into account scaling of nside with last in the filter, up to lmax;
!  currently nside is set computing the actual degree of freedom present in the map.
!  We wanna change that to match the Needlet ILC procedure, where the map are computed at nside
!  such that lmax in the bin is le 2*nside. Actually we have the boost parameter so that we'll set
!  lmac_bin le nside_boost*nside
!
! Correction for pwf and w8r added
!
! factor computed by the code: set to 1./sqrt(npix) (=> needlets are normalised by multiplying to norm)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program calculates ths "psi" function necessary to build
! wavelets basis
!
!  B_jk = 4*Pi/B**j * sum_l b(l/B**j) sum_m a_lm*Y_lm(Omega_k)
!
! Note that this code computes bsqr, i.e. b(l/B**j)**2; when building
!   needlets, beta_jk, remember the normalization factor (actually pixel counting)!  
!
!  Omega_k = ( theta_k, phi_k ), polar angles
!
! Transormation is done using Healpix routine: I can give a_lm * psi_l
! as input of synfast: in this way I have one map for each j
!
!
!1)        f(t) = exp( -1/(1-t**2) ), -1 <= t <= 1
!             = 0                   otherwise
!
!2)        psi(u) = int_-1^u f(t)dt / int_-1^1 f(t)dt
!
!3)        phi(t) = 1                            , 0 < t < 1/B
!                 = psi( 1 - 2B/(B-1)(t - 1/B) ) , 1/B <= t <= 1
!                 = 0                            , t > 1
!4)        bsqr(xi) = phi(xi/B) - phi(xi)
!

        MODULE needlet_variables
          
          USE healpix_types

          implicit none

          REAL(dp)     :: B=2, oneob=0.5, deltaj=1.
          INTEGER(i4b) :: n_resol, lmax=512, jmax

        END MODULE needlet_variables

! **********************************************************************

        MODULE needlets_mod

          USE healpix_types
          USE alm_tools
          USE fitstools
          USE utilities
          USE pix_tools
          USE udgrade_nr, ONLY: udgrade_ring
          USE rngmod, ONLY: rand_init, planck_rng

          USE needlet_variables

          IMPLICIT NONE

          PRIVATE

          PUBLIC synneed_sub, ananeed_sub, &
	parsing_hpx, SYNCODE, ANACODE, MYVERSION, &
	bl2, io_B, io_nresol, set_input, &
	set_needlet_environment, which_l, bj_of_l, nside, npix, &
	lst_l, jflag, nside_boost, set_nside, build_needlet, dof, i_ns

          CHARACTER(LEN=10), PARAMETER               :: SYNCODE = "SYNNEED"
          CHARACTER(LEN=10), PARAMETER               :: ANACODE = "ANANEED"
          CHARACTER(LEN=10), PARAMETER               :: MYVERSION = "1.0 - Jun 2010"
          INTEGER(i4b)                               :: nside, npix, mapnside, mapnpix, nside_boost=2, in_order
          INTEGER(i4b), PARAMETER                    :: i_lw=1, i_ave=2, i_up=3, i_dl=4, i_dof=5, i_ns=6 
          INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE  :: dof
          REAL(dp)                                   :: need_norm
          REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: bl2
          REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: need_map, need_mask
          COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: alm

          REAL(dp), ALLOCATABLE, DIMENSION(:)    :: masked_map, lst_l
          
          LOGICAL :: mask_applied = .FALSE.
          LOGICAL :: beam_applied = .FALSE.
          LOGICAL :: do_needlets  = .TRUE.
          LOGICAL, DIMENSION(:), ALLOCATABLE :: jflag

          CHARACTER(LEN=filenamelen)   :: healpix_dir
          CHARACTER(LEN=filenamelen)   :: clsfile
          CHARACTER(LEN=filenamelen)   :: almsfile
          CHARACTER(LEN=filenamelen)   :: mapfile
          CHARACTER(LEN=filenamelen)   :: maskfile
          CHARACTER(LEN=filenamelen)   :: beamfile
          CHARACTER(LEN=filenamelen)   :: need_maskfile
          CHARACTER(LEN=filenamelen)   :: bl2_root
          CHARACTER(LEN=filenamelen)   :: bl2_infile
          CHARACTER(LEN=filenamelen)   :: need_root
          CHARACTER(LEN=filenamelen)   :: need_resol_scheme = 'constant_nside'

          CHARACTER(LEN=5) :: io_B
          CHARACTER(LEN=5) :: io_dj
          CHARACTER(LEN=3) :: io_nresol
          CHARACTER(LEN=4) :: io_nside
 
        CONTAINS

! **********************************************************************

          subroutine set_input

            implicit none

            WRITE(*,*) " Type Needlets parameter B (1.8):"
            READ(*,*) B

            WRITE(*,*) " Type the lmax at which you'd like to analyze maps (500)"
            READ(*,*) lmax

            WRITE(*,*) " Type the nside_boost (2)"
            READ(*,*) nside_boost

            return
          end subroutine set_input

! **********************************************************************

          subroutine set_needlet_environment(voice)

            implicit none

            character(len=10) :: code
            
            logical, optional, intent(inout) :: voice

           if ( .not. present(voice) ) voice = .true.

            oneob = 1./B
            write(io_B,'(f4.2)') B

! ------ deltaj .ne. 1. is not fully implemented for the reconstruction
            deltaj = 1.
            write(io_dj,'(1f5.2)') deltaj
! ---
            jmax = 1 + INT( LOG( REAL(lmax) ) / LOG(B) / deltaj )
            n_resol = jmax
            WRITE(*,*) " n_resol = ", n_resol
            WRITE(io_nresol,"(i3.3)") n_resol

            ALLOCATE( bl2(n_resol, 0:lmax) )
            bl2 = 0.d0
            ALLOCATE( jflag(n_resol) )
            jflag = .FALSE.
            ALLOCATE( dof(n_resol, 6) )
            dof = 0.

            CALL bj_of_l

            if (voice) CALL which_l

            CALL get_dof(voice)

            code = SYNCODE
            need_resol_scheme = 'constant_nside'
            if (code .eq. SYNCODE) then

               if ( (trim(adjustl(need_resol_scheme)) .eq. 'constant_nside') .or. (trim(adjustl(need_resol_scheme)) .eq. 'scaling_boost') ) nside = set_nside(lmax, voice=voice)

               if ( trim(adjustl(need_resol_scheme)) .eq. 'scaling_dof') nside = maxval( dof(:,i_ns) )
               npix = 12. * nside**2
            endif

            return



          end subroutine set_needlet_environment

! **********************************************************************

          subroutine get_dof(voice)

            implicit none

            integer(i4b) :: lw, up, ave, dl, d, ns, hns, j, i

            logical, optional, intent(inout) :: voice

            print*, " ***************************"
            print*, " *    Subroutine get_dof   *"
            print*, " ***************************"


            if (.not. present(voice) ) voice = .true.

            do j=1,jmax
               lw = int( B**(j-1) )
               up = int( B**(j+1) )
               ave = int( ( up+lw )/2. )
               dl = int( B**(j+1)-B**(j-1) )
               d = int( dl * (ave+3./2.) )
               ns = int( sqrt( d / 12. ) )
               hns = 2
               do i=1, 10
                  if (ns .gt. 2**i .and. ns .le. 2**(i+1)) hns = int( 2**(i+1) )
               enddo
               dof(j,:) = (/lw, ave, up, dl, d, hns/)
               if (j .eq. jmax) dof(j,i_ns) = dof(j-1,i_ns)
            enddo

            if (voice) then
               print*, " Computation of actual d.o.f. for each j:"
               print*, ""
               write(*,'(7a15)'), "j", "lw", "ave", "up", "dl", "dof", "hns"
               print*, ""
               do j=1,jmax
                  write(*,'(7i15)'), j, dof(j,:)
               enddo
            endif

            return

          end subroutine get_dof

! **********************************************************************

          SUBROUTINE write_bl2(in_root)

            USE head_fits, only: add_card, write_minimal_header
            USE fitstools, only: write_asctab

            IMPLICIT NONE

            CHARACTER(LEN=filenamelen), INTENT(IN) :: in_root

            INTEGER(I4B) :: j, il
            CHARACTER(LEN=80), DIMENSION(80) :: bl2_header
            CHARACTER(LEN=filenamelen) :: tempfile
            CHARACTER(LEN=3) :: iores

            bl2_header = ''
            call write_minimal_header(bl2_header, 'CL', nlmax=lmax)
 !!$            CALL add_card(bl2_header)
            CALL add_card(bl2_header, 'filecont', 'MYCODE', SYNCODE)
            CALL add_card(bl2_header, 'filecont', 'version', MYVERSION)
            CALL add_card(bl2_header, 'filecont', 'B', io_b)
            CALL add_card(bl2_header, 'filecont', 'NUM js', io_nresol)
            CALL add_card(bl2_header, 'filecont', 'bl^2', 'filter functions')

 !!$            print*, bl2_header(1:5)

            print*, " ***************************"
            print*, " *   Subroutine write_bl2  *"
            print*, " ***************************"

            tempfile = TRIM(in_root)//'_B'//TRIM(io_B)//'_Nj'//TRIM(io_nresol)//'.fits'
            print*, ' ...writing on file_root: ', TRIM(tempfile)

            call write_asctab(transpose(bl2), lmax, n_resol, bl2_header, 80, tempfile)   

            if (.false.) then 
               DO j = 1, n_resol

                  IF (j .LT. 10) WRITE(iores,'(1i1)') j
                  IF ( (j .LT. 100).AND.(j.GT.9) ) WRITE(iores,'(1i2)') j
                  IF ( (j .LT. 1000).AND.(j.GT.99) ) WRITE(iores,'(1i3)') j
!               print*, iores
                  tempfile = TRIM(in_root)//'_B'//TRIM(io_B)//'_Nj'//TRIM(io_nresol)//'_'//TRIM(iores)//'.dat'
!!$               print*, ' ...writing on file: ', TRIM(tempfile)
                  IF (tempfile(1:1) .EQ. "!") OPEN(unit=23,file=tempfile(2:),form='formatted',status='unknown',action='write')
                  IF (tempfile(1:1) .NE. "!") OPEN(unit=23,file=tempfile,form='formatted',status='new',action='write')
                  DO il = 0, lmax
                     WRITE(23,'(1i10,f15.6)') il, bl2(j,il)
                  ENDDO
                  CLOSE(23)
               ENDDO

            endif

            RETURN

          END SUBROUTINE write_bl2

! **********************************************************************

          subroutine parsing_hpx(paramfile, code)

            use paramfile_io
            use extension

            implicit none

            type(paramfile_handle)         :: handle
            character(len=256), intent(in) :: paramfile
            character(len=10), intent(in)  :: code


            handle = parse_init(paramfile)

            healpix_dir       = parse_string(handle, 'healpix_dir', default='/Users/dpietrob/work/cmbcodes/Healpix_2.13a', descr='Healpix package path')
            lmax              = parse_int(handle, 'l_max', default=500, descr='Maximum ell at which the analysis is carried out')
            B                 = parse_double(handle, 'B', default=2._dp, descr='Filter range parameter')
!!$            deltaj            = parse_double(handle, 'delta_j', default=1._dp)
!!$            write(io_dj,'(1f5.2)') deltaj
            nside_boost       = parse_int(handle, 'nside_boost', default=2, descr='Sets Nside of the needlet coefficient map (lmax < nside_boost*nside)', vmin=1, vmax=4)
!!$            need_resol_scheme = parse_string(handle, 'need_resol_scheme', default='constant_nside')
            bl2_infile        = parse_string(handle, 'bl2_infile', default="''", descr='Fits file containing square window functions') 
! --- synneed parameters
            if (trim(adjustl(code)) == SYNCODE) then
               do_needlets       = parse_lgt(handle, 'compute_needlets', default=.true.)
               mapfile           = parse_string(handle, 'mapfile', default='input/lcdm_map_lmax500.fits', descr='fits file containing the map to be decomposed onto needlets')
               mapnside          = parse_int(handle, 'mapnside', default=256, descr='Nside of the input map')
               maskfile          = parse_string(handle, 'maskfile', default='')
               if (trim(adjustl(maskfile)) .ne. '') mask_applied = .true.
               beamfile          = parse_string(handle, 'beamfile', default='')
               if (trim(adjustl(beamfile)) .ne. '') beam_applied = .true.
               bl2_root          = parse_string(handle, 'bl2_root', default='!test_bl2', descr='Tag for the bl2 filter files')
               need_root         = parse_string(handle, 'need_root', default='!test_needlet_coefficients', descr='Tag for the needlet coefficient file')
            endif

! --- ananeed parameters
            if (trim(adjustl(code)) == ANACODE) then
               mapfile           = parse_string(handle, 'mapfile', default='test_needlet_coefficients_2.00_Nj009.fits', descr='fits file containing needlet coefficients')
               mapnside          = parse_int(handle, 'mapnside', default=256, descr='Nside of the input needlet coefficients')
               need_maskfile     = parse_string(handle, 'need_maskfile', default="''")
               need_root         = parse_string(handle, 'need_root', default='!test_recmap', descr='Tag for the reconstructed map file')
            endif

            call parse_summarize(handle, code=code)

            call parse_check_unused(handle, code=code)

            call parse_finish(handle)

            return

          end subroutine parsing_hpx

! **********************************************************************

          FUNCTION set_nside(l_max, voice)

            IMPLICIT NONE

            INTEGER(i4b), INTENT(IN)  :: l_max
            INTEGER(i4b)              :: set_nside

            INTEGER(i4b) :: i
            
            LOGICAL, OPTIONAL :: voice

            print*, " ***************************"
            print*, " *   Subroutine set_nside  *"
            print*, " ***************************"


            if ( .not. present(voice) ) voice= .true.

            if (voice) WRITE(*,*) " Setting nside: mode '", trim(adjustl(need_resol_scheme)), "'"

            i = 0

            DO WHILE (l_max > nside_boost*2**i )
               i = i+1               
            ENDDO

            set_nside = 2**i

            if (voice) WRITE(*,'(A14,I8)') " nside     = ", set_nside

            RETURN

          END FUNCTION set_nside

! **********************************************************************

          subroutine synneed_sub

            USE head_fits, only: add_card, write_minimal_header

            IMPLICIT NONE

            REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: mask, noise, w5map, w8r
            REAL(sp), DIMENSION(:,:), ALLOCATABLE    :: ioneedmap
            LOGICAL                               :: anynull = .FALSE.
            REAL(dp)                              :: nullval = -1.63750e+30
! ----------------
! Variables required by map2alm
            INTEGER(i4b), PARAMETER               :: p = 1, nlheader=80
            REAL(dp), DIMENSION(0:lmax,1:p)       :: pwf
            REAL(dp), DIMENSION(1:2)              :: zbounds
            REAL(dp), DIMENSION(:), ALLOCATABLE   :: multipoles_fit
!            REAL(dp), DIMENSION(1:nside_boost*mapnside,1:p) :: w8r
! ----------------
! Variables required by alm2cl
            REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: cl
! ----------------
! Variables necessary to call "create_alms" - "write_bintab"
!              INTEGER(I4B), PARAMETER          :: polar = 0
!              TYPE(planck_rng)                 :: rng_handle
!              REAL(sp), PARAMETER              :: fwhm_arcmin = 0.
            CHARACTER(LEN=80), DIMENSION(80)  :: header_PS
            CHARACTER(LEN=80), DIMENSION(120) :: need_header
            CHARACTER(LEN=80), DIMENSION(80)  :: cl_header
!              CHARACTER(LEN=100)               :: hd_comment
!        CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL :: window_file
!              CHARACTER(LEN=80), DIMENSION(1)  :: units !,  INTENT(OUT),    OPTIONAL
!        CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL :: beam_file
! ----------------
! Variables necessary to call "udgrade_ring"
!              REAL(sp), PARAMETER :: bad_data = 0.
              LOGICAL :: speak = .true.
! ----------------        
!              REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: sim_need!, sim_need2
!              REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: tot_need, tot_need2, var
!              REAL(sp), DIMENSION(0:npix-1)         :: map

            INTEGER(i4b) :: imap, i_resol, j, jpr, counter, iw, ist, iter, l, ordering, multipole_remov_deg, nmap
            INTEGER(i4b), PARAMETER :: itemp = 1, lmin = 2
            INTEGER(i4b) :: mmax, nbins, nmaps, in_nside
            INTEGER(i8b) :: in_npix

            REAL(sp)     :: t1, t2, lower, higher

            CHARACTER(LEN=1) :: choice
            CHARACTER(LEN=2) :: io_j
            CHARACTER(LEN=4) :: io_hns

            CHARACTER(LEN=filenamelen) :: tempfile

            CALL CPU_TIME(t1)

            WRITE(*,*) " ...computing needlets coefficients:"

            CALL set_needlet_environment( speak )

            CALL write_bl2(bl2_root)

            need_required:IF (do_needlets) THEN

               iter = 2
               mmax = lmax
                 
               ALLOCATE(alm(1:p,0:lmax,0:mmax))

               ordering = 1 ! 1:RING, 2:NEST
               multipole_remov_deg = 2 ! 0:none, 1:monopole, 2:monopole and dipole
               allocate( multipoles_fit(0:multipole_remov_deg*multipole_remov_deg-1) )
               multipoles_fit = 0.
               zbounds = 0.

! ------ Reading map
               IF (mapfile /= "''") THEN
                  WRITE(*,*) " ...reading input map: "//TRIM(mapfile)
                  in_npix = getsize_fits(mapfile, ordering=in_order, nside=in_nside, nmaps=nmap)
                  if (in_nside .ne. mapnside) then
                     print*, " WARNING - The nside found in the header does not match your entry!"
                     print*, " The header value is assumed the correct one."
                     mapnside = in_nside
                  endif
                  mapnpix = in_npix
                  if (in_order .eq. 0) print*, " WARNING - ordering not found; RING assumed!"

                  WRITE(*,*) "    map nside: ", mapnside
                  WRITE(*,*) "    map npix:  ", mapnpix
                  WRITE(*,*) "    ordering:  ", in_order
                  WRITE(*,*) "    map nmap:  ", nmap
                  print*, ""
                  WRITE(*,*) "    mask:      ", mask_applied
                  WRITE(*,*) "    beam:      ", beam_applied

                  allocate( w5map(0:mapnpix-1,1:nmap) )
                  w5map = 0._dp
!               allocate( noise(0:mapnpix-1,1:nmap) )
                  allocate( mask(0:mapnpix-1,1:nmap) )
                  mask = 0._dp
                  CALL read_bintab(mapfile, w5map, mapnpix, nmap, nullval, anynull)
                  if (in_order .eq. 2) then
                     print*, " Input map in nest order: reordering..."
                     call convert_nest2ring(mapnside, w5map)
                     print*," Map converted to RING order."
                  endif
               ENDIF

! ------ Reading mask
!                 print*, " w5map size: ",size(mask)
                 mask = 1._dp
                 IF (mask_applied) THEN
                    WRITE(*,*) " ...reading " //TRIM(maskfile)
                    in_npix = getsize_fits(maskfile, ordering=in_order, nside=in_nside, nmaps=nmap)
                    CALL read_bintab(maskfile, mask, mapnpix, nmap, nullval, anynull)
                    if (in_nside .ne. mapnside) stop " ERROR - Mismatch between map and mask nside!" 
                    if (in_order .eq. 2) then
                       print*, " WARNING - ordering of the provided mask is NEST: converting into RING"
                       call convert_nest2ring(mapnside,mask)
                    endif
                    print*, "read mask"
                 ENDIF
! TEST
! --- removing residual monopole and dipole from all needlets
                 CALL remove_dipole(mapnside, w5map(:,1), ordering, multipole_remov_deg, multipoles_fit, zbounds, mask=mask(:,1))

! ------ Extracting alms
! including possibility to choose all anafast parameters

                 WRITE(*,*) " ...extracting alm from masked map:"

! ------ Correcting for the w8ring
                 WRITE(io_nside,'(1i4.4)') mapnside
                 allocate( w8r(1:nside_boost*mapnside,1:p) )
                 w8r = 0.
                 tempfile = TRIM(adjustl(healpix_dir))//"/data/weight_ring_n0"//io_nside//".fits"
                 print*, " Reading w8r correction..."
                 CALL input_map(tempfile, w8r, 2*mapnside, p)

!                 print*, " !!! WARNING - w8r removed in synneed !!!"
                 w8r = 1._dp + w8r
!                 CALL fits2cl(tempfile, w8r, 2*mapnside, p, cl_header)
!                 print*, w8r(1:10,1)

                 IF (.NOT. mask_applied) CALL map2alm_iterative(mapnside, lmax, mmax, iter, w5map, alm, w8ring=w8r)
                 IF (mask_applied) CALL map2alm_iterative(mapnside, lmax, mmax, iter, w5map, alm, w8ring=w8r, mask=mask)

! ------ Correcting for the pwf
! more precise the reconstruction if not used
!!$                 WRITE(io_nside,'(1i4.4)') mapnside
!!$                 tempfile = TRIM(healpix_dir)//"/data/pixel_window_n"//io_nside//".fits"
!!$                 print*, " Reading pwf correction..."
!!$                 CALL fits2cl(tempfile,pwf,lmax,p,cl_header)
!!$                 DO l = 0, lmax 
!!$                    alm(:,l,0:mmax) = alm(:,l,0:mmax) !/ pwf(l,1)
!!$                 ENDDO

! --- beam reading
                 if (beam_applied) then 
                    tempfile = TRIM(beamfile)
                    print*, " Reading beam_file..."
                    pwf = 1.
                    CALL fits2cl(tempfile,pwf,lmax,p,cl_header)
                    DO l = 0, lmax 
                       alm(:,l,0:mmax) = alm(:,l,0:mmax) * pwf(l,1)
                    ENDDO
                 endif

! ------ Building needlets
                 WRITE(*,*) " ...building needlets:"
                 WRITE(*,*) "    needlets nside: ", nside
                 WRITE(*,*) "    needlets npix:  ", npix

                 ALLOCATE(need_map(0:npix-1,0:n_resol-1))
                 need_map = hpx_sbadval
                 ALLOCATE(need_mask(0:npix-1,0:n_resol-1))
                 need_mask = hpx_sbadval

                 WRITE(*,*) " ...building needlets and masks for each resolution:"
                 DO i_resol=0,n_resol-1

                    IF ( jflag(i_resol+1) ) then
                       if (need_resol_scheme .eq. 'constant_nside') CALL Build_Needlet(i_resol+1, alm, p, lmax, mmax, nside, need_map(:,i_resol))

                       if (need_resol_scheme .eq. 'scaling_dof') then
                          npix = 12.*dof(i_resol+1,i_ns)**2
                          call Build_Needlet(i_resol+1, alm, p, lmax, mmax, dof(i_resol+1,i_ns), need_map(0:npix-1,i_resol)) 
                       endif

                       if (need_resol_scheme .eq. 'scaling_boost') then
                          npix = 12.*set_nside( dof(i_resol+1,i_up) )**2
                          CALL Build_Needlet(i_resol+1, alm, p, lmax, mmax, set_nside(dof(i_resol+1,i_up)), need_map(0:npix-1,i_resol))
                       endif
                    endif
! to be modified if the masking procedure depends on n_resol
! In principle should depend on n_resol to take into account different needlets
!   responses to the mask due to higher or lower localization
!!                 need_mask(:,i_resol) = mask(:,1)

                    need_mask(0:npix-1,i_resol) = 1.!mask(:,1)

                 ENDDO
        
 ! ------ masking needlets
 !!$                 WRITE(*,*) " ...masking needlets:"
 !!$                 need_map = need_map * need_mask / need_norm
                 print*, " Normalization included in the needlet construction"
 !!$                 print*, " normalization: ", need_norm
 !!$                 need_map = need_map * need_norm

                 need_header = ''
                 call write_minimal_header(need_header, 'MAP', nside=nside, ordering='Ring', coordsys='G', creator=SYNCODE, version=MYVERSION, nlmax=lmax)
 !!$                 CALL add_card(need_header)
 !!$                 CALL add_card(need_header, 'MYCODE:', SYNCODE)
 !!$                 CALL add_card(need_header, 'version', MYVERSION)
                 CALL add_card(need_header, 'filecont', 'B', io_B)
                 CALL add_card(need_header, 'filecont', 'NUM js', io_nresol)
 !!$                 CALL add_card(need_header, 'nside', nside)
 !!$                 CALL add_card(need_header, 'ordering', 'ring')
                 CALL add_card(need_header, 'filecont', 'beta_jk', 'Needlets coefficients')
                 CALL add_card(need_header, 'filecont', 'nr_schem', trim(adjustl(need_resol_scheme)))
                 do j=1,jmax
                    write(io_j,'(1i2.2)') j
                    write(io_hns,'(1i4.4)') dof(j,i_ns)
!                    print*, io_j
                    CALL add_card(need_header,  'filecont', 'ns '//io_j, io_hns)
                 enddo

                 tempfile = TRIM(need_root)//'_'//TRIM(io_B)//'_Nj'//TRIM(io_nresol)//'.fits'

                 WRITE(*,*) " ...writing needlets on file: ", TRIM(tempfile)

                 CALL output_map(need_map, need_header, TRIM(tempfile))

                 DEALLOCATE( need_mask, need_map, alm, jflag, bl2 )

              ENDIF need_required

              CALL CPU_TIME(t2)

              WRITE(*,*) "elapsed time: ", INT((t2-t1)/60),"min", MOD((t2-t1),60.),"sec"

! ------ deallocating variables
!!$              DEALLOCATE( need_mask, need_map, masked_map, cl, alm )

              RETURN

            END SUBROUTINE synneed_sub

! **********************************************************************

            SUBROUTINE ananeed_sub

              USE head_fits, only: add_card, write_minimal_header
              USE fitstools, only: fits2cl

              IMPLICIT NONE

              CHARACTER(LEN=filenamelen)                 :: tempfile
!!$              REAL(sp), DIMENSION(:), ALLOCATABLE        :: temp_mask!temp_map, 
              INTEGER(i4b)                               :: imap, mmax, l, nl, nm, iter, ordering, multipole_remov_deg, nmap, jnside, jnpix
              INTEGER(i4b), PARAMETER                    :: p = 1!, nlheader=80
              integer(i8b) :: in_npix

              REAL(dp), DIMENSION(1:2)                   :: zbounds
              REAL(dp), DIMENSION(0:lmax,1:p)            :: pwf
!              REAL(dp), DIMENSION(1:nside_boost*nside,1:p)         :: w8r
              REAL(dp), DIMENSION(:), ALLOCATABLE        :: multipoles_fit
              REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: mask, temp_map
!              REAL(dp), DIMENSION(0:npix-1,1:p)          :: temp_map
!              REAL(dp), DIMENSION(0:mapnpix-1,1:nmap)    :: w5map
              REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: w5map, w8r
              REAL(sp) :: t1, t2

              COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:) :: temp_alm

              CHARACTER(LEN=80), DIMENSION(120)          :: need_header

              logical :: speak = .true.

              CALL CPU_TIME(t1)

              iter = 2 ! number of iteration for map2alm_iterative

              WRITE(*,*) " ...computing bl2:"

              CALL set_needlet_environment( speak)
             
              mmax = lmax
                 
              nl = lmax
              nm = nl
              ALLOCATE( alm(1:p, 0:nl, 0:nm) )
              alm = 0.
              ALLOCATE( temp_alm(1:p, 0:lmax, 0:mmax) )
              temp_alm = 0.

              tempfile = mapfile
              WRITE(*,*) " ...reading needlet coefficients from file: ", TRIM(tempfile)
              print*, " retriving information from: ", trim(adjustl(mapfile))
              in_npix = getsize_fits(mapfile, ordering=in_order, nside=nside, nmaps=nmap)
              npix = in_npix
! --- removed from set_needlet_environment"
              WRITE(*,*) "    needlet nside:  ", nside
              WRITE(*,*) "    needlet npix :  ", npix
              WRITE(*,*) "    needlet nmap :  ", nmap
              WRITE(*,*) "    needlet order:  ", in_order
              WRITE(*,*) "need_resol_scheme:  ", trim(adjustl(need_resol_scheme))
              ALLOCATE( temp_map(0:npix-1, 1:p) )
              temp_map = 0.
              if (nmap .ne. n_resol) then
                 print*, " WARNING - Mismatch between number of resolutions and maps found in the fits file:" 
                 print*, n_resol, " VS ", nmap
                 print*, " n_resol set equal to nmap. Recovered map may be not accurate."
                 n_resol = nmap
              endif
              WRITE(*,*) "    mask:      ", mask_applied
              print*, ""
              
              ALLOCATE(need_map(0:npix-1,0:n_resol-1))
              need_map = 0.
              ALLOCATE( mask(0:npix-1,0:n_resol-1) )
              mask = 0.
              CALL input_map(tempfile, need_map, npix, n_resol)
              if (in_order .eq. 0) print*, " Undefined ordering: assumed RING."
              if (in_order .eq. 2) then
                 print*, " WARNING - Order of the input needlet coefficients is NEST."
                 print*, " Converted to RING."
                 call convert_nest2ring(nside, need_map)
              endif

              WRITE(*,*) " ...analysing needlet coefficients"


              IF (mask_applied) THEN
!!$                 ALLOCATE( temp_mask(0:npix-1) )
                 tempfile = TRIM(need_maskfile)
                 WRITE(*,*) "...reading needlet mask from file: ", TRIM(tempfile)
                 in_npix = getsize_fits(need_maskfile, ordering=in_order, nmaps=nmap)
                 print*, " You provided ", nmap, " masks: be sure this is consistent with the input needlet coefficient maps."
                 mask = 1.
                 if (in_npix .ne. npix) print*, " WARNING - The mask provided does not have the same resolution as needlet coefficients: skipped!" 
                 if (in_npix .eq. npix) then
                    CALL input_map(tempfile, mask, npix, nmap)
                    if (in_order .eq. 0) print*, " Undefined ordering: assumed RING."
                    if (in_order .eq. 2) then
                       print*, " WARNING - Order of the input needlet coefficients is NEST."
                       print*, " Converted to RING."
                       call convert_nest2ring(nside, mask)
                    endif
                    if ( ( (nmap .eq. 1) .or. (nmap .ne. n_resol) ) .and. (need_resol_scheme .ne. 'constant_nside') ) then
                       print*, " One mask provided but needlet resolution scales with j. Setting mask = 1."
                       mask = 1.
                    endif
                 endif
!!$                 temp_mask = mask(:,1)
              ENDIF

! The following has to be coherent with synneed_sub
!!$              zbounds = 0._dp
!!$              w8r = 1.

!!$              WRITE(io_nside,'(1i4.4)') mapnside

! --- This assumes the same resolution for each needlet set.
! move inside the loop if that is not satisfied
!!$              WRITE(io_nside,'(1i4.4)') nside
!!$              tempfile = TRIM(healpix_dir)//"/data/pixel_window_n"//io_nside//".fits"
!!$              CALL fits2cl(tempfile, pwf, lmax, p, need_header)
!!$!              print*, pwf(0:10,1)
!!$
!!$              tempfile = TRIM(healpix_dir)//"/data/weight_ring_n0"//io_nside//".fits"
!!$              allocate( w8r(1:nside_boost*nside, 1:p) )
!!$              w8r = 0.
!!$              CALL input_map(tempfile, w8r, 2*nside, p)
!!$              w8r = 1._dp + w8r
!!$!              CALL fits2cl(tempfile, w8r, 2*nside, p, need_header)
!!$!              print*, w8r(1:10,1)

              ordering = 1 ! 1:RING, 2:NEST
              multipole_remov_deg = 2 ! 0:none, 1:monopole, 2:monopole and dipole
              allocate( multipoles_fit(0:multipole_remov_deg*multipole_remov_deg-1) )
              multipoles_fit = 0.
              zbounds = 0.
!print*, size(temp_map)
              DO imap = 0, n_resol-1
                 IF ( jflag(imap+1) ) THEN 
                    print*, " j ", imap+1
                    temp_map = 0.
                    temp_alm = 0.
                    if (need_resol_scheme .eq. 'constant_nside') jnside = nside

                    if (need_resol_scheme .eq. 'scaling_dof')    jnside = dof(imap+1,i_ns)

                    if (need_resol_scheme .eq. 'scaling_boost')  jnside = set_nside( dof(imap+1,i_up) )
                    jnpix = 12 * jnside**2
                    need_norm = 1./sqrt(real(jnpix))
!print*, jnside, jnpix, need_norm
                    temp_map(0:jnpix-1,1) = need_map(0:jnpix-1,imap)
!print*, "here"
                    WRITE(io_nside,'(1i4.4)') jnside
                    tempfile = TRIM(healpix_dir)//"/data/pixel_window_n"//io_nside//".fits"
                    CALL fits2cl(tempfile, pwf, lmax, p, need_header)
!              print*, pwf(0:10,1)

                    tempfile = TRIM(healpix_dir)//"/data/weight_ring_n0"//io_nside//".fits"
                    allocate( w8r(1:2*jnside, 1:p) )
                    w8r = 0.
                    CALL input_map(tempfile, w8r, 2*jnside, p)
                    w8r = 1._dp + w8r

                    IF (.NOT. mask_applied) CALL remove_dipole(jnside, temp_map(0:jnpix-1,1), ordering, multipole_remov_deg, multipoles_fit, zbounds)
                    IF (mask_applied) CALL remove_dipole(jnside, temp_map(0:jnpix-1,1), ordering, multipole_remov_deg, multipoles_fit, zbounds, mask=mask(0:jnpix-1,imap-1))

                    CALL map2alm_iterative(jnside, lmax, mmax, iter, temp_map(0:jnpix-1,1:p), temp_alm, w8ring=w8r)
                    deallocate(w8r)
!!$                    DO l = 0, lmax
! TEST
! --- keeping monopole and dipole
                    DO l = 2, lmax
! ---
!!$                       alm(:,l,0:mmax) = alm(:,l,0:mmax) + temp_alm(:,l,0:mmax) * sqrt(bl2(imap+1,l)) / need_norm
                       alm(:,l,0:mmax) = alm(:,l,0:mmax) + temp_alm(:,l,0:mmax) * sqrt(bl2(imap+1,l)) / need_norm !/ pwf(l,1)
                    ENDDO
                 ENDIF
              ENDDO

              WRITE(*,*) " ...synthesizing map"
              mapnside = jnside
              mapnpix = jnpix
              print*, " Output map nside: ", mapnside
              print*, " Output map npix : ", mapnpix
              allocate( w5map(0:mapnpix-1,1:p) )
              w5map = 0.
              CALL alm2map( mapnside, nl, nm, alm, w5map(:,1) )

 !!$              need_header = ''
              call write_minimal_header(need_header, 'MAP', nside=mapnside, ordering='Ring', coordsys='G', creator=ANACODE, version=MYVERSION, nlmax=lmax)
!!$              CALL add_card(need_header)
!              CALL add_card(need_header, 'MYCODE:', ANACODE)
!              CALL add_card(need_header, 'version', MYVERSION)
!!$              CALL add_card(need_header, 'B', B)
!!$              CALL add_card(need_header, 'NUM js', n_resol)
!              CALL add_card(need_header, 'nside', mapnside)
!              CALL add_card(need_header, 'ordering', 'ring')
              CALL add_card(need_header, 'filecont', 'Reconstructed map')

              tempfile = TRIM(need_root)//'_'//TRIM(io_B)//'_Nj'//TRIM(io_nresol)//'_reconstructed_map.fits'

              WRITE(*,*) " ...writing map on file: ", TRIM(tempfile)

!              CALL output_map(w5map(0:,1:1), need_header, TRIM(tempfile))
              CALL write_bintab(w5map, mapnpix, p,  need_header, 120, TRIM(tempfile))

              DEALLOCATE( need_map, temp_alm, alm, jflag, bl2 )

              CALL CPU_TIME(t2)

              WRITE(*,*) "elapsed time: ", INT((t2-t1)/60),"min", MOD((t2-t1),60.),"sec"

            END SUBROUTINE Ananeed_sub

!***********************************************************************
          subroutine bj_of_l

            implicit none

            REAL(dp):: Bx

! ---
            INTEGER(i4b), PARAMETER :: lmin = 2
            integer(i4b) :: lstr, lstp, lb, lbstr, lbstp, indx, j, l!, jmax
            integer(i4b) :: err, m, nside, npx, pj, bin_lines, bin, i_line, cen_l
            real(dp) :: cc, re, im, bj, needlet, aj, area
            real(dp) :: x
            real(dp) :: cen_x, low_x, high_x, top_hat, top_hat_c, top_hat_r
            integer(i4b) :: choice, bin_mode


            print*, " ***************************"
            print*, " *    Subroutine bj_of_l   *"
            print*, " ***************************"

            choice   = 1
            bin_mode = 1
 
            CALL intg(oneob,B,area,bsqr,1.d-5)

            write(*,'(1a27,1f15.6)') "  - bsqr has area equal to:", area
            top_hat_r = area

            cen_x = ( B + oneob ) / 2.d0
            low_x = cen_x - top_hat_r / 2.d0
            high_x = cen_x + top_hat_r / 2.d0
 
            do j = 1, jmax
               do l = 0, lmax
                  top_hat = 0.d0
!                  x = l / B**(j*deltaj)
                  x = l / B**j
                  needlet = bsqr(x)
                  if ( (x.ge.low_x) .and. (x.le.high_x) ) top_hat = 1.d0
                  if ( needlet .le. 1.e-30 ) needlet = 0.d0
                  if ( top_hat .le. 1.e-30 ) top_hat = 0.d0

                  bl2(j,l) = needlet
 
               enddo
            enddo

!!$ test
            print*, ""
            write(*,'(a22,f9.2)') " CHECK: sum of bl2: ", sum(bl2)
            IF ( (sum(bl2) - lmax-2+1.) .GT. 1.e-4) then 
               print*, " WARNING - If not equal to the sum ell=2,...,lmax you re missing power at low ell; the reconstruction may fail"
!               pause
            endif
            print*, ""
!!$            do j=1,jmax
!!$               print*,bl2(j,2)
!!$            enddo
!!$ = 499 -> correct, sum from ell=2,...,500
!!$              

            RETURN

        END subroutine bj_of_l

! **********************************************************************

        subroutine which_l
          
          IMPLICIT NONE

          INTEGER(i4b) :: j, l, fstl, lstl
          INTEGER(i4b), DIMENSION(0:lmax) :: lrange, temp
          INTEGER(i4b), DIMENSION(n_resol,2) :: lr

          lr = 0.
          temp = 0.
          DO l = 0, lmax
             lrange(l) = l
          ENDDO
          WRITE(*,*) " l-range per each j:"

          DO j = 1, n_resol
             IF (B**(j+1) .LT. 2.) THEN
                lr(j,:) = 0
             ELSE
                temp = 0.
                DO l = 0, lmax
                   IF (bl2(j,l) .GT. 0.) temp(l) = lrange(l)
                ENDDO
                DO l = 0, lmax
                   IF (temp(l) .NE. 0_i8b) THEN
                      lr(j,1) = temp(l)
                      EXIT
                   ENDIF
                ENDDO
                DO l = 0, lmax
                   IF (temp(lmax-l) .NE. 0_i8b) THEN
                      lr(j,2) = temp(lmax-l)
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
             IF (SUM(lr(j,:)) .GT. 0.) jflag(j) = .TRUE.
             WRITE(*,*) " j = ", j, lr(j,:), jflag(j)
          ENDDO

          allocate( lst_l(n_resol) )
          lst_l = lr(:,2)

          RETURN
          
        end subroutine which_l

! **********************************************************************
!
! New definition of psi(x) following the step-by-step file of Domenico
!
!1)        f(t) = exp( -1/(1-t**2) ), -1 <= t <= 1
!             = 0                   otherwise
!
!2)        psi(u) = int_-1^u f(t)dt / int_-1^1 f(t)dt
!
!3)        phi(t) = 1                            , 0 < t < 1/B
!                 = psi( 1 - 2B/(B-1)(t - 1/B) ) , 1/B <= t <= 1
!                 = 0                            , t > 1
!4)        bsqr(xi) = phi(xi/B) - phi(xi)
!
! We do not need to scale the function any more.

        function bsqr(xi)

          implicit none
          real(dp) :: bsqr
          real(dp) :: xi
 
          bsqr = phi(xi/B) - phi(xi)

          RETURN

        end function bsqr

!!!   ******************************************************8***********

        function f_sx(x)

          implicit none
          
          real(dp) :: x, arg
          real(dp) :: f_sx

          f_sx = 0.d0

          if ( x.lt.-1.d0 .and. x.gt.1.d0 ) f_sx = 0.d0
          if ( x.ge.-1.d0 .and. x.le.1.d0 ) then
             arg = 1.d0 - x**2
             if (arg.lt.1d-6) f_sx = 0.d0 
             if (arg.ge.1d-6) f_sx = dexp( -1.d0/arg )
          endif

          RETURN

        end function f_sx

!!!   ******************************************************************

        function psi(u)

          implicit none

          real(dp) :: u
          real(dp) :: psi, norm
          real(dp), parameter :: prec = 1.d-5

          call intg(-1.d0,1.d0,norm,f_sx,prec)
          call intg(-1.d0,u,psi,f_sx,prec)

          psi = psi / norm

          RETURN

        end function psi

!!!   ******************************************************************

        function phi(t)

          implicit none

          real(dp) :: phi
          real(dp) :: t, arg

          phi = 0.d0
          if (t.ge.0 .and. t.le.oneob)             phi = 1.d0

          arg = ( 1 - 2.d0*B/(B-1) * (t-oneob) )

          if ( (t.ge.oneob) .and. (t.le.1) )       phi = psi(arg)
          if (t.gt.1)                              phi = 0.d0

        end function phi


!!!   ******************************************************************

! --- :D --- Eventually to substitute with 'rombint'

          subroutine intg(xi, xf, y, f, tol)

            implicit none

            real(dp),intent(in) :: xi, xf
            real(dp) :: dx, f, yi, yf
            real(dp),intent(out) :: y
            real(dp), intent(in) :: tol
            integer(i8b), parameter :: ninterv = 10
            integer(i8b) :: n, i, cyc

            external f

            n = ninterv
            yi = 0.
            yf = 1.
            y = 0.d0
            cyc = 0
            do while ( abs((yi-yf))/yf .gt. tol )
               dx = abs(xf - xi)/n
               yi = y
               y = 0.d0
               do i = 0, n-1
                  y = y + dx/2.d0 * ( f(xi + i*dx) + f(xi + (i+1)*dx) )
               end do
               n = 2 * n
               yf = y
            enddo

            RETURN

            end subroutine intg


! **********************************************************************

            SUBROUTINE Build_Needlet(j, alms, nmap, nlmax, nmmax, need_nside, needlet)  

              USE alm_tools, only: alm2map

              IMPLICIT NONE

              REAL(dp), DIMENSION(0:),INTENT(OUT)         :: needlet
              INTEGER(i4b),INTENT(IN)                     :: j, nlmax, nmmax, need_nside, nmap
              COMPLEX(dpc), DIMENSION(1:nmap, 0:nlmax,0:nmmax), intent(in) :: alms

              COMPLEX(dpc), ALLOCATABLE, DIMENSION(:,:,:) :: alms_temp
              INTEGER(i4b) :: l, n_pix

              n_pix = 12.*need_nside**2
              ALLOCATE( alms_temp(1:1,0:nlmax,0:nmmax) )

!!$              DO l = 0, lmax
! TEST
              alms_temp(1,0:1,:) = 0. !alms(1,0:1,:)
! --- keeping monopole and dipole
              DO l = 2, nlmax
! ---
                 alms_temp(1,l,:) = alms(1,l,:) * sqrt(bl2(j,l)) / sqrt( real(n_pix) )
              ENDDO
              
              CALL alm2map(need_nside, nlmax, nmmax, alms_temp, needlet(0:n_pix-1))

!              needlet = needlet / sqrt(real(n_pix))

              DEALLOCATE(alms_temp)

              RETURN

            END SUBROUTINE Build_Needlet

! ----------------------------------------------------------------------

          END MODULE needlets_mod

!!! --------------------------------------------------------------------
!!! --------------------------------------------------------------------
!!! ********************************************************************
!!! ____________________________________________________________________

!!$          MODULE direct_space
!!$
!!$            USE healpix_types
!!$
!!$            private
!!$
!!$            public need_dir, corr_jjpkkp
!!$
!!$            REAL(dp), DIMENSION(0:5000)            :: plll
!!$
!!$          CONTAINS
!!$
!!$
!!$            function need_dir(j, costheta)
!!$
!!$!-----------------------------------------------------------------------
!!$!-   JPL - D.Pietrobon April 2010
!!$!-   Computes the covariance matrix according to -The normalization should be 1/sqrt(Npix) 
!!$!-----------------------------------------------------------------------
!!$
!!$              USE needlets_mod, only: bl2, lst_l
!!$
!!$              IMPLICIT NONE
!!$
!!$
!!$              RETURN
!!$              
!!$            END function need_dir
!!$
!!$! **********************************************************************
!!$
!!$            function corr_jjpkkp(j, jp, costheta, cls, nlmax)
!!$
!!$!-----------------------------------------------------------------------
!!$!-   JPL - D.Pietrobon April 2010
!!$!-   Computes the covariance matrix according to -The normalization should be 1/sqrt(Npix) 
!!$!-----------------------------------------------------------------------
!!$
!!$!              USE needlet_variables, only: lmax
!!$              USE needlets_mod, only: bl2, lst_l
!!$
!!$              IMPLICIT NONE
!!$
!!$              RETURN
!!$              
!!$            END function corr_jjpkkp
!!$
!!$! **********************************************************************
!!$            FUNCTION pl(l,x)
!!$
!!$              IMPLICIT NONE
!!$
!!$
!!$              RETURN
!!$
!!$            END FUNCTION pl
!!$
!!$!  ------------------------------------------------
!!$
!!$          END MODULE direct_space
