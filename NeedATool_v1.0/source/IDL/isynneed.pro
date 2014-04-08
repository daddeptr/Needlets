; Davide Pietrobon June 2010
;
; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2008  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.jpl.nasa.gov
;
; -----------------------------------------------------------------------------
pro isynneed, map_in, B $
              , need_bl2=need_bl2
              , needlet_coef=needlet_coef $ 
              , bl2_root=bl2_root $ 
              , need_root=need_root $ 
;              , apply_windows=apply_windows $
              , beam_file=beam_file $
              , binpath=binpath $
;              , double=double $   
;              , fwhm_arcmin=fwhm_arcmin $
              , help=help $
;              , iseed=iseed $
              , keep_tmp_files=keep_tmp_files $
              , lmax=lmax ,       nlmax=nlmax $
;              , nside=nside ,     nsmax=nsmax $
;              , plmfile=plmfile $
;              , simul_type=simul_type $
              , silent=silent $
              , tmpdir=tmpdir $
;              , windowfile=windowfile $
;              , winfiledir=winfiledir 

;+
; NAME:
;    isynneed
;
; PURPOSE:
;    interface to 'synneed' F90 facility
;      can be used to    map -> needlet coefficients [ & needlet filters]
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    isynneed, map_in, B [,
;      , need_bl2=, needlet_coef=, bl2_root=, need_root=, beamfile=, help=,
;      keep_tmp_files=, lmax=, nside=, silent=, binpath=]
;
;
; INPUTS:
;   map_in: input map, can be a FITS file, or a memory array containing the
;        temperature map to be decomposed onto needlets.
;
;   B     : needlet harmonic filter parameter.
;
; OPTIONAL OUTPUTS:
;    map_out : map synthetised from the power spectrum or from constraining alm
;
; KEYWORD PARAMETERS:
;
;   beamfile=: beam window function, either a FITS file or an array
;
;   binpath=: full path to back-end routine [default: $HEXE/synfast, then $HEALPIX/bin/synfast]
;                 a binpath starting with / (or \), ~ or $ is interpreted as absolute
;                 a binpath starting with ./ is interpreted as relative to current directory
;                 all other binpaths are relative to $HEALPIX
;
;
;   /help:      if set, prints extended help
;
;   /keep_tmp_files: if set, temporary files are not discarded at the end of the
;                     run
;
;   lmax=,nlmax= : maximum multipole of simulation [default: 2*Nside]
;
;    /silent:    if set, works silently
;
;    tmpdir=:      directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;
; COMMON BLOCKS:
;    hxp_xface_com
;
; SIDE EFFECTS:
;    writes temporary and data files, spawn external processes
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;
;
; EXAMPLE:
;    isynneed, '$HEALPIX/test/map.fits', 2., lmax=500, /silent
;    mollview, map, 1, title='I'
;    mollview, map, 2, title='Q'
;
; will synthetize and plot I needlet coefficients of WMAP-1yr map.
;
; MODIFICATION HISTORY:
;
; June 2010 Davide Pietrobon v 1.0
;-
local = {routine: 'isynneed', exe: 'synneed', double: keyword_set(double)}
syntax = [local.routine+', map_in, B [,
      , need_bl2=, needlet_coef=, bl2_root=, need_root=, beamfile=, help=,
      keep_tmp_files=, lmax=, nside=, silent=, binpath=]']

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (undefined(cl_in) && undefined(alm_in)) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some input : cl_in or alm_in'
    return
endif
if (~arg_present(map_out) && undefined(map_out) && undefined(alm_out)) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some output: map_out or alm_out'
    return
endif

; check that no keyword is duplicated
;solve_kw_conflict,'alm_in', 'almsfile',     k1=alm_in,  k2=almsfile, kout=alm_in, /defined
solve_kw_conflict,'lmax',   'nlmax',         k1=lmax,    k2=nlmax,    kout=lmax, /defined
solve_kw_conflict,'nside',  'nsmax',         k1=nside,   k2=nsmax,    kout=nside, /defined

;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, tmpdir=tmpdir

NoFile = " '' "
;if (~arg_present(map_out)) then map_out = NoFile

; deal with online data
x_cl_in = set_parameter(cl_in,     NoFile, /ifempty, /ifzero)
tmp_cl_in     = hpx_mem2file(x_cl_in, /cl,   /in)
tmp_beam_file = hpx_mem2file(set_parameter(beam_file, NoFile, /ifempty), /beam, /in)
tmp_map_out   = hpx_mem2file((arg_present(map_out) || defined(map_out)) ? map_out : NoFile,  /out)

no_cl = undefined(x_cl_in) || (size(/tname,x_cl_in) eq 'STRING' && x_cl_in eq NoFile)
if (no_cl && undefined(alm_in)) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some input : cl_in or alm_in'
    return
endif

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '
tmp_nside                     = set_parameter(nside,       128)
printf,lunit,hpx_add_parameter('nsmax',       tmp_nside)
printf,lunit,hpx_add_parameter('simul_type', simul_type,  def=2,           /ifempty)
printf,lunit,hpx_add_parameter('nlmax',       lmax,        def=2*tmp_nside, /ifempty)
printf,lunit,hpx_add_parameter('iseed',       iseed,       def=0,           /ifempty)
printf,lunit,hpx_add_parameter('fwhm_arcmin',fwhm_arcmin, def=0,           /ifempty)
printf,lunit,hpx_add_parameter('apply_windows',apply_windows, def='F',    /ifempty)
;
printf,lunit,hpx_add_parameter('winfiledir', winfiledir, /expand, /skip_if_not_set)
printf,lunit,hpx_add_parameter('windowfile', windowfile,          /skip_if_not_set)
;
printf,lunit,hpx_add_parameter('infile',     tmp_cl_in,      /expand)
printf,lunit,hpx_add_parameter('beam_file',  tmp_beam_file, /expand)
             
printf,lunit,hpx_add_parameter('almsfile',  alm_in,   def=NoFile, /ifempty, /expand)
printf,lunit,hpx_add_parameter('plmfile',   plmfile,  def=NoFile, /ifempty, /expand)
;            
printf,lunit,hpx_add_parameter('outfile',      tmp_map_out,             /expand, /force)
printf,lunit,hpx_add_parameter('outfile_alms', alm_out,    def=NoFile, /expand, /force, /ifempty)
free_lun, lunit

; execute command
hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent

; deal with online data
if (arg_present(map_out)) then hpx_file2mem, tmp_map_out, map_out,/map ; 2008-08-27

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)

return
end
