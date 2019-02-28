	program mc_hrs_single

C modified for 25cm Tritium target. multiple scattering???
C CHANGES FOR OPTICS TESTING:
C 1. Remove target multiple scattering/energy loss

C+______________________________________________________________________________
!
! Monte-Carlo of HRS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
!		which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
C-_____________________________________________________________________________

	implicit none

	include 'hrsl/struct_hrsl.inc'
	include 'spectrometers.inc'
	include 'g_dump_all_events.inc'
	include 'constants.inc'

C HBOOK/NTUPLE common block and parameters.
	integer*4	pawc_size
	parameter	(pawc_size = 80000)
	common		/pawc/ hbdata(pawc_size)
	integer*4	hbdata,nhut
	parameter	(nhut = 27)
	character*8	hut_nt_names(nhut)/
     >			'hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >			'hsytari', 'hsdeltai', 'hsyptari', 'hsxptari',
     >			'hsytar', 'hsdelta', 'hsyptar', 'hsxptar','fry',
     >                  'frx', 'ok_spec','hsztari','hsztar',
     >                  'stopwhen', 'x_stop','y_stop',
     >                  'p_init','e_init','p_recon','e_recon',
     >                  'elossi','elossf',"eloss0"/


	real*4		hut(nhut)

C Local declarations.
	integer*4	i,
     >			chanin	/1/,
     >			chanout	/2/,
     >			n_trials,trial,
     >			tmp_int

	logical*4	iss,extended

C Event limits, topdrawer limits, physics quantities
	real*8 gen_lim(8)			!M.C. phase space limits.

	real*8 gen_lim_up(3)
	real*8 gen_lim_down(3)

	real*8 cut_dpp,cut_dth,cut_dph,cut_z	!cuts on reconstructed quantities
	real*8 xoff,yoff,zoff                   !Beam offsets (z target offset)
	real*8 spec_xoff,spec_yoff,spec_zoff    !Spectrometer offsets
	real*8 spec_xpoff, spec_ypoff           !Spectrometer angle offsets
	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 th_ev,cos_ev,sin_ev		!cos and sin of event angle
	real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4	!temporaries
	real*8 musc_targ_len			!target length for multiple scattering
	real*8 m2				!particle mass squared.
	real*8 mass                             !particle mass  

	real*8 rad_len_cm			!conversion r.l. to cm for target
	real*8 pathlen				!path length through spectrometer.
C DJG Variables used for target calcs
	real*8 t,atmp,btmp,ctmp
	real*8 side_path,costmp,th_can,s_Al,zz
	real*8 y_coff,x_beam,arg1,aa1,aa2 !used when reconstructing z vertex


C from barak-------for energy loss
	real*8 targ_Z,targ_A  
	real*8 targ_rho                         !target density in g/cm^3 

	real*8 s_air, s_mylar, s_target0

    	real*8 theta_rec,cos_rec,sin_rec        !cos and sin of reconstructed angle

	real*8 momentumi,momentums,momentumf    !particle momentum (vertex,spectrometer,reconstructed at spec.)
	real*8 energyi,energys,energyf          !particle energy (vertex,spectrometer,reconstructed at spec.)
	integer*4 typeflag	                !1=generate eloss, 2=min, 3=max, 4=most probable
	real*8 Elossi,Elossf,Etemp,eloss0              !energy loss (real,reconstructed)
	real*8 Eloss_target,Eloss_Al,Eloss_air,Eloss_mylar !temporary elosses
	logical*4 elossi_flag
	logical*4 elossf_flag /.false./


C shujie    -------for beer can target multiple scattering-------------------
	real*8 s_wall, s_end ! thickness of the target side wall and end cap
	real*8 target_angle, target_r, s_target
	real*8 forward_path
C	------------------------------------------

C Miscellaneous
	logical*4 ok_spec			!indicates whether event makes it in MC
c	integer*4 hit_calo                      !flag for hitting the calorimeter

C Initial and reconstructed track quantities.
	real*8 dpp_init,dth_init,dph_init,xtar_init,ytar_init,ztar_init
	real*8 dpp_recon,dth_recon,dph_recon,ytar_recon,ztar_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		!at focal plane
	real*8 frx,fry,fr1,fr2
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 resmult
	real*8 good_evt

C Control flags (from input file)
	integer*4 p_flag			!particle identification
	logical*4 ms_flag
	logical*4 wcs_flag
	integer*4 col_flag
	logical*4 gen_evts_file_flag

C Hardwired control flags.
	logical*4 hut_ntuple	/.true./
	logical*4 decay_flag	/.false./

	real*8	dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
	real*4	stime,etime,zero

	character*132	str_line

C Random Number Generation
	integer*4 rnd_seed_flag
	character*80 random_state_file
	logical restorerndstate

C Function definitions.

	integer*4	last_char
	logical*4	rd_int,rd_real
	real*8          grnd,gauss1

        character*80 rawname, filename
	character*5  id
	integer*4 tarid
	integer   runnum,runid,stat
	real*4  secnds

        integer iquest
        common/quest/iquest(100)

	save		!Remember it all!

C ================================ Executable Code =============================

C Initialize
C xiaochao:
C using SIMC unstructured version
C
	lSTOP_trials	= 0
	lSTOP_col_entr  = 0
	lSTOP_col_exit  = 0
	lSTOP_spec_entr = 0
	lSTOP_Q1_in	= 0
	lSTOP_Q1_mid	= 0
	lSTOP_Q1_out	= 0
	lSTOP_Q2_in	= 0
	lSTOP_Q2_mid	= 0
	lSTOP_Q2_out	= 0
	lSTOP_Q3_in	= 0
	lSTOP_Q3_mid	= 0
	lSTOP_Q3_out	= 0
	lSTOP_D1_in	= 0
	lSTOP_D1_out	= 0
	lSTOP_hut	= 0
	lSTOP_dc1	= 0
	lSTOP_dc2	= 0
	lSTOP_s0	= 0
	lSTOP_cer       = 0
	lSTOP_s2	= 0
	lSTOP_prl1	= 0
	lSTOP_prl2      = 0
	lSTOP_successes	= 0

C Open setup file.

C	write(*,*)'Enter input filename (assumed to be in infiles dir)'
C	read(*,1968) rawname

	call getarg(1,rawname)
 1968	format(a)
	filename = 'infiles/'//rawname(1:last_char(rawname))//'.inp'
	write(6,*) filename,'opened'
	open(unit=chanin,status='old',file=filename)


!       read runplan by runid
	call getarg(2,id)
	read(id,'(i5)') runid
	filename = 'infiles/runplan.inp'
	write(6,*) filename,'opened'

	open(unit=12,status='old',file=filename)
	read(12,*) str_line
	write(6,*) str_line
	do
	   read(12,*,iostat=stat) runnum,p_spec,p_spec,th_spec
	   if(runnum .EQ. runid) exit
	   if(stat/=0) then            
	      write(6,*) "Can't find run #",
     &           id, " info in infiles/runplan.inp. Stopped"
	      stop
	   endif
	enddo
	write(6,*) "--found run #",id
	p_spec = p_spec*1000 !gev to mev
	write(6,*)  "ep =",p_spec,"angle =",th_spec
	

!       Read the flag for vertex generating location 
!      ( 1/2/3: D2,H3,He3 gas
!      /4/5: entrance window, exit window)
!       6: single carbon foil

	call getarg(3,rawname)
	iss = rd_int(rawname,tarid)

	if(tarid.GT.6.or.tarid.LT.0)  then
	   tarid = 2
	   write(6,*) 'Invalid input.Use default vertex flag (tarid = 2, H3 gas)'
	endif
	write(6,*) '----generate event at tarid = ', tarid



C Initialize HBOOK/NTUPLE if used.
	if (hut_ntuple) then
	  call hlimit(pawc_size)
c	  filename = 'worksim/'//rawname(1:last_char(rawname))//'.rzdat'
c	  filename = 'worksim/'//rawname(1:last_char(rawname))
c	  write(6,*) filename  
	  write(filename,'("worksim/Ar_",a5,"_",i1,".rzdat")') id,tarid
!	  call hropen(30,'HUT',filename,'N',1024,i)
	  iquest(10) = 256000
	  iquest(10) = 510000
c	  write(6,*) filename
! see for example
!   http://wwwasd.web.cern.ch/wwwasd/cgi-bin/listpawfaqs.pl/7
! the file size is limited to ~260M no matter how I change iquest !
	  call hropen(30,'HUT',filename,'NQ',4096,i) !CERNLIB
 
	  if (i.ne.0) then
	    write(6,*),'HROPEN error: istat = ',i
	    stop
	  endif
	  call hbookn(1,'HUT NTUPLE',nhut,'HUT',10000,hut_nt_names)
	endif	   

C Open Output file.
c	write(filename,'("outfiles/src_",a5,"_",i1,".hist")') id,tarid
	write(filename,'("outfiles/Ar_",a5,"_",i1,".hist")') id,tarid

c	filename = "outfiles/"//runid//"_"//tarid//".hist"
	open (unit=chanout,status='unknown',file=filename)

C Read in real*8's from setup file

	str_line = '!'

C Strip off header

	do while (str_line(1:1).eq.'!')
	  write(6,*) str_line(1:last_char(str_line))
	  read (chanin,1001) str_line
	enddo

! Read data lines.

! N_TRIALS:
c	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_int(str_line,n_trials)
	if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer momentum: !!now read from runplan.inp
	read (chanin,1001) str_line
!	write(6,*) str_line(1:last_char(str_line))
!	iss = rd_real(str_line,p_spec)
!	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! Spectrometer angle:  !!now read from runplan.inp
	read (chanin,1001) str_line
!	write(6,*) str_line(1:last_char(str_line))
!	iss = rd_real(str_line,th_spec)
!	if (.not.iss) stop 'ERROR (Spec theta) in setup!'
	th_spec = abs(th_spec) / degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
	do i=1,3
	  read (chanin,1001) str_line
	  write(6,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_down(i) = gen_lim(i)
	  read (chanin,1001) str_line
	  write(6,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_up(i) = gen_lim(i)
	enddo
!	do i=1,3
!	   write(*,*)'gen_lim_up/down = ',gen_lim_up(i),' ',gen_lim_down(i)
!	enddo

	do i = 4,6
	  read (chanin,1001) str_line
	  write(6,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	enddo
	extended = .false.
	if(gen_lim(6)>3) extended = .true.
! Raster size

	do i=7,8
	   read (chanin,1001) str_line
	   write(6,*) str_line(1:last_char(str_line))
	   iss = rd_real(str_line,gen_lim(i))
	   if (.not.iss) stop 'ERROR (Fast Raster) in setup'
	enddo

! Cuts on reconstructed quantities
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dpp)) stop 'ERROR (CUT_DPP) in setup!'

	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dth)) stop 'ERROR (CUT_DTH) in setup!'

	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dph)) stop 'ERROR (CUT_DPH) in setup!'

	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_z)) stop 'ERROR (CUT_Z) in setup!'

! Read in radiation length of target material in cm
	read (chanin,1001) str_line
!	write(6,*) str_line(1:last_char(str_line))
!	if (.not.rd_real(str_line,rad_len_cm)) stop 'ERROR (RAD_LEN_CM) in setup!'
	if(tarid==0) rad_len_cm = 999!dummy
	if(tarid==1) rad_len_cm = 122.6 !D2
	if(tarid==2) rad_len_cm = 183.6 !H3
	if(tarid==3) rad_len_cm = 71.07 !He3
	if(tarid==4.or.tarid==5) rad_len_cm = 8.89 !Al
	if(tarid==6) rad_len_cm = 19.32 !C`
	
! Beam and target offsets
	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,xoff)
	if(.not.iss) stop 'ERROR (xoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,yoff)
	if(.not.iss) stop 'ERROR (yoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,zoff)
	if(.not.iss) stop 'ERROR (zoff) in setup!'

! Spectrometer offsets
	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xoff)
	if(.not.iss) stop 'ERROR (spect. xoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_yoff)
	if(.not.iss) stop 'ERROR (spect. yoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_zoff)
	if(.not.iss) stop 'ERROR (spect. zoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xpoff)
	if(.not.iss) stop 'ERROR (spect. xpoff) in setup!'

	read (chanin, 1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_ypoff)
	if(.not.iss) stop 'ERROR (spect. ypoff) in setup!'

! read in flag for particle type.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) stop 'ERROR: p_flag in setup file!'

! Read in flag for aerogel usage in the HMS.
!	read (chanin,1001) str_line
!	write(6,*) str_line(1:last_char(str_line))
!	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: use_aer in setup file!'
!	if (tmp_int.eq.1) use_aer = .true.

! Read in flag for multiple scattering.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: ms_flag in setup file!'
	if (tmp_int.eq.1) ms_flag = .true.

! Read in flag for wire chamber smearing.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: wcs_flag in setup file!'
	if (tmp_int.eq.1) wcs_flag = .true.

! Read in flag for dumping ALL events into the HUT NTUPLE (note that
! the ...recon quantities will be ill defined but the FAIL_ID could be
! used to tell...
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:dump_all_in_ntuple in setup file!'
	if (tmp_int.eq.1) dump_all_in_ntuple = .true.

! Read in flag that sets the collimator option.
	read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,col_flag)) stop 'ERROR:col_flag in setup file!'

! Read in flag that sets random number generation option
	read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,rnd_seed_flag)) stop 'ERROR:rnd_seed_flag in setup file!'

! Read in flag for using generated events file.
        read (chanin,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: gen_evts_file_flag in setup file!'
	if (tmp_int.eq.1) gen_evts_file_flag = .true.

 
! Read in flag for ionization energy loss.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: elossi_flag in setup file!'
	if (tmp_int.eq.1) elossi_flag = .true.
	write(6,*)  'elossi_flag =', elossi_flag
! Read in flag for ionization energy loss.
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: elosf_flag in setup file!'
	if (tmp_int.eq.1) elossf_flag = .true.

! Read in density of target in g/cm^3
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,targ_rho)) stop 'ERROR (TARG_RHO) in setup!'

! Read in target atomic number (Z)
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,targ_Z)) stop 'ERROR (TARG_Z) in setup!'

! Read in target standard atomic weight (A)
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,targ_A)) stop 'ERROR (TARG_A) in setup!'


!       Read target width(radius)
	read (chanin,1001) str_line
	write(6,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,target_r)
	if(.not.iss) stop 'please provide target radius in setup!'
	


C       For solid target
	s_wall = 0 
	s_end  = 0
C	target_r = 0.26*2.54
C       For extended target
	if(extended) then
	
!       Read the end cap thickness for beer can target
	   read (chanin,1001) str_line
	   write(6,*) str_line(1:last_char(str_line))
	   iss = rd_real(str_line,s_end)
	   if(.not.iss) stop 'please provide target window thickness in setup!'
	   
	   
!       Read the side wall thickness for beer can target
	   read (chanin,1001) str_line
	   write(6,*) str_line(1:last_char(str_line))
	   iss = rd_real(str_line,s_wall)
	   if(.not.iss) stop 'please provide target sidewall thickness in setup!'
	   
	endif
	
	write(6,*)  "endcap_thickness = ",s_end,"side_wall =",s_wall


C =============================== Format Statements ============================

1001	format(a)

	end
