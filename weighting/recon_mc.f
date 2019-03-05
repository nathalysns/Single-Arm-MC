	program recon_mc
c       messed up with xsection table exttab
c       discard
	implicit none

        INCLUDE 'imodel.cmn'
        INCLUDE 'rad.cmn'
        INCLUDE 'logicals.cmn'
        INCLUDE 'flags.cmn'
       
        integer	nwpawc
        integer stat, runid,tarid
c	integer*4 wall
        character*5 id
	character*3 tarname
c	character*16 mergef(3) /
c     >            'mc_1.rzdat','mc_2.rzdat','mc_3.rzdat'/
        real    memor
        integer lrecl,nevt,ievt,ierr,istat
        real*4 ntuple_contents(27),xfoc,dydz,ypcor
        real*4 hse,hsp,p_rec,hsec,hsev,thetactemp,ebeam,radcon,thetac
        real*4 thetacrad,mp,mp2,nu,sin2,q2,w2,rcic,rcib,rcec
        real*4 rceb,trad,targetdata(6),eff_cal,eff_cer,emc
        real*4 pie,bcm1charge,bcm2charge,gbcm1charge,hpre,edt
        real*4 bcmavecharge,bmcur1,bmcur2,bmcur,eltime,cltime2
        real*4 positron_weight,posal,born,rci,rce,delup,deldown
        real*4 delini,ytarini,yptarini,xptarini,thetaini,ztari
        real*4 delrec,ytarrec,yptarrec,xptarrec,zrec,thetarec
        real*4 w,normfac,charge,ydata,lumfract,sigmacent,sigmac
        real*4 lumdata
        real*4 fract,sigave,sigtot,dxp,dyp,dep
        real*4 prescale,cltime,trackeff,trigeff
	real*4 xpup,ypup,dtup,hstheta
        real*4 dt,phasespcor,denscor,delcor,phase_space
        real*4 hmsprlo,hmstof,hms34,rate,rate_dat,poscs,t1
        real*8 xb
	

        real*4 ex,z,emean,de

	character*80 infile,outfile,title,ishit,new,name
        character*80 directory,icycle
        integer nentries,target,tartemp,fail_id,runnum,maxev
	integer	i, m, NtupleSize, bank,ngen,ngentot,ntrecl
        integer cs_flag,rc_flag
        logical firstr,goodfit,newrc,docs,dorc,use_rcmod

        parameter (nwpawc=500000)
        parameter (nentries = 27)        
        parameter(bank = 1000)
        parameter(title = 'RECONTUPLE')
        character*80 NtupleTag(nentries)
        real*4 ntu(nentries)
        common /pawc/ memor(nwpawc)
	logical*4	rd_int,iss

        ntrecl = 4096

	rate=0

        call getarg(1,id)        
        read(id,'(i5)') runnum

	call getarg(2,tarname)
       read(tarname,'(i1)') tarid 

c	iss=rd_int(tarname,wall)
c	write(6,*) runnum,wall,tarname

        
        firstr = .true.
	use_rcmod = .true.
        dorc = .true.        !   if true then include radiative contributions         !

         ngen = 0        

        mp = .9382723
        mp2 = mp*mp

  
c	inquire(file="rate.txt", exist=iss)
c	if (iss) then
c	   open(12, file="rate.txt", status="old", 
c    <                       position="append", action="write")
c	else
c	   open(12, file="rate.txt", status="new", action="write")
c	   write(12,*) '    runid       tarid    rate'
c	end if

c
 
        open(unit=18, file='input.dat',status='old') 
        read(18,*) name,ngentot,maxev,dxp,dyp,delup,deldown,cs_flag,
     &       rc_flag

        if(cs_flag.EQ.0) docs = .false.
        if(rc_flag.EQ.0) dorc = .false.



	
        radcon = 180./3.141592654
        
        sigave = 0.0

        

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCC      Read run info from database for data           CCCCCCC

        open(unit=16,file='runplan.inp',status='old')
c        open(unit=16,file='reconnmc.in',status='old')

	read(16,*)
        do 
          read(16,*,iostat=stat) runid,ebeam,hsec,thetac,
     &            prescale,
     &            bmcur1,bmcur2,
     &            bcm1charge,bcm2charge,cltime,eltime,trackeff,
     &            trigeff,rate_dat
            if(runnum .EQ. runid) exit
          if(stat /= 0) then 
            write(6,*) "Can't find run #",
     &                  id, " info in runplan.inp.  Stop"
            stop
          endif
        enddo
        write(6,*) "--found run #",id," info in reconmc.in"
        bcmavecharge = (bcm1charge+bcm2charge)/2. !!! Average over BCMs !!!
        bmcur = (bmcur1+bmcur2)/2.
        close(16)

        write(6,*) 'Run#        = ', id
        write(6,*) 'Target#     = ', tarid
        write(6,*) 'Ebeam       = ', ebeam  
        write(6,*) 'P0          = ', hsec 
        write(6,*) 'Theta       = ', thetac
        write(6,*) 'Prescale    = ', prescale 
 	write(6,*) 'Comp. Ltime = ', cltime
        write(6,*) 'Elect. Ltime = ', eltime        
        write(6,*) 'Track Eff.  = ', trackeff
        write(6,*) 'Trig Eff.   = ', trigeff
        write(6,*) 'Charge BCM1 = ', bcm1charge
        write(6,*) 'Charge BCM2 = ', bcm2charge
        write(6,*) 'Charge      = ', bcmavecharge
        write(6,*) 'Current     = ', bmcur
        write(6,*) 'Rate (kHz)  = ', rate_dat
        write(6,*)
        write(6,*) 'Including  Radiative Contributions? ', dorc
        write(6,*) 'Including Charge Symmetric Contributions?', docs
        write(6,*) 'Using Model from RC file?', use_rcmod
        thetacrad = thetac/radcon
        dep = (delup-deldown)/100.*hsec

c	if(target==1) tarname="d2"
c	if(target==2) tarname="h3"
c	if(target==3) tarname="he3"

	write(outfile,'("mc_",i1,".rzdat")') tarid
	write(infile,'(i1,".rzdat")') tarid

c	tarid = target
c	if(wall==1) tarid = 4
c	if(wall==3) tarid = 5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        
CCCCCC              Initialize radcor arrays               CCCCCCC
 
      
         firstr = .true. 
         trad = thetac/radcon
         call rc_mod(firstr,trad,trad,hsec,tarid,rcic,t1,id)
         write(6,*) firstr,trad,hsec,rcic
         firstr = .false.
        


CCCCCC        Read in target data from file                CCCCCCC
        open(unit=17, file='targetdata.dat',status='old') 
        do i=1,20
          read(17,*) targetdata
          
          if(i==tarid) then
c             write(6,*) "8888888888888888888888888"
c             write(6,*) tarid
c             write(6,*)    targetdata
c             write(6,*) "-----",targetdata(6),"-------",targetdata(4)
             lumdata = targetdata(6)*6.022137e-10/targetdata(4)  !!g/cm2 * 6.02e23/mol*1e-33 cm2/nb/(g/mol)
     &            *bcmavecharge/1.602177e-13          !!           * 1e-6 C/s /1.602e-19 C/e
           endif
       enddo        

	
	phase_space = 4.0*dxp/1000.*dyp/1000.*dep*1000    !rad*rad*Mev for XEMC(nb/MeV/sr)!!

	
c       scal factor.
	fract = lumdata*phase_space/ngentot
	fract = fract/1000.                        !if use Erics model (ub)
c	write(6,*) 'fract = ',fract
!!      density correction for gas target
	denscor = 1.0 - bmcur/25*0.2
	if(tarid.GT.3) denscor = 1.0
	
	eff_cer = 0.996		!!!!  Cer Efficiency  !!!!
	eff_cal = 0.999

	fract = fract*trackeff*trigeff*cltime*eltime
     &            *eff_cer*eff_cal/prescale*denscor

c	write(6,*) 'fract =',fract,trackeff,trigeff,cltime,eltime,
c     &           eff_cer,eff_cal,'/',prescale,denscor
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        call hlimit(nwpawc)
       
        lrecl = 0
        ishit = 'U'
        new = 'N'
        istat = 0
	
c	infile = trim(name)//'.rzdat'
c	if(wall.gt.3) then 
c	   infile = "combined.rzdat"
c	   call hmerge(3,mergef,infile)
c	   goto 777
c	endif
c	infile = '1.rzdat'
	call hropen (1, 'MCntuple', infile, ' ',ntrecl,istat) 
	call hgnpar (1, 'readdat') 
        call hnoent (1, nevt)  

c        call HCDIR(directory,'R') 
        write(6,*) '                    Old ntuple is: ', infile
        call HLDIR(' ',' ')


                 
 	m = 0
        m = m+1
        NtupleTag(m) = 'xfoc'
        m = m+1
        NtupleTag(m) = 'yfoc'
        m = m+1 
        NtupleTag(m) = 'xpfoc'
        m = m+1
        NtupleTag(m) = 'ypfoc'
        m = m+1
        NtupleTag(m) = 'ytari'
        m = m+1
        NtupleTag(m) = 'deltai'
        m = m+1
        NtupleTag(m) = 'yptari'
        m = m+1
        NtupleTag(m) = 'xptari'
        m = m+1
        NtupleTag(m) = 'ytar'
        m = m+1
        NtupleTag(m) = 'delta'  
        m = m+1
        NtupleTag(m) = 'yptar' 
        m = m+1
        NtupleTag(m) = 'xptar'        !dphr !m=12
         m = m+1          
        NtupleTag(m) = 'ztari'      !16
	m = m+1          
        NtupleTag(m) = 'ztar'
	m =m+1
        NtupleTag(m) =  'stopwhen'
        m = m+1
        NtupleTag(m) =  'x_stop'
        m = m+1
        NtupleTag(m) =  'y_stop'   !20
	m = m+1          
	NtupleTag(m) = 'born'
        m = m+1
        NtupleTag(m) = 'rci'
        m = m+1
        NtupleTag(m) = 'p_spec'
        m = m+1
        NtupleTag(m) = 'th_spec'
        m = m+1
        NtupleTag(m) = 'q2'
        m = m+1
        NtupleTag(m) = 'w2'
	m=m+1
	NtupleTag(m) = 'xbj'
	m = m+1
        NtupleTag(m) = 'yield'
	m = m+1
        NtupleTag(m) = 'rate'

	m=m+1
	NtupleTag(m) = 'tarid'
	
        ntuplesize = m
	
        call hropen (2, 'reconmc', outfile,'N',ntrecl,istat) 
        call hrin (2, 99999,0)

        call HBOOKN(9040,'reconmc',NtupleSize,'reconmc',1000,NtupleTag)
        call hcdir('//reconmc',' ')
        write(6,*) '                    New ntuple is:  ', outfile


        write(6,*)
        write(6,*) 'Number of events in master file =',nevt  
        if(maxev.GT.nevt) maxev = nevt
        write(6,*)
        write(6,*) 'Number of events analyzing = ',maxev
        write(6,*)

c	maxev = 800
        do ievt = 1, maxev      

    
         if(ievt.EQ.10000) write(6,*) ' analyzed 10000 events'
         if(ievt.EQ.50000) write(6,*) ' analyzed 50000 events'
         if(mod(ievt,10000).EQ.0.) then
	    write(6,*) ' analyzed',ievt,' events'
	 endif

         call hcdir('//MCntuple',' ')
	call hgnf(1, ievt, ntuple_contents, ierr)
c	write(6,*) ntuple_contents(17)
	if (ierr .ne. 0) then
	   WRITE(6,*) ievt
          stop 'hgnf err'
         endif 

         xfoc = ntuple_contents(1) !??
         dydz = ntuple_contents(4)  !??
         delini = ntuple_contents(6)          
         yptarini = ntuple_contents(7)*1000.
         xptarini = ntuple_contents(8)*1000.  !convert to mrad
         delrec = ntuple_contents(10)
         yptarrec = ntuple_contents(11)*1000. 
         xptarrec = ntuple_contents(12)*1000.
         fail_id = ntuple_contents(15)
         ytarini = ntuple_contents(5)
         ytarrec = ntuple_contents(9) 
	 ztari = ntuple_contents(16) 
	 hse = ntuple_contents(21) /1000.  !p_init in MeV to GeV
	 p_rec = ntuple_contents(23) /1000.  !p_recon in MeV to GeV

c         hse = hsec*(1.+delini/100.)
         
         hsev = hse

c         write(6,*) hse,emean,de

         hsp = hse

c         yptarrec = yptarrec - 0.90    !!!!!!   Test
c         xptarrec = xptarrec - 0.69 
c         yptarini = yptarini - 0.90    !!!!!!   Test
c         xptarini = xptarini - 0.69

         call yp_optcor(xfoc,dydz*1000.,ypcor)

         ypcor = 0.0

c         write(6,*) xfoc,dydz,yptarrec,ypcor

         yptarrec = yptarrec-ypcor 

         thetaini = acos(cos(thetacrad+yptarini/1000.)
     &              *cos(xptarini/1000.))             
      
         hstheta = acos(cos(thetacrad+yptarrec/1000.)
     &              *cos(xptarrec/1000.))
c          zrec  = ytarrec/cos(thetacrad+yptarrec/1000.)

CCCCCCC    Calculate the vertex kinematics for the event    CCCCCCCCC


         sin2 = sin(thetaini/2.)*sin(thetaini/2.)
         nu = ebeam - hsev
         q2 = 4.*hsev*ebeam*sin2
         xb = q2/2./mp/nu        
	 
	 w2 = mp2 + 2.*mp*nu-q2

c	 if((w2.LT.0).OR.(nu.LT.0)) then
c	    born = 0
c	    rci = 1
c
c	    goto 777
c	 endif

c	 w = sqrt(w2) 
        
CCCCCCC           Get Model Cross section in nb/SR/GeV, and born/rad factor         CCCCCCCCC


         call rc_mod(firstr,thetaini,thetacrad,
     &                                    hsev,tarid,rci,born,id)
c       

c	 if(born<0) born=0
         dt = thetaini - trad
         phasespcor = 1./cos(dt)/cos(dt)/cos(dt)

!!!!!!!     Changed on 6/19/01  !!!!!!!!
c          phasespcor = 1.

c          write(6,*) dt,phasespcor
 
c         write(6,*) ebeam,hse,thetaini,target,born

c         write(6,*) xb,AA,emc

         if(abs(xptarini).LT.dxp.AND.abs(yptarini).LT.dyp.AND.
     &      delini.GT.deldown.AND.delini.LT.delup.AND.born.GE.0.) then
c          sigave = sigave + born  
	    ngen = ngen + 1
         endif


C        the edge falls into (naive)w2<0 in grid.C, so the rc_mod can't interpolate correctly
c        make it run for now. fix later?
         if(rci.GT.100.OR.rci.LT.0.0001) rci = 1
	 if(born.LT.0) then
c	    write(6,*) "--------bad interpolation!!!--------------"
c	    write(6,*)  "--------thetaini,thetacrad,hsev,tarid,rci,born="
c	    write(6,*)  thetaini,thetacrad,hsev,tarid,rci,born
c	    write(6,*) "------------------!!!--------------"
	    born = 0
c	    fract = 0
	 endif


	 delcor =  1.009+.4260E-02*delrec-.8603E-03*delrec*
     &       delrec-.10942E-03*delrec*delrec*delrec+.12697E-04*
     &       delrec*delrec*delrec*delrec+.1094E-07*delrec*
     &       delrec*delrec*delrec*delrec

c          delcor =  1.0077+.4707E-02*delrec-.84067E-03*delrec*
c     &       delrec-.16119E-03*delrec*delrec*delrec+.13636E-04*
c     &       delrec*delrec*delrec*delrec+.64174E-06*delrec*
c     &       delrec*delrec*delrec*delrec



          if(delrec.GT.-10.0.AND.delrec.LT.-9.0) delcor=delcor/1.04
c          if(delrec.GT.-10.0.AND.delrec.LT.-9.0) delcor=1.00
          if(delrec.GT.-9.0.AND.delrec.LT.-8.0) delcor=delcor/0.99
c          if(delrec.GT.-9.0.AND.delrec.LT.-8.0) delcor=1.00
c          if(delrec.GT.-6.0.AND.delrec.LT.-5.0) delcor=delcor/1.012
          if(delrec.GT.9.0.AND.delrec.LT.10.0) delcor=delcor/0.99

CCCCCCC   Now Get RC corrections for event    CCCCCCCC

c          call rc_mod(firstr,thetaini,thetacrad,hsev,tarid,rci,t1,id)
c
     
          
c          if(w2.LT.1.18) then 
c	     rce = rce+1
c             rci = 1.
          
c          endif
                             
          if(.not.dorc) then
             rci = 1.
             rce = 1.
          endif



c          	  write(6,*) rci,rce

CCCCCCC          Calculate the efficiencies for Calorimeter and Cerenkov          CCCCCCCCC
 
c         eff_cal = 1.*exp(-0.0192/hse)
c         eff_cal = .999*exp(-0.0017114*hse**(-2.73))
c         eff_cer = .99622     


c          eff_cer = 0.998                                      !!!!  Cer Efficiency  !!!!
c         eff_cal = .999*exp(-0.0039258*hse**(-1.4930))        !!!!  Cal Efficiency  !!!!
cc          eff_cal = 0.96503+0.75590E-01*hse-0.65283E-01*hse**2+
c     &       0.26938E-01*hse**3-0.53013E-02*hse**4+0.39896E-03*hse**5  !!!! For Jan05 !!!!
c          eff_cal = eff_cal*(1.-0.000018*rate)


         positron_weight = 1.                                !!!!  reset just in case   !!!!
CCCCCCC            Fill new Ntuple                   CCCCCCCCC



 777	 do j = 1, 12
          ntu(j) = ntuple_contents(j)
         enddo
	 do j = 13, 17
	    ntu(j) = ntuple_contents(j+3)
	 enddo
         

         j = 18
         ntu(j) = born *phasespcor*delcor   !nb/MeV/sr
         j = j+1
         ntu(j) = rci
	 j = j+1
         ntu(j) = p_rec !use reconstructed p
         j = j+1
         ntu(j) = hstheta
	 j = j+1
         ntu(j) = q2
         j = j+1
         ntu(j) = w2
	 j = j+1
	 ntu(j) = xb
	 j = j+1
	 ntu(j) = ntu(18)/rci*fract
	 j = j+1
	 ntu(j) = ntu(j-1)*bmcur1/bcm1charge
c	 j = j+1
	 rate = rate + ntu(j)
c	 ntu(j) = rate
c
c	 if(rate>1000000.or.rate<0) then
c	 write(6,*) 'rate =',rate
c
cc	    write(6,*) born,phasespcor,delcor,rci,fract,ntu(j-1),rate
c
c        call rc_mod(firstr,thetaini,thetacrad,
c     &                                    hsev,tarid,rci,born,id)
c	 write(6,*) thetaini,thetacrad,born,rci
	 
c	 stop "wrong rate"
c	 endif
	j=j+1
	 ntu(j) = tarid


   
         call hcdir('//reconmc',' ') 
         call HFN(9040,ntu)
        
c     endif 
         
c     write(6,*)  sigave,trackeff,prescale,cltime,
c     &              eff_cer,eff_cal     
         
	enddo
	call hcdir('//reconmc',' ')
	call HROUT(9040,ICYCLE,' ')
	call HREND('reconmc')	!CERNLIB close file   
	close (1)
	close (2)


	write(6,*) "                 *****************"
	write(6,*) id,tarid, '    DONE!'
	write(6,*) 'rate =',rate
	write(6,*) 'fract =',fract
	write(6,*) 'lumdata =',lumdata
	write(6,*) 'phase_space = ', phase_space
        write(6,*) "                 *****************"

c	write(12,*) id,tarid,rate
c	close(12)
	close(16)
	close(17)
	close(18)


c 777	if(wall.gt.3) then

c	   write(6,*) " generated   ",infile

c	endif
	end
