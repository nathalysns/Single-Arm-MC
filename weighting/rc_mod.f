      subroutine rc_mod(firstr,theta,thspect,xin,tarid,rci,sig,runid)
      implicit none
      INCLUDE 'rad.cmn'
      character*80 infile,filename
      integer*4 i,j,k,a,z
      integer aa(6) /2,3,3,27,27,12/
      integer zz(6) /1,1,2,13,13,6/
      real*8 m(6)  /0,0,0,0,0,0/
      integer temp
      character*5  runid
      real*4 xin,rci,theta,thspect,thcentdeg,sig
      real*8 x(6),rc(6)
      real*8 radtab_temp(7),thetadeg,thcent,thetalow,thetahigh
      real*8 thetatab,diffxL1,diffxL2,diffxH1,diffxH2
      real*4 deglow,deghigh     !to check if deg in table range      
      logical  xl1,xl2,xh1,xh2  ! check the p range
      integer*4 tarnum,tarid,tar(20),eof,tdiff,tdiff_min,tdiff_max
      real*8 mp,mp2,radcon,thetarad
      real*8 xtab,xtab_next,xtab_pre
      logical firstr,endof,extrap1,extrap2,extrap_x_hi,extrap_x_lo


      xl1=.false.
      xl2=.false.
      xh1=.false.
      xh2=.false.

c      write(filename,"(A,A5,A,I1)") "table.out.",runid,".",tarid
      write(filename,"(A,A5,A,I1)")"XStables/table.out.",runid,".",tarid


      infile = trim(filename)
c      write(6,*) filename, infile
      
c      infile = 'table.out.'//runid//'.'//tarid
c      infile = 'rc94.dat'
      radcon = 180./3.141593
      thetadeg = theta*radcon       !!! convert rad to deg !!!
c      thcentdeg = thspect*radcon    !!! convert rad to deg !!!

      do i=1,6
        rc(i) = 0.
        m(i) = 0.
      enddo

CCCCCC              read in radcor table              CCCCCC

c      write(6,*)"here",tarid,tar(tarid),firstr
      
      if (firstr) then 
         open(unit=34,file=infile,status='old')    
         i = 1
         eentries = 0 
         endof = .false.
         temp = tarid
 
         do while(.not.endof)
            read(34,*,END=1001) radtab_temp
c           write(6,*) radtab_temp
            a = radtab_temp(1)
            z = radtab_temp(2)
c            write(6,*) a,z,aa(temp),zz(temp)
            if(aa(temp)==a.and.zz(temp)==z) then
c            write(6,*) runid,tarid,a,z,eentries
               do j=1,5
                  exttab(i,j) = radtab_temp(j+2)
c                  write(6,*) exttab(i,j)
               enddo 
               eentries = eentries + 1
            endif
            i = i + 1 
         enddo 
      endif
 1001 endof = .true.
 
       if (firstr) then 
         write(6,*) "Nentries in radcor table is:  ",eentries
      endif

 

      close(34) 

c      write(6,*) "Nentries in radcor table is:  ",eentries

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    
CCCCCC       Calculate radiative correction and model by doing    CCCCCC
CCCCCC       linear interpolation in theta and xin                CCCCCC

        tdiff_min = 10.0
        tdiff_max = 10.0

        do j=1,eentries
          thetatab = exttab(j,3)
          tdiff = abs(thetatab - thetadeg)
          if(thetatab.LT.thetadeg) then
             deglow = 1 ! the table has a point lower than thetadeg
            if(tdiff.LT.tdiff_min) then 
              tdiff_min = tdiff 
              thetalow = thetatab
            endif
          else
             deghigh = 1        ! the table has a point higher than thetadeg
            if(tdiff.LT.tdiff_max) then 
              tdiff_max = tdiff
              thetahigh = thetatab
            endif 
          endif
        enddo
        
        if((deglow*deghigh).LT.0.1) then ! if thetadeg is out of table coverage
          write(6,*) "theta=",thetadeg,"out of range,check you table!"
          write(6,*) "set born=0, rci = 1 for this point"
          sig = 0
          rci = 1
          goto 777
       endif



c        thetalow = int(thetadeg)-tdiff_min     !!! find integer angle below !!! 
c        thetahigh = int(thetadeg)+tdiff_max    !!! find integer angle above !!!



CCCCCC     do search for rcs to interpolate in theta and xin.     CCCCCC 
CCCCCC     thetahigh is the integer theta above the               CCCCCC
CCCCCC     central theta.                                         CCCCCC
 

        extrap_x_lo = .true.
        extrap_x_hi = .true.
        diffxh1 = 1000.
        diffxL1  = 1000.
        diffxh2 = 1000.
        diffxL2  = 1000.
        xl1 = 0
        xl2 = 0
        xh1 = 0
        xh2 = 0
     

        do j=1,eentries
          thetatab = exttab(j,3)
      
  
c            extrap_th_lo = .false.
       if(abs(thetatab-thetalow)<0.0001) then  
         
          if(exttab(j,2).LE.xin) then    
             xl1 = .true.
              if(abs(exttab(j,2)-xin).LE.diffxL1) then
               
                diffxL1 = abs(exttab(j,2)-xin)
                x(1) = exttab(j,2)
                m(1) = exttab(j,4)
                rc(1) = exttab(j,5)
                   endif
            elseif(exttab(j,2).GT.xin) then
               xh1 = .true.
              if(abs(exttab(j,2)-xin).LE.diffxH1) then  
                diffxH1 = abs(exttab(j,2)-xin)       
                x(2) = exttab(j,2)
                m(2) = exttab(j,4)
                rc(2) = exttab(j,5)
              endif
            endif 
          endif

          if(abs(thetatab-thetahigh)<0.0001) then
          
c            extrap_th_hi = .false.
            if(exttab(j,2).LE.xin) then
               xl2 = .true.
              if(abs(exttab(j,2)-xin).LE.diffxL2) then
                diffxL2 = abs(exttab(j,2)-xin)
                x(3) = exttab(j,2)
                m(3) = exttab(j,4)
                rc(3) = exttab(j,5)
              endif
            elseif(exttab(j,2).GT.xin) then
               xh2 = .true.
              if(abs(exttab(j,2)-xin).LE.diffxH2) then       
                diffxH2 = abs(exttab(j,2)-xin)       
                x(4) = exttab(j,2)
                m(4) = exttab(j,4)
                rc(4) = exttab(j,5)
              endif
            endif 
          endif
     
        enddo 

c        if((m(1)*m(2)*m(3)*m(4))<0.1) then
        if(.not.(xl1.and.xl2.and.xh1.and.xh2)) then
c           write(6,*) m(1),m(2),m(3),m(4)
           write(6,*) "theta, Ep=",thetadeg, xin,"out of range"
           write(6,*) "set born = 0, rci = 1. check your table!"
          sig = 0
          rci = 1
          goto 777
       endif
 


        m(5) = (m(2)*(xin-x(1))+m(1)*(x(2)-xin))/(x(2)-x(1))
        rc(5) = (rc(2)*(xin-x(1))+rc(1)*(x(2)-xin))/(x(2)-x(1))
        m(6) = (m(4)*(xin-x(3))+m(3)*(x(4)-xin))/(x(4)-x(3))
        rc(6) =(rc(4)*(xin-x(3))+rc(3)*(x(4)-xin))/(x(4)-x(3))


c 


        sig = m(6)*(thetadeg-thetalow)+m(5)*(thetahigh-thetadeg)
        sig = sig/(thetahigh-thetalow)      
        rci = rc(6)*(thetadeg-thetalow)+rc(5)*(thetahigh-thetadeg)
        rci = rci/(thetahigh-thetalow)

c        if(firstr)  write(6,*)"here: ",m(3),m(6),sig,rce,firstr

CCCCCC                          End search                            CCCCCC

   
CCCCCC             Now do interpolation in theta                      CCCCCC


c        sig = (m(3)-m(6))*(thetadeg-thetahigh)/
c     &        (thetalow-thetahigh)+rc(6)
c        rce = (rcr(3)-rcr(6))*(thetadeg-thetahigh)/
c     &        (thetalow-thetahigh)+rcr(6)

c        if(rce.LE.0) rce = 0.000001
c        if(sig.LE.0) sig = 0.000001

c        write(6,*) thetadeg,xin,sig,rce

 8000 format(a80) 
      close(34)
 777  return

      end





















