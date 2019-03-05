c       Program yptest

CCCCC     xfp in cm, ypfp and ypcor in mrad,   CCCCC

	Subroutine yp_optcor(xfp,ypfp,ypcor)
	implicit none

* Original version made on 01/07/06 by E. Christy
* to correct the MC_HMS Recon MEs (hms_cosey_refit_1.57.dat).  
* Interpolate in Xfp.

	real*4 xfp,ypfp,ypcor,coef(4,5),xfpb(5),xfplo,xfphi,dxfp
        real*4 ypcorlo,ypcorhi,slope
        integer i,j
        logical llo,hhi

        data coef/0.40767E-01,-0.12160E-01,-0.63210E-04,-0.79303E-05, 
     >           0.76542E-02, -0.44591E-02, -0.39353E-04, -0.28240E-04,
     >           0.88832E-02, -0.56952E-02, -0.20058E-03, -0.11670E-03,         
     >           -0.10000E-01, -0.50000E-01,  0.0,             0.0,
     >            0.20000E-01,  -0.64286E-01,  0.0,             0.0 /

        data xfpb/-20.,-10.,0.,10.,20./  

        llo = .false.
        hhi = .true. 
        xfphi = 1000.
        i = 1   

        if(abs(xfp).GT.40..OR.abs(ypfp).GT.40.) then
         ypcor = 0.0
         return
        endif  


c        write(6,*) "Enter xfp,ypfp"
c        read(5,*) xfp,ypfp

        if(xfp.LT.xfpb(1)) then 
          xfphi = xfpb(2)
          xfplo = xfpb(1)
          dxfp = xfp - xfpb(1)  
          i=2
          llo = .true.
          hhi = .false.
        elseif(xfp.GT.xfpb(5)) then
          xfphi = xfpb(5)
          xfplo = xfpb(4)
          dxfp = xfp - xfpb(5)
          i=5
          llo = .false.
          hhi = .true.
        else 
          dowhile(xfphi.GE.100)
            i = i+1
            if(xfpb(i).GT.xfp) then
              xfphi = xfpb(i)
              xfplo = xfpb(i-1)
              dxfp = xfp-xfplo
              hhi = .false.
            endif
          enddo
        endif

        ypcorlo = 0.0
        ypcorhi = 0.0

        do j=1,4
          ypcorlo = ypcorlo+coef(j,i-1)*ypfp**float(j-1)
          ypcorhi = ypcorhi+coef(j,i)*ypfp**float(j-1)
        enddo

        slope = (ypcorhi-ypcorlo)/(xfphi-xfplo)
        ypcor = ypcorlo + slope*dxfp
        if(hhi) ypcor = ypcorhi + slope*dxfp

c        write(6,*) xfp,ypfp,xfplo,xfphi,i,ypcorlo,ypcorhi,ypcor

	return
	end




