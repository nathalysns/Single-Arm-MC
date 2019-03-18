C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C       SUBROUTINE HE4_ELASTIC
C
C       AUTHOR: D.W.Higinbotham
C       DATE:   February 2000
C       PURPOSE:
C               Calculate 4He cross sections for elastic scattering.
C               The form factor  calculations use the fit of  C. R. 
C               Ottermann et  al.  Nucl. Phys. A436 (1985) 688-698.
C               The Coulomb correction is made by using Q effective.  
C
C----------------------------------------------------------------------
C
      SUBROUTINE HE4_ELASTIC(E0,PF_E,TSCAT,POL_BEAM,HE4_SIG,HE4_ASYM,
     #                       PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
      IMPLICIT NONE
C
      DOUBLE PRECISION E0,PF_E,TSCAT,POL_BEAM,HE4_SIG,HE4_ASYM
      DOUBLE PRECISION PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD
      DOUBLE PRECISION ALPHA,AMOTT,COSTH2,SINTH2,TANTH2
      DOUBLE PRECISION EP,ETA,FACC,FACM,FACS,FCH,FMAG,HC,Q,Q4SQ
      DOUBLE PRECISION SIG,TAU,XMT,XMU,XMUA,Z,A,Qeff
C
C----------------------------------------------------------------------
C     Constants Used In The Calculation
C----------------------------------------------------------------------
C
      HC=197.327             !Conversion constant in MeV fm
      ALPHA=1./137.035989    !Fine Struction Constant 
      XMT=3727.40841         !4He Nuclear Mass of 3He           
      XMU=0.0                !4He Nuclear Magnetic Moment
      A=4.                   !4He Atomic Number
      Z=2.                   !4He Charge 
C
C----------------------------------------------------------------------
C     Calculation of Kinematic Quantities
C----------------------------------------------------------------------
C
      ETA= 1./(1.+ 2.*E0/XMT*SIN(TSCAT/2)**2)
      EP= ETA*E0
      Q4SQ= 4.*E0*EP*SIN(TSCAT/2)**2
      Q=SQRT(Q4SQ)/HC           !q4 in fm**-1
      SINTH2= SIN(TSCAT/2)**2
      COSTH2= COS(TSCAT/2)**2
      TANTH2= SINTH2/COSTH2
      AMOTT= 10000.*(HC*ALPHA/2./E0)**2* COSTH2/SINTH2**2
C
C----------------------------------------------------------------------
C     Simple Calculation Of The Effective Q  
C----------------------------------------------------------------------
C
      Qeff=Q*
     >     ( 1. + 1.5*Z*ALPHA*HC/E0/
     >        (1.12*A**0.33333)  )
C
C----------------------------------------------------------------------
C     Calculation Of The Elastic Form Factors  
C----------------------------------------------------------------------
C
      FCH=0.
      FMAG=0.
      Call HE4_FACTORS(FMAG,FCH,Qeff)
C
C----------------------------------------------------------------------
C     Approxiamte Calculation The DWBA Cross Section  
C----------------------------------------------------------------------
C
      TAU= Q4SQ/(4.*XMT*XMT)
      FACC= 1./((1.+TAU))
      XMUA= 3./Z*XMU
      FACM= TAU*XMUA**2* (FACC+ 2.*TANTH2)
      FACS= FCH*FCH*FACC+ FMAG*FMAG*FACM
      SIG= Z*Z*AMOTT*ETA*FACS ! microbarns/sr'
      HE4_SIG =SIG/10000.      !He4 cross section in fm^2/sr
C
C----------------------------------------------------------------------
C     Asymmetry Set To Zero 
C----------------------------------------------------------------------
C
      HE4_ASYM = 0.
C
C----------------------------------------------------------------------
C     Recoil Polarizations P_n Set To Zero
C----------------------------------------------------------------------
C
      PN_HI = 0.
      PT_HI = 0.
      PL_HI = 0.
      PN_HD = 0.
      PT_HD = 0.
      PL_HD = 0.
C      
      RETURN
      END
C
C----------------------------------------------------------------------
C    Subroutine to generate the form factors of 4He using 
C    C.R.Ottermann et al., Nucl. Phys A436(1985)688-698.   
C----------------------------------------------------------------------
C
      SUBROUTINE HE4_FACTORS(FMAG,FCH,Q)

      IMPLICIT NONE

      DOUBLE PRECISION A,B,Q,FMAG,FCH

      A=0.316
      B=0.675
 
      FMAG=0.0 
      FCH=(1-(A**2*Q**2)**6)*exp(-1*B**2*Q**2)

      RETURN
      End 
