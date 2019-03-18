C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C       SUBROUTINE he3_elastic
C
C       AUTHOR: A. Deur
C       DATE:   May 1999
C       PURPOSE:
C               Calculate cross section and polarizations for 3He
C               elastic scattering.
C
C       MODIFICATIONS:
C
C       AUTHOR: D.W.Higinbotham
C       DATE:   January 2000
C       PURPOSE:
C               Update the form factor calculations using the global
C               fit of A. Amroun et al., Nucl. Phys. A579 (1994) 596
C               and to use Qeff to approximate a DWBA calculation of
C               the form factors. Future improvement would be to use
C               a code like HADES to do an exact DWBA calculation.  
C
C----------------------------------------------------------------------
C
      SUBROUTINE HE3_ELASTIC(E0,PF_E,TSCAT,POL_BEAM,HE3_SIG,HE3_ASYM,
     #                       PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
      IMPLICIT NONE
C
      DOUBLE PRECISION E0,PF_E,TSCAT,POL_BEAM,HE3_SIG,HE3_ASYM
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
      XMT=2808.413           !3He Nuclear Mass of 3He           
      XMU=-2.12762485        !3He Nuclear Magnetic Moment
      A=3.                   !3He Atomic Number
      Z=2.                   !3He Charge 
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
C     Hesenberg and Blok, Ann. Rev. Nucl. Part. Sci. 33 (1983) 574.
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
      call HE3_FACTORS(FMAG,FCH,Qeff)
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
      HE3_SIG =SIG/10000.      !he3 cross section in fm^2/sr
C
C----------------------------------------------------------------------
C     Asymmetry Calculation Still Needs To Be Added
C----------------------------------------------------------------------
C
      He3_ASYM = 0.
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
C    Subroutine to generate the form factors of 3He
C    written by D.W.Higinbotham using the parameterization
C    A. Amroun et al., Nucl. Phys. A579 (1994) 596-626,
C    with the correct form of equation (1) from
C    I. Sick, Nucl. Phys. A218 (1974) 509.
C----------------------------------------------------------------------
C
      SUBROUTINE HE3_FACTORS(FMAG,FCH,Q)

      DOUBLE PRECISION Q,GAMMA,FM,FC,FMAG,FCH,HE3_CFF
      DOUBLE PRECISION HE3_MFF,QM(12),R(12),QC(12)

      DATA QM/0.059785,0.138368,0.281326,0.000037,0.289808,
     _0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.0/
      DATA QC/0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,
     _0.074310,0.053970,0.023689,0.017502,0.002034,0.004338/
      DATA R/0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.0,4.6,5.2/

        GAMMA=0.65320

        FM  = 0.0
        FC  = 0.0
        HE3_MFF = 0.0
        HE3_CFF = 0.0

        FM = 0.0
        do 100 i=1,12
         A= QM(i)/(1.+2.*(R(i)**2)/GAMMA**2)
         B = cos(Q*R(i))
         C = (2.*R(i)**2)/GAMMA**2
         D = (sin(Q*R(i)))/(Q*R(i))
         FM = A*(B+C*D) +FM
100     continue

        HE3_MFF=FM*exp(-0.25*Q**2*GAMMA**2)
        HE3_MFF=sqrt(HE3_MFF*HE3_MFF)
        FMAG=HE3_MFF

        do 200 i=1,12
         A= qc(i)/(1.+2.*(R(i)**2)/GAMMA**2)
         B = cos(Q*R(i))
         C = (2.*R(i)**2)/GAMMA**2
         D = (sin(Q*R(i)))/(Q*R(i))
         FC = A*(B+C*D) +FC
200     continue

        HE3_CFF=FC*exp(-0.25*Q**2*GAMMA**2)
        HE3_CFF=sqrt(HE3_CFF*HE3_CFF)
        FCH=HE3_CFF
      RETURN
      END
