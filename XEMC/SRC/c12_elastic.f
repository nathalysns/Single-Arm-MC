C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C       SUBROUTINE C12_ELASTIC
C
C       AUTHOR: D.W.Higinbotham
C       DATE:   February 2000
C       PURPOSE:
C               To calculate the elastic cross section from C12  using 
C               the  Fourier  Bessel  coefficents  from  the  paper by
C               E.A.J.M. Offermann et al., Phys. Rev. C 44 (1991) 1096.
C               The Coulomb correction is made by using an effective Q.
C
C----------------------------------------------------------------------
C
      SUBROUTINE C12_ELASTIC(E0,PF_E,TSCAT,POL_BEAM,C12_SIG,C12_ASYM,
     #                       PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
      IMPLICIT NONE
C
      DOUBLE PRECISION E0,PF_E,TSCAT,POL_BEAM,C12_SIG,C12_ASYM
      DOUBLE PRECISION PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD
      DOUBLE PRECISION ALPHA,AMOTT,COSTH2,SINTH2,TANTH2
      DOUBLE PRECISION EP,ETA,FF,HC,Q,Q4SQ
      DOUBLE PRECISION SIG,XMT,Z,A,Qeff
      DOUBLE PRECISION Q2MEV,XOM,Q3,PI 
C
      PI = 4*atan(1.0)
C
C----------------------------------------------------------------------
C     Constants Used In The Calculation
C----------------------------------------------------------------------
C
      HC=197.327                   ! Conversion constant in MeV fm
      ALPHA=1/137.035989           ! Fine Struction Constant 
      XMT=11177.928                ! 12C Nuclear Mass            
      A=12.                        ! 12C Atomic Number
      Z=6.                         ! 12C Charge 
C
C----------------------------------------------------------------------
C     Calculation of Kinematic Quantities
C----------------------------------------------------------------------
C

      ETA= 1./(1.+ 2.*E0/XMT*SIN(TSCAT/2)**2)
      EP= ETA*E0
      Q4SQ= 4.*E0*EP*SIN(TSCAT/2)**2
      Q=SQRT(Q4SQ)/HC           !q4 in fm**-1
      Q2MEV= (Q*HC)**2
      XOM= Q2MEV/(2.*XMT)             !in MeV
      Q3=sqrt(Q2MEV+XOM*XOM)
      SINTH2= SIN(TSCAT/2)**2
      COSTH2= COS(TSCAT/2)**2
      TANTH2= SINTH2/COSTH2
      AMOTT= 10000.*(HC*ALPHA/2./E0)**2* COSTH2/SINTH2**2
C      AMOTT= 10000.*(HC*ALPHA/2./E0)**2*cos(tscat/2)**2/
C     _sin(tscat/2)**4
C
C----------------------------------------------------------------------
C     Simple Calculation Of The Effective Q
C----------------------------------------------------------------------
C
            Qeff=Q3/HC*
     >     ( 1. + 1.5*Z*ALPHA*HC/E0/
     >        (1.29*2.4776)  )
C
C----------------------------------------------------------------------
C     Calculation Of The Form Factor  
C----------------------------------------------------------------------
C
      FF=0.
      call C12_FACTOR(FF,Qeff)
C
C----------------------------------------------------------------------
C     Calculation Of The Cross Section  
C----------------------------------------------------------------------
C
      SIG= AMOTT*ETA*FF**2 ! microbarns/sr'
      C12_SIG =SIG/10000.      !he3 cross section in fm^2/sr
C
C----------------------------------------------------------------------
C     Asymmetry Polarizations Set To Zero 
C----------------------------------------------------------------------
C
      C12_ASYM = 0.
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
C Formalism for Form-factor parmaeterization can be found in
C  B. Dreher et al., Nucl. Phys. A235 (1974) 219-248, eq 2.18
C    checked by R.J.Feuerbach, and made to better match the formalism.
C
C Useful reference, but not what was used for below.
C    Subroutine to generate the form factors of 3He
C    written by D.W.Higinbotham using the parameterization
C    A. Amroun et al., Nucl. Phys. A579 (1994) 596-626,
C    with the correct form of equation (1) from
C    I. Sick, Nucl. Phys. A218 (1974) 509.
C
C    
C----------------------------------------------------------------------
C
      SUBROUTINE C12_FACTOR(FF,Q)

      IMPLICIT NONE

      DOUBLE PRECISION Q,FBB(18),FF,R,PI
      INTEGER n

      DATA FBB/
     _1.5709E-2,
     _3.8610E-2,
     _3.6418E-2,
     _1.4293E-2,
     _-4.4628E-3,
     _-9.8420E-3,
     _-6.6518E-3,
     _-2.7066E-3,
     _-5.6697E-4,
     _-2.7453E-4,
     _-1.7093E-4,
     _1.2433E-4,
     _-4.8496E-5,
     _1.5675E-5,
     _-4.5194E-6,
     _1.1920E-6,
     _-2.9065E-7,
     _6.5845E-8/

      R  = 8.                  ! The cutoff radius in fm
      PI = 4*atan(1.0) 
      FF = 0.0

      do 100 n=1,18
         FF=FF+FBB(n)*(-1)**n/((Q*R)**2-(n*PI)**2)*sin(Q*R)
 100  continue
      FF=FF*4*pi*R**2/Q
      RETURN
      END



