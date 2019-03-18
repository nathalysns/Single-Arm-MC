C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C       SUBROUTINE HYD_ELASTIC
C
C       AUTHOR: P.E. Ulmer
C       DATE:   6-SEP-1991
C
C       PURPOSE:
C               Calculate cross section and polarizations for proton
C               elastic scattering.
C
C       MODIFICATIONS:
C               09-SEP-1991   P.E. Ulmer/J.M. Finn
C               Added calculation of parity violating asymmetry.
C
C----------------------------------------------------------------------
C
      SUBROUTINE HYD_ELASTIC(E0,PF_E,TSCAT,POL_BEAM,HYD_SIG,HYD_ASYM,
     #                       PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
      IMPLICIT NONE
C
      COMMON /KINVAR/ KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      COMMON /FF_MOD/ FFTYPE
C
      DOUBLE PRECISION KI(3),Q(3),PF(3),PR(3),Q2,QMU2,EP,THETA_EP,CTH_PQ
      DOUBLE PRECISION E0,PF_E,TSCAT,POL_BEAM,MP,ETA
      DOUBLE PRECISION TTH2,TAU,QMU,EPSILON
      DOUBLE PRECISION FF,GE2,GM2,F1P,F2P,F1N,F2N
      DOUBLE PRECISION G_EP,G_EN,G_MP,G_MN,G_AP,G_AN,MU_N
      DOUBLE PRECISION S_MOTT,REC,HYD_SIG,HYD_ASYM
      DOUBLE PRECISION FTAU,FACT0,D_LT,D_LL
      DOUBLE PRECISION PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD
      DOUBLE PRECISION ASYM_FACT, A_P_NUM, A_P_DENOM, ASYM
      DOUBLE PRECISION PI,ALPHA,SINTHW,GV,K,KA,G_F
      CHARACTER*7 FFTYPE
C
      PARAMETER (MP     = 938.2796D0)        !Proton mass in MeV
      PARAMETER (MU_N   = -1.913)            !neutron magnetic moment
      PARAMETER (SINTHW = 0.230D0)           !sin^2(theta_W)
      PARAMETER (GV     = 0.080D0)           !1-4.0*sinthw
      PARAMETER (K      = 0.020D0)           !GV/4
      PARAMETER (ALPHA  = 7.29735D-3)        !fine structure constant
      PARAMETER (PI     = 3.141592654D0)
      PARAMETER (G_F    = 1.16637E-11)       !Fermi constant
C
C----------------------------------------------------------------------
C     Calculate some kinematical quantities.  EPSILON is the
C     longitudinal polarization of the virtual photon.
C----------------------------------------------------------------------
C
C     added 2018.4
      ETA= 1./(1.+ 2.*E0/MP*SIN(TSCAT/2)**2)
      PF_E= ETA*E0 ! EP
      QMU2= 4.*E0*EP*SIN(TSCAT/2)**2

     
      TTH2 = TAN(TSCAT/2.D0)
      TAU = QMU2/(4.D0*MP*MP)
      QMU = SQRT(QMU2)
      EPSILON = 1.D0/(1.D0 + 2.D0 * (1.D0 + TAU) * TAN(TSCAT/2.D0)**2)
C
C----------------------------------------------------------------------
C     Get form factors.
C----------------------------------------------------------------------
C
      FFTYPE = 'HALLA4'         ! temp choice, should not matter so much
      CALL FORM_FACT(FFTYPE,QMU,1.D0,G_EP,G_MP,G_EN,G_MN,G_AP,G_AN,
     #                       F1P,F2P,F1N,F2N)
C
C----------------------------------------------------------------------
C     Get Mott cross section and calculate unpolarized 1H elastic
C     cross section.
C----------------------------------------------------------------------
C
      CALL SMOTT(TSCAT,PF_E,QMU2,S_MOTT) !sigma_mott in fm^2/sr
      REC = PF_E/E0                      !recoil factor
      GE2 = G_EP**2
      GM2 = G_MP**2
      FF  = 2.D0*TAU*GM2*TTH2**2 + (GE2+TAU*GM2)/(1.D0+TAU)
      HYD_SIG = S_MOTT*FF*REC            !1H cross section in fm^2/sr
C
C----------------------------------------------------------------------
C     Calculate parity-violating asymmetry.
C----------------------------------------------------------------------
C
      ASYM_FACT = -(G_F * QMU2)/(PI * ALPHA * SQRT(2.0D0))
      KA        =  -0.5D0*GV*SQRT((1.0D0-EPSILON**2)*TAU*(1.D0 + TAU))
      A_P_NUM   = EPSILON*G_EP*(K*G_EP-G_EN/4.D0)+TAU*G_MP*(K*G_MP-G_MN
     #              /4.D0)
      A_P_NUM   = A_P_NUM + KA * G_AP * G_MP
      A_P_DENOM = EPSILON*G_EP**2+TAU*G_MP**2
      ASYM      = ASYM_FACT*A_P_NUM/A_P_DENOM
C
C----------------------------------------------------------------------
C     Calculate the part of the cross section which arises from a
C     parity violating asymmetry (see J. Napolitano, Phys. Rev. C 43,
C     1473 (1991).)
C----------------------------------------------------------------------
C
      HYD_ASYM = POL_BEAM*HYD_SIG*ASYM
C
C----------------------------------------------------------------------
C     Calculate longitudinal and transverse proton polarization
C     components.
C
C     Note:  FACT0 = I_0 from Arnold, Carlson and Gross,
C            Phys. Rev. C 23, 363 (1981).  It seems a bit odd that
C            they pull out this factor which is (1+tau)FF rather than
C            FF.  In any case the cross section is given by:
C
C       sigma_pol = sig_unpol x (1 + h*A + pL + pT)
C                 = sig_mott x rec x FF x (1 + h*A + pL + pT)
C                 = sig_mott x rec x FACT0/(1+tau) x (1 + h*A + pL + pT)
C
C            where h is the beam polarization and A is the
C            parity-violating asymmetry.
C----------------------------------------------------------------------
C
      FTAU = SQRT(TAU*(1.D0+TAU))
      FACT0 = (1.D0+TAU)*FF
      D_LT = -2.D0*FTAU*SQRT(GE2*GM2)*TTH2/FACT0    !pol. transf. coeff.
      D_LL = ((E0+PF_E)/MP)*FTAU*GM2*TTH2**2/FACT0  !pol. transf. coeff.
C
      PN_HD = 0.D0
      PN_HI = 0.D0                   !no normal polarization for 1H elastic
      PT_HD = POL_BEAM*D_LT*HYD_SIG  !cross section for transverse pol.
      PT_HI = 0.D0
      PL_HD = POL_BEAM*D_LL*HYD_SIG  !cross section for longitudinal pol.
      PL_HI = 0.D0
C
      RETURN
      END
