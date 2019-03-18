C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C       SUBROUTINE deut_elastic
C
C       AUTHOR: S. Van Verst
C       DATE:   DEC-1991
C
C       PURPOSE:
C               Calculate cross section and polarizations for deuteron
C               elastic scattering.
C
C       MODIFICATIONS:
C
C
C
C----------------------------------------------------------------------
C
      SUBROUTINE DEUT_ELASTIC(E0,PF_E,TSCAT,POL_BEAM,DEUT_SIG,DEUT_ASYM,
     #                       PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
      IMPLICIT NONE
C
      COMMON /KINVAR/  KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
c      COMMON /DEUT_FF/ FC,FM,FQ,Q_2
      COMMON /N_DEUT_FF/ N_D_PTS
C
      DOUBLE PRECISION KI(3),Q(3),PF(3),PR(3),Q2,QMU2,EP,THETA_EP,CTH_PQ
      DOUBLE PRECISION E0,PF_E,TSCAT,POL_BEAM,MD,XMT
      DOUBLE PRECISION TTH2,STH2,ETA,QMU2_FM,RINTERPQ,CTH2
      DOUBLE PRECISION FF,G_C,G_M,G_Q,GC2,GM2,GQ2,A,B
      DOUBLE PRECISION S_MOTT,REC,DEUT_SIG,DEUT_ASYM
C
      DOUBLE PRECISION FETA,D_LT,D_LL
      DOUBLE PRECISION PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD
      DOUBLE PRECISION PI,ALPHA,SINTHW,GV,K,G_F,HBARC
C

      CHARACTER*200    FC_FILE,FM_FILE,FQ_FILE
      DOUBLE PRECISION FC(100),FQ(100),FM(100),Q_2(100)

      INTEGER N_D_PTS,J,ILUN
C
      PARAMETER (MD     = 1875.630)        !Deuteron mass in MeV
      PARAMETER (SINTHW = 0.230)           !sin^2(theta_W)
      PARAMETER (GV     = 0.080)           !1-4.0*sinthw
      PARAMETER (K      = 0.020)           !GV/4
      PARAMETER (ALPHA  = 7.29735E-3)      !fine structure constant
      PARAMETER (HBARC   = 197.3286)
      PARAMETER (PI     = 3.141592654)
      PARAMETER (G_F    = 1.16637E-11)     !Fermi constant

C        READ parameters 
C        see meceep physics.f for details
           DO J=1,100
              FC(J)  = 0.
              FM(J)  = 0.
              FQ(J)  = 0.
              Q_2(J) = 0.
           ENDDO
C
C ------------------------------------------------------------------------
C       Ask for deuteron form factor model (from Tjon), then read them
C       in from files in DAT$D
C ------------------------------------------------------------------------
C
c$$$           WRITE(6,'(a)')' Model: ia,iamec,rsc,rscmec (1,2,3,4) >'
c$$$           READ(5,*)IMODEL
c$$$           IF(IMODEL .EQ. 1)THEN
c$$$              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/ia.fc'
c$$$              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/ia.fq'
c$$$              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/ia.fm'
c$$$           ELSEIF(IMODEL .EQ. 2) THEN
c$$$              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/iamec.fc'
c$$$              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/iamec.fq'
c$$$              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/iamec.fm'
c$$$           ELSEIF(IMODEL .EQ. 3) THEN
c$$$              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rsc.fc'
c$$$              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rsc.fq'
c$$$              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rsc.fm'
c$$$           ELSE
              FC_FILE =  'SRC/rscmec.fc'
              FQ_FILE =  'SRC/rscmec.fq'
              FM_FILE =  'SRC/rscmec.fm'
c          ENDIF
C
           OPEN(UNIT=11,FILE=FM_FILE,STATUS='OLD')
           OPEN(UNIT=12,FILE=FC_FILE,STATUS='OLD')
           OPEN(UNIT=13,FILE=FQ_FILE,STATUS='OLD')
C
           DO J=1,100
              READ(11,*,END=9) Q_2(J),FM(J)
              READ(12,*) Q_2(J),FC(J)
              READ(13,*) Q_2(J),FQ(J)
              N_D_PTS = J
           ENDDO
C
9          DO ILUN = 11,13
              CLOSE(UNIT=ILUN)
           ENDDO
C
C----------------------------------------------------------------------
C     Calculate some kinematical quantities
C----------------------------------------------------------------------
C
      XMT=1875.612928 
      ETA= 1./(1.+ 2.*E0/XMT*SIN(TSCAT/2)**2)
      PF_E= ETA*E0 ! actually Ep
     
      QMU2= 4.*E0*PF_E*SIN(TSCAT/2)**2
      write(6,*) E0, PF_E, ETA,QMU2
      CTH2    = COS(TSCAT/2.D0)
      TTH2    = TAN(TSCAT/2.D0)
      STH2    = (SIN(TSCAT/2.D0))**2
      ETA     = QMU2/(4.D0 * MD*MD) ! actually tau
C
C----------------------------------------------------------------------
C     Get deuteron form factors
C----------------------------------------------------------------------
C
      QMU2_FM  = QMU2/HBARC**2
      G_C = RINTERPQ(Q_2,FC,QMU2_FM,N_D_PTS)
      G_Q = RINTERPQ(Q_2,FQ,QMU2_FM,N_D_PTS)*(MD/HBARC)**2
      G_M = RINTERPQ(Q_2,FM,QMU2_FM,N_D_PTS)

C      write(6,*) G_C,G_Q,G_M
C
C----------------------------------------------------------------------
C     Get Mott cross section and calculate unpolarized Deuteron
C     elastic cross section.
C----------------------------------------------------------------------
C
      CALL SMOTT(TSCAT,PF_E,QMU2,S_MOTT) !sigma_mott in fm^2/sr
      REC = PF_E/E0                      !recoil factor
      GC2 = G_C**2
      GM2 = G_M**2
      GQ2 = G_Q**2
      A   = GC2 + (2.D0*ETA/3.D0) * GM2 + (8.D0*ETA*ETA/9.D0) * GQ2
      B   = (4.D0*ETA/3.D0) * (1.D0+ETA) * GM2
      FF  = A + B*TTH2*TTH2
      DEUT_SIG = S_MOTT * FF * REC       !deuteron cross section in fm^2/sr
      write(6,*) S_MOTT, FF, REC, A, B, DEUT_SIG 

C
C----------------------------------------------------------------------
C     Asymmerty         set to zero
C
C     Polarizations     P_n set to zero
C----------------------------------------------------------------------
C
      DEUT_ASYM = 0.
C
      FETA = SQRT(ETA*(1.D0+ETA))
C
      D_LT = -(4.D0*FETA/3.D0) * G_M *
     #         (G_C + (ETA/3.D0) * G_Q) * TTH2 / FF
      D_LL = (2.D0*ETA/3.D0) *
     #         SQRT( (1.D0 + ETA) * (1.D0 + ETA*STH2) ) * GM2 / FF *
     #         TTH2/CTH2
C
      PN_HD = 0.D0
      PN_HI = 0.D0                   !no normal polarization for 1H elastic
      PT_HD = POL_BEAM*D_LT*DEUT_SIG  !cross section for transverse pol.
      PT_HI = 0.D0
      PL_HD = POL_BEAM*D_LL*DEUT_SIG  !cross section for longitudinal pol.
      PL_HI = 0.D0
C
      RETURN
      END

