C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       AUTHOR: B.H.Cottman
C       DATE:   V1.0    23-SEP-83
C       MODIFICATIONS:
C               09-SEP-1991   P.E. Ulmer
C                   Added calculation of proton and neutron axial form
C                   factors in dipole model.
C               14-MAR-2000   S. Dieterich/S. Strauch
C                   GEP (for FFTYPE='DIPOLE') now employs fit to
C                   JLAB-Hall A GEP/GMP ratio in E93027 (GMP is taken
C                   to be standard dipole form factor times MUP). 
C               20-SEP-2001   P. Ulmer
C                   GEP for above Hall A fit (E93027) now corresponds to
C                   FFTYPE='HALLA1' (GMP is still taken
C                   to be standard dipole form factor times MUP).
C                   The FFTYPE='DIPOLE' is now the standard dipole + Galster
C                   GEn.  Another Hall A fit (from Ed Brash)
C                   was added and corresponds to FFTYPE='HALLA2'.
C                   Other parameterizations added as well.
C       PURPOSE:
C               Calculate dipole, MAINZ fit or FIT 8.2 nucleon formfactors.
C       VARIABLES:
C       INPUT:
C               FFTYPE:  Prescription for form factors (dipole, Mainz, 8.2)
C               Q_MU:    4-momentum transfer (MeV/c)
C               M_EFF:   Factor to scale nucleon mass.
C
C       OUTPUT:
C               GEP:    Electric formfactor of proton
C               GMP:    Magnetic formfactor of proton
C               GEN:    Electric formfactor of neutron
C               GMN:    Magnetic formfactor of neutron
C               GAP:    Axial formfactor of proton (dipole model)
C               GAN:    Axial formfactor of neutron (dipole model)
C               F1P,N:  Dirac formfactor for proton and neutron
C               F2P,N:  Pauli formactor for proton and neutron
C       ROUTINES CALLED:
C                       None
C
C       Note that the rather hideous-looking expressions for the
C       dipole nucleon form factors reduce to the "more familiar"
C       forms (for M_EFF=1) after a little algebra:
C                GEP =          GDP
C                GMP =      MUP*GDP
C                GEN = -TAU*MUN*GDP/(1+5.6TAU)
C                GMN =      MUN*GDP
C       where GDP is the dipole form factor, MUP and MUN are the proton
C       and neutron magnetic moments and TAU = Q_mu^2/4M^2.  Note that
C       F1N and F2N in the dipole model are chosen so that the above
C       expressions for GEN and GMN are obtained.
C------------------------------------------------------------------------------
C
      SUBROUTINE FORM_FACT(FFTYPE,Q_MU,M_EFF,GEP,GMP,GEN,GMN,GAP,GAN,
     #                       F1P,F2P,F1N,F2N)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      CHARACTER*7 FFTYPE
C
      DOUBLE PRECISION F(2,2),G(2,2)
C
      PARAMETER (PI = 3.1415927)
      PARAMETER (MP = 938.2796)
      PARAMETER (MA   =  1030.D0)
      PARAMETER (HBARC = 197.3286)
C
      GDP(QQ)=1./(1.+QQ/0.71E+6)**2                           ! DIPOLE FIT
C
      F1RHO(T)=(.955+.09/(1.-T/.355)**2)/2./(1.-T/.536)       ! RHO MESON...
      F2RHO(T)=(5.335+.962/(1.-T/.268))/2./(1.-T/.603)        ! CONTRIBUTIONS
      F1V(T)=.05/(1.46-T)-.52/(6.-T)+.28/(8.7-T)
      F2V(T)=-1.99/(1.46-T)+.2/(6.-T)+.19/(8.7-T)
      F1S(T)=.71/(0.782**2-T)-.64/(1.02**2-T)-.13/(3.24-T)
      F2S(T)=-.11/(0.782**2-T)+.13/(1.02**2-T)-.02/(3.24-T)
C
      G_E_M(T)=0.312/(1+T/6.0)+                               !MAINZ FIT
     &               1.312/(1+T/15.02)-
     &               0.709/(1+T/44.08)+
     &               0.085/(1+T/154.2)
      G_M_M(T)=( 0.694/(1+T/8.5)+
     &                 0.719/(1+T/15.02)-
     &                 0.418/(1+T/44.08)+
     &                 0.005/(1+T/355.4) )
C
      Q_MU_GEV = Q_MU/1000.   !in GeV
      Q2_GEV = Q_MU_GEV**2    !in (GeV/c)^2
      Q2 = Q_MU**2
      Q2_F=Q2/HBARC**2
      MN = M_EFF*MP           !Nucleon mass
      MUN=-1.913              !Neutron magnetic moment
      MUP=2.793               !Proton magnetic moment
      KAPPAN=MUN
      KAPPAP=MUP-1.
      TAU=Q2/4./MP**2
      X=1.+TAU
      XX=  Q2/4.0/MN/MN
C
C       Dipole fit
C
      IF(FFTYPE .EQ. 'DIPOLE ') THEN
        F1P=(1.+TAU*MUP)*GDP(Q2)/X                     ! DIPOLE CHOSEN
        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        F2P=GDP(Q2)/X
        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        GEN=F1N-XX*KAPPAN*F2N
        GEP=F1P-XX*KAPPAP*F2P
        GMN=F1N+KAPPAN*F2N
        GMP=F1P+KAPPAP*F2P
C ..................................
      ELSEIF(FFTYPE .EQ. 'MAINZ  ') THEN
        F1P=(G_E_M(Q2_F)+(X-1.)*MUP*G_M_M(Q2_F))/X     ! MAINZ FIT
        F1N=(Q2*MUN/2./MP**2)*G_M_M(Q2_F)/X
        F2P=(MUP*G_M_M(Q2_F)-G_E_M(Q2_F))/(KAPPAP*X)
        F2N=(2.-X)*G_M_M(Q2_F)/X
        GEN=F1N-XX*KAPPAN*F2N
        GEP=F1P-XX*KAPPAP*F2P
        GMN=F1N+KAPPAN*F2N
        GMP=F1P+KAPPAP*F2P
C ..................................
      ELSEIF(FFTYPE .EQ. 'FIT 8.2') THEN
        TT=-Q2/1.E+6                                   ! FIT 8.2 CHOSEN
        F1P=F1S(TT)+F1V(TT)+F1RHO(TT)
        F2P=(F2S(TT)+F2V(TT)+F2RHO(TT))/(MUP-1.)
        F1N=F1S(TT)-F1V(TT)-F1RHO(TT)
        F2N=(F2S(TT)-F2V(TT)-F2RHO(TT))/MUN
        GEN=F1N-XX*KAPPAN*F2N
        GEP=F1P-XX*KAPPAP*F2P
        GMN=F1N+KAPPAN*F2N
        GMP=F1P+KAPPAP*F2P
C ..................................
      ELSEIF(FFTYPE .EQ. 'HALLA1 ') THEN               ! Hall A Fit 1
                                             ! Take GMP from dipole
                                             ! and take Hall A data for GEP
                                             ! and take Galster neutron ff's
        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
C S&S
        A = -0.45253
        B =  3.3915
        GEGM_RATIO = 1. + A*(Q2/1.E6)**3/(1.+B*(Q2/1.E6)**2)
        F1P=(GEGM_RATIO+TAU*MUP)*GDP(Q2)/X
        F2P=((MUP-GEGM_RATIO)/KAPPAP)*GDP(Q2)/X
        GEN=F1N-XX*KAPPAN*F2N
        GEP=F1P-XX*KAPPAP*F2P
        GMN=F1N+KAPPAN*F2N
        GMP=F1P+KAPPAP*F2P
C ..................................
      ELSEIF(FFTYPE .EQ. 'HALLA2 ') THEN               ! Hall A Fit 2
                                             ! GMP from E. Brash fit/extraction
                                             ! and take Hall A data for GEP
                                             ! and take Galster neutron ff's
        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        GEN=F1N-XX*KAPPAN*F2N
        GMN=F1N+KAPPAN*F2N
        GMP=MUP/( 1. + 0.35681*Q_MU_GEV    + 2.0126 *Q_MU_GEV**2
     &               + 1.3052 *Q_MU_GEV**3 + 0.46952*Q_MU_GEV**4
     &               + 0.43753*Q_MU_GEV**5 )
        A = -0.45253
        B =  3.3915
        GEGM_RATIO = 1. + A*(Q2/1.E6)**3/(1.+B*(Q2/1.E6)**2)
        GEP=GEGM_RATIO*GMP/MUP
        F1P=(GEP+XX*GMP)/X
        F2P=(GMP-GEP)/(KAPPAP*X)
C ..................................
c$$$      ELSEIF(FFTYPE .EQ. 'HALLA3 ') THEN               ! Hall A Fit 3
c$$$                                             ! Take GMP from MMD alt. 3-pole
c$$$                                             ! and take Hall A data for GEP
c$$$                                             ! and take Galster neutron ff's
c$$$        CALL MMD(Q2_GEV,2,F)   ! isoscalar/isovector Dirac ff's
c$$$        CALL FTOG_PEU(Q2_GEV,F,GEP,GEN,GMP,GMN) ! Get Sachs ff's
c$$$        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
c$$$        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
c$$$        GEN=F1N-XX*KAPPAN*F2N
c$$$        GMN=F1N+KAPPAN*F2N
c$$$        A = -0.45253
c$$$        B =  3.3915
c$$$        GEGM_RATIO = 1. + A*(Q2/1.E6)**3/(1.+B*(Q2/1.E6)**2)
c$$$        GEP=GEGM_RATIO*GMP/MUP
c$$$        F1P=(GEP+XX*GMP)/X
c$$$        F2P=(GMP-GEP)/(KAPPAP*X)
C ..................................
      ELSEIF(FFTYPE .EQ. 'HALLA4 ') THEN               ! Hall A Fit 4
                                             ! Take GMP from Simon (Mainz)
                                             ! and take Hall A data for GEP
                                             ! and take Galster neutron ff's
        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
        GEN=F1N-XX*KAPPAN*F2N
        GMN=F1N+KAPPAN*F2N
        GMP = MUP*G_M_M(Q2_F)
        A = -0.45253
        B =  3.3915
        GEGM_RATIO = 1. + A*(Q2/1.E6)**3/(1.+B*(Q2/1.E6)**2)
        GEP=GEGM_RATIO*GMP/MUP
        F1P=(GEP+XX*GMP)/X
        F2P=(GMP-GEP)/(KAPPAP*X)
C ..................................
c$$$      ELSEIF(FFTYPE .EQ. 'MMD    ') THEN               ! MMD alt. 3-pole
c$$$        CALL MMD(Q2_GEV,2,F)   ! isoscalar/isovector Dirac ff's
c$$$        CALL FTOG_PEU(Q2_GEV,F,GEP,GEN,GMP,GMN) ! Get Sachs ff's
c$$$        F1N=(GEN+XX*GMN)/X
c$$$        F2N=(GMN-GEN)/(KAPPAN*X)
c$$$        F1P=(GEP+XX*GMP)/X
c$$$        F2P=(GMP-GEP)/(KAPPAP*X)
c$$$C ..................................
c$$$      ELSEIF(FFTYPE .EQ. 'JRA660 ') THEN               ! JohnA sigma fit
c$$$        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
c$$$        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
c$$$        GEN=F1N-XX*KAPPAN*F2N
c$$$        GMN=F1N+KAPPAN*F2N
c$$$        GMP=MUP/(1.D0
c$$$     &            + Q_MU_GEV    * (-0.61313)
c$$$     &            + Q_MU_GEV**2 * (7.45299)
c$$$     &            + Q_MU_GEV**3 * (-8.98875)
c$$$     &            + Q_MU_GEV**4 * (9.10233)
c$$$     &            + Q_MU_GEV**5 * (-2.79047)
c$$$     &            + Q_MU_GEV**6 * (0.44317)  )
c$$$        GEP=1.D0/(1.D0
c$$$     &            + Q_MU_GEV    * (1.25084)
c$$$     &            + Q_MU_GEV**2 * (-3.94200)
c$$$     &            + Q_MU_GEV**3 * (14.37442)
c$$$     &            + Q_MU_GEV**4 * (-12.35122)
c$$$     &            + Q_MU_GEV**5 * (6.81231)
c$$$     &            + Q_MU_GEV**6 * (-1.19725)  )
c$$$        F1P=(GEP+XX*GMP)/X
c$$$        F2P=(GMP-GEP)/(KAPPAP*X)
c$$$C ..................................
c$$$      ELSEIF(FFTYPE .EQ. 'JRA661 ') THEN          ! JohnA sigma/poltran fit
c$$$        F1N=5.6*TAU*TAU*MUN*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
c$$$        F2N=(1.+6.6*TAU)*GDP(Q2)/((1.+TAU)*(1.+5.6*TAU))
c$$$        GEN=F1N-XX*KAPPAN*F2N
c$$$        GMN=F1N+KAPPAN*F2N
c$$$        GMP=MUP/(1.D0
c$$$     &            + Q_MU_GEV    * (-0.46545)
c$$$     &            + Q_MU_GEV**2 * (6.37617)
c$$$     &            + Q_MU_GEV**3 * (-6.67047)
c$$$     &            + Q_MU_GEV**4 * (6.97610)
c$$$     &            + Q_MU_GEV**5 * (-1.95558)
c$$$     &            + Q_MU_GEV**6 * (0.31844)  )
c$$$        GEP=1.D0/(1.D0
c$$$     &            + Q_MU_GEV    * (-0.22313)
c$$$     &            + Q_MU_GEV**2 * (6.69319)
c$$$     &            + Q_MU_GEV**3 * (-13.76028)
c$$$     &            + Q_MU_GEV**4 * (22.74632)
c$$$     &            + Q_MU_GEV**5 * (-14.29633)
c$$$     &            + Q_MU_GEV**6 * (3.98071)  )
c$$$        F1P=(GEP+XX*GMP)/X
c$$$        F2P=(GMP-GEP)/(KAPPAP*X)
c$$$C ..................................
c$$$      ELSEIF(FFTYPE .EQ. 'LOMON  ') THEN               ! Lomon fit
c$$$        CALL LOMONFF(4,Q2_GEV,F1P,F2P,F1N,F2N,GEP,GMP,GEN,GMN)
c$$$C
      ENDIF
C
C ---------------------------------------------------------------------
C     Calculate neutron and proton axial vector form factors
C     in dipole model.
C ---------------------------------------------------------------------
C
      GAP = -(1.262D0/2.0D0)/(1.0D0 + Q2/MA**2)**2
      GAN = - GAP    !Assume isovector character for axial vector FF
C
      RETURN
      END
c
c ---------------------------------------------------------------------
c     SUBROUTINE ftog_peu
c
c     Converts form factors from Isoscalar/Isovector Dirac to
c     Neutron/Proton Sachs.  This assumes that the anomalous moments
c     have already been factored into the input f's.
c ---------------------------------------------------------------------
c
      SUBROUTINE ftog_peu(q2,f,gep,gen,gmp,gmn)
c
      implicit none
c
      double precision f(2,2),fn(0:1,2),gep,gen,gmp,gmn
      double precision mn,q2,tau
      integer          i
c
      parameter (mn=0.9389)            ! nucleon mass (GeV)
c
      tau = q2/(4.d0*mn*mn)
c
      do i=1,2
         fn(0,i)=0.5d0*(f(2,i)+f(1,i))   ! proton  F1 (i=1) and F2 (i=2)
         fn(1,i)=0.5d0*(f(2,i)-f(1,i))   ! neutron F1 (i=1) and F2 (i=2)
      enddo
c
      gep = fn(0,1)-tau*fn(0,2)          ! proton  electric
      gen = fn(1,1)-tau*fn(1,2)          ! neutron electric
      gmp = fn(0,1)+fn(0,2)              ! proton  magnetic
      gmn = fn(1,1)+fn(1,2)              ! neutron magnetic
c
      return
      end



