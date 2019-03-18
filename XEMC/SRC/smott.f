C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE SMOTT
C
C       AUTHOR: P.E. Ulmer
C       DATE:   V1.0    27-AUG-1990
C       MODIFICATIONS:
C               none
C       PURPOSE:
C               Calculates Mott cross section (in fm^2/sr).
C------------------------------------------------------------------------------
C
      SUBROUTINE SMOTT(SCAT_ANG,EFIN,QMU2,S_MOTT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (HBARC = 197.3286D0)
      PARAMETER (ALPHA = 7.297D-3)            !Fine structure constant
C
      S_MOTT = 4.D0*HBARC**2*((ALPHA*COS(SCAT_ANG/2.D0)*EFIN)/QMU2)**2
C
      RETURN
      END
