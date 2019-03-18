C     
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Linear Interpolation Routine:
C     
C     Calculates Y(X0) given array Y(X) assuming a linear
C     dependence:  Y = AX + B
C     
C     This routine assumes the X values are arranged in
C     ascending order but does not assume a uniform X spacing
C     of the array.
C------------------------------------------------------------------------------
C     
      DOUBLE PRECISION FUNCTION RINTERP(X,Y,X0,NBINS)
      IMPLICIT NONE
C     
C     INCLUDE 'summary.cmn'
C     
      DOUBLE PRECISION X(100),Y(100),X0,DET,A,B
      INTEGER NBINS,I1,I2,I

     
  
C     
C      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
CXXX  #     TYPE *,' Extrapolating outside range:  X = ',X0
C     #     I_EXTRAP = I_EXTRAP + 1
C     
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-1
         I2 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C     
 1    DET = X(I1)-X(I2)
      A = (Y(I1)-Y(I2))/DET
      B = (X(I1)*Y(I2)-X(I2)*Y(I1))/DET
C     
      RINTERP = A*X0 + B
C     
      RETURN
      END
C     
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Quadratic Interpolation Routine:
C     
C     Calculates Y(X0) given array Y(X) assuming a quadratic
C     dependence:  Y = AX^2 + BX + C
C     
C     This routine assumes the X values are arranged in
C     ascending order but does not assume a uniform X spacing
C     of the array.
C------------------------------------------------------------------------------
C     
      DOUBLE PRECISION FUNCTION RINTERPQ(X,Y,X0,NBINS)
      IMPLICIT NONE
C     
C     INCLUDE 'summary.cmn'
C     
      DOUBLE PRECISION X(100),Y(100),X0,DET,A,B,C
      INTEGER NBINS,I1,I2,I3,I
C     
C      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
CXXX  #     TYPE *,' Extrapolating outside range:  X = ',X0
C     #     I_EXTRAP = I_EXTRAP + 1
C     

      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C     
 1    DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y(I1)*(X(I2)-X(I3)) - X(I1)*(Y(I2)-Y(I3)) + Y(I2)*X(I3)
     #- X(I2)*Y(I3) )/DET
      B = ( X(I1)**2*(Y(I2)-Y(I3)) - Y(I1)*(X(I2)**2-X(I3)**2)
     #+ X(I2)**2*Y(I3) - X(I3)**2*Y(I2) )/DET
      C = ( X(I1)**2*(X(I2)*Y(I3)-X(I3)*Y(I2))
     #- X(I1)*(X(I2)**2*Y(I3)-X(I3)**2*Y(I2))
     #+ Y(I1)*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C     
      RINTERPQ = A*X0**2 + B*X0 + C
C     
      RETURN
      END
C     
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Exponential Interpolation Routine:
C     
C     Calculates Y(X0) given array Y(X) assuming the exponential
C     form:  Y = EXP(AX^2 + BX + C)
C     
C     This routine assumes the X values are arranged in
C     ascending order but does not assume a uniform X spacing
C     of the array.
C------------------------------------------------------------------------------
C     
      DOUBLE PRECISION FUNCTION RINTERPEXP(X,Y,X0,NBINS)
      IMPLICIT NONE
C     
C     INCLUDE 'summary.cmn'
C     
      DOUBLE PRECISION X(100),Y(100),X0,DET,A,B,C,Y1,Y2,Y3
      INTEGER NBINS,I1,I2,I3,I
C     
C      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
C     XXX #     TYPE *,' Extrapolating outside range:  X = ',X0
C     #     I_EXTRAP = I_EXTRAP + 1
C     
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C     
C     ----------------------------------------------------------------------
C     If all three Y-values are > 0, perform quadratic interpolation
C     on their logarithms; otherwise return zero as the result.
C     ----------------------------------------------------------------------
C     
 1    IF(Y(I1).GT.0.D0 .AND. Y(I2).GT.0.D0 .AND. Y(I3).GT.0.D0) THEN
         Y1 = LOG(Y(I1))
         Y2 = LOG(Y(I2))
         Y3 = LOG(Y(I3))
      ELSE
CXXX  TYPE *,' RINTERPEXP:  non-positive y-value; setting to zero '
C     I_INTERPEXP_NP = I_INTERPEXP_NP + 1
         RINTERPEXP = 0.D0
         RETURN
      ENDIF
C     
      DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y1*(X(I2)-X(I3)) - X(I1)*(Y2-Y3) + Y2*X(I3)
     #- X(I2)*Y3 )/DET
      B = ( X(I1)**2*(Y2-Y3) - Y1*(X(I2)**2-X(I3)**2)
     #+ X(I2)**2*Y3 - X(I3)**2*Y2 )/DET
      C = ( X(I1)**2*(X(I2)*Y3-X(I3)*Y2)
     #- X(I1)*(X(I2)**2*Y3-X(I3)**2*Y2)
     #+ Y1*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C     
      RINTERPEXP = EXP(A*X0**2 + B*X0 + C)
C     
      RETURN
      END
C     
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     2-Dimensional Linear Interpolation Routine:
C     
C     Calculates F(X0,Y0) given array F(X,Y) by two successive
C     interpolations, first in X and then in Y.
C     
C     Assumes uniform spacing of array in X and Y.
C------------------------------------------------------------------------------
C     
      DOUBLE PRECISION FUNCTION RINTERP2D(X,Y,F,X0,Y0,NX,NY,DELX,DELY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
C     INCLUDE 'summary.cmn'
C     
      DIMENSION X(200),Y(200),F(200,200)
C     
      I = DINT((X0-X(1))/DELX) + 1
      J = DINT((Y0-Y(1))/DELY) + 1
      IF((I+1.GT.NX).OR.(I.LT.1).OR.(J+1.GT.NY).OR.(J.LT.1))THEN
C     I_EXTRAP = I_EXTRAP + 1
         RINTERP2D = 0.D0
         RETURN
      ENDIF
C     
      RINTX1 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J)-F(I,J))
     #+ F(I,J)
      RINTX2 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J+1)-F(I,J+1))
     #+ F(I,J+1)
      RINTERP2D = ((Y0-Y(J))/(Y(J+1)-Y(J)))*(RINTX2-RINTX1) + RINTX1
C     
      RETURN
      END
