C***********************************************************************BEGIN-HDR
C  Solver written by Richard West rwest@mit.edu November 2010
C  Based in part off the examples in 
C   SLATEC's dnsqe.f 
C  and CHEMKIN's 
C   /chemkin/4.1/prebuild/samples/sample_apps_cpp/conp/cpsolve.f
C***********************************************************************END-HDR
 
      SUBROUTINE DNSQSOLVE(N, PPRES, TTEMP, X, 
     1                   IINITRO, FFIXEDMF, 
     2                   LRW, LIW, RRWORK, IIWORK,
     3                   LOUT, ISTATE )
C
C      DNSQSOLVE(iSpeciesCount, dPres, dTemp, dMoleFractions,
C              iNitro, dFixedMoleFractions, iCKsizeD, iCKsizeI, dCKwork, iCKwork,
C              iOutputfileUnit, iState);
C
C     Interface from the C++ driver program
C
C     Input and output variables
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      INTEGER LRW, LIW
      DIMENSION RRWORK(LRW), IIWORK(LIW)
      DOUBLE PRECISION FFIXEDMF(N)
C     Other variables
      DOUBLE PRECISION WA((3*N**2 + 13*N)/2 + 1)
      DOUBLE PRECISION TOL,FNORM
      DOUBLE PRECISION X(N),FVEC(N), XLOG(N)
      DOUBLE PRECISION DENORM,D1MACH
      INTEGER J,K,IOPT,NPRINT,INFO,LWA,NWRITE
      EXTERNAL FCN
      
C ****** START OF COMMON BLOCK INITIALIZATION
      COMMON /RCKBLK/ RWORK(500000),FIFIXEDMF(128),PRES,TEMP
      COMMON /ICBBLK/ IWORK(500000),IFIXEDMF(128),INITRO,LLOUT
C     Copy the chemkin work arrays into my common block
      IF (LRW .GT. 500000 .OR. LIW .GT. 500000) THEN
        WRITE (LOUT,*) ' Chemkin work arrays are too large!'
        WRITE (LOUT,*) ' Edit dnsqsolve.f, recompile, and try again.'
        ISTATE = -99
        RETURN
      END IF
      DO K = 1, LRW
        RWORK(K) = RRWORK(K)
      END DO
      DO K = 1, LIW
        IWORK(K) = IIWORK(K)
      END DO
C     Copy other variables into my common block
      PRES = PPRES
      TEMP = TTEMP
      INITRO = IINITRO
      LLOUT = LOUT
C Store the fixed species (AT MOST 128)
C FIRST SET THEM ALL TO ZERO
      DO K = 1, 128
        IFIXEDMF(K) = 0
        FIFIXEDMF(K) = 0.0
      END DO
C THEN STORE INDICES OF THOSE SPECIES THAT SHOULD BE FIXED IN IFIXEDMF
C AND THEIR FIXED VALUES IN FIFXEDMF
      J = 1
      DO K = 1, N
        IF (FFIXEDMF(K) .NE. 0.0) THEN
            IFIXEDMF(J) = K
            FIFIXEDMF(J) = FFIXEDMF(K)
            J = J + 1
        END IF
C   Also fix species whose initial mole fractions (at this stage) are zero!!
C   These should only be inert species like Ar and He, but this assumes
C   that the initial guess is well formed and has no other zeros in.
        IF (X(K) .EQ. 0.0) THEN
            IFIXEDMF(J) = K
            FIFIXEDMF(J) = 0.0
            J = J + 1
        END IF
      END DO
      WRITE(LOUT,*) 'Stored ',J-1,' fixed mole fractions'
C ****** END OF COMMON BLOCK INITIALIZATION

      WRITE(LOUT,*) 'If INITRO=',INITRO,' then X(N2) = ',X(INITRO)
      
      DO K = 1, N
        XLOG(K) = LOG10(X(K))
      END DO

      ISTATE = 0
C       WA is a work array of length LWA.
C       LWA is a positive integer input variable not less than
C         (3*N**2+13*N))/2.
      LWA = (3*N**2 + 13*N)/2 + 1
C
      IOPT = 2
      NPRINT = 1
C
C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C     THIS IS THE RECOMMENDED SETTING.
C
      TOL = SQRT(D1MACH(4))
      TOL = TOL * 2.0
C
      CALL DNSQE(FCN,JAC,IOPT,N,XLOG,FVEC,TOL,NPRINT,INFO,WA,LWA)
      FNORM = DENORM(N,FVEC)
      
      DO K = 1, N
        X(K) = 10.0 ** XLOG(K)
      END DO
      
      WRITE (LOUT,1000) FNORM,INFO,(X(J),J=1,N)
      WRITE(LOUT,1001) (FVEC(J),J=1,N)
      
 1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
     *        5X,' EXIT PARAMETER',16X,I10 //
     *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
 1001 FORMAT (5X,' FINAL RESIDUALS' // (5X,3E15.7))
     
      RETURN
      END
C **************************
C       THE SUBROUTINE THAT WE'RE TRYING TO SOLVE:
C       
      SUBROUTINE FCN(N,XLOG,FVEC,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      COMMON /RCKBLK/ RWORK(500000),FIFIXEDMF(128),PRES,TEMP
      COMMON /ICBBLK/ IWORK(500000),IFIXEDMF(128),INITRO,LLOUT
      INTEGER N,IFLAG,K
      DOUBLE PRECISION XLOG(N),FVEC(N), X(N),SUM
      LOUT = LLOUT

      DO K = 1, N
        X(K) = 10.0 ** XLOG(K)
      END DO
C      IF (IFLAG .EQ. 0) THEN
C         WRITE(LLOUT, 1002) (X(J),J=1,N)
C 1002 FORMAT (5X,' CURRENT SOLUTION' // (4X,4E15.7))
C      END IF
C     We want the solution FVEC=0
C
C     Returns the molar production rates of the species given pressure,
C     temperature(s) and mole fractions. Result returned in FVEC.
      CALL CKWXP (PRES, TEMP, X, IWORK, RWORK, FVEC)
      
      DO K = 1, N
        FVEC(K) = FVEC(K) * 1E6
      END DO
      
C For species which have a nonzero FIXEDMF we set the residual to the 
C difference between its value and its FIXED value.
      DO 200 K = 1, 128
         J = IFIXEDMF(K)
         IF (J .NE. 0) THEN
           FVEC(J) = LOG10(X(J) / FIFIXEDMF(K))
C          Fix for things being set to zero
           IF (FIFIXEDMF(K) .EQ. 0.0) FVEC(J) = X(J) 
         ELSE 
C           WRITE(LOUT,*) 'Fixed ', K-1,' residuals.'
           GOTO 201
         END IF
 200  CONTINUE
 201  CONTINUE
 
C For Nitrogen, the equation we solve is that the sum of everything equals 1
      SUM = 0.0
      DO K = 1, N
        SUM = SUM + X(K)
      END DO
      FVEC(INITRO) = LOG10(1.0 / SUM)
C      WRITE(LOUT,*) '1.0 - SUM = ',FVEC(INITRO),'  N2 = ',X(INITRO)

      IF (IFLAG .EQ. 0) THEN
        WRITE(LOUT,*) ' CURRENT NORM OF THE RESIDUAL = ', DENORM(N,FVEC)
C        WRITE(LOUT,*) ' CURRENT X(N2) = ', X(INITRO)
      END IF
      
      RETURN
      END