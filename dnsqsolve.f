      PROGRAM TEST
C
C     DRIVER FOR DNSQE EXAMPLE.
C
      INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE
      DOUBLE PRECISION TOL,FNORM
      DOUBLE PRECISION X(9),FVEC(9),WA(180)
      DOUBLE PRECISION DENORM,D1MACH
      EXTERNAL FCN
      DATA NWRITE /6/
C
      IOPT = 2
      N = 9
C
C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
C
      DO 10 J = 1, 9
         X(J) = -1.E0
   10    CONTINUE

      LWA = 180
      NPRINT = 0
C
C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C     THIS IS THE RECOMMENDED SETTING.
C
      TOL = SQRT(D1MACH(4))
C
      CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
      FNORM = DENORM(N,FVEC)
      WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
      STOP
 1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
     *        5X,' EXIT PARAMETER',16X,I10 //
     *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
      END
      SUBROUTINE FCN(N,X,FVEC,IFLAG)
      INTEGER N,IFLAG
      DOUBLE PRECISION X(N),FVEC(N)
      INTEGER K
      DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
      DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C
      DO 10 K = 1, N
         TEMP = (THREE - TWO*X(K))*X(K)
         TEMP1 = ZERO
         IF (K .NE. 1) TEMP1 = X(K-1)
         TEMP2 = ZERO
         IF (K .NE. N) TEMP2 = X(K+1)
         FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
   10    CONTINUE
      RETURN
      END