C***********************************************************************BEGIN-HDR
C
C  Copyright (c) 1997-2005 Reaction Design.  All rights reserved.
C
C  $Header: /chemkin/4.1/prebuild/samples/sample_apps_cpp/conp/cpsolve.f 3     12/21/05 8:43p Frupley $
C
C  $NoKeywords: $
C***********************************************************************END-HDR
      SUBROUTINE CPSOLVE(NKK, NNP, NNWT, NNH, NNWDOT, NEQ, Z, TT1,
     1                   TT2, ITOL, RTOL, ATOL, ITASK, IOPT, RVODE,
     2                   LRW, IVODE, LIW, MF, RWORK, IWORK,
     3                   LOUT, ISTATE)
C
C     Interface to DVODE and FUN from the C++ driver program
C
C*****EXCLUDE-IF precision_quad precision_single
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END EXCLUDE-IF precision_quad precision_single
C*****INCLUDE-IF precision_single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END INCLUDE-IF precision_single
      DIMENSION RVODE(LRW),IVODE(LIW),Z(NEQ),RWORK(*),IWORK(*)
      COMMON /ICONS/ KK, NP, NWT, NH, NWDOT
      EXTERNAL FUN
      KK   = NKK
      NP   = NNP
      NWT  = NNWT
      NH   = NNH
      NWDOT= NNWDOT
C*****INCLUDE-IF precision_single
C      CALL SVODE
C*****END INCLUDE-IF precision_single
C*****EXCLUDE-IF precision_quad precision_single
      CALL DVODE
C*****END EXCLUDE-IF precision_quad precision_single
     *           (FUN, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RVODE, LRW, IVODE,
     2            LIW, JAC, MF, RWORK, IWORK)
C
      IF (ISTATE .LE. -2) WRITE (LOUT,*) ' ISTATE=',ISTATE
      RETURN
      END

      SUBROUTINE FUN (N, TIME, Z, ZP, RPAR, IPAR)
C
C*****EXCLUDE-IF precision_quad precision_single
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
C*****END EXCLUDE-IF precision_quad precision_single
C*****INCLUDE-IF precision_single
C      IMPLICIT REAL (A-H,O-Z), INTEGER(I-N)
C*****END INCLUDE-IF precision_single
C
      COMMON /ICONS/ KK, NP, NWT, NH, NWDOT
      DIMENSION Z(*), ZP(*), RPAR(*), IPAR(*)
C
C     Variables in Z are:  Z(1)   = T
C                          Z(K+1) = Y(K)
C
C     Call CHEMKIN subroutines
C
      CALL CKRHOY (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RHO)
      CALL CKCPBS (Z(1), Z(2), IPAR, RPAR, CPB)
      CALL CKWYP  (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RPAR(NWDOT))
      CALL CKHMS  (Z(1), IPAR, RPAR, RPAR(NH))
C
C     Form governing equation
C
      SUM = 0.0
      DO 100 K = 1, KK
         H    = RPAR(NH    + K - 1)
         WDOT = RPAR(NWDOT + K - 1)
         WT   = RPAR(NWT   + K - 1)
         ZP(K+1) = WDOT * WT / RHO
         SUM = SUM + H * WDOT * WT
 100  CONTINUE
      ZP(1) = -SUM / (RHO*CPB)
C
C     end of SUBROUTINE FUN
      RETURN
      END
