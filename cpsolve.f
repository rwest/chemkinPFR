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
     3                   LOUT, ISTATE, FFIXEDMF)
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
      DIMENSION FFIXEDMF(KK)
      COMMON /ICONS/ IFIXEDMF(128), KK, NP, NWT, NH, NWDOT
      EXTERNAL FUN
      KK   = NKK
      NP   = NNP
      NWT  = NNWT
      NH   = NNH
      NWDOT= NNWDOT
      
C Store the fixed species
C FIRST SET THEM ALL TO ZERO
      DO 10 K = 1, 128
        IFIXEDMF(K) = 0
 10   CONTINUE
C THEN STORE 
      J = 1
      DO 20 K = 1, KK
        IF (FFIXEDMF(K) .NE. 0.0) THEN
            IFIXEDMF(J) = K
            J = J + 1
        END IF
 20   CONTINUE
      
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
      COMMON /ICONS/ IFIXEDMF(128), KK, NP, NWT, NH, NWDOT
      DIMENSION Z(*), ZP(*), RPAR(*), IPAR(*)
C
C     Variables in Z are:  Z(1)   = T
C                          Z(K+1) = Y(K)
C
C     Call CHEMKIN subroutines
C
C Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions. RHO
      CALL CKRHOY (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RHO)
C Returns the mean specific heat at constant pressure. CPB
      CALL CKCPBS (Z(1), Z(2), IPAR, RPAR, CPB)
C Returns the molar production rates of the species given pressure, temperature(s) and mass fractions. WDOT(*)
      CALL CKWYP  (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RPAR(NWDOT))
C Returns the enthalpies in mass units. HMS(*)
C isothermal so following line commented out:
C      CALL CKHMS  (Z(1), IPAR, RPAR, RPAR(NH))
C
C     Form governing equation
C
      SUM = 0.0
      DO 100 K = 1, KK
C         H    = RPAR(NH    + K - 1)
         WDOT = RPAR(NWDOT + K - 1)
         WT   = RPAR(NWT   + K - 1)
         ZP(K+1) = WDOT * WT / RHO
C isothermal so following line commented out:
C         SUM = SUM + H * WDOT * WT
 100  CONTINUE

C For species which have a nonzero FIXEDMF we set the creation rate to zero
      DO 200 K = 1, 128
         IF (IFIXEDMF(K) .NE. 0) ZP(IFIXEDMF(K)+1) = 0.0
 200  CONTINUE
      
C      ZP(1) = -SUM / (RHO*CPB)
      ZP(1) = 0.0
C
C     end of SUBROUTINE FUN
      RETURN
      END
