/***********************************************************************BEGIN-HDR
*
*  Copyright (c) 1997-2006 Reaction Design.  All rights reserved.
*
*  $Header: /chemkin/4.1/prebuild/samples/sample_apps_cpp/conp/conp.hpp 7     12/21/05 8:43p Frupley $
*
*  $NoKeywords: $
************************************************************************END-HDR
*/
#include <string.h>
#include <ctype.h>
#include <chemkinbase.h>

#if !defined( __CONP_H__ )
#define __CONP_H__

// C++ function prototypes
int resetFixedMoleFractions( int iSpeciesCount, 
							double *dMoleFractions, double *dFixedMoleFractions );
int resetN2 (char *sOutputfileName, int iOutputfileUnit,
			 int iSpeciesCount, int iStringLength,
			 char *sSpeciesNames, double *dMoleFractions);
int cpsize (char *sOutputfileName, int iOutputfileUnit,
            char *sGasLinkfileName, int iGasLinkfileUnit,
            int &iCKsizeI, int &iCKsizeD, int &iCKsizeS,
            int &iSpeciesCount, int &iStringLength);
int cpinp (char *sOutputfileName, int iOutputfileUnit,
           char *sInputfileName,
           int *iCKwork, double *dCKwork,
           int iSpeciesCount, int iStringLength,
           char *sSpeciesNames, double *dMoleFractions, double *dFixedMoleFractions,
           double &dPres, double &dTemp, double &dTend, double &dTdelta);
int cpout (char *sOutputfileName, int iOutputfileUnit,
           int *iCKwork, double *dCKwork, double *dSolution,
           int iSpeciesCount, char *sSpeciesNames, double *dMoleFractions,
           int iStringLength, double dTemp, double dTime);

// Fortran function prototypes
#if !defined( FORTRAN_FROM_CPLUSPLUS )
#include <fortran.h>    // include the cppf77 headers
#endif // !defined( FORTRAN_FROM_CPLUSPLUS )

#if !defined( F77_STUB_REQUIRED ) //This is for Win32 only
SUBROUTINE CPSOLVE(INTEGER& NKK, INTEGER& NNP, INTEGER& NNWT, INTEGER& NNH, INTEGER& NNWDOT,
                   INTEGER& NEQ, DOUBLE_PRECISION* Z, DOUBLE_PRECISION& TT1,
                   DOUBLE_PRECISION& TT2, INTEGER& ITOL, DOUBLE_PRECISION& RTOL,
                   DOUBLE_PRECISION& ATOL, INTEGER& ITASK,
                   INTEGER& IOPT, DOUBLE_PRECISION* RVODE, INTEGER& LRW,
                   INTEGER* IVODE, INTEGER& LIW, INTEGER& MF,
                   DOUBLE_PRECISION* RWORK, INTEGER* IWORK, INTEGER& LOUT,
                   INTEGER& ISTATE);
#else // F77_STUB_REQUIRED for unix platforms
#if defined ( F77_APPEND_UNDERSCORE ) //some unix platforms require appended underscore
SUBROUTINE_F77 cpsolve_ (int&, int&, int&, int&, int&,
                         int&, double*, double&,
                         double&, int&, double&,
                         double&, int&,
                         int&, double*, int&,
                         int*, int&, int&,
                         double*, int*, int&,
                         int&);
SUBROUTINE CPSOLVE(INTEGER& NKK, INTEGER& NNP, INTEGER& NNWT, INTEGER& NNH, INTEGER& NNWDOT,
                   INTEGER& NEQ, DOUBLE_PRECISION* Z, DOUBLE_PRECISION& TT1,
                   DOUBLE_PRECISION& TT2, INTEGER& ITOL, DOUBLE_PRECISION& RTOL,
                   DOUBLE_PRECISION& ATOL, INTEGER& ITASK,
                   INTEGER& IOPT, DOUBLE_PRECISION* RVODE, INTEGER& LRW,
                   INTEGER* IVODE, INTEGER& LIW, INTEGER& MF,
                   DOUBLE_PRECISION* RWORK, INTEGER* IWORK, INTEGER& LOUT,
                   INTEGER& ISTATE)
{   cpsolve_ (NKK, NNP, NNWT, NNH, NNWDOT, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL, ITASK,
              IOPT, RVODE, LRW, IVODE, LIW, MF,
              RWORK, IWORK, LOUT, ISTATE); }

#else // !defined( F77_APPEND_UNDERSCORE ); some unix platforms don't use underscore
SUBROUTINE_F77 cpsolve (int&, int&, int&, int&, int&,
                         int&, double*, double&,
                         double&, int&, double&,
                         double&, int&,
                         int&, double*, int&,
                         int*, int&, int&,
                         double*, int*, int&,
                         int&);
SUBROUTINE CPSOLVE(INTEGER& NKK, INTEGER& NNP, INTEGER& NNWT, INTEGER& NNH, INTEGER& NNWDOT,
                   INTEGER& NEQ, DOUBLE_PRECISION* Z, DOUBLE_PRECISION& TT1,
                   DOUBLE_PRECISION& TT2, INTEGER& ITOL, DOUBLE_PRECISION& RTOL,
                   DOUBLE_PRECISION& ATOL, INTEGER& ITASK,
                   INTEGER& IOPT, DOUBLE_PRECISION* RVODE, INTEGER& LRW,
                   INTEGER* IVODE, INTEGER& LIW, INTEGER& MF,
                   DOUBLE_PRECISION* RWORK, INTEGER* IWORK, INTEGER& LOUT,
                   INTEGER& ISTATE)
{   cpsolve (NKK, NNP, NNWT, NNH, NNWDOT, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL, ITASK,
              IOPT, RVODE, LRW, IVODE, LIW, MF,
              RWORK, IWORK, LOUT, ISTATE); }
#endif
#endif // !defined( F77_STUB_REQUIRED )
#endif // !defined( __CONP_H__ )

