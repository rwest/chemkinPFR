/***********************************************************************BEGIN-HDR
*
*  Copyright (c) 1997-2006 Reaction Design.  All rights reserved.
*
*  $Header: /chemkin/4.1/prebuild/samples/sample_apps_cpp/conp/conp.cpp 9     12/21/05 8:44p Frupley $
*
*  $NoKeywords: $
************************************************************************END-HDR
*/
#include "conp.hpp"

// function prototypes

int main()
{   char * sInputfileName  = (char *)"conp.inp";   // input file for the problem
    char * sOutputfileName = (char *)"conp.out";   // output file for the problem
    char * sGasLinkfileName= (char *)"chem.asc";   // gas-phase mechanism linkfile
    int iOutputfileUnit    = 9;            // Fortran output file unit number
    int iGasLinkfileUnit   = 25;           // Fortran linkfile unit number

    int iFlag=0;
    if (!strstr(sOutputfileName,"stdout") && iOutputfileUnit!=6) {
       // Open diagnostics file for Fortran WRITE (requires integer Fortran unit)
       CCOPEN (sOutputfileName,
               (char *)"FORMATTED",
               (char *)"UNKNOWN", iOutputfileUnit, iFlag);
       if (iFlag) return iFlag;
    }

    // print program identification
    int nlines=1;
    BANNER (iOutputfileUnit,
            (char *)"CONP",
            (char *)"Gas Phase Kinetics Sample Program", nlines);

    // Open gas-phase linkfile for Fortran READ (requires integer Fortran unit)
    CCOPEN (sGasLinkfileName,
            (char *)"FORMATTED",
            (char *)"OLD", iGasLinkfileUnit, iFlag);
    if (iFlag) {
        CFMESS(sOutputfileName, (char *)"Error opening linkfile...");
        return iFlag;
    }

    // get CHEMKIN memory requirements, species count, string length
    int iCKsizeI=0, iCKsizeD=0, iCKsizeS=0, iSpeciesCount=0,iStringLength=0;
    if ( iFlag = cpsize (sOutputfileName, iOutputfileUnit,
                         sGasLinkfileName, iGasLinkfileUnit,
                         iCKsizeI, iCKsizeD, iCKsizeS,
                         iSpeciesCount, iStringLength) )
    {  CCCLOS (sOutputfileName);  return iFlag; }

    // extend real workspace size for function processes
    int iPPtr        = iCKsizeI + 1;
    int iWtPtr       = iPPtr    + 1;
    int iHPtr        = iWtPtr   + iSpeciesCount;
    int iWdotPtr     = iHPtr    + iSpeciesCount;
    int iCKworkTotal = iWdotPtr - 1;

    // allocate workspaces
    int    *iCKwork = new int[ iCKsizeI ];
    double *dCKwork = new double[ iCKworkTotal ];
    int nchar = iCKsizeS*iStringLength + 1;
    char *sCKwork = new char [ nchar ];
    // initialize CHEMKIN work arrays
    CKINIT (iCKsizeI, iCKsizeD, iCKsizeS,
            iGasLinkfileUnit, iOutputfileUnit,
            iCKwork, dCKwork, sCKwork, iFlag);
    CCCLOS (sGasLinkfileName);   // done with linkfile
    if (iFlag) {
       CFMESS(sOutputfileName,(char *)"Error initializing linkfile data...");
       delete [] iCKwork; delete [] dCKwork;  delete [] sCKwork;
       return iFlag;
    }

    // initialize species symbols
    nchar = iSpeciesCount*iStringLength + 1;
    char *sSpeciesNames = new char [ nchar ];
    for (int ichar=0; ichar<nchar-1; ichar++) sSpeciesNames[ichar] = ' ';
    sSpeciesNames[nchar-1] = '\0';
    CKSYMS (sCKwork, iOutputfileUnit, sSpeciesNames, iFlag);
    if (iFlag) {
       CFMESS(sOutputfileName,(char *)"Error initializing species symbols...");
       delete [] iCKwork; delete [] dCKwork; delete [] sCKwork;
       delete [] sSpeciesNames;
       return iFlag;
    }

    // get molecular weights and store into extension of work array
    // for availability to Fortran function
    CKWT (iCKwork, dCKwork, &dCKwork[iWtPtr-1]);

    // problem parameters
    double dPres=0.0;  double dTemp=0.0;
    double *dMoleFractions = new double[ iSpeciesCount ];
    for (int i=0; i<iSpeciesCount; i++) dMoleFractions[i] = 0.0;
    double dTend=0.0; double dTdelta=0.0;
    // initialization from input file
    if ( iFlag = cpinp (sOutputfileName, iOutputfileUnit,
                        sInputfileName, iCKwork, dCKwork,
                        iSpeciesCount, iStringLength,
                        sSpeciesNames, dMoleFractions,
                        dPres, dTemp, dTend, dTdelta) )
    {  CCCLOS (sOutputfileName);
       delete [] iCKwork; delete [] dCKwork;  delete [] sCKwork;
       delete [] sSpeciesNames; delete [] dMoleFractions;
       return iFlag;
    }

    // initial conditions
    double dTstart = 0.0; double dTlast = 0.0;             // times
    int     iEquationCount = iSpeciesCount + 1;            // equations to solve
    double *dSolution      = new double[ iEquationCount ]; // solution array
    dSolution[0] = dTemp;                                  // temperature
    dCKwork[iPPtr-1] = dPres;                              // pressure
    CKXTY(dMoleFractions, iCKwork, dCKwork, dSolution+1);  // mass fractions
    int bVprint=true;                                      // for VODE output
    LUNSAV(iOutputfileUnit, bVprint);

    // ODE solver work arrays
    int     iODEsizeI = iEquationCount + 30;
    int     iODEsizeD = 22 + 9*iEquationCount +
                        (2*iEquationCount)*(2*iEquationCount);
    int    *iODEwork  = new int[ iODEsizeI ];
    double *dODEwork  = new double[ iODEsizeD ];

    // Integration control parameters for VODE
    int iMethodFlag=22, iState=1, iItol=1, iOpt=0, iTask=1;
    double dAtol = 1.0e-15; double dRtol = 1.0e-6;

    // print initial conditions to the Fortran output file
    iFlag = cpout (sOutputfileName, iOutputfileUnit, iCKwork, dCKwork,
                   dSolution, iSpeciesCount, sSpeciesNames, dMoleFractions,
                   iStringLength, dTemp, dTlast);
    while ( iState>-2 && dTlast<dTend && iFlag==0) {
       // Call the differential equation solver
       dTlast =  min(dTlast + dTdelta, dTend);
       CPSOLVE(iSpeciesCount, iPPtr, iWtPtr, iHPtr, iWdotPtr,
               iEquationCount, dSolution, dTstart, dTlast, iItol, dRtol,
               dAtol, iTask, iOpt, dODEwork, iODEsizeD,
               iODEwork, iODEsizeI, iMethodFlag, dCKwork, iCKwork,
               iOutputfileUnit, iState);
       // append the current solution to the Fortran output file
       iFlag = cpout (sOutputfileName, iOutputfileUnit, iCKwork, dCKwork,
                      dSolution, iSpeciesCount, sSpeciesNames, dMoleFractions,
                      iStringLength, dTemp, dTlast);
	   // reset the amount of N2 in the dMoleFractions array.
	   iFlag = resetN2 (sOutputfileName, iOutputfileUnit, iCKwork, dCKwork, dSolution,
						iSpeciesCount, iStringLength, sSpeciesNames, dMoleFractions);
	   //CKXTY(dMoleFractions, iCKwork, dCKwork, dSolution+1);
    }

    CFMESS (sOutputfileName,(char *)"END OF INTEGRATION...");
    CCCLOS (sOutputfileName);
    CKDONE (iCKwork, dCKwork);
    //delete [] iCKwork; delete [] dCKwork; delete [] sCKwork;
    //delete [] sSpeciesNames; delete [] dMoleFractions; delete [] dSolution;
    //delete [] iODEwork; delete [] dODEwork;
    if (iState<=-2) return iState;
    else return iFlag;
}

int cpsize (char *sOutputfileName, int iOutputfileUnit,
            char *sGasLinkfileName, int iGasLinkfileUnit,
            int &iCKsizeI, int &iCKsizeD, int &iCKsizeS,
            int &iSpeciesCount, int &iStringLength)
{  // Get information about the mechanism, and array requirements
   int iFileNameLength=0, iReactionStringLength=0;
   CKCLEN (iStringLength, iFileNameLength, iReactionStringLength );
   int iElementCount=0, iReactionCount=0, iMaxReactionSpecies,
       iMaxSpeciesFitTemperatures=0, iMaxReactionThirdBodies=0,
       iMaxReactionOrderChanges=0, iIonSpeciesCount=0, iFlag=0;
   CKLEN2 (iGasLinkfileUnit, iOutputfileUnit,
           iCKsizeI, iCKsizeD, iCKsizeS,
           iElementCount, iSpeciesCount, iReactionCount,
           iMaxReactionSpecies, iMaxSpeciesFitTemperatures,
           iMaxReactionThirdBodies, iMaxReactionOrderChanges,
           iIonSpeciesCount, iFlag);
   if (iFlag) {
      CFMESS(sOutputfileName, (char *)"Error reading linkfile...");
      CCCLOS (sGasLinkfileName);
   }
   return iFlag;
}

int resetN2 (char *sOutputfileName, int iOutputfileUnit,
			 int *iCKwork, double *dCKwork, double *dSolution,
			 int iSpeciesCount, int iStringLength,
			 char *sSpeciesNames, double *dMoleFractions)
{
	/*  
	 This function takes the vector dMoleFractions (which may or may not sum to 1.0)
	 and replaces the N2 entry with (1-sum(others)). It does not update the chemkin
	 simulation work arrays - you must do this yourself. e.g. call
	   CKYTX(dSolution+1, iCKwork, dCKwork, dMoleFractions);
	 before calling this function and 
	   CKXTY(dMoleFractions, iCKwork, dCKwork, dSolution+1);
	 after calling this function.
	 
	 */
	int iFlag=0;
	FILE *fpOutfile = stdout;
	if (!strstr(sOutputfileName,"stdout") && iOutputfileUnit != 6) {
		// switch from Fortran to C++ formatted output file
		CCCLOS(sOutputfileName);
		fpOutfile = fopen(sOutputfileName, "a");
	}
	
	// convert current mass fractions to mole fractions
	//CKYTX(dSolution+1, iCKwork, dCKwork, dMoleFractions);  // get mole fractions
	
	char *sReactant= new char [ iStringLength + 1 ];
	char *sName    = new char [ iStringLength + 1 ];
	int i,iFound; double dTotal;
	
	dTotal = 0;
	sReactant = "N2";
	
	iFound=-1;
	for (i=0; i < iSpeciesCount; i++) {
		sscanf(sSpeciesNames+(i*iStringLength), "%s", sName);
		if (strcmp(sReactant, sName)==0) {
			iFound = i;
		}
		else {
			dTotal = dTotal + dMoleFractions[i];
		}
		i++;
	}
	if (iFound < 0) {
		iFlag = 1; fprintf(fpOutfile, "%s\n", "Couldn't find N2 in species list.");
	}
	else {
		dMoleFractions[iFound] = 1.0 - dTotal;
	}
	fprintf(fpOutfile, "%s %10.3e\n", "N2 mole fraction should equal ",(1.0 - dTotal));
	
	// CKXTY(dMoleFractions, iCKwork, dCKwork, dSolution+1);  // get mass fractions, i.e set mole fractions
	
	delete [] sName;
	
	if (fpOutfile != stdout) {                   // switch to Fortran output
		fflush(fpOutfile); fclose(fpOutfile);
		int ifFlag=0;
		CCOPEN (sOutputfileName,
				(char *)"FORMATTED",
				(char *)"UNKNOWN", iOutputfileUnit, ifFlag);
		if (ifFlag==0) ifFlag = CFEND (sOutputfileName); // go to end of output file
		return max(iFlag,ifFlag);
	}
	return iFlag;
}

int cpinp (char *sOutputfileName, int iOutputfileUnit,
           char *sInputfileName,
           int *iCKwork, double *dCKwork,
           int iSpeciesCount, int iStringLength,
           char *sSpeciesNames, double *dMoleFractions,
           double &dPres, double &dTemp, double &dTend, double &dTdelta)
{  int iFlag=0;
   FILE *fpInfile  = stdin;
   FILE *fpOutfile = stdout;

   if (!strstr(sOutputfileName,"stdin")) {
       // open an input file
       fpInfile = fopen(sInputfileName, "r");
       if (!fpInfile) {
          CFMESS(sOutputfileName, (char *)"Error opening input file...");
          return 1;
      }
   }
   if (!strstr(sOutputfileName,"stdout") && iOutputfileUnit != 6) {
      // switch from Fortran to C++ formatted output file
      CCCLOS(sOutputfileName);
      fpOutfile = fopen(sOutputfileName, "a");
   }

   // pressure (atm), temperature (K)
   fprintf(fpOutfile, "%s\n\n%s\n",
                      "ADIABATIC FIXED PRESSURE PROBLEM",
                      "INPUT PRESSURE(ATM) AND TEMPERATURE(K)");
   fscanf(fpInfile, "%lf%lf", &dPres, &dTemp);
   fprintf(fpOutfile, "%10.3e %10.3e\n", dPres, dTemp);

   // CHEMKIN gas and pressure constants
   double dRu=0.0, dRuc=0.0, dPatm=0.0;
   CKRP (iCKwork, dCKwork, dRu, dRuc, dPatm);
   dPres = dPres * dPatm;

   // Initial non-zero moles
   char *sReactant= new char [ iStringLength + 1 ];
   char *sName    = new char [ iStringLength + 1 ];
   int i,iFound; double dValue;

   while (strncmp(sReactant,"END",3) !=0 ) {
      fprintf(fpOutfile, "\n%s\n", "INPUT MOLES OF NEXT SPECIES");
      fscanf(fpInfile, "%s", sReactant);
      fprintf(fpOutfile, "%s", sReactant);

      if (strncmp(sReactant,"END",3)!=0) {
         fscanf(fpInfile, "%lf", &dValue);
         fprintf(fpOutfile, " %10.3e\n", dValue);
         // CKCOMP(sReactant, sSpeciesNames, iSpeciesCount, iSpeciesIndex);
         i=0; iFound=-1;
         while (iFound<0 && i<iSpeciesCount) {
             sscanf(sSpeciesNames+(i*iStringLength), "%s", sName);
             if (strcmp(sReactant, sName)==0) {
                iFound = i;
                dMoleFractions[iFound] = dValue;
             }
             i++;
         }
         if (iFound < 0) {
            iFlag = 1; fprintf(fpOutfile, "%s\n", "Error reading moles...");
         }
      }
   }
   fprintf(fpOutfile, "\n");
   delete [] sName;

   // Normalize the mole fractions
   CKNORM(dMoleFractions, iSpeciesCount);

   // Final time and print interval
   fprintf(fpOutfile, "\n%s\n", "INPUT FINAL TIME AND DT");
   fscanf(fpInfile,  "%lf%lf", &dTend, &dTdelta);
   fprintf(fpOutfile, "%10.3e %10.3e\n\n", dTend, dTdelta);
   if (fpInfile != stdin) fclose(fpInfile);     // done with input file
   if (fpOutfile != stdout) {                   // switch to Fortran output
      fflush(fpOutfile); fclose(fpOutfile);
      int ifFlag=0;
      CCOPEN (sOutputfileName,
              (char *)"FORMATTED",
              (char *)"UNKNOWN", iOutputfileUnit, ifFlag);
      if (ifFlag==0) ifFlag = CFEND (sOutputfileName); // go to end of output file
      return max(iFlag,ifFlag);
   }
   return iFlag;
}

int cpout (char *sOutputfileName, int iOutputfileUnit,
           int *iCKwork, double *dCKwork, double *dSolution,
           int iSpeciesCount, char *sSpeciesNames, double *dMoleFractions,
           int iStringLength, double dTemp, double dTime)
{
   FILE *fpOutfile = stdout;
   if (iOutputfileUnit!=6 ) {  // switch from Fortran to C++ file
      CCCLOS(sOutputfileName);
      fpOutfile = fopen(sOutputfileName, "a");
   }

   int i, iNcol=7, iCol = 3, iLine=1;
   if (dTime == 0.0) {           // print column headings
      fprintf(fpOutfile, "\n%8s  %11s   ", "T(SEC)", "TMP(K)");
      char *sName    = new char [ iStringLength + 1 ];
      for (i=0; i < iSpeciesCount; i++) {
         sscanf(sSpeciesNames+(i*iStringLength), "%10s  ", sName);
         if (iCol==3 && iLine>1) fprintf(fpOutfile, "%24s", " ");
         fprintf(fpOutfile, "%10s  ", sName);
         if (++iCol > iNcol || i+1 == iSpeciesCount)
             { iCol=3; iLine++; fprintf(fpOutfile,"\n"); }
      }
      delete [] sName;
   }

   // convert current mass fractions to mole fractions
   CKYTX(dSolution+1, iCKwork, dCKwork, dMoleFractions);
   // print time, temperature, mole fractions
   dTemp = dSolution[0];
   fprintf(fpOutfile, "%10.3e  %10.3e  ", dTime, dTemp);
   iCol = 3; iLine = 1;
   for (i=0; i < iSpeciesCount; i++) {
       if (iCol==3 && iLine>1) fprintf(fpOutfile, "%24s", " ");
       fprintf(fpOutfile, "%10.3e  ", dMoleFractions[i]);
       if (++iCol > iNcol || i+1 == iSpeciesCount)
          { iCol=3; iLine++; fprintf(fpOutfile,"\n"); }
   }
   fflush(fpOutfile);

   int iFlag = 0;
   if (fpOutfile != stdout) {    // switch back to Fortran output
      CCOPEN (sOutputfileName,
              (char *)"FORMATTED",
              (char *)"UNKNOWN", iOutputfileUnit, iFlag);
      if (iFlag==0) {
         iFlag = CFEND (sOutputfileName);   // go to end of file
         if (iFlag != 0) fprintf(fpOutfile,"\n%s\n","Error re-reading output file"); }
      else
         fprintf(fpOutfile,"\n%s\n","Error re-opening output file");
      fclose (fpOutfile);
   }
   return iFlag;
}


