#***********************************************************************BEGIN-HDR
#
#  Copyright (c) 1997-2006 Reaction Design.  All rights reserved.
#
#  $Header: /chemkin/4.1/prebuild/samples/sample_apps_cpp/conp/conp_make.inc 9     12/21/05 8:45p Frupley $
#
#  $NoKeywords: $
#***********************************************************************END-HDR
#
# This make file is used to build the Chemkin CONP sample.

#
# dependencies
#

SLATEC_OBJS = \
 slatec/src/dqform$(OBJECT)     \
 slatec/src/dqrfac$(OBJECT)    \
 slatec/src/i1mach$(OBJECT)    \
 slatec/src/d1mach$(OBJECT)    \
 slatec/src/d1mpyq$(OBJECT)    \
 slatec/src/d1updt$(OBJECT)    \
 slatec/src/ddoglg$(OBJECT)    \
 slatec/src/denorm$(OBJECT)    \
 slatec/src/dfdjc1$(OBJECT)    \
 slatec/src/dnsq$(OBJECT)      \
 slatec/src/dnsqe$(OBJECT)     
# already defined in chemkin_v4.1.1/distrib/@sys/lib/chemkin_public.a :
# slatec/src/j4save$(OBJECT)    \
# slatec/src/fdump$(OBJECT)     \
# slatec/src/xercnt$(OBJECT)    \
# slatec/src/xerhlt$(OBJECT)    \
# slatec/src/xermsg$(OBJECT)    \
# slatec/src/xerprn$(OBJECT)    \
# slatec/src/xersve$(OBJECT)    \
# slatec/src/xgetua$(OBJECT)    \

CONP_OBJS = cpsolve$(OBJECT) $(SLATEC_OBJS) dnsqsolve$(OBJECT)

conp$(OBJECT) : conp.cpp
	$(CXX) $(CPPFLAGS) -c conp.cpp $(CPPF77FLAGS)

CHEMKIN_PKG = $(CHEMKIN_LIBRARY) $(CHEMKIN_PUB_LIBRARY) \
		$(CHEMKIN_USER_LIBRARY) $(CPPF77_SYSTEM_LIBS)

conp$(EXE)     :  conp$(OBJECT) $(CONP_OBJS)
	$(CK_DISPLAY) Creating $@
	$(CPPF77LINK) conp$(OBJECT) $(CPP2F77_LIBRARY) $(CPPF77_LIBFLAGS) $(CONP_OBJS) $(CHEMKIN_PKG)

#
# run the programs
#

#*************************************
#   Run gas-phase Chemkin interpreter
#     the input file is 'chem.inp'
#     the output file is 'chem.out'
#*************************************
chem.asc     : $(CHEMKINBIN)chem$(EXE) chem.inp
	$(CK_DISPLAY) Running chem
	$(REMOVE) $(CHEM_OUTPUTS)
	$(CHEMKINBIN)chem -i chem.inp -o chem.out

#*************************************
#   Run CONP
#     the input file is 'conp.inp'
#     the output file is 'conp.out'
#*************************************
CONP_OUTPUTS = conp.out
conp.out     : conp$(EXE) chem.asc conp.inp
	$(CK_DISPLAY) Running conp
	.$(SLASH)conp$(EXE)

CLEAN       :
	$(REMOVE) $(EXTRA_CLEAN) $(CONP_OUTPUTS) $(ALL_EXES) $(ALL_LIBS)

# ----------------------------------------------------------------------------

