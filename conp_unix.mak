#***********************************************************************BEGIN-HDR
#
#  Copyright (c) 1997-2006 Reaction Design.  All rights reserved.
#
# $Header: /chemkin/4.1/prebuild/samples/sample_apps_cpp/conp/conp_unix.mak 5     12/21/05 8:44p Frupley $
#
#  $NoKeywords: $
#***********************************************************************END-HDR
#
# This make file is used to build the Chemkin CONP sample.

ALL : conp.out

# relative path to the CHEMKIN install directory
CKROOT  = /afs/athena.mit.edu/software/chemkin_v4.1.1/distrib/@sys

# global CHEMKIN make directives
include $(CKROOT)/include/chemkin_make_unix.inc

# make directives for this problem
include conp_make.inc

# ----------------------------------------------------------------------------

