#ifndef CAPSFUNIDS_H
#define CAPSFUNIDS_H
/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Function ID & Journal Header
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


/* base-level object functions */

#define CAPS_REVISION         1
#define CAPS_INFO             2
#define CAPS_SIZE             3
#define CAPS_CHILDBYINDEX     4
#define CAPS_CHILDBYNAME      5
#define CAPS_BODYBYINDEX      6
#define CAPS_GETHISTORY       7
#define CAPS_OWNERINFO        8
#define CAPS_DELETE           9
#define CAPS_ERRORINFO       10
#define CAPS_FREEERROR       11
#define CAPS_FREEVALUE       12


/* problem & I/O functions */

#define CAPS_OPEN             0
#define CAPS_PHASESTATE      20
#define CAPS_CLOSE           21
#define CAPS_OUTLEVEL        22
#define CAPS_WRITEPARAMETERS 23
#define CAPS_READPARAMETERS  24
#define CAPS_WRITEGEOMETRY   25
#define CAPS_GETROOTPATH     26
#define CAPS_INTENTPHRASE    27


/* attribute functions */

#define CAPS_ATTRBYNAME      30
#define CAPS_ATTRBYINDEX     31
#define CAPS_SETATTR         32
#define CAPS_DELETEATTR      33


/* analysis functions */

#define CAPS_QUERYANALYSIS   40
#define CAPS_GETBODIES       41
#define CAPS_GETINPUT        42
#define CAPS_GETOUTPUT       43
#define CAPS_MAKEANALYSIS    44
#define CAPS_DUPANALYSIS     45
#define CAPS_DIRTYANALYSIS   46
#define CAPS_ANALYSISINFO    47
#define CAPS_PREANALYSIS     48
#define CAPS_RUNANALYSIS     49
#define CAPS_CHECKANALYSIS   50
#define CAPS_POSTANALYSIS    51
#define CAPS_RESETANALYSIS   52
#define CAPS_GETTESSELS      53
#define CAPS_SYSTEM          54
#define CAPS_AIMBACKDOOR     59
  
  
/* bound, vertexset and dataset functions */

#define CAPS_MAKEBOUND       60
#define CAPS_BOUNDINFO       61
#define CAPS_CLOSEBOUND      62
#define CAPS_MAKEVERTEXSET   63
#define CAPS_VERTEXSETINFO   64
#define CAPS_OUTPUTVERTEXSET 65
#define CAPS_FILLUNVERTEXSET 66
#define CAPS_MAKEDATASET     67
#define CAPS_LINKDATASET     68
#define CAPS_INITDATASET     69
#define CAPS_SETDATA         70
#define CAPS_GETDATA         71
#define CAPS_GETDATASETS     72
#define CAPS_GETTRIANGLES    73
#define CAPS_DATASETINFO     74

  
/* value functions */

#define CAPS_GETVALUE        80
#define CAPS_MAKEVALUE       81
#define CAPS_SETVALUE        82
#define CAPS_GETLIMITS       83
#define CAPS_SETLIMITS       84
#define CAPS_GETVALUEPROPS   85
#define CAPS_SETVALUEPROPS   86
#define CAPS_CONVERT         87
#define CAPS_CONVERTVALUE    88
#define CAPS_TRANSFERVALUES  89
#define CAPS_LINKVALUE       90
#define CAPS_HASDERIV        91
#define CAPS_GETDERIV        92

/* string conversion for the funIDs */

#define CAPS_NFUNID          93

extern const char* caps_funID[CAPS_NFUNID];

#endif
