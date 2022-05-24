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

#include "capsFunIDs.h"

const char *caps_funID[CAPS_NFUNID];

void
caps_initFunIDs()
{

  /* base-level object functions */
/*@-observertrans@*/
  caps_funID[CAPS_REVISION    ] = "revision";
  caps_funID[CAPS_INFO        ] = "info";
  caps_funID[CAPS_SIZE        ] = "size";
  caps_funID[CAPS_CHILDBYINDEX] = "childByIndex";
  caps_funID[CAPS_CHILDBYNAME ] = "childByName";
  caps_funID[CAPS_BODYBYINDEX ] = "bodyByIndex";
  caps_funID[CAPS_GETHISTORY  ] = "getHistory";
  caps_funID[CAPS_OWNERINFO   ] = "ownerInfo";
  caps_funID[CAPS_DELETE      ] = "delete";
  caps_funID[CAPS_ERRORINFO   ] = "errorInfo";
  caps_funID[CAPS_FREEERROR   ] = "freeError";
  caps_funID[CAPS_FREEVALUE   ] = "freeValue";

  /* problem & I/O functions */

  caps_funID[CAPS_OPEN           ] = "open";
  caps_funID[CAPS_PHASESTATE     ] = "phaseState";
  caps_funID[CAPS_CLOSE          ] = "close";
  caps_funID[CAPS_OUTLEVEL       ] = "outLevel";
  caps_funID[CAPS_WRITEPARAMETERS] = "writeParameters";
  caps_funID[CAPS_READPARAMETERS ] = "readParameters";
  caps_funID[CAPS_WRITEGEOMETRY  ] = "writeGeometry";
  caps_funID[CAPS_GETROOTPATH    ] = "getRootPath";
  caps_funID[CAPS_INTENTPHRASE   ] = "intentPhrase";


  /* attribute functions */

  caps_funID[CAPS_ATTRBYNAME ] = "attrByName";
  caps_funID[CAPS_ATTRBYINDEX] = "attrByIndex";
  caps_funID[CAPS_SETATTR    ] = "setAttr";
  caps_funID[CAPS_DELETEATTR ] = "deleteAttr";


  /* analysis functions */

  caps_funID[CAPS_QUERYANALYSIS] = "queryAnalysis";
  caps_funID[CAPS_GETBODIES    ] = "getBodies";
  caps_funID[CAPS_GETINPUT     ] = "getInput";
  caps_funID[CAPS_GETOUTPUT    ] = "getOutput";
  caps_funID[CAPS_MAKEANALYSIS ] = "makeAnalysis";
  caps_funID[CAPS_DUPANALYSIS  ] = "dupAnalysis";
  caps_funID[CAPS_DIRTYANALYSIS] = "dirtyAnalysis";
  caps_funID[CAPS_ANALYSISINFO ] = "analysisInfo";
  caps_funID[CAPS_PREANALYSIS  ] = "preAnalysis";
  caps_funID[CAPS_RUNANALYSIS  ] = "runAnalysis";
  caps_funID[CAPS_CHECKANALYSIS] = "checkAnalysis";
  caps_funID[CAPS_POSTANALYSIS ] = "postAnalysis";
  caps_funID[CAPS_RESETANALYSIS] = "resetAnalysis";
  caps_funID[CAPS_GETTESSELS   ] = "getTessels";
  caps_funID[CAPS_SYSTEM       ] = "system";
  caps_funID[CAPS_AIMBACKDOOR  ] = "aimBackdoor";


  /* bound, vertexset and dataset functions */

  caps_funID[CAPS_MAKEBOUND      ] = "makeBound";
  caps_funID[CAPS_BOUNDINFO      ] = "boundInfo";
  caps_funID[CAPS_CLOSEBOUND     ] = "closeCound";
  caps_funID[CAPS_MAKEVERTEXSET  ] = "makeVertexset";
  caps_funID[CAPS_VERTEXSETINFO  ] = "vertexSetInfo";
  caps_funID[CAPS_OUTPUTVERTEXSET] = "outputVertexSet";
  caps_funID[CAPS_FILLUNVERTEXSET] = "fillUnVertexSet";
  caps_funID[CAPS_MAKEDATASET    ] = "makeDataSet";
  caps_funID[CAPS_LINKDATASET    ] = "linkDataSet";
  caps_funID[CAPS_INITDATASET    ] = "initDataSet";
  caps_funID[CAPS_SETDATA        ] = "setData";
  caps_funID[CAPS_GETDATA        ] = "getData";
  caps_funID[CAPS_GETDATASETS    ] = "getDataSets";
  caps_funID[CAPS_GETTRIANGLES   ] = "getTriangles";
  caps_funID[CAPS_DATASETINFO    ] = "dataSetInfo";


  /* value functions */

  caps_funID[CAPS_GETVALUE      ] = "getValue";
  caps_funID[CAPS_MAKEVALUE     ] = "makeValue";
  caps_funID[CAPS_SETVALUE      ] = "setValue";
  caps_funID[CAPS_GETLIMITS     ] = "getLimits";
  caps_funID[CAPS_SETLIMITS     ] = "setLimits";
  caps_funID[CAPS_GETVALUEPROPS ] = "getValueProps";
  caps_funID[CAPS_SETVALUEPROPS ] = "setValueProps";
  caps_funID[CAPS_CONVERT       ] = "convert";
  caps_funID[CAPS_CONVERTVALUE  ] = "convertValue";
  caps_funID[CAPS_TRANSFERVALUES] = "transferValues";
  caps_funID[CAPS_LINKVALUE     ] = "linkValue";
  caps_funID[CAPS_HASDERIV      ] = "hasDeriv";
  caps_funID[CAPS_GETDERIV      ] = "getDeriv";
/*@+observertrans@*/

}
