/*
 *
 */

#ifndef __IOINFO_H_
#define __IOINFO_H_

/*
 * $Id: IOinfo.h,v 1.1 2004/03/22 18:53:14 smurman Exp $
 */

#include "c3d_global.h"
#include "basicTypes.h"

/* -------------------------------------------------------------- */
typedef struct IOinfoStructure  tsIOinfo;
typedef tsIOinfo *p_tsIOinfo;

/*                                     ...04.03.04 new IOinfoStructure
                                          getting rid of discrepancy
                                          between SM and MA versions */
struct IOinfoStructure {
  char   InputCntlFile[FILENAME_LEN];
  char  OutputSolnFile[FILENAME_LEN];
  char     SurfTriFile[FILENAME_LEN];
  char    MeshInfoFile[FILENAME_LEN];
  char   InputMeshFile[FILENAME_LEN];
  char     PreSpecFile[FILENAME_LEN];
  char  OutputMeshFile[FILENAME_LEN];
  char     RestartFile[FILENAME_LEN];
  char   ConfigXMLFile[FILENAME_LEN];
  char ScenarioXMLFile[FILENAME_LEN];
  char     HistoryFile[FILENAME_LEN];  FILE *p_Historystrm;
  char       ForceFile[FILENAME_LEN];  FILE *p_Forcestrm; 
  char      MomentFile[FILENAME_LEN];  FILE *p_Momentstrm; 
  bool binaryIO;
};


#endif /* __IOINFO_H_ */
