
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: c3dio_lib.h,v 1.18 2022/11/07 18:45:42 mnemec Exp $
 */

/* open source */

#ifndef __C3DIO_LIB_H_
#define __C3DIO_LIB_H_

#include <stdio.h>

#include "basicTypes.h"
#include "stateVector.h"
#include "geomStructures.h"
#include "cartCells.h"
#include "gridStructures.h"
#include "memory_util.h"

#include "SolverInfo.h"
#include "GridInfo.h"
#include "CaseInfo.h"
#include "MGinfo.h"
#include "PostProcInfo.h"
#include "IOinfo.h"
#include "ConvInfo.h"

#include "dualtime.h"          /* get the time-dependent structure tsTDinfo */

#define IO_VERBOSE 1

/* -------------- datatypes --------------------- */

typedef struct controlFileStructure tsControl;
typedef tsControl *p_tsControl;

struct controlFileStructure {
  tsSinfo            SolverInfo;   /* ...local copies of InfoStructures   */
  tsGinfo              GridInfo;
  tsCinfo              CaseInfo;
  p_tsHinfo      p_convergeInfo;   /* ...Point to main Converg  history   */
  p_tsIOinfo         p_fileInfo;   /* ...Point to file name information   */
  p_tsMGinfo           p_mgInfo;   /* ...Point to multigrid scheme info   */
  p_tsPPinfo           p_ppInfo;   /* ...Point to postProcessing info     */
  int               nSubDomains;   /*    (Highest level grid for the M-G  */
  int                  PartType;   /* ...type of dom-decomp being used    */
};

/* -------------- define macros --------------------- */

/*          ...integer bitmask used for options for TRIX reader and writer */
#define TRIX_VERBOSE    1      /* bit 1: verbose option for TRIX reader    */
#define IS_TRIX         2      /* bit 2: set to TRUE if valid TRIX file    */
#define TRIX_DP_VERTS   4      /* bit 3: read/write double precision verts */

#define VALID_TRIX_FILE 1

/*                            ...define file extensions for tecplot files */
#define TECPLOT_ASCII_FILE_EXTENSION  "dat"
#define TECPLOT_BINARY_FILE_EXTENSION "plt"

/**
 *  ==================================================================
 *
 *  1. Top-Level application interface io_()
 *  2. Mid-Level application interface trix_()
 *  3. Low-Level application interface c3d_()
 *
 *  ==================================================================
 */

/* ------- 2.  Top-Level Application Interface --------------------- */

#ifdef __cplusplus
extern "C" {
#endif
  void          io_setVerboseIO(bool);
  int           io_readCNTLfile(p_tsControl);
  void          io_readGinfo(p_tsControl);
  void          io_newC3Dcntl(p_tsControl);
  void          io_freeC3Dcntl(p_tsControl p_c3dInfo);
  void          io_initializeInfoStructs(p_tsControl);
  /* negative value on return represents error condition */
  int           io_readSurfTri(const char*, p_tsTriangulation);
  int           io_writeSurfTri(char*, const p_tsTriangulation, bool is_ascii);
  int           io_readFineMesh(char*, const p_tsBasicGrid p_g);
  void          io_freeBasicGrid(p_tsBasicGrid p_g);
  FILE*         io_strmMetaInfo(p_tsControl, FILE *);
  int           io_writeMesh(char*, const p_tsBasicGrid);

  /*  NOTE: io_read/writeCkpt() ...if p_tsTDinfo is NULL then we read/write a
	                                 steady-state checkpoint file, else we
	                                 read/write the time-dependent  */
  void          io_readCkPt(char *, p_tsBasicGrid , p_tsMGinfo , p_tsHinfo,
                            const int, p_tsTDinfo);
  int           io_writeCkPt(char *, p_tsBasicGrid, p_tsMGinfo,  p_tsHinfo,
                             p_tsTDinfo);

  /**
   * io_readSurfTrix(): read VTK unstructured grid ASCII file.
   *                    Function returns the number of components,
   *                    nComps, and a triangulation for each
   *                    component. nComps must be initialized to zero
   *                    by the caller to read in the first component
   *                    set. If p_compName is specified then only that
   *                    component is retrieved. For each component,
   *                    you may also request specific vert and tri
   *                    data. Initializing p_compName,
   *                    p_vertDataNames, p_triDataNames to ALL
   *                    retrieves everything from the file. see
   *                    documentation of bitwise trix 'options'
   *                    variable above
   */
  int io_readSurfTrix(const char *const p_fileName, p_tsTriangulation *pp_config,
                      int *const p_nComps, const char *const p_compName,
                      const char *const p_vertDataNames,
                      const char *const p_triDataNames,
                      const int options);

  /**
   *  io_readTrixIntersectDims(): basic function that returns
   *                              dimensions of triangulation in
   *                              p_fileName. Passing pointers and
   *                              returning void so it is xcallable
   *                              from fortran, i.e. Intersect. If any
   *                              returned value is -1, this indicates
   *                              an error condition.
   */
  void io_readTrixIntersectDims(const char *const p_fileName,
                                int *p_nVerts, int *p_nTris,
                                int *p_nVertScalars, int *p_nTriScalars,
                                int *p_nVertData, int *p_nTriData,
                                const int options);
  /**
   * io_writeSurfTrix(): write configuration in VTK XML
   *                     UnstructuredGrid ASCII format.
   */
  int io_writeSurfTrix(const p_tsTriangulation p_config, const int nComps,
                       const char *const p_fileName, const int options);

  /**
   *  io_writeTecplot(): writes a tecplot file of p_surf with triData
   *                     written as cell-centered data. Data in
   *                     a_Tris[].Comp
   */
  int io_writeTecplot(FILE *p_strm, const p_tsTriangulation p_surf);


  /* -------  Mid Level Application Interface ------------------ */

  /**
   *  trix_initLibXML2(): mid-level-wrapper for LibXML2 parser
   *                      initialization only needed by apps that
   *                      handle parser init and shutdown externally
   */
  void trix_initLibXML2(void);

  /**
   *  trix_stopLibXML2(): mid-level-wrapper forfor xmlCleanupParser()
   *                      needed by apps directly calling
   *                      trix_initLibXML2()
   */
  void trix_stopLibXML2(void);

  /**
   *  trix_readSurf(): mid-level trix surface reader. Application Must
   *                   call trix_initLibXML2() before calling this and
   *                   clean up using call to trix_stopLibXML2() after
   *                   all reading is completed.  see documentation of
   *                   trix 'options' above
   */
  int trix_readSurf(const char *const p_fileName, p_tsTriangulation *pp_config,
                    int *const p_nComps, const char *const p_compName,
                    const char *const p_vertDataNames,
                    const char *const p_triDataNames, const int options);

  /* -------  Low Level Application Interface ------------------ */

  /**
   * c3d_isASCII()  Return 1 if file is ASCII, using unix system calls
   */
  int c3d_isASCII(const char *const p_name);

#ifdef __cplusplus
}
#endif

#endif /* __C3DIO_LIB_H_ */
