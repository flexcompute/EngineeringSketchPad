
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
 * $Id: test_c3dio.c,v 1.1 2022/11/07 18:47:10 mnemec Exp $
 */

/* open source */

/**
 * Test reading and writing tri files. Usage:
 * test_c3dio <input_tri_filename> <out_tri_filename>
 * If input is a VTK file, output will be VTK.
 */

#include <string.h>
#include "c3dio_lib.h"

int main(const int argc, char* p_argv[]){

  p_tsTriangulation p_surf = NULL;
  int nComps = 0;
  int  rc, trixOpts=0, i, j;
  char in_fname[FILENAME_LEN], out_fname[FILENAME_LEN];
  bool is_ascii = TRUE;
  bool is_trix = TRUE;

  if( argc != 3 ) {
    ERR("Needs 2 arguments: test_c3dio <input_tri_filename> <out_tri_filename>\n");
    exit(1);
  }

  strncpy(in_fname, p_argv[1], FILENAME_LEN);
  strncpy(out_fname,p_argv[2], FILENAME_LEN);

  if (!strlen(in_fname)) {
    ERR("Check input triangulation name\n");
    exit(1);
  }

  io_setVerboseIO(TRUE);
  trixOpts |= TRIX_VERBOSE;

  rc = io_readSurfTrix(in_fname, &p_surf, &nComps, "ALL", "ALL", "ALL", trixOpts);
  if (rc != 0) {
    is_trix = FALSE;
    nComps = 1;
    c3d_newTriangulation(&p_surf, 0, nComps);
    rc = io_readSurfTri(in_fname, p_surf);
    if (rc != 0) {
      ERR("Could not read file \"%s\", check name and format\n",in_fname);
      exit(FILE_ERROR);
    }
    strncpy(p_surf->geomName, in_fname, FILENAME_LEN);
  }

  NOTE("Parsed %d components\n", nComps);

  for (i=0; i<nComps; ++i) {
    printf("       Comp %d Name: \"%s\" nVerts %d nTris %d nVertData %d nTriData %d\n",
           i, p_surf[i].geomName, p_surf[i].nVerts, p_surf[i].nTris,
           p_surf[i].nVertData, p_surf[i].nTriData);

    for (j=0; j<p_surf[i].nTriData; ++j) {
      printf("         triData: %s\n",p_surf[i].p_triData[j].name);
    }

    for (j=0; j<p_surf[i].nVertData; ++j) {
      printf("         vertData: %s\n",p_surf[i].p_vertData[j].name);
    }
  }

  if (strlen(out_fname)) {
    if (is_trix) {
      rc = io_writeSurfTrix(p_surf, nComps, out_fname, trixOpts);
      if (rc != 0) {
        WARN("trix_writeSurfTrix failed for \"%s\"\n", out_fname);
        exit(FILE_ERROR);
      }
      NOTE("Wrote %s in extended VTK format\n", out_fname);
    }
    else {
      rc = io_writeSurfTri(out_fname, p_surf, is_ascii);
      if (0 != rc){
        ERR("writing file failed\n");
      }
      NOTE("Wrote %s in traditional format\n", out_fname);
    }
  }

  exit(rc);
}
