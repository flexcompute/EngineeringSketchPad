
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
 * $Id: c3dio_lib.c,v 1.79 2022/11/07 18:45:42 mnemec Exp $
 */

/* open source */

/**
 * I/O library for Cart3D
 */

/* ----- include files ------ */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h> /* for fsync and access */

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include "c3dio_lib.h"      /* <- Lib datatypes, macros, and prototypes */

#define BADINDEX_ERROR   -7

#define REPORTSTATUS(A,F_){                                           \
    if (verbose_io){              /*   ...print a progress report */  \
      printf("              ...Reading %s: ",F_);                     \
      printf((A));                                                    \
      printf("               \r"); fflush(stdout);                    \
    }                                                                 \
  }

#define OUTREPORTSTATUS(A,F_){   /*   ...print a progress report */   \
    if (verbose_io){                                                  \
      printf("              ...Writing %s: ",F_);                     \
      printf((A));                                                    \
      printf("               \r"); fflush(stdout);                    \
    }                                                                 \
  }

/* ------------- private prototypes ------------- */

int   ioPrivate_tokenizeLine(char *, float * );
void  triform_(char *, int *, int *, int *, int *, int *,int *);
int   ioPrivate_intLog2(const double);
void  read_tri_file_(char*, int*, int*, int*, int*, float*, float*, float*,
                     int*, int  *, int  *, int  *, int*);
void  read_triq_file_(char*,  int*, int*, int*, int*, int*, float*,float*,
                      float*,int*, int*, int*, int*, float*, int*);
void  ioPrivate_copyHinfo(const p_tsHinfo p_orig, const p_tsHinfo p_dest);

/* ------ vars public to library ---- */
bool verbose_io;

/**
 * -------------readCNTLfile()----....parse control file, parsing is done by
 *                                    case and keyword testing, using a
 *                                    while(!EOF){} control structure.
 *
 *                     NOTES:    1. Sections can appear in any order
 *                               2. keywords within sections can be in
 *                                  any order
 *                               3. comments begin with a "#"
 *                               4. Sections start with a "$"
 *
 *                     CURRENT SECTION TITLES:
 *                                 $__Case_Information:
 *                                 $__Solver_Control_Information:
 *                                 $__File_Name_Information:
 *                                 $__Convergence_History_reporting:
 *                                 $__Partiion_Information:
 *                                 $__Post_Processing:
 *                                 $__Boundary_Conditions:
 */
#define KEY_MATCHES(S) (0 == strcmp(name,S))

int io_readCNTLfile(p_tsControl p_cntl){

  char *inputCNTLfile = p_cntl->p_fileInfo->InputCntlFile;
  FILE *p_strm        = NULL;
  char  line[STRING_LEN];
  int   j,istage = 0, nOrders = UNSET;
  int   rc;
  char  name[100], value[100], value2[100], value3[100], value4[100];
  char  section_name[100];
  bool  in_unrecognized_section = FALSE;
  int   td_Method = UNSET;
  tsTDinfo td;
                                                      /*...parse categories  */
  p_tsCinfo  p_C   =&p_cntl->CaseInfo;
  p_tsSinfo  p_S   =&p_cntl->SolverInfo;
  p_tsIOinfo p_F   = p_cntl->p_fileInfo;
  p_tsMGinfo p_MG  = p_cntl->p_mgInfo;
  p_tsHinfo  p_H   = p_cntl->p_convergeInfo;
  p_tsPPinfo p_PP  = p_cntl->p_ppInfo;
  p_tsTDinfo p_TD  = &td;

  p_PP->nSlices[X] = p_PP->nSlices[Y] = p_PP->nSlices[Z] = 0;

  /* Set defaults on reference state for backwards compatibility with
     old input files that specify rhoinf and ainf */
  p_C->rhoinf = 1.0;
  p_C->ainf = 1.0;
  p_C->Minf = 0.0;
  p_C->alpha = 0.0;
  p_C->beta = 0.0;
  p_C->roll = 0.0;

                                           /* ...open file and begin parsing */
  if ((p_strm = fopen(inputCNTLfile,"r")) == NULL) {
    ERR("error opening parameter file %s\n",inputCNTLfile);
    exit(FILE_ERROR);
  }
  /* initialize line: */
  line[0]='\0';
  /* ...main parse loop */
  while (fgets(line, STRING_LEN, p_strm) != NULL) {

    /* Skip any lines starting with '#' or ' ' */
    if (line[0]=='#' || line[0]==' ' || line[0]=='\n') {
      line[0]='\0';
      continue;
    }
    /* ...read section title */
    if (line[0]=='$') {
      if (sscanf(line,"%s\n",section_name)!=1) {
        printf("Error reading section name\n");
        return PARSE_ERROR;
      }
      in_unrecognized_section = FALSE;
      continue;
    }

    /* Wait until we encounter a section that we recognize */
    if (in_unrecognized_section)
      continue;

    /*..."RK" line setup  */
    if (strcmp(section_name,"$__Solver_Control_Information:")==0 &&
        line[0]=='R' && line[1]=='K') {
      if (sscanf(line,"%s %s %s\n",name, value, value2)!=3) {
        printf("Error reading parameter file %s\n",inputCNTLfile);
        printf("line='%s', line[0]='%c'\n",line,line[0]);
        continue;
      }                                         /* ..."Dir_Lo_Hi" line setup*/
    } else if (strcmp(section_name,"$__Boundary_Conditions:")==0 &&
               'D'==line[0] && 'i'==line[1] && 'r'==line[2] &&
               'L'==line[4] && 'H'==line[7]){
      if (sscanf(line,"%s %s %s %s\n",name, value, value2, value3)!=4) {
        printf("Error reading parameter file %s\n",inputCNTLfile);
        printf("line='%s', line[0]='%c'\n",line,line[0]);
        continue;
      }                                       /* ..."Froude" line setup */
    } else if (strcmp(section_name,"$__Case_Information:")==0 &&
               'F'==line[0] && 'r'==line[1] && 'o'==line[2] &&
               'u'==line[3] && 'd'==line[4] && 'e'==line[5]) {
      if (sscanf(line,"%s %s %s %s %s\n",name, value, value2, value3, value4)!=5) {
        printf("Error reading parameter file %s\n",inputCNTLfile);
        printf("line='%s', line[0]='%c'\n",line,line[0]);
        continue;
      }                              /* ...All other lines are name, value pairs */
    } else if (sscanf(line,"%s %s\n",name, value)!=2) {
      printf("Error reading parameter file %s\n",inputCNTLfile);
      printf("line='%s', line[0]='%c'\n",line,line[0]);
      continue;
    }
                                          /*... $__Case_Information: Parsing */
    if (strcmp(section_name,"$__Case_Information:")==0) {
      rc = 1;
      if (KEY_MATCHES("Restart")) {
        if (sscanf(value,"%d",  &p_C->restart )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("Mach")) {
        if (sscanf(value,"%lf", &p_C->Minf    )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("alpha")) {
        if (sscanf(value,"%lf", &p_C->alpha   )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("beta")) {
        if (sscanf(value,"%lf", &p_C->beta    )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("roll")) {
        if (sscanf(value,"%lf", &p_C->roll    )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("Rho_inf")) {
        if (sscanf(value,"%lf", &p_C->rhoinf  )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("gamma")) {
        if (sscanf(value,"%lf", &p_C->gamma   )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("Froude")) { // ...bf: handle Froude & gravity
        if (sscanf(value, "%lf",&p_C->Froude)!=1) return PARSE_ERROR;
        if (sscanf(value2,"%lf",&p_C->gravity[X])!=1) return PARSE_ERROR;
        if (sscanf(value3,"%lf",&p_C->gravity[Y])!=1) return PARSE_ERROR;
        if (sscanf(value4,"%lf",&p_C->gravity[Z])!=1) return PARSE_ERROR;
      } else {
        printf("Warning:  Variable '%s' is not recognized by parser\n",name);
      }
    }                                      /* ...$__File_Name_Information: */
    else if (strcmp(section_name,"$__File_Name_Information:")==0) {
      if (KEY_MATCHES("MeshInfo")) {
        strcpy( p_F->MeshInfoFile  ,value);
      } else if (KEY_MATCHES("MeshFile")) {
        strcpy( p_F->InputMeshFile ,value);
      } else if (KEY_MATCHES("HistoryFile")) {
        strcpy( p_F->HistoryFile   ,value);
      } else {
        printf("Warning:  Variable '%s' is not recognized by parser\n",name);
      }
    }                                    /* ...$__Solver_Control_Information:*/
    else if (strcmp(section_name,"$__Solver_Control_Information:")==0) {
      rc = 1;
      if        (KEY_MATCHES("CFL")) {
        if (sscanf(value,"%lf", &p_S->cfl          )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("Limiter")) {
        if (sscanf(value,"%d",  &p_S->limiter      )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("FluxFun")) {
        if (sscanf(value,"%d",  &p_S->fluxFunction )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("Precon")) {
        if (sscanf(value,"%d",  &p_S->pc           )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("wallBCtype")) {
        if (sscanf(value,"%d",  &p_S->doSubcellSurf)!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("nMGlev")) {
        if (sscanf(value,"%d",  &p_MG->nMGlev      )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("maxCycles")) {
        if (sscanf(value,"%d",  &p_MG->maxCycles   )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("MG_cycleType")) {
        if (sscanf(value,"%d",  &p_MG->cycleType   )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("MG_nPre")) {
        if (sscanf(value,"%d",  &p_MG->nPre        )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("MG_nPost")) {
        if (sscanf(value,"%d",  &p_MG->nPost       )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("TD_nSteps")) {
        if (sscanf(value,"%d",  &p_TD->nTDSteps    )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("TD_dt_phys")) {
        if (sscanf(value,"%lf", &p_TD->dt_p        )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("TD_method")) {
        if (sscanf(value,"%d",  &td_Method         )!=1) return PARSE_ERROR;
      }  else
        if (KEY_MATCHES("RK")) {
          if (sscanf(value,"%lf",&p_S->a_stageCoef[istage])!=1)
            return PARSE_ERROR;
          if (sscanf(value2,"%d",&p_S->a_gradEval [istage])!=1)
            return PARSE_ERROR;
          istage++;                            /* <- running stage counter */
        } else {
          printf("Warning:  Variable '%s' is not recognized by parser\n",name);
        }
      if (rc != 1)
        return PARSE_ERROR;
      /* ...$__Convergence_History_reporting:*/
    } else if (strcmp(section_name,"$__Convergence_History_reporting:")==0) {
      if (KEY_MATCHES("iForce")) {
        if (sscanf(value,"%d", &p_H->iForce)!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("iHist")) {
        if (sscanf(value,"%d", &p_H->iHist )!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("nOrders")) {
        if (sscanf(value,"%d",&nOrders)!=1) return PARSE_ERROR;
        p_H->initResidual   = 1.;
        p_H->targetResidual = pow(10,(-1.*(double)nOrders));
      } else {
        printf("Warning:  Variable '%s' is not recognized by parser\n",name);
      }
    }                                       /* ... $__Partition_Information: */
    else if (strcmp(section_name,"$__Partition_Information:")==0) {
      if (KEY_MATCHES("nPart")) {
        if (sscanf(value,"%d",&p_cntl->nSubDomains)!=1) return PARSE_ERROR;
      } else if (KEY_MATCHES("type")) {
        if (sscanf(value,"%d",&p_cntl->PartType)!=1)    return PARSE_ERROR;
      } else {
        printf("Warning:  Variable '%s' is not recognized by parser\n",name);
      }
    }                                      /* ... $__Post_Processing:   */
    else if (strcmp(section_name,"$__Post_Processing:")==0) {
      if        (KEY_MATCHES("Xslices")) {
        if ((p_PP->nSlices[X] = ioPrivate_tokenizeLine(line,p_PP->a_slicesX)) <= 0)
          return PARSE_ERROR;
      } else if (KEY_MATCHES("Yslices")) {
        if ((p_PP->nSlices[Y] = ioPrivate_tokenizeLine(line,p_PP->a_slicesY)) <= 0)
          return PARSE_ERROR;
      } else if (KEY_MATCHES("Zslices")) {
        if ((p_PP->nSlices[Z] = ioPrivate_tokenizeLine(line,p_PP->a_slicesZ)) <= 0)
          return PARSE_ERROR;
      } else {
        /* Suppress messages about known variables that are simply not parsed
           in the library. */
        if (!KEY_MATCHES("lineSensor") &&
            !KEY_MATCHES("eaSensor") &&
            !KEY_MATCHES("pointSensor"))
          {
            printf("Warning:  Variable '%s' is not recognized by parser\n",name);
          }
      }
      if (MAX(MAX(p_PP->nSlices[X],p_PP->nSlices[Y]),p_PP->nSlices[Z]) > MAXSLICES) {
        ERR("Too many slices requested in X/Y/Z (requested %d/%d/%d)\n",
            p_PP->nSlices[X],p_PP->nSlices[Y],p_PP->nSlices[Z]);
        CONT("Max num slices in each direction is limited to %d\n",MAXSLICES);
        CONT("Please contact the development team if this is too restrictive\n");
        return PARSE_ERROR;
      }
    }                                      /* ... $__Boundary_Conditions:   */
    else if (strcmp(section_name,"$__Boundary_Conditions:")==0) {
      int idir,Lo,Hi;
      if (KEY_MATCHES("Dir_Lo_Hi")) {     /* ...Simple BC on cartesian box  */
        if (sscanf(value ,"%d",&idir)!=1) return PARSE_ERROR;
        if (sscanf(value2,"%d",&Lo  )!=1) return PARSE_ERROR;
        if (sscanf(value3,"%d",&Hi  )!=1) return PARSE_ERROR;
        if (DIM <= idir || 0 > idir){
          ERR("Bad entry for idir in %s\n",name);
          return PARSE_ERROR;
        }
        p_S->bboxBCs[idir*2  ] = Lo;     /* [LoX,HiX, LoY,HiY, LoZ,HiZ]    */
        p_S->bboxBCs[idir*2+1] = Hi;
      } else {
        /* Suppress warning for sections that are recognized but not parsed */
        if (!KEY_MATCHES("SurfBC") &&
            !KEY_MATCHES("InletBC") &&
            !KEY_MATCHES("ExitBC"))
          {
            printf("Warning:  Variable '%s' is not recognized by parser\n",name);
          }
      }
    } else {
      in_unrecognized_section = TRUE;
      /* Suppress warning for sections that are recognized but not parsed */
      if ((0 == strcmp(section_name,"$__Force_Moment_Processing:")) ||
          (0 == strcmp(section_name,"$__Design_Info:")) ||
          (0 == strcmp(section_name,"$__Steering_Information"))){
        printf("Skipping %s section\n",section_name);
      } else {
        printf("Warning:  Section name '%s' not recognized\n", section_name);
      }

    }

    line[0]='\0';                                  /* ...re-initialize line */
  } /* -- end while loop-- */

                       /* Done parsing, check for any errors or missing info */
  if (UNSET == nOrders){
    WARN(" # of orders for convergence criteria not set, defaulting to %f\n",
         log10(p_H->targetResidual));
  }

                              /* Set number of RK stages based on input file */
  p_S->nStage = istage;
  for(j=0;j<p_S->nStage;j++){
    if (0 != p_S->a_gradEval[j]){
      p_S->first_order = FALSE;
      break;
    }
  }

  if (p_strm)
    (void)fclose(p_strm);

  return 0;
}
#undef KEY_MATCHES

/**
 * ------ ioPrivate_tokenizeLine()---...takes a string and keeps parsing
 *                                      off floats till it hits something
 *                                      else --------------------------------
 */
int ioPrivate_tokenizeLine(char * string, float * a_floatArray){
  int     nTokens = 0;
  int     status  = 1;
  char *p_space   = &string[0];

  /* ...stop if sscanf() hits a non-float, or if    */
  /*    you try to get beyond the end of the string */
  while(status == 1 && p_space < (&string[0]+STRING_LEN)){
    while(p_space[0] != ' ') p_space++;     /* <--skip the current token    */
    while(p_space[0] == ' ') p_space++;     /* <--skip to end of whitespace */
    status = sscanf(p_space,"%f",&a_floatArray[nTokens]);
    nTokens++;
  }
  nTokens--;              /* <- last one was "error" so uncount it */
  return nTokens;
}

/**
 * --------strmMetaInfo----------------...stream to file descriptive info
 *                                        about the case we're working on
 *                                        -------------------------------
 */
#define DEG(A) ((A)*180./pi)
FILE * io_strmMetaInfo(p_tsControl p_cntl, FILE *p_strm){
  int j;
  p_tsCinfo  p_C    =&p_cntl->CaseInfo;
  p_tsSinfo  p_S    =&p_cntl->SolverInfo;
  p_tsIOinfo p_io   = p_cntl->p_fileInfo;
  p_tsMGinfo p_MG   = p_cntl->p_mgInfo;
  double    pi = atan(1.)*4.;

  fprintf(p_strm,"# cmd: %s\n",p_C->cmdLine);
  fprintf(p_strm,"#\n");
  fprintf(p_strm,"# Surface Triangulation file: %s\n",p_io->SurfTriFile);
  fprintf(p_strm,"#                  Mesh file: %s\n",p_io->InputMeshFile);
  fprintf(p_strm,"#         Input control file: %s\n",p_io->InputCntlFile);
  fprintf(p_strm,"#\n");
  fprintf(p_strm,
          "# Mach No = %6.4f,  alpha = %5.3f,  beta = %6.4f,  gamma = %4.3g\n",
          p_C->Minf, DEG(p_C->alpha), DEG(p_C->beta),p_C->gamma);
  fprintf(p_strm,"#     cfl = %4.2f,     FluxFun = %d,   Limiter = %d, \n" ,
          p_S->cfl, p_S->fluxFunction, p_S->limiter);
  fprintf(p_strm,
          "#  nMGlev = %d,    MGcycleType = %d,   MG_nPre = %d, MG_nPost = %d\n",
          p_MG->nMGlev,p_MG->cycleType,p_MG->nPre,p_MG->nPost);
  fprintf(p_strm,"#\n");
  for(j=0;j<p_S->nStage;j++)
    fprintf(p_strm,
            "#      RK   %7.5f         %d   # stageCoef GradEval \n",
            p_S->a_stageCoef[j],p_S->a_gradEval[j]);
  fprintf(p_strm,"#\n");


  return p_strm; /* ...send back the file pointer */
}
#undef DEG

/**
 * --------io_newC3Dcntl(): mallocs memory for pointer in cntl struct
 */
void io_newC3Dcntl(p_tsControl p_c3dInfo){
  NEW(p_c3dInfo->p_convergeInfo, tsHinfo);
  NEW(p_c3dInfo->p_fileInfo,     tsIOinfo);
  NEW(p_c3dInfo->p_mgInfo,       tsMGinfo);
  NEW(p_c3dInfo->p_ppInfo,       tsPPinfo);

  strcpy(p_c3dInfo->p_fileInfo->InputCntlFile , "input.cntl"        );
  strcpy(p_c3dInfo->p_fileInfo->InputMeshFile , "Mesh.c3d"          );
  strcpy(p_c3dInfo->p_fileInfo->RestartFile ,   "preAdapt.ckpt"     );
  strcpy(p_c3dInfo->p_fileInfo->OutputSolnFile, "postAdapt.ckpt"    );
  strcpy(p_c3dInfo->p_fileInfo->MeshInfoFile  , "Mesh.c3d.Info"     );
  strcpy(p_c3dInfo->p_fileInfo->OutputMeshFile, "adaptedMesh.R.c3d" );
  strcpy(p_c3dInfo->p_fileInfo->SurfTriFile,    "Components.i.tri"  );
  strcpy(p_c3dInfo->p_fileInfo->PreSpecFile,    "UNSET"             );

  return;
}

/**
 * ----------io_freeC3Dcntl(): free memory in cntl struct (including
 *                             triangulation)
 */
void io_freeC3Dcntl(p_tsControl p_c3dInfo) {
  if (p_c3dInfo->p_convergeInfo) {
    free(p_c3dInfo->p_convergeInfo);
    p_c3dInfo->p_convergeInfo=NULL;
  }
  if (p_c3dInfo->p_fileInfo) {
    free(p_c3dInfo->p_fileInfo);
    p_c3dInfo->p_fileInfo = NULL;
  }
  if (p_c3dInfo->p_mgInfo) {
    free(p_c3dInfo->p_mgInfo);
    p_c3dInfo->p_mgInfo = NULL;
  }
  if (p_c3dInfo->p_ppInfo) {
    free(p_c3dInfo->p_ppInfo);
    p_c3dInfo->p_ppInfo = NULL;
  }
  if (p_c3dInfo->GridInfo.p_surf) {
    c3d_freeTriangulation(p_c3dInfo->GridInfo.p_surf, 0);
    free(p_c3dInfo->GridInfo.p_surf);
    p_c3dInfo->GridInfo.p_surf=NULL;
  }
  return;
}

/**
 * ------io_freeBasicGrid(): free memory of basic grid struct,
 */
void io_freeBasicGrid(p_tsBasicGrid p_g){
  if (p_g->tPolys.p_tpEntry)      free(p_g->tPolys.p_tpEntry);
  if (p_g->tPolys.p_tpIntTriList) free(p_g->tPolys.p_tpIntTriList);
  if (p_g->tPolys.p_tpCentroids)  free(p_g->tPolys.p_tpCentroids);
  if (p_g->tPolys.p_tpAreas)      free(p_g->tPolys.p_tpAreas);

  if (p_g->a_Cells)  free(p_g->a_Cells);
  if (p_g->a_Faces)  free(p_g->a_Faces);

  if (p_g->a_cCells) free(p_g->a_cCells);
  if (p_g->a_cFaces) free(p_g->a_cFaces);

  return;
}

/**
 * ---------------- initializeInfoStructs() --------------------------------
 */
void io_initializeInfoStructs(p_tsControl p_cntl){
  p_tsMGinfo p_mgInfo = p_cntl->p_mgInfo;
  p_tsCinfo  p_C      =&p_cntl->CaseInfo;
  p_tsSinfo  p_S      =&p_cntl->SolverInfo;
  p_tsHinfo  p_H      = p_cntl->p_convergeInfo;
  p_tsMGinfo p_MG     = p_mgInfo;
  int        i;
  /*                                  ...Case Info defaults                  */
  p_C->restartMP = p_C->restart = FALSE;
  p_C->Minf      = p_C->alpha  = p_C->beta = p_C->qinf = 0.0;
  p_C->gamma     = 1.4;   p_C->ainf  = p_C->rhoinf = 1.0;
  p_C->pinf      = 1./p_C->gamma;
  /*           ...File Info defaults          */
  /* should be done elsewhere-usually in main */
  /* to allow cmd line parsing of filenames   */

                                                  /* ...Solver Info defaults */
  for(i=0;i<MAXNUMSTAGES;i++)
    {p_S->a_stageCoef[i]=UNSET;p_S->a_gradEval[i]=FALSE;}
  p_S->nStage        = p_S->cfl = UNSET;
  p_S->pc            = (pretype)0;
  p_S->limiter       = NOLIMITER;
  p_S->fluxFunction  = (fftype)0;
  p_S->doSubcellSurf = FALSE;
  p_S->first_order   = TRUE;
  p_S->bboxBCs[0]    =  p_S->bboxBCs[1] = p_S->bboxBCs[2] =
    p_S->bboxBCs[3]    =  p_S->bboxBCs[4] = p_S->bboxBCs[5] = 0;

                                            /* ...MultiGrid info defaults    */
  p_MG->nMGlev    = 1;
  p_MG->nPre      = p_MG->nPost = p_MG->nCycles = p_MG->maxCycles = 0;
  p_MG->maxCycles = p_MG->cycleType = 0;

                                            /* ...Converg Hist Info Defaults */
  p_H->initResidual  = 1.;
  p_H->targetResidual = 1.E-8;
  p_H->iForce = 1; p_H->iHist = 1;


  p_cntl->nSubDomains = 1;  p_cntl->PartType = 1;

  return ;
}

/**
 *  ------------------- readGInfo() ---------....gets meta info about
 *                                               mesh from Mesh.c3d.Info
 *                                               file generated by cubes
 *                                               ---------------------------
 */
#define BUFFER_LEN 256                      /*...set consistent with cubes! */

void io_readGinfo(p_tsControl p_cntl) {

  const     p_tsGinfo p_GridInfo = &p_cntl->GridInfo;
  int       idim, minAllow, j,  minMaxAllowRef;
  int       maxAR, nBufLayers, size, minSurfRefine, NumRefBits;
  int       meshInternal, weightNorms;
  FILE     *p_strm = NULL;
  char      buffer[BUFFER_LEN];
  char      s1[BUFFER_LEN], s2[BUFFER_LEN];
  char     *p_MeshInfoFile = p_cntl->p_fileInfo->MeshInfoFile;
  char     *p_surfTriFile  = p_cntl->p_fileInfo->SurfTriFile;

  if (NULL ==(p_strm= fopen(p_MeshInfoFile,"r")) ){ /* open for reading*/
    ERR( "Could not find input data file '%s'\n",p_MeshInfoFile);
    exit(FILE_ERROR);
  }

  for (j=0;j<4;j++)                              /* skip first four lines */
    fgets(buffer, BUFFER_LEN, p_strm);

                          /* these are the input.c3d vars that made the mesh */
  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%s ",p_surfTriFile);

  for(j=0;j<DIM;j++) {                        /* read bbox for each dir */
    fgets(buffer, BUFFER_LEN, p_strm);
    sscanf(buffer,"%s %s ",s1,s2);
    p_GridInfo->minBound[j] = atof(s1);
    p_GridInfo->maxBound[j] = atof(s2);
  }

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&p_GridInfo->nBits);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&maxAR);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&p_GridInfo->finestCellLevel);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&minSurfRefine);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&meshInternal);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&weightNorms);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d ",&nBufLayers);

  fgets(buffer, BUFFER_LEN, p_strm);
  sscanf(buffer,"%d %d  %d ",
         &p_GridInfo->nDiv[X],&p_GridInfo->nDiv[Y],&p_GridInfo->nDiv[Z]);

  if (p_strm) (void)fclose(p_strm);

                             /* file read.  set some internal variables next */
  if (1 == maxAR)
    p_GridInfo->isotropic = TRUE;
  else
    p_GridInfo->isotropic = FALSE;

  p_GridInfo->meshInternal = (1 == meshInternal) ? TRUE  : FALSE;

  NumRefBits = p_GridInfo->nBits;
  for (idim=0;idim<DIM;idim++){
    p_GridInfo->maxAllowRef[idim] = ioPrivate_intLog2( (double)(
                                                                ((float)TWO_TO_THE(p_GridInfo->nBits)-1.)
                                                                /(p_GridInfo->nDiv[idim]-1.)) );
    if (p_GridInfo->maxAllowRef[idim] < NumRefBits)
      NumRefBits = p_GridInfo->maxAllowRef[idim];
  }

  minAllow = 100000;
  for (idim=0;idim<DIM;idim++)
    minAllow = MIN(minAllow , p_GridInfo->maxAllowRef[idim]);

  if (p_GridInfo->isotropic) { /*  insure all directions have same limits */
    for (idim=0;idim<DIM;idim++) p_GridInfo->maxAllowRef[idim] = minAllow;
  }

  for (j=0;j<DIM;j++)
    p_GridInfo->coarseIntDelta[j] = TWO_TO_THE(p_GridInfo->maxAllowRef[j]);

  for(j=0;j<DIM;j++) {

    p_GridInfo->M[j] = TWO_TO_THE(p_GridInfo->maxAllowRef[j])
      * (p_GridInfo->nDiv[j] - 1) + 1;
    p_GridInfo->fine_spacing[j] = (((double)p_GridInfo->maxBound[j])-
                                   ((double)p_GridInfo->minBound[j]))/
      ((double)p_GridInfo->M[j] - 1);
  }

  minMaxAllowRef = MIN(p_GridInfo->maxAllowRef[X],p_GridInfo->maxAllowRef[Y]);
  minMaxAllowRef = MIN(minMaxAllowRef,p_GridInfo->maxAllowRef[Z]);
  for (j=0;j<=minMaxAllowRef;j++){
    for (idim=0;idim<DIM;idim++){
      size = (int)TWO_TO_THE(p_GridInfo->maxAllowRef[idim] - j);
      p_GridInfo->h[j][idim] = size * p_GridInfo->fine_spacing[idim];
    }
  }

  return;
}
#undef BUFFER_LEN

/* ---------ioPrivate_intLog2()--------...return int(log2(x) ------------ */

int ioPrivate_intLog2(const double d){
  int i=0;
  while( (double)(TWO_TO_THE(i)) < d ) i++;
  return (i-1);
}

/**
 * --c3d_is_tri_binary(): reads the first line of a tri file to figure
 *                        out if it is ascii or binary, tri or
 *                        triq. Returns the length of the record
 *                        marker in bytes if the file is binary,
 *                        either 4 or 8, returns 0 if it is ascii
 *                        tri-file, 2 if it is ascii triq-file,
 *                        negative means error
 */
int c3d_is_tri_binary(const char * fileName, const int io_opts) {
  FILE *p_strm;
  int rlen = -1;
  int rec41, rec42, nVerts, nTris;
  long rec81, rec82;
  char ascii_line[STRING_LEN];
  char *buff;

  if ((p_strm = fopen(fileName, "rb")) == NULL) {
    ERR("error opening file %s\n", fileName);
    return FILE_ERROR;
  }
                               /* check for 4 byte markers or 8 byte markers */
  fread(&rec41, 4, 1, p_strm);
  if (rec41 != 8) {
    if (io_opts & IO_VERBOSE) {
      NOTE("Did not find valid 4 byte marker for binary file, trying 8 byte ...\n");
    }
    rewind(p_strm);
    fread(&rec81, 8, 1, p_strm);
    if (rec81 != 8) {
      if (io_opts & IO_VERBOSE) {
        NOTE("Did not find valid 8 byte marker for binary file, trying ascii ...\n");
      }
    }
    else {
      rlen = 8;
    }
  }
  else {
    rlen = 4;
  }

  if (rlen > 0){                          /* if found a valid binary marker */
    fread(&nVerts, 4, 1, p_strm);
    fread(&nTris, 4, 1, p_strm);

    if (rlen == 4) fread(&rec42, rlen, 1, p_strm);
    else fread(&rec82, rlen, 1, p_strm);

    if (rlen == 4 && rec41 != rec42) { /* check for matching record numbers */
      if (io_opts & IO_VERBOSE) {
        NOTE("%d byte record markers did not match, %d and %d, trying ascii ...\n",
             rlen, rec41, rec42);
      }
    }
    else if (rlen == 8 && rec81 != rec82) {
      if (io_opts & IO_VERBOSE) {
        NOTE("%d byte record markers did not match, %ld and %ld, trying ascii ...\n",
             rlen, rec81, rec82);
      }
    }
    else if (nVerts < 4) { /* check for valid nVerts */
      if (io_opts & IO_VERBOSE) {
        NOTE("nVets = %d not valid for triangulation, trying ascii ...\n",
             nVerts);
      }
    }
    else if (nTris < nVerts) { /* check for valid nTris */
      if (io_opts & IO_VERBOSE) {
        NOTE("nTris = %d not valid for nVerts = %d, trying ascii ...\n",
             nTris, nVerts);
      }
    }
    else { /* seems like a good binary file */
      fclose(p_strm);
      return rlen;
    }
  }

  fclose(p_strm);
  /* checking for valid ascii file */
  if ((p_strm = fopen(fileName, "r")) == NULL) {
    ERR("error opening file %s\n", fileName);
    return FILE_ERROR;
  }

  fscanf(p_strm, "%d %d", &nVerts, &nTris);

  if (nVerts < 4) {
    ERR("nVets = %d not valid for triangulation, expected > 3\n",
        nVerts);
    fclose(p_strm);
    return PARSE_ERROR;
  }
  else if (nTris < nVerts) {
    ERR("nTris = %d not valid for nVerts = %d, expected nVerts > nTris\n",
        nTris, nVerts);
    fclose(p_strm);
    return PARSE_ERROR;
  }

  rewind(p_strm);
  fgets(ascii_line, STRING_LEN, p_strm); /* checking for triq file */
  fclose(p_strm);

  buff = strtok(ascii_line, " ");
  buff = strtok(NULL, " ");
  if (buff == NULL) {
    ERR("Unknown file format, found one entry in first line.\n");
    return PARSE_ERROR;
  }
  buff = strtok(NULL, " ");
  if (buff == NULL) { /* there were two entries in first line */
    /* seems like a good ascii file */
    return 0;
  }
  if (*buff == '\n') {
    /* two entries on first line, terminating with white space and newline */
    /* seems like a good ascii file, tri2comp generates these files */
    return 0;
  }
  buff = strtok(NULL, " ");
  if (buff == NULL) {              /* there were three entries in first line */
    /* try triq file */
    WARN("Three entries found on first line, parsing triq file\n");
    return 2;
  }

  ERR("Unkown file format, found more than three entries on first line\n");
  return PARSE_ERROR;
} /* end c3d_is_tri_binary */

/**
 * c3d_read_tri_ascii(): read tri files in ascii format into p_surf
 */
int c3d_read_tri_ascii(const char *fileName, p_tsTriangulation p_surf,
                       const int io_opts) {
  FILE *p_strm;
  int i, comp, v1, v2, v3;

  if ((p_strm = fopen(fileName, "r")) == NULL) {
    ERR("error opening file %s\n", fileName);
    return FILE_ERROR;
  }
                          /* expect the next two ints to be nVerts and nTris */
  if ( 2 != fscanf(p_strm, "%d %d", &p_surf->nVerts, &p_surf->nTris) ) {
    ERR("c3d_read_tri_ascii: did not find nVerts nTris\n");
    return PARSE_ERROR;
  }

  c3d_allocTriangulation(&p_surf);

  for (i=0; i < p_surf->nVerts; ++i){ /* skip error checking for speed */
    fscanf(p_strm, "%f %f %f", p_surf->a_Verts[i].x + X,
           p_surf->a_Verts[i].x + Y, p_surf->a_Verts[i].x + Z);
  }

  for (i=0; i < p_surf->nTris; ++i){
    fscanf(p_strm, "%d %d %d", &v1, &v2, &v3);
    p_surf->a_Tris[i].vtx[0] = v1 - 1;
    p_surf->a_Tris[i].vtx[1] = v2 - 1;
    p_surf->a_Tris[i].vtx[2] = v3 - 1;
  }

  for (i=0; i < p_surf->nTris; ++i){
    if ( 1 != fscanf(p_strm, "%d", &comp) ) {
      p_surf->a_Tris[i].Comp = 0; /* initialize Comp field */
    }
    else {
      p_surf->a_Tris[i].Comp = comp - 1;
    }
  }

  fclose(p_strm);
  return 0;
}

/**
 * c3d_read_tri_binary(): read tri files in binary format into p_surf
 */
int c3d_read_tri_binary(const char *fileName, p_tsTriangulation p_surf,
                        const int rlen, const int io_opts) {
  FILE *p_strm;
  int rec41, rec42, i, comp;
  long rec81, rec82;

  if ((p_strm = fopen(fileName, "rb")) == NULL) {
    ERR("error opening file %s\n", fileName);
    return FILE_ERROR;
  }

  if (rlen == 4) {
    fread(&rec41, rlen, 1, p_strm);
    if (rec41 != 8) {
      ERR("Unexpected record marker value for header: %d, expected: 8\n",
          rec41);
      fclose(p_strm);
      return PARSE_ERROR;
    }
  }
  else {
    fread(&rec81, rlen, 1, p_strm);
    if (rec81 != 8){
      ERR("Unexpected record marker value for header: %ld, expected: 8\n",
          rec81);
      fclose(p_strm);
      return PARSE_ERROR;
    }
  }
                          /* expect the next two ints to be nVerts and nTris */
  fread(&p_surf->nVerts, 4, 1, p_strm);
  fread(&p_surf->nTris, 4, 1, p_strm);
                                                      /* two markers up next */
  if (rlen == 4){
                                                 /* ending marker for header */
    fread(&rec42, rlen, 1, p_strm);
    if (rec42 != rec41){
      ERR("Unexpected record marker values: %d and %d\n", rec41, rec41);
      fclose(p_strm);
      return PARSE_ERROR;
    }
                                     /* starting record marker for verticies */
    fread(&rec41, rlen, 1, p_strm);
    if (rec41 != p_surf->nVerts*12) {
      ERR("Unexpected record marker value for verticies: %d, expected: %d\n",
          rec41, p_surf->nVerts*12);
      fclose(p_strm);
      return PARSE_ERROR;
    }
  }
  else{
                                                 /* ending marker for header */
    fread(&rec82, rlen, 1, p_strm);
    if (rec82 != rec81) {
      ERR("Unexpected record marker values: %ld and %ld\n", rec81, rec81);
      fclose(p_strm);
      return PARSE_ERROR;
    }
                                     /* starting record marker for verticies */
    fread(&rec81, rlen, 1, p_strm);
    if (rec81 != p_surf->nVerts*12) {
      ERR("Unexpected record marker value for verticies: %ld, expected: %d\n",
          rec81, p_surf->nVerts*12);
      fclose(p_strm);
      return PARSE_ERROR;
    }
  }

  c3d_allocTriangulation(&p_surf);

  for (i = 0; i < p_surf->nVerts; ++i) {                   /* read verticies */
    fread(p_surf->a_Verts[i].x, 4, 3, p_strm);
  }
                                                         /* two more markers */
  if (rlen == 4){
                                                  /* ending marker for verts */
    fread(&rec42, rlen, 1, p_strm);
    if (rec42 != rec41){
      ERR("Unexpected record marker values: %d and %d\n", rec41, rec41);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
                                     /* starting record marker for triangles */
    fread(&rec41, rlen, 1, p_strm);
    if (rec41 != p_surf->nTris*12) {
      ERR("Unexpected record marker value for triangles: %d, expected: %d\n",
          rec41, p_surf->nTris*12);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
  }
  else {
                                                  /* ending marker for verts */
    fread(&rec82, rlen, 1, p_strm);
    if (rec82 != rec81){
      ERR("Unexpected record marker values: %ld and %ld\n", rec81, rec81);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
                                     /* starting record marker for triangles */
    fread(&rec81, rlen, 1, p_strm);
    if (rec81 != p_surf->nVerts*12) {
      ERR("Unexpected record marker value for triangles: %ld, expected: %d\n",
          rec81, p_surf->nTris*12);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
  }

  for (i = 0; i < p_surf->nTris; ++i) {                    /* read triangles */
    fread(p_surf->a_Tris[i].vtx, 4, 3, p_strm);
    p_surf->a_Tris[i].vtx[0]--;
    p_surf->a_Tris[i].vtx[1]--;
    p_surf->a_Tris[i].vtx[2]--;
  }
                                                         /* two more markers */
  if (rlen == 4) {
                                              /* ending marker for triangles */
    fread(&rec42, rlen, 1, p_strm);
    if (rec42 != rec41){
      ERR("Unexpected record marker values: %d and %d\n", rec41, rec41);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }

                                    /* starting record marker for components */
    fread(&rec41, rlen, 1, p_strm);
    if (rec41 != p_surf->nTris*4) {
      ERR("Unexpected record marker value for components: %d, expected: %d\n",
          rec41, p_surf->nTris*4);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
  }
  else {
                                              /* ending marker for triangles */
    fread(&rec82, rlen, 1, p_strm);
    if (rec82 != rec81) {
      ERR("Unexpected record marker values: %ld and %ld\n",
          rec81, rec81);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
                                    /* starting record marker for components */
    fread(&rec81, rlen, 1, p_strm);
    if (rec81 != p_surf->nVerts*12) {
      ERR("Unexpected record marker value for components: %ld, expected: %d\n",
          rec81, p_surf->nTris*4);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
  }

  for (i=0; i < p_surf->nTris; ++i){
    if (fread(&comp, 4, 1, p_strm) != 1) p_surf->a_Tris[i].Comp = 0;
    else p_surf->a_Tris[i].Comp = comp - 1;
  }
                                                 /* expect the ending marker */
  if (rlen == 4){
                                             /* ending marker for components */
    fread(&rec42, rlen, 1, p_strm);
    if (rec42 != rec41){
      ERR("Unexpected record marker values: %d and %d\n", rec41, rec41);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
  }
  else {
                                             /* ending marker for components */
    fread(&rec82, rlen, 1, p_strm);
    if (rec82 != rec81){
      ERR("Unexpected record marker values: %ld and %ld\n", rec81, rec81);
      fclose(p_strm);
      c3d_freeTriangulation(p_surf, 0);
      return PARSE_ERROR;
    }
  }

  fclose(p_strm);
  return 0;
} /* end c3d_read_tri_binary */

int c3d_read_triq_ascii(const char *fileName, p_tsTriangulation p_surf,
                        const int io_opts) {
  FILE *p_strm;
  int v1, v2, v3, i, comp;

  if ((p_strm = fopen(fileName, "r")) == NULL) {
    ERR("error opening parameter file %s\n", fileName);
    return FILE_ERROR;
  }

  if ( 3 != fscanf(p_strm, "%d %d %d", &p_surf->nVerts, &p_surf->nTris,
                   &p_surf->nVertScalars) ) {
    ERR("Wrong triq header in %s\n", fileName);
    return PARSE_ERROR;
  }

  if (p_surf->nVerts < 4){ /* check for valid nVerts */
    ERR("nVets = %d not valid for triangulation\n", p_surf->nVerts);
    fclose(p_strm);
    return PARSE_ERROR;
  }
  else if (p_surf->nTris < p_surf->nVerts){ /* check for valid nTris */
    ERR("nTris = %d not valid for nVerts = %d\n", p_surf->nTris, p_surf->nVerts);
    fclose(p_strm);
    return PARSE_ERROR;
  }
  else if (p_surf->nVertScalars != 6){ /* check for valid nScal */
    ERR("nScal = %d not valid, expecting '6' for a triq file\n",
        p_surf->nVertScalars);
    fclose(p_strm);
    return PARSE_ERROR;
  }

  c3d_allocVertData(&p_surf, 4);

  strcpy(p_surf->p_vertData[0].name,"Cp");
  p_surf->p_vertData[0].dim = 1;
  p_surf->p_vertData[0].offset = 0;
  p_surf->p_vertData[0].type = VTK_Float64;
  p_surf->p_vertData[0].info = TRIX_flowVariable;

  strcpy(p_surf->p_vertData[1].name,"Density");
  p_surf->p_vertData[1].dim = 1;
  p_surf->p_vertData[1].offset = 1;
  p_surf->p_vertData[1].type = VTK_Float64;
  p_surf->p_vertData[1].info = TRIX_flowVariable;

  strcpy(p_surf->p_vertData[2].name,"Velocity");
  p_surf->p_vertData[2].dim = 3;
  p_surf->p_vertData[2].offset = 2;
  p_surf->p_vertData[2].type = VTK_Float64;
  p_surf->p_vertData[2].info = TRIX_flowVariable;

  strcpy(p_surf->p_vertData[3].name,"Pressure");
  p_surf->p_vertData[3].dim = 1;
  p_surf->p_vertData[3].offset = 5;
  p_surf->p_vertData[3].type = VTK_Float64;
  p_surf->p_vertData[3].info = TRIX_flowVariable;

  c3d_allocTriangulation(&p_surf);

  for (i=0; i < p_surf->nVerts; ++i) {
    fscanf(p_strm, "%f %f %f", p_surf->a_Verts[i].x + X,
           p_surf->a_Verts[i].x + Y,
           p_surf->a_Verts[i].x + Z);
  }

  for (i=0; i < p_surf->nTris; ++i) {
    fscanf(p_strm, "%d %d %d", &v1, &v2, &v3);
    p_surf->a_Tris[i].vtx[0] = v1 - 1;
    p_surf->a_Tris[i].vtx[1] = v2 - 1;
    p_surf->a_Tris[i].vtx[2] = v3 - 1;
  }

  for (i=0; i < p_surf->nTris; ++i) {
    fscanf(p_strm, "%d", &comp);
    p_surf->a_Tris[i].Comp = comp - 1;
  }

  for (i=0; i < p_surf->nVerts; ++i) {
    fscanf(p_strm, "%lf %lf %lf %lf %lf %lf", &p_surf->a_scalar0[i],
           &p_surf->a_scalar0[i + p_surf->nVerts],
           &p_surf->a_scalar0[3*i + 2*p_surf->nVerts],
           &p_surf->a_scalar0[3*i+1 + 2*p_surf->nVerts],
           &p_surf->a_scalar0[3*i+2 + 2*p_surf->nVerts],
           &p_surf->a_scalar0[i + 5*p_surf->nVerts]);
  }

  fclose(p_strm);
  return 0;
}

/**
 * io_readSurfTri(): read triangulation into p_surf. Negative value on
 *                   return represents error condition. Rewritten
 *                   Sept-Dec 2020 to avoid Fortran calls
 */
int io_readSurfTri(const char* surfTriFile, p_tsTriangulation p_surf) {
  int rlen, rc;
  int io_opts = 0;

  if (verbose_io) {
    io_opts |= IO_VERBOSE;
  }

  /* use access() from unistd.h to see if file exists to avoid a bug where
   * an empty tri-file is created if it does not exist -- GRA
   */
  if (-1 == access(surfTriFile, F_OK)) {
    ERR("%s does not exist.\n",surfTriFile);
    return FILE_ERROR;
  }
  if (-1 == access(surfTriFile, R_OK)) {
    ERR("%s does not have read access.\n",surfTriFile);
    return FILE_ERROR;
  }

  rlen = c3d_is_tri_binary(surfTriFile, io_opts);
  if ( 0 > rlen ) {
    ERR("c3d_is_tri_binary failed to determine file type for %s\n",surfTriFile);
    rc = rlen;
  }
  else if (0 == rlen) { /* should have ascii tri file */
    rc = c3d_read_tri_ascii(surfTriFile, p_surf, io_opts);
  }
  else if (2 == rlen) {
    rc = c3d_read_triq_ascii(surfTriFile, p_surf, io_opts);
  }
  else {
    rc = c3d_read_tri_binary(surfTriFile, p_surf, rlen, io_opts);
  }

  return rc;
}

/**
 * io_writeSurfTri(): write ascii tri file.  Negative value on return
 *                    represents error condition. Rewritten Sept-Dec
 *                    2020 to avoid Fortran calls. is_ascii flag is
 *                    not used.
 */
int io_writeSurfTri(char *fileName, const p_tsTriangulation p_surf,
                    bool is_ascii) {
  FILE *p_strm = NULL;
  int i;
  int cp_offset = -1, v_offset = -1, rho_offset = -1, p_offset = -1;

  if ( ! is_ascii ) {
    WARN("Binary output not supported, writing ASCII\n");
  }

  if ((p_strm = fopen(fileName, "w")) == NULL) {
    ERR("error opening file %s\n", fileName);
    return FILE_ERROR;
  }

  if ( 0 >= p_surf->nVertScalars ) {        /* write a traditional tri-file */
    fprintf(p_strm, "%12d %12d\n", p_surf->nVerts, p_surf->nTris);

    for (i=0; i<p_surf->nVerts; ++i) {
      /* verts are single precision */
      fprintf(p_strm, "%16.9e %16.9e %16.9e\n", p_surf->a_Verts[i].x[X],
              p_surf->a_Verts[i].x[Y],
              p_surf->a_Verts[i].x[Z]);
    }

    for (i=0; i<p_surf->nTris; ++i) {
      fprintf(p_strm, "%12d %12d %12d\n", p_surf->a_Tris[i].vtx[0] + 1,
              p_surf->a_Tris[i].vtx[1] + 1,
              p_surf->a_Tris[i].vtx[2] + 1);
    }

    for (i=0; i<p_surf->nTris; ++i){
      fprintf(p_strm, "%12d\n", (int)p_surf->a_Tris[i].Comp + 1);
    }
  }
  else { /* writing triq-file */

    /* need to make sure the appropriate scalar values to write out exist */
    for (i = 0; i < p_surf->nVertData; ++i) {
      if (0 == strcmp(p_surf->p_vertData[i].name, "Cp")) {
        cp_offset = p_surf->p_vertData[i].offset;
        if (p_surf->p_vertData[i].dim != 1) {
          WARN("Invalid dimension for Cp, %d, cannot write triq file\n",
               p_surf->p_vertData[i].dim);
          return BADINDEX_ERROR;
        }
      }
      else if (0 == strcmp(p_surf->p_vertData[i].name, "Density")) {
        rho_offset = p_surf->p_vertData[i].offset;
        if (p_surf->p_vertData[i].dim != 1) {
          WARN("Invalid dimension for Density, %d, cannot write triq file\n",
               p_surf->p_vertData[i].dim);
          return BADINDEX_ERROR;
        }
      }
      else if (0 == strcmp(p_surf->p_vertData[i].name, "Velocity")) {
        v_offset = p_surf->p_vertData[i].offset;
        if (p_surf->p_vertData[i].dim != 3) {
          WARN("Invalid dimension for Velocity, %d, cannot write triq file\n",
               p_surf->p_vertData[i].dim);
          return BADINDEX_ERROR;
        }
      }
      else if (0 == strcmp(p_surf->p_vertData[i].name, "Pressure")) {
        p_offset = p_surf->p_vertData[i].offset;
        if (p_surf->p_vertData[i].dim != 1) {
          WARN("Invalid dimension for Pressure, %d, cannot write triq file\n",
               p_surf->p_vertData[i].dim);
          return BADINDEX_ERROR;
        }
      }
    } /* end for */

    if (-1 == cp_offset) {
      WARN("Could not find 'Cp' scalar field, cannot write .triq file\n");
      return BADINDEX_ERROR;
    }
    else if (-1 == v_offset) {
      WARN("Could not find 'Velocity' scalar field, cannot write .triq file\n");
      return BADINDEX_ERROR;
    }
    else if (-1 == rho_offset) {
      WARN("Could not find 'Density' scalar field, cannot write .triq file\n");
      return BADINDEX_ERROR;
    }
    else if (-1 == p_offset) {
      WARN("Could not find 'Pressure' scalar field, cannot write .triq file\n");
      return BADINDEX_ERROR;
    }

    if ((p_strm = fopen(fileName, "w")) == NULL) {
      ERR("error opening file %s\n", fileName);
      return FILE_ERROR;
    }

    fprintf(p_strm, "%12d %12d %12d\n", p_surf->nVerts, p_surf->nTris,
            p_surf->nVertScalars);

    for (i=0; i<p_surf->nVerts; ++i){
      /* verts are single precision */
      fprintf(p_strm, "%16.9e %16.9e %16.9e\n", p_surf->a_Verts[i].x[X],
              p_surf->a_Verts[i].x[Y],
              p_surf->a_Verts[i].x[Z]);
    }

    for (i=0; i<p_surf->nTris; ++i){
      fprintf(p_strm, "%12d %12d %12d\n", p_surf->a_Tris[i].vtx[0]+1,
              p_surf->a_Tris[i].vtx[1]+1,
              p_surf->a_Tris[i].vtx[2]+1);
    }

    for (i=0; i<p_surf->nTris; ++i){
      fprintf(p_strm, "%12d\n", (int)p_surf->a_Tris[i].Comp + 1);
    }

    for (i=0; i<p_surf->nVerts; ++i){
      fprintf(p_strm, "%22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
              p_surf->a_scalar0[i     + cp_offset*p_surf->nVerts],
              p_surf->a_scalar0[i     + rho_offset*p_surf->nVerts],
              p_surf->a_scalar0[3*i   + v_offset*p_surf->nVerts],
              p_surf->a_scalar0[3*i+1 + v_offset*p_surf->nVerts],
              p_surf->a_scalar0[3*i+2 + v_offset*p_surf->nVerts],
              p_surf->a_scalar0[i     + p_offset*p_surf->nVerts]);
    }
  }

  fclose(p_strm);

  return 0;
} /* end io_writeSurfTri */

/**
 *  ----------------- io_readFineMesh () -----...alloc's and reads in the
 *                                              cells of the finest
 *                                              mesh. From a Mesh.*.c3d
 *                                              file, There may be more
 *                                              than 1 level of mesh in
 *                                              these files, but the rest
 *                                              is ignored ----------------
 */

int io_readFineMesh(char* InputMeshFile, const p_tsBasicGrid p_g) {

  int        i, j, nTinyHexes, totalTri=0, nBoundCells, itri;
  int        nFacesXYZ[DIM];
  int       *p_tpEntry, *p_tpIntTriList;
  p_dpoint3  p_tpCentroids;
  double    *p_tpAreas;
  FILE      *p_Inputstrm;
  double thisArea, thisCentroid[DIM];

  REPORTSTATUS("Opening", InputMeshFile);
  if (NULL ==( p_Inputstrm = fopen(InputMeshFile,"r")) ){
    ERR( "Could not find input Mesh file '%s'\n", InputMeshFile);
    exit(FILE_ERROR);
  }
  /*                           read directly into this grids structures */

  REPORTSTATUS("Header", InputMeshFile);          /* 1. number of items */
  /* (void) fseek(p_Inputstrm,0L,SEEK_SET);               */
  fread(&p_g->nVolHexes,  sizeof(int),  1, p_Inputstrm);
  fread(&p_g->nCutHexes,  sizeof(int),  1, p_Inputstrm);
  fread(&p_g->nSplitCells,sizeof(int),  1, p_Inputstrm);
  fread(&p_g->nFaces, sizeof(int),  1, p_Inputstrm);
  fread(nFacesXYZ,   sizeof(int),DIM, p_Inputstrm);
  fread(&p_g->nCutFaces,  sizeof(int),  1, p_Inputstrm);

  REPORTSTATUS("HexInformation", InputMeshFile);    /* 2. tiny hex info */
  nTinyHexes = p_g->nVolHexes + p_g->nCutHexes;
  NEW_ARRAY(p_g->a_Cells, tsTinyHex, nTinyHexes);
  for(j=0;j<nTinyHexes;j++){
    fread(&p_g->a_Cells[j].name,   sizeof(INT64),  1, p_Inputstrm);
    fread(&p_g->a_Cells[j].ref[0], sizeof(char), DIM, p_Inputstrm);
    fread(&p_g->a_Cells[j].flagByte, sizeof(byte), 1, p_Inputstrm);
  }

  /*                                         3. cut and split cell info */
  nBoundCells = p_g->nCutHexes + p_g->nSplitCells;
  p_g->nCells = p_g->nCutHexes + p_g->nSplitCells + p_g->nVolHexes;

  p_g->a_cCells = NULL;
  p_tpEntry     = NULL;
  if ( nBoundCells > 0 ) {
    NEW_ARRAY(p_g->a_cCells, tsCutCell,nBoundCells);
    NEW_ARRAY(p_tpEntry    , int      ,nBoundCells);
  }


  REPORTSTATUS("CutCellInfo", InputMeshFile);
  for(j=0;j<nBoundCells;j++){
    fread(&p_tpEntry[j],               sizeof(int     ),1,p_Inputstrm);
    fread(&p_g->a_cCells[j].nIntTri,   sizeof(int     ),1,p_Inputstrm);
    fread( p_g->a_cCells[j].normal,    sizeof(dpoint3 ),1,p_Inputstrm);
    fread( p_g->a_cCells[j].centroid,  sizeof(dpoint3 ),1,p_Inputstrm);
    fread(&p_g->a_cCells[j].volume,    sizeof(double  ),1,p_Inputstrm);
    fread(&p_g->a_cCells[j].splitIndex,sizeof(int     ),1,p_Inputstrm);
  }
  /*                             ...conventions in the Mesh file ensure */
  /*                                    tPoly lists get built like this */

  if (0 < nBoundCells){ /* totalTri initialized to 0 above */
    totalTri = (p_g->tPolys.nTotTriPolys = p_tpEntry[nBoundCells - 1] +
                p_g->a_cCells[nBoundCells - 1].nIntTri);
  }

  p_tpIntTriList = NULL;
  p_tpCentroids  = NULL;
  p_tpAreas      = NULL;
  if ( totalTri > 0 ) {
    NEW_ARRAY(p_tpIntTriList,   int     , totalTri);/* ...Get space for ALL */
    NEW_ARRAY(p_tpCentroids,    dpoint3 , totalTri);/*    surf tPoly's info */
    NEW_ARRAY(p_tpAreas,        double  , totalTri);/*      ("tpoly lists") */
  }

  p_g->tPolys.p_tpEntry      = p_tpEntry;     /*   ...hook poly lists */
  p_g->tPolys.p_tpIntTriList = p_tpIntTriList;/*        into the grid */
  p_g->tPolys.p_tpCentroids  = p_tpCentroids;
  p_g->tPolys.p_tpAreas      = p_tpAreas;

  REPORTSTATUS("tPolyInfo", InputMeshFile);
  for(j=0;j<totalTri;j++){         /* ...read all the lists of tPoly info */
    fread(&p_tpIntTriList[j], sizeof(int)     , 1, p_Inputstrm);
    fread( p_tpCentroids[j] , sizeof(dpoint3) , 1, p_Inputstrm);
    fread(&p_tpAreas[j]     , sizeof(double)  , 1, p_Inputstrm);
  }

  for(j=0;j<nBoundCells;j++){/*...set pointers to hook info in tPoly lists */
    p_g->a_cCells[j].p_IntTriList =         &p_tpIntTriList[ p_tpEntry[j] ];
    p_g->a_cCells[j].p_centroid   = (dpoint3*)p_tpCentroids[ p_tpEntry[j] ];
    p_g->a_cCells[j].p_area       =              &p_tpAreas[ p_tpEntry[j] ];
    thisArea = thisCentroid[X] = thisCentroid[Y]= thisCentroid[Z] = 0.0;
    for(itri=0;itri<p_g->a_cCells[j].nIntTri;itri++){
      thisArea += p_g->a_cCells[j].p_area[itri];        /* ...accumulation */
      for(i=0;i<DIM;i++) {
        thisCentroid[i] += p_g->a_cCells[j].p_area[itri]
          * p_g->a_cCells[j].p_centroid[itri][i];
      }
    }
    for(i=0;i<DIM;i++){                    /*...save the surf Centroid  */
      thisCentroid[i] /= thisArea;
      p_g->a_cCells[j].surfCentroid[i] = thisCentroid[i];
    }
  }

  p_g->a_Faces = NULL;
  if ( p_g->nFaces > 0 ) {
    NEW_ARRAY(p_g->a_Faces, tsFace, p_g->nFaces);
  }
  REPORTSTATUS("CartFaceLists", InputMeshFile);
  {
    const p_tsFace p_Faces = p_g->a_Faces;
    for(j=0;j<p_g->nFaces;j++){
      fread(&p_Faces[j].adjCell[0], sizeof(int), 1, p_Inputstrm);
      fread(&p_Faces[j].adjCell[1], sizeof(int), 1, p_Inputstrm);
      fread(&p_Faces[j].dir,          sizeof(byte), 1, p_Inputstrm);
    }
  }

  p_g->a_cFaces = NULL;
  if ( p_g->nCutFaces > 0 ) {
    NEW_ARRAY(p_g->a_cFaces, tsCutFace, p_g->nCutFaces);
  }
  REPORTSTATUS("CutFaces", InputMeshFile);
  {
    const p_tsCutFace p_cFaces = p_g->a_cFaces;
    for(j=0;j<p_g->nCutFaces;j++){
      fread(&p_cFaces[j].adjCell[0], sizeof(int),    1,p_Inputstrm);
      fread(&p_cFaces[j].adjCell[1], sizeof(int),    1,p_Inputstrm);
      fread( p_cFaces[j].centroid  , sizeof(dpoint3),1,p_Inputstrm);
      fread(&p_cFaces[j].area,       sizeof(double), 1,p_Inputstrm);
      fread(&p_cFaces[j].dir,        sizeof(char),   1,p_Inputstrm);
      if (IS_NAN(p_cFaces[j].area)       || IS_NAN(p_cFaces[j].centroid[X]) ||
          IS_NAN(p_cFaces[j].centroid[Y])|| IS_NAN(p_cFaces[j].centroid[Z]) ){
        WARN("nan found on incomming face list, cut Face # %d\n",j);
      }
    }
  }
  REPORTSTATUS("                                           ",InputMeshFile);

  fclose(p_Inputstrm);

  return (totalTri);
}

/**
 *  ------ io_setVerboseIO() ---------...set the c3d_io_lib verbose mode T/F
 */
void io_setVerboseIO(bool desiredVerbosity){
  verbose_io = desiredVerbosity;
  return;
}

/**
 *  ------ ioPrivate_copyHinfo()----- ....safe (deep) copy of convergence
 *                                        history info for restarts ---------
 */
void ioPrivate_copyHinfo(const p_tsHinfo p_orig, const p_tsHinfo p_dest) {
  int i;
  p_dest->initResidual   = p_orig->initResidual;
  p_dest->targetResidual = p_orig->targetResidual;
  p_dest->elapsedCPUtime = p_orig->elapsedCPUtime;
  p_dest->iHist          = p_orig->iHist;
  p_dest->iForce         = p_orig->iForce;
  p_dest->refArea        = p_orig->refArea;
  p_dest->refLength      = p_orig->refLength;
  p_dest->xStop          = p_orig->xStop;
  p_dest->xStart         = p_orig->xStart;

  for (i=0;i<3;i++) {
    p_dest->timeArray[i] = p_orig->timeArray[i];
    p_dest->netForce[i]  = p_orig->netForce[i];
    p_dest->momentCtr[i] = p_orig->momentCtr[i];
    p_dest->netMoment[i] = p_orig->netMoment[i];
  }
  return;
}

/**
 *  ----------io_readCkpt()----...read the state vector into a BasicGrid
 *                                 from a 'check.#####' restart file ------
 */
void io_readCkPt(char *inputCkPtFile,   p_tsBasicGrid p_g,
                 p_tsMGinfo p_mgInfo,   p_tsHinfo p_Hinfo,
                 const int nOlapFlowHexes, p_tsTDinfo p_TDinfo){
  FILE      *p_strm;
  p_tsState  p_start;
  tsHinfo    dummy_his;

  REPORTSTATUS("OpeningStateVector",inputCkPtFile);
  if (NULL ==( p_strm= fopen(inputCkPtFile,"r")) ){
    ERR( "error opening checkpoint file '%s' \n",inputCkPtFile);
    exit(FILE_ERROR);
    return;
  }
                                     /* 1. get the mgInfo from restart file  */
  fread(p_mgInfo, sizeof(tsMGinfo),1,p_strm);

                                    /* 2. Read  the  convergence  info next  */
  fread(&dummy_his, sizeof(tsHinfo),1,p_strm);
  /* careful with sensor info: dont copy the */
  ioPrivate_copyHinfo(&dummy_his, p_Hinfo);             /* sensor pointer */

  /* remember, the soltn is on the */
  /*           "current" best grid */

  REPORTSTATUS("ReadingStateVector",inputCkPtFile);

  /*                                3. Read stagevec values into a_U      */
  p_start = p_g->a_U;            /*    3.a Read through  nVolHexes        */
  fread(p_start, sizeof(tsState), p_g->nVolHexes, p_strm);

  /*                                   3.b Read into nCutHexes            */
  p_start = &p_g->a_U[p_g->nVolHexes+nOlapFlowHexes]; /* skip olap   */
  fread(p_start, sizeof(tsState), p_g->nCutHexes, p_strm);

  /*                                   3.c Read into nSplitCells          */
  p_start = p_start + (size_t)p_g->nCutHexes;
  fread(p_start, sizeof(tsState), p_g->nSplitCells, p_strm);


  /*                                4. If Unsteady read the TDinfo struct */
  if (NULL != p_TDinfo) {          /*  4.a Store the TDinfo structs       */
    fread(p_TDinfo, sizeof(tsTDinfo),1,p_strm);
    /*                              4.b. Read Unm1 into "other" state a_Uo */
    p_start = p_g->a_Uo;           /*   a. Read through  nVolHexes        */
    fread(p_start, sizeof(tsState), p_g->nVolHexes, p_strm);

    /*                                  b. Read into nCutHexes            */
    p_start = &p_g->a_Uo[p_g->nVolHexes+nOlapFlowHexes]; /*   skip olap   */
    fread(p_start, sizeof(tsState), p_g->nCutHexes, p_strm);

    /*                                  c. Read into nSplitCells          */
    p_start = p_start + (size_t)p_g->nCutHexes;
    fread(p_start, sizeof(tsState), p_g->nSplitCells, p_strm);
  }


  REPORTSTATUS("                                              ",inputCkPtFile);
  if (verbose_io) printf("\n");

  (void)fclose(p_strm);
  return;
}
/* --------------------------------------------------------------------------
 * NOTES:
 *  o remember the index space of a_U is composed of:
 *              [nVolHexes, nCutHexes, nOlapFlowHexes nSplitCells]
 *              !! p_g->a_U MUST BE EXTERNALLY ALLOC'ED !!
 *     this is because we dont know how much olap padding it may
 *     require at the end.
 *  o nOlapFlowHexes is an input param which you can use to leave
 *     a gap between the CutHexes and the SplitCells
 * -------------------------------------------------------------------------- 
 */

/**
 *  ---io_writeMesh() ---------------------...  write the Mesh to file .-----
 */

#define NOT_SORTED -111
#define ONE_MB     1048576

int io_writeMesh(char* OutputMeshFile, const p_tsBasicGrid p_g){

  int         j, i, nBoundCells, nTotalTri=0;
  size_t      nBytes = 0;
  FILE       *p_strm;
  int         dummyXYZ[DIM];
  p_tsCutCell p_cC;

  NOTE("Writing mesh: ..........\r"); fflush(stdout);

  p_strm = fopen(OutputMeshFile,"wb");       /*   ..open the OutputMeshFile */
  if (NULL == p_strm) {
    ERR("FILE_ERROR writing %s\n",OutputMeshFile);
    exit(FILE_ERROR);
    return(1);
  }

  dummyXYZ[X] = NOT_SORTED; dummyXYZ[Y] = NOT_SORTED; dummyXYZ[Z] = NOT_SORTED;

  OUTREPORTSTATUS("Header",OutputMeshFile);  /* 1.  write HEADER Information */
  fwrite(&p_g->nVolHexes,  sizeof(int),  1, p_strm);
  fwrite(&p_g->nCutHexes,  sizeof(int),  1, p_strm);
  fwrite(&p_g->nSplitCells,sizeof(int),  1, p_strm);
  fwrite(&p_g->nFaces,     sizeof(int),  1, p_strm);
  fwrite( dummyXYZ,        sizeof(int),DIM, p_strm);
  fwrite(&p_g->nCutFaces,  sizeof(int),  1, p_strm);
  nBytes += 8*sizeof(int);

  /* 3. Write CELL info           */
  OUTREPORTSTATUS("HexInformation",OutputMeshFile);/*HexInfo of all Volume */
  for(j=0;j<p_g->nVolHexes + p_g->nCutHexes ;j++){/*   hexes  and cutHexes */
    fwrite(&p_g->a_Cells[j].name,    sizeof(INT64),   1,p_strm);
    fwrite(&p_g->a_Cells[j].ref[0],  sizeof(char),  DIM,p_strm);
    fwrite(&p_g->a_Cells[j].flagByte,sizeof(char),    1,p_strm);
  }
  nBytes += (p_g->nVolHexes + p_g->nCutHexes)*(sizeof(INT64) + 4*sizeof(char));

  nBoundCells = p_g->nCutHexes + p_g->nSplitCells;
  OUTREPORTSTATUS("CutCellInfo",OutputMeshFile);
  for(j=0;j<nBoundCells;j++){              /*..Write cut cell info (slow!) */
    fwrite(&nTotalTri,                   sizeof(int)    , 1, p_strm);
    fwrite(&p_g->a_cCells[j].nIntTri,    sizeof(int)    , 1, p_strm);
    fwrite( p_g->a_cCells[j].normal,     sizeof(dpoint3), 1, p_strm);
    fwrite( p_g->a_cCells[j].centroid,   sizeof(dpoint3), 1, p_strm);
    fwrite(&p_g->a_cCells[j].volume,     sizeof(double ), 1, p_strm);
    fwrite(&p_g->a_cCells[j].splitIndex, sizeof(int)    , 1, p_strm);
    nTotalTri+=p_g->a_cCells[j].nIntTri;
  }
  nBytes+= nBoundCells*(2*sizeof(dpoint3) + 3*sizeof(int) + sizeof(double));

  OUTREPORTSTATUS("tPolyInfo",OutputMeshFile);
  for(j=0;j<nBoundCells;j++){  /* ...write triangle list for each cut Cell */
    p_cC = &p_g->a_cCells[j];
    for(i=0;i<p_cC->nIntTri;i++){
      fwrite(&p_cC->p_IntTriList[i], sizeof(int)     , 1, p_strm);
      fwrite( p_cC->p_centroid[i]  , sizeof(dpoint3) , 1, p_strm);
      fwrite(&p_cC->p_area[i]      , sizeof(double)  , 1, p_strm);
    }
    nBytes+= p_cC->nIntTri*(sizeof(int) + sizeof(dpoint3) + sizeof(double));
  }

  OUTREPORTSTATUS("CartFaceLists",OutputMeshFile);/* 4.Write FACE information*/
  for(j=0;j<p_g->nFaces;j++){
    fwrite(&p_g->a_Faces[j].adjCell[0],   sizeof(int), 1,p_strm);
    fwrite(&p_g->a_Faces[j].adjCell[1],   sizeof(int), 1,p_strm);
    fwrite(&p_g->a_Faces[j].dir,          sizeof(byte),1,p_strm);
  }
  nBytes+= p_g->nFaces*(2*sizeof(int) + sizeof(byte));

  OUTREPORTSTATUS("cutFaces",OutputMeshFile); /* write tsCutFaces elm-by-elm */
  for(j=0;j<p_g->nCutFaces;j++){
    fwrite( p_g->a_cFaces[j].adjCell, sizeof(int    ), 2, p_strm);
    fwrite( p_g->a_cFaces[j].centroid,sizeof(dpoint3), 1, p_strm);
    fwrite(&p_g->a_cFaces[j].area    ,sizeof(double ), 1, p_strm);
    fwrite(&p_g->a_cFaces[j].dir     ,sizeof(char   ), 1, p_strm);
  }
  nBytes+= p_g->nCutFaces *
    (sizeof(int) + sizeof(dpoint3) + sizeof(double) + sizeof(char));

  OUTREPORTSTATUS("....written",OutputMeshFile);
  if (verbose_io) printf("\n");

  if (p_strm) fflush(p_strm);

  /* force to disk: experimental to see if it reduces intermitent failures in
   * downstream executables
   */
  if (p_strm) {
    int fd, rc;
    fd = fileno(p_strm);
    rc = fsync(fd);
    if ( 0 != rc ) {
      WARN("fsync on mesh io failed\n");
    }
  }

  if (p_strm) (void)fclose(p_strm);

  printf("              ...Wrote %3.2fMb file %s with mesh\n",
         ((float)nBytes)/ONE_MB, OutputMeshFile);

  fflush(stdout);

  return(0);
}

#undef NOT_SORTED
#undef ONE_MB

/* --------------------------------------------------------------------------
 * NOTES:
 * 1. dummyXYZ[DIM]: this is a relic. We dont keep the faces sorted in xyz
 *    anymore, so just dump a flag in the file and move on
 *--------------------------------------------------------------------------
 */

/* ----------io_writeCkpt()----...read the state vector into a BasicGrid
 *                                from a 'check.#####' restart file -------
 */
int io_writeCkPt(char *outputCkPtFile,  p_tsBasicGrid p_g,
                 p_tsMGinfo p_mgInfo,    p_tsHinfo p_Hinfo,
                 p_tsTDinfo p_TDinfo){
  FILE      *p_strm;
  p_tsState  p_start;

  OUTREPORTSTATUS("OpeningStateVector",outputCkPtFile);
  if (NULL ==( p_strm= fopen(outputCkPtFile,"wb")) ){
    ERR( "error opening checkpoint file '%s' \n",outputCkPtFile);
    exit(FILE_ERROR);
    return(1);
  }
  /*                               1. Write the mgInfo from restart file  */
  fwrite(p_mgInfo, sizeof(tsMGinfo),1,p_strm);
  /*                               2. Write the  convergence  info next   */
  fwrite(p_Hinfo, sizeof(tsHinfo),1,p_strm);


  OUTREPORTSTATUS("WritingStateVector",outputCkPtFile);

  p_start = p_g->a_U;  /*          3. Mesh is sorted. Write all at once.  */
  fwrite(p_start, sizeof(tsState), p_g->nCells, p_strm);

  /*                               4.  If Unsteady read the TDinfo struct */
  if (NULL != p_TDinfo) {
    fwrite(p_TDinfo, sizeof(tsTDinfo),1,p_strm);
    OUTREPORTSTATUS("Writing_Unm1_StateVector",outputCkPtFile);
    p_start = p_g->a_Uo;                  /* ...a_Uo holds U_n-1 solution */
    fwrite(p_start, sizeof(tsState), p_g->nCells, p_strm);
  }


  OUTREPORTSTATUS("                                          ",outputCkPtFile);
  if (verbose_io) printf("\n");

  fflush(p_strm);

  /* force to disk: experimental to see if it reduces intermittent failures in
   * downstream executables
   */
  {
    int fd, rc;
    fd = fileno(p_strm);
    rc = fsync(fd);
    if ( 0 != rc ) {
      WARN("fsync on mesh io failed\n");
    }
  }

  (void)fclose(p_strm);

  fflush(stdout);
  return(0);
}

/* unused code, may come in handy later */

/**
 *---c3d_isASCII()-------------------------(currently unused)---------
 *                               ...Check if file format is ASCII or
 *                               binary by using the system 'file'
 *                               command and checking for 'ASCII' in its
 *                               output. Return 1 if file is ASCII, zero
 *                               otherwise.
 */
int
c3d_isASCII(const char *const p_name)
{
  int istatus = 0;
  char cTestFile[STRING_LEN], line[STRING_LEN];
  char *p_char=NULL;
  FILE *p_strm;

  strcpy(cTestFile,"file -i ");
  strcat(cTestFile, p_name);
  strcat(cTestFile, "  > __tmp_comp2trix.ascii");

  if ( 0 != system(cTestFile) ) {
    ERR("Unable to execute %s\n",cTestFile);
    exit(FILE_ERROR);
  }

  if ((p_strm = fopen("__tmp_comp2trix.ascii","r")) == NULL) {
    ERR("error opening parameter file __tmp_comp2trix.ascii\n");
    exit(FILE_ERROR);
  }

  line[0]='\0';
  while (fgets(line, STRING_LEN, p_strm) != NULL) {
    p_char = strstr(line,"ASCII");
    if ( p_char ) {
      istatus = 1;
      break;
    }

    p_char = strstr(line,"ascii");
    if ( p_char ) {
      istatus = 1;
      break;
    }

    p_char = strstr(line,"Ascii");
    if ( p_char ) {
      istatus = 1;
      break;
    }
  }

  fclose(p_strm);

  strcpy(cTestFile,"/bin/rm -f __tmp_comp2trix.ascii");
  if ( 0 != system(cTestFile) ) {
    ERR("Unable to execute %s\n",cTestFile);
    exit(FILE_ERROR);
  }

  return istatus;
}
/**
 *---readASCIItriFile() ------------------(currently unused)---------
 */
int readASCIItriFile(const char *const p_name, p_tsTriangulation p_surf,
                     const bool verbose)
{
  FILE *p_strm;
  char line[STRING_LEN], value1[100], value2[100];
  int i, nV, nT;

  const int nVerts            = p_surf->nVerts;
  const int nTris             = p_surf->nTris;
  const p_tsVertex p_V        = p_surf->a_Verts;
  const p_tsTri p_T           = p_surf->a_Tris;

  ASSERT(p_surf);

  if ((p_strm = fopen(p_surf->geomName,"r")) == NULL) {
    ERR("error opening trix file %s\n",p_surf->geomName);
    exit(FILE_ERROR);
  }
  /* ...comment skipper */
  while  (fgets(line, STRING_LEN, p_strm) != NULL) {
    if (line[0]=='#' || line[0]=='\n') /* ...if 1st char is #  or */
      { line[0]='\0'; continue;}       /*     '\n' then skip line */
    break;                    /* assume it is the first data line */
  }

  if (sscanf(line,"%s %s\n",value1,value2) !=2){
    ERR("Error reading tri file header parameters\n");
    exit(PARSE_ERROR);
  }

  if (sscanf(value1,"%d",&nV) !=1) exit(PARSE_ERROR);
  if (sscanf(value2,"%d",&nT) !=1) exit(PARSE_ERROR);

  if ( nV != nVerts || nT != nTris ){
    ERR("nVerts or nTris do not agree\n");
  }

  for (i=0; i<nVerts; ++i){
    fscanf(p_strm,"%e %e %e",&p_V[i].x[X], &p_V[i].x[Y], &p_V[i].x[Z]);
  }

  for (i=0; i<nTris; ++i){
    fscanf(p_strm,"%d %d %d",&p_T[i].vtx[0], &p_T[i].vtx[1], &p_T[i].vtx[2]);
    --p_T[i].vtx[0]; --p_T[i].vtx[1]; --p_T[i].vtx[2]; /* reduce to zero offset */
  }

  fclose(p_strm);
  return 0;
}

int c3d_readASCIItriFileHeader(const char *const p_name, int *const p_nVerts,
                               int *const p_nTris,
                               const bool verbose)
{
  FILE *p_strm;
  char  line[STRING_LEN], value1[100], value2[100];

  if ((p_strm = fopen(p_name,"r")) == NULL) {
    ERR("error opening tri file %s\n",p_name);
    exit(FILE_ERROR);
  }

  while  (fgets(line, STRING_LEN, p_strm) != NULL) {
    if (line[0]=='#' || line[0]=='\n') /* ...if 1st char is #  or */
      { line[0]='\0'; continue;}       /*     '\n' then skip line */
    break;  /* ...assume it's the first data line */
  }

  if (sscanf(line,"%s %s\n",value1,value2) !=2){
    ERR("Error reading tri file header parameters\n");
    exit(PARSE_ERROR);
  }

  if (sscanf(value1,"%d",p_nVerts) !=1) exit(PARSE_ERROR);
  if (sscanf(value2,"%d",p_nTris)  !=1) exit(PARSE_ERROR);

  if (verbose){
    NOTE("nVerts=%d nTris=%d\n", *p_nVerts, *p_nTris);
  }

  fclose(p_strm);
  return 0;
}

/*
 * ----------------------------------------------------------------------
 */

#undef REPORTSTATUS
#undef OUTREPORTSTATUS
#undef BADINDEX_ERROR
