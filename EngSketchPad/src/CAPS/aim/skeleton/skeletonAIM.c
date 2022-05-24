/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Skeleton AIM Example Code
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

/*!\mainpage Introduction
 *
 * \section overviewSkeleton Skeleton AIM Overview
 * This is an example skeleten AIM intended to provude examples of best practices when writing an AIM
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"

#ifdef WIN32
#define strcasecmp stricmp
#endif

#define DEBUG

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

enum aimInputs
{
  InputVariable = 1,           /* index is 1-based */
  SkeletonAIMin,
  Mach,
  Mesh_Format,
  Table,
  NUMINPUT = Table             /* Total number of inputs */
};

enum aimOutputs
{
  SkeletonAIMout = 1,          /* index is 1-based */
  NUMOUT = SkeletonAIMout      /* Total number of outputs */
};

typedef struct {
  int nBody;
  ego *tess;
} aimStorage;

/****************** exposed AIM entry points -- Analysis **********************/

/* aimInitialize: Initialization Information for the AIM */
int
aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
              /*@unused@*/ void **instStore, /*@unused@*/ int *major,
              /*@unused@*/ int *minor, int *nIn, int *nOut,
              int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int status = CAPS_SUCCESS;
  int  *ints=NULL, i;
  char **strs = NULL;
  aimStorage *aimStore = NULL;

#ifdef DEBUG
/*@-nullpass@*/
  printf("\n skeletonAIM/aimInitialize  instance = %d  unitSys = %s!\n",
         inst, unitSys);
/*@+nullpass@*/
#endif

  /* specify the number of analysis inputs  defined in aimInputs
   *     and the number of analysis outputs defined in aimOutputs */
  *nIn    = NUMINPUT;
  *nOut   = NUMOUT;

  /* return if "query" only */
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 2;

  /* specify the name of each field variable */
  AIM_ALLOC(strs, *nFields, char *, aimInfo, status);

  strs[0]  = EG_strdup("scalar");
  strs[1]  = EG_strdup("coordinates");
  for (i = 0; i < *nFields; i++)
    if (strs[i] == NULL) { status = EGADS_MALLOC; goto cleanup; }
  *fnames  = strs;

  /* specify the dimension of each field variable */
  AIM_ALLOC(ints, *nFields, int, aimInfo, status);
  ints[0]  = 1;
  ints[1]  = 3;
  *franks  = ints;
  ints = NULL;

  /* specify if a field is an input field or output field */
  AIM_ALLOC(ints, *nFields, int, aimInfo, status);

  ints[0]  = FieldIn;
  ints[1]  = FieldOut;
  *fInOut  = ints;
  ints = NULL;

  /* setup our AIM specific state */
  AIM_ALLOC(aimStore, 1, aimStorage, aimInfo, status);
  *instStore = aimStore;

  aimStore->nBody = 0;
  aimStore->tess  = NULL;

  status = CAPS_SUCCESS;

cleanup:
  if (status != CAPS_SUCCESS) {
    /* release all possibly allocated memory on error */
    if (*fnames != NULL)
      for (i = 0; i < *nFields; i++) AIM_FREE((*fnames)[i]);
    AIM_FREE(*franks);
    AIM_FREE(*fInOut);
    AIM_FREE(*fnames);
    AIM_FREE(*instStore);
    *nFields = 0;
  }
  
  return status;
}


/* aimInputs: Input Information for the AIM */
int
aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, int index,
          char **ainame, capsValue *defval)
{
  int i, status = CAPS_SUCCESS;
  capsTuple *tuple=NULL;
  
#ifdef DEBUG
  printf(" skeletonAIM/aimInputs  instance = %d  index  = %d!\n",
         aim_getInstance(aimInfo), index);
#endif

  /* fill in the required members based on the index */
 if (index == InputVariable) {
    *ainame              = AIM_NAME(InputVariable);
    defval->type         = Boolean;
    defval->vals.integer = false;

    /*! \page aimInputsSkeleton
     * - <B> InputVariable = false</B> <br>
     * A boolean input variable.
     */

  } else if (index == SkeletonAIMin) {
    *ainame           = AIM_NAME(SkeletonAIMin);
    defval->type      = Double;
    defval->vals.real = 5.0;
    defval->units     = EG_strdup("cm");

    /*! \page aimInputsSkeleton
     * - <B> SkeletonAIMin = 5 cm</B> <br>
     * A real input with units
     */

  } else if (index == Mach) {
    *ainame              = AIM_NAME(Mach); // Mach number
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = NULL;
    defval->lfixed       = Change;
    defval->dim          = Scalar;

    /*! \page aimInputsSkeleton
     * - <B> Mach = NULL</B> <br>
     * Mach number
     */

  } else if (index == Mesh_Format) {
    *ainame              = AIM_NAME(Mesh_Format);
    defval->type         = String;
    defval->lfixed       = Change;
    AIM_STRDUP(defval->vals.string, "AFLR3", aimInfo, status);

    /*! \page aimInputsSkeleton
     * - <B> Mesh_Format = AFLR3 </B> <br>
     * String mesh format
     */

  } else if (index == Table) {

    /* an example of filling in a Tuple */
    *ainame = AIM_NAME(Table);

    AIM_ALLOC( tuple, 3, capsTuple, aimInfo, status);
    for (i = 0; i < 3; i++) {
      tuple[i].name = NULL;
      tuple[i].value = NULL;
    }
    AIM_STRDUP(tuple[0].name , "Entry1", aimInfo, status);
    AIM_STRDUP(tuple[1].name , "Entry2", aimInfo, status);
    AIM_STRDUP(tuple[2].name , "Entry3", aimInfo, status);
    AIM_STRDUP(tuple[0].value, "Value1", aimInfo, status);
    AIM_STRDUP(tuple[1].value, "Value2", aimInfo, status);
    AIM_STRDUP(tuple[2].value, "Value3", aimInfo, status);

    defval->type       = Tuple;
    defval->dim        = Vector;
    defval->nrow       = 1;
    defval->ncol       = 3;
    defval->vals.tuple = tuple;

    /*! \page aimInputsSkeleton
     * - <B> Table = {Entry1:Value1, Entry2:Value2, Entry3:Value3} </B> <br>
     * An example of a tuple input
     */
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  return status;
}


/* aimPreAnalysis: Parse Inputs, Generate Input File(s) & Possibly Tessellate */
int
aimPreAnalysis(void *instStore, void *aimInfo, capsValue *inputs)
{
  int        i, n, status = CAPS_SUCCESS, nBody;
  const char *name, *intents;
  double     size, box[6], params[3], mach;
  ego        *bodies;
  capsValue  *val;
  aimStorage *aimStore;

#ifdef DEBUG
  printf("\n skeletonAIM/aimPreAnalysis instance = %d!\n",
         aim_getInstance(aimInfo));
#endif
  aimStore = (aimStorage *) instStore;

  /* look at the CSM design parameters */
  printf("   GeometryIn:\n");
  n = aim_getIndex(aimInfo, NULL, GEOMETRYIN);
  for (i = 0; i < n; i++) {
    status = aim_getName(aimInfo, i+1, GEOMETRYIN, &name);
    if (status != CAPS_SUCCESS) continue;
    status = aim_getValue(aimInfo, i+1, GEOMETRYIN, &val);
    if (status == CAPS_SUCCESS)
      printf("       %d: %s  %d  (%d,%d)\n", i+1, name,
             val->type, val->nrow, val->ncol);
  }

  /* look at the CSM output parameters */
  printf("\n   GeometryOut:\n");
  n = aim_getIndex(aimInfo, NULL, GEOMETRYOUT);
  for (i = 0; i < n; i++) {
    status = aim_getName(aimInfo, i+1, GEOMETRYOUT, &name);
    if (status != CAPS_SUCCESS) continue;
    status = aim_getValue(aimInfo, i+1, GEOMETRYOUT, &val);
    if (status == CAPS_SUCCESS)
      printf("       %d: %s  %d  (%d,%d)\n", i+1, name,
             val->type, val->nrow, val->ncol);
  }

  /* write out input list of values */
  if (inputs != NULL) {
    printf("\n   AnalysisIn:\n");
    for (i = 0; i < NUMINPUT; i++) {
      status = aim_getName(aimInfo, i+1, ANALYSISIN, &name);
      AIM_STATUS(aimInfo, status);
      printf("       %d: %s  %d  (%d,%d) %s\n", i+1, name,
             inputs[i].type, inputs[i].nrow, inputs[i].ncol, inputs[i].units);
    }
    mach = inputs[Mach-1].vals.real;
    if (mach < 0) {
      AIM_ANALYSISIN_ERROR(aimInfo, Mach, "Mach number must be >= 0\n");
      AIM_ADDLINE(aimInfo, "Negative mach = %f is non-physical\n", mach);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  }

  /* create the tessellation */
  
  status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
  AIM_STATUS(aimInfo, status, "aim_getBodies");
  if ((bodies == NULL) || (nBody == 0)) return CAPS_SUCCESS;

  aimStore->nBody = nBody;
  AIM_REALL(aimStore->tess, aimStore->nBody, ego, aimInfo, status);

  for (i = 0; i < nBody; i++) {
    /* tessellate with EGADS tessellation in this example */
    status = EG_getBoundingBox(bodies[i], box);
    AIM_STATUS(aimInfo, status);

    size = box[3]-box[0];
    if (size < box[4]-box[1]) size = box[4]-box[1];
    if (size < box[5]-box[2]) size = box[5]-box[2];
    params[0] =  0.025*size;
    params[1] =  0.001*size;
    params[2] = 15.0;
    status = EG_makeTessBody(bodies[i], params, &aimStore->tess[i]);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" skeletonAIM/aimPreAnalysis tess %d!\n", i);
#endif
    
    /* store the tessellation in CAPS */
    status = aim_newTess(aimInfo, aimStore->tess[i]);
    AIM_STATUS(aimInfo, status);
  }
  printf("\n");
  
cleanup:
  return status;
}


/* aimExecute: runs the Analysis & specifies the AIM does the execution */
int
aimExecute(/*@unused@*/ void *instStor, /*@unused@*/ void *aimInfo, int *state)
{
  
  *state = 0;
  return CAPS_SUCCESS;
}


/* aimPostAnalysis: Perform any processing after the Analysis is run */
int
aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS;

  return status;
}


/* aimOutputs: Output Information for the AIM */
int
aimOutputs(/*@unused@*/ void *instStore, void *aimInfo, /*@unused@*/ int index,
           char **aoname, capsValue *form)
{
  int status = CAPS_SUCCESS;
#ifdef DEBUG
  printf(" skeletonAIM/aimOutputs instance = %d  index  = %d!\n",
         aim_getInstance(aimInfo), index);
#endif
  
  if (index == SkeletonAIMout) {
    *aoname = AIM_NAME(SkeletonAIMout);
    form->type = Double;
  }

  return status;
}


/* aimCalcOutput: Calculate/Retrieve Output Information */
int
aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              /*@unused@*/ int index, capsValue *val)
{
#ifdef DEBUG
  int        status = CAPS_SUCCESS;
  const char *name;

  status = aim_getName(aimInfo, index, ANALYSISOUT, &name);
  printf(" skeletonAIM/aimCalcOutput instance = %d  index = %d %s %d!\n",
         aim_getInstance(aimInfo), index, name, status);
#endif

  if (index == SkeletonAIMout) {
    val->vals.real = 3.1415926;
  }

  return status;
}


/* aimCleanup: Free up the AIMs storage */
void aimCleanup(/*@unused@*/ void *instStore)
{
  aimStorage *aimStore=NULL;
#ifdef DEBUG
  printf(" skeletonAIM/aimCleanup!\n");
#endif

  /* clean up any allocated data 
   * tessellation objects are deleted by CAPS 
   */

  aimStore = (aimStorage*)instStore;
  aimStore->nBody = 0;
  AIM_FREE(aimStore->tess);
  AIM_FREE(aimStore);
}


/************ exposed AIM entry points -- Discretization Structure ************/

/* aimFreeDiscrPtr: Frees up the pointer in Discrete Structure -- Optional */
void
aimFreeDiscrPtr(void *ptr)
{
  AIM_FREE(ptr);
}


/* aimDiscr: Fill-in the Discrete data for a Bound Object -- Optional */
int
aimDiscr(char *tname, capsDiscr *discr)
{
  int           ibody, iface, i, vID=0, itri, ielem, status, found;
  int           nFace, atype, alen, tlen, nBodyDisc;
  int           ntris, state, nGlobal, global;
  int           *vid = NULL;
  const int     *ints, *ptype, *pindex, *tris, *nei;
  const double  *reals, *xyz, *uv;
  const char    *string;
  ego           body, *faces=NULL, *tess = NULL;
  capsBodyDiscr *discBody;
  aimStorage    *aimStore=NULL;

#ifdef DEBUG
  printf(" skeletonAIM/aimDiscr: tname = %s, instance = %d!\n",
         tname, aim_getInstance(discr->aInfo));
#endif
  aimStore = (aimStorage*)discr->instStore;
  tess = aimStore->tess;

  /* find any faces with our boundary marker */
  for (nBodyDisc = ibody = 0; ibody < aimStore->nBody; ibody++) {
    status = EG_statusTessBody(tess[ibody], &body, &state, &nGlobal);
    AIM_STATUS(discr->aInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &faces);
    AIM_STATUS(discr->aInfo, status, "getBodyTopos (Face) for Body %d",
               ibody+1);
    AIM_NOTNULL(faces, discr->aInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;
#ifdef DEBUG
      printf(" skeletonAIM/aimDiscr: Body %d/Face %d matches %s!\n",
             ibody+1, iface+1, tname);
#endif
      /* count the number of Bodys with capsBound */
      nBodyDisc++;
      break;
    }
    AIM_FREE(faces);
  }
  if (nBodyDisc == 0) {
    printf(" skeletonAIM/aimDiscr: No Faces match %s!\n", tname);
    return CAPS_SUCCESS;
  }

  /* specify our single triangle element type */
  discr->nTypes = 1;
  AIM_ALLOC(discr->types, discr->nTypes, capsEleType, discr->aInfo, status);

  discr->types[0].nref  = 3;
  discr->types[0].ndata = 0;         /* data at geom reference positions
                                        (i.e. vertex centered/iso-parametric) */
  discr->types[0].ntri  = 1;
  discr->types[0].nseg  = 0;
  discr->types[0].nmat  = 0;         /* match points at geom ref positions */
  discr->types[0].tris  = NULL;
  discr->types[0].segs  = NULL;
  discr->types[0].gst   = NULL;
  discr->types[0].dst   = NULL;
  discr->types[0].matst = NULL;

  /* specify the numbering for the points on the triangle */
  AIM_ALLOC(discr->types[0].tris, discr->types[0].nref, int, discr->aInfo,
            status);

  discr->types[0].tris[0] = 1;
  discr->types[0].tris[1] = 2;
  discr->types[0].tris[2] = 3;

  /* specify the reference coordinates for each point on the triangle */
  AIM_ALLOC(discr->types[0].gst, 2*discr->types[0].nref, double, discr->aInfo,
            status);

  discr->types[0].gst[0] = 0.0;   /* s = 0, t = 0 */
  discr->types[0].gst[1] = 0.0;
  discr->types[0].gst[2] = 1.0;   /* s = 1, t = 0 */
  discr->types[0].gst[3] = 0.0;
  discr->types[0].gst[4] = 0.0;   /* s = 0, t = 1 */
  discr->types[0].gst[5] = 1.0;

  /* allocate the body discretizations */
  discr->nBodys = nBodyDisc;
  AIM_ALLOC(discr->bodys, discr->nBodys, capsBodyDiscr, discr->aInfo, status);

  /* get the tessellation and
   make up a linear continuous triangle discretization */
  vID = nBodyDisc = 0;
  for (ibody = 0; ibody < aimStore->nBody; ibody++) {

    status = EG_statusTessBody(tess[ibody], &body, &state, &nGlobal);
    AIM_STATUS(discr->aInfo, status);

    ntris = 0;
    AIM_FREE(faces);
    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &faces);
    AIM_STATUS(discr->aInfo, status, "getBodyTopos (Face) for Body %d", ibody+1);
    AIM_NOTNULL(faces, discr->aInfo, status);

    found = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;

      status = EG_getTessFace(tess[ibody], iface+1, &alen, &xyz, &uv,
                              &ptype, &pindex, &tlen, &tris, &nei);
      if (status != EGADS_SUCCESS) {
        printf(" skeletonAIM: EG_getTessFace %d = %d for Body %d!\n",
               iface+1, status, ibody+1);
        continue;
      }
      ntris += tlen;
      found = 1;
    }
    if (found == 0) continue;
#ifdef DEBUG
    printf(" skeletonAIM/aimDiscr: ntris = %d!\n", ntris);
#endif
    if (ntris == 0) {
      AIM_ERROR(discr->aInfo, "No faces with capsBound = %s", tname);
      status = CAPS_SOURCEERR;
      goto cleanup;
    }

    discBody = &discr->bodys[nBodyDisc];
    aim_initBodyDiscr(discBody);

    discBody->nElems = ntris;
    discBody->tess   = tess[ibody];

    AIM_ALLOC(discBody->elems,      ntris, capsElement, discr->aInfo, status);
    AIM_ALLOC(discBody->gIndices, 6*ntris, int,         discr->aInfo, status);

    AIM_FREE(vid);
    AIM_ALLOC(vid, nGlobal, int, discr->aInfo, status);
    for (i = 0; i < nGlobal; i++) vid[i] = 0;

    ielem = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;

      status = EG_getTessFace(tess[ibody], iface+1, &alen, &xyz, &uv,
                              &ptype, &pindex, &tlen, &tris, &nei);
      AIM_STATUS(discr->aInfo, status);

      /* construct global vertex indices */
      for (i = 0; i < alen; i++) {
        status = EG_localToGlobal(tess[ibody], iface+1, i+1, &global);
        AIM_STATUS(discr->aInfo, status);
        if (vid[global-1] != 0) continue;
        vid[global-1] = vID+1;
        vID++;
      }

      /* fill elements */
      for (itri = 0; itri < tlen; itri++, ielem++) {
        discBody->elems[ielem].tIndex      = 1;
        discBody->elems[ielem].eIndex      = iface+1;
/*@-immediatetrans@*/
        discBody->elems[ielem].gIndices    = &discBody->gIndices[6*ielem];
/*@+immediatetrans@*/
        discBody->elems[ielem].dIndices    = NULL;
        discBody->elems[ielem].eTris.tq[0] = itri+1;

        status = EG_localToGlobal(tess[ibody], iface+1, tris[3*itri  ], &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[0] = vid[global-1];
        discBody->elems[ielem].gIndices[1] = tris[3*itri  ];

        status = EG_localToGlobal(tess[ibody], iface+1, tris[3*itri+1], &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[2] = vid[global-1];
        discBody->elems[ielem].gIndices[3] = tris[3*itri+1];

        status = EG_localToGlobal(tess[ibody], iface+1, tris[3*itri+2], &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[4] = vid[global-1];
        discBody->elems[ielem].gIndices[5] = tris[3*itri+2];
      }
    }
    nBodyDisc++;
  }

  /* set the total number of points */
  discr->nPoints = vID;

  status = CAPS_SUCCESS;

cleanup:

  /* free up our stuff */
  AIM_FREE(faces);
  AIM_FREE(vid);

  return status;
}


/* aimTransfer: Data Transfer using the Discrete Structure -- Optional */
int
aimTransfer(capsDiscr *discr, const char *fname, int npts, int rank,
            double *data, /*@unused@*/ char **units)
{
  int        status = CAPS_SUCCESS;
  int        i, j, bIndex, ptype, pindex, global;
  double     xyz[3];
  capsBodyDiscr *discBody;

#ifdef DEBUG
  printf(" skeletonAIM/aimTransfer name = %s  instance = %d  npts = %d/%d!\n",
         fname, aim_getInstance(discr->aInfo), npts, rank);
#endif

  if (strcmp(fname, "scalar") == 0 ) {

    /* fill in a simple scalar field */
    for (i = 0; i < npts; i++) {
      data[i] = i;
    }

  } else if (strcmp(fname, "coordinates") == 0 ) {

    /* fill in with our coordinates -- as an example */
    for (i = 0; i < npts; i++) {
      bIndex = discr->tessGlobal[2*i  ];
      global = discr->tessGlobal[2*i+1];

      discBody = &discr->bodys[bIndex-1];
      status = EG_getGlobal(discBody->tess, global, &ptype, &pindex, xyz);
      AIM_STATUS(discr->aInfo, status);
      for (j = 0; j < rank; j++) data[rank*i+j] = xyz[j];
    }
  }

cleanup:
  return status;
}


/********* exposed AIM entry points -- Interpolation Functions ***************/

/* aimLocateElement: Element in the Mesh -- Optional */
int
aimLocateElement(capsDiscr *discr, double *params, double *param,
                 int *bIndex, int *eIndex, double *bary)
{
  int    i, ib, in[3], stat, ibsmall = 0, ismall = 0;
  double we[3], w, smallw = -1.e300;
  capsBodyDiscr *discBody=NULL;

  if (discr == NULL) return CAPS_NULLOBJ;

  /* locate the element that contains param and return the bary centric
     coordinates of that element */

  for (ib = 0; ib < discr->nBodys; ib++) {
    discBody = &discr->bodys[ib];
    for (i = 0; i < discBody->nElems; i++) {
      in[0] = discBody->elems[i].gIndices[0] - 1;
      in[1] = discBody->elems[i].gIndices[2] - 1;
      in[2] = discBody->elems[i].gIndices[4] - 1;
      stat  = EG_inTriExact(&params[2*in[0]], &params[2*in[1]], &params[2*in[2]],
                            param, we);
      if (stat == EGADS_SUCCESS) {
        *bIndex = ib+1;
        *eIndex = i+1;
        bary[0] = we[1];
        bary[1] = we[2];
        return CAPS_SUCCESS;
      }
      w = we[0];
      if (we[1] < w) w = we[1];
      if (we[2] < w) w = we[2];
      if (w > smallw) {
        ibsmall = ib+1;
        ismall = i+1;
        smallw = w;
      }
    }
  }

  /* must extrapolate! */
  if (ismall == 0 || discBody == NULL) return CAPS_NOTFOUND;
  in[0] = discBody->elems[ismall-1].gIndices[0] - 1;
  in[1] = discBody->elems[ismall-1].gIndices[2] - 1;
  in[2] = discBody->elems[ismall-1].gIndices[4] - 1;
  EG_inTriExact(&params[2*in[0]], &params[2*in[1]], &params[2*in[2]],
                param, we);
  *bIndex = ibsmall;
  *eIndex = ismall;
  bary[0] = we[1];
  bary[1] = we[2];

  return CAPS_SUCCESS;
}


/* aimInterpolation: Interpolation on the Bound -- Optional */
int
aimInterpolation(capsDiscr *discr, /*@unused@*/ const char *name,
                 int bIndex, int eIndex, double *bary, int rank,
                 double *data, double *result)
{
  int    in[3], i;
  double we[3];
  capsBodyDiscr *discBody=NULL;

  /* interpolate data to barycentric coordinates of the element */

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" skeletonAIM/Interpolation: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_RANGEERR;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" skeletonAIM/Interpolation: eIndex = %d [1-%d]!\n",
           eIndex, discBody->nElems);
    return CAPS_RANGEERR;
  }

  we[0] = 1.0 - bary[0] - bary[1];
  we[1] = bary[0];
  we[2] = bary[1];
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  for (i = 0; i < rank; i++)
    result[i] = data[rank*in[0]+i]*we[0] + data[rank*in[1]+i]*we[1] +
                data[rank*in[2]+i]*we[2];

  return CAPS_SUCCESS;
}


/* aimInterpolateBar: Interpolation on the Bound -- Optional */
int
aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name,
                  int bIndex, int eIndex, double *bary, int rank,
                  double *r_bar, double *d_bar)
{
  int    in[3], i;
  double we[3];
  capsBodyDiscr *discBody=NULL;

  /* reverse differentiation of aimInterpolation */

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" skeletonAIM/InterpolateBar: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_RANGEERR;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" skeletonAIM/InterpolateBar: eIndex = %d [1-%d]!\n",
           eIndex, discBody->nElems);
    return CAPS_RANGEERR;
  }

  we[0] = 1.0 - bary[0] - bary[1];
  we[1] = bary[0];
  we[2] = bary[1];
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  for (i = 0; i < rank; i++) {
/*  result[i] = data[rank*in[0]+i]*we[0] + data[rank*in[1]+i]*we[1] +
                data[rank*in[2]+i]*we[2];  */
    d_bar[rank*in[0]+i] += we[0]*r_bar[i];
    d_bar[rank*in[1]+i] += we[1]*r_bar[i];
    d_bar[rank*in[2]+i] += we[2]*r_bar[i];
  }

  return CAPS_SUCCESS;
}


/* aimIntegration: Element Integration on the Bound -- Optional */
int
aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name,
               int bIndex, int eIndex, int rank, double *data,
               double *result)
{
  int        status = CAPS_SUCCESS;
  int        i, in[3], ptype, pindex, global[3];
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area;
  capsBodyDiscr *discBody=NULL;

  /* computes the integral of the data on a given element */

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" skeletonAIM/Integration: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_RANGEERR;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems)) {
    printf(" skeletonAIM/Integration: eIndex = %d [1-%d]!\n",
           eIndex, discBody->nElems);
    return CAPS_RANGEERR;
  }

  /* element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++)
    global[i] = discr->tessGlobal[2*in[i]+1];

  /* get coordinates */
  status = EG_getGlobal(discr->bodys[bIndex-1].tess, global[0], &ptype, &pindex,
                        xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discr->bodys[bIndex-1].tess, global[1], &ptype, &pindex,
                        xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discr->bodys[bIndex-1].tess, global[2], &ptype, &pindex,
                        xyz3);
  AIM_STATUS(discr->aInfo, status);
  
  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];
  CROSS(x3, x1, x2);
  area  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */
  if (data == NULL) {
    *result = 3.0*area;
    return CAPS_SUCCESS;
  }

  for (i = 0; i < rank; i++)
    result[i] = (data[rank*in[0]+i] + data[rank*in[1]+i] +
                 data[rank*in[2]+i])*area;

cleanup:
  return status;
}


/* aimIntegration: Element Integration on the Bound -- Optional */
int
aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name,
                int bIndex, int eIndex, int rank,
                double *r_bar, double *d_bar)
{
  int        status = CAPS_SUCCESS;
  int        i, in[3], ptype, pindex, global[3];
  double     x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area;
  capsBodyDiscr *discBody=NULL;

  /* reverse differentiation of aimIntegration */
  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" skeletonAIM/Integration: name = %s, bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_RANGEERR;
  }

  discBody = &discr->bodys[bIndex-1];
  if ((eIndex <= 0) || (eIndex > discBody->nElems))
    printf(" skeletonAIM/aimIntegrateBar: eIndex = %d [1-%d]!\n",
           eIndex, discBody->nElems);

  /* element indices */
  in[0] = discBody->elems[eIndex-1].gIndices[0] - 1;
  in[1] = discBody->elems[eIndex-1].gIndices[2] - 1;
  in[2] = discBody->elems[eIndex-1].gIndices[4] - 1;
  
  /* convert to global indices */
  for (i = 0; i < 3; i++)
    global[i] = discr->tessGlobal[2*in[i]+1];
    
  /* get coordinates */
  status = EG_getGlobal(discr->bodys[bIndex-1].tess, global[0], &ptype, &pindex,
                        xyz1);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discr->bodys[bIndex-1].tess, global[1], &ptype, &pindex,
                        xyz2);
  AIM_STATUS(discr->aInfo, status);
  status = EG_getGlobal(discr->bodys[bIndex-1].tess, global[2], &ptype, &pindex,
                        xyz3);
  AIM_STATUS(discr->aInfo, status);

  x1[0] = xyz2[0] - xyz1[0];
  x2[0] = xyz3[0] - xyz1[0];
  x1[1] = xyz2[1] - xyz1[1];
  x2[1] = xyz3[1] - xyz1[1];
  x1[2] = xyz2[2] - xyz1[2];
  x2[2] = xyz3[2] - xyz1[2];
  CROSS(x3, x1, x2);
  area  = sqrt(DOT(x3, x3))/6.0;      /* 1/2 for area and then 1/3 for sum */

  for (i = 0; i < rank; i++) {
/*  result[i] = (data[rank*in[0]+i] + data[rank*in[1]+i] +
                 data[rank*in[2]+i])*area;  */
    d_bar[rank*in[0]+i] += area*r_bar[i];
    d_bar[rank*in[1]+i] += area*r_bar[i];
    d_bar[rank*in[2]+i] += area*r_bar[i];
  }

cleanup:
  return status;
}
