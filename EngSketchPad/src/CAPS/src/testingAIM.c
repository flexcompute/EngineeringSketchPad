/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Testing AIM Example Code
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"
#include "aimMesh.h"

#define DEBUG
//#define SENSITIVITY

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


typedef struct {
  int        instance;
  int        nBody;
  ego        *tess;
  aimMeshRef *mesh;
} aimStorage;



void
aimCleanup(/*@only@*/ void *instStore)
{
  aimStorage *aimStore;

  aimStore = (aimStorage *) instStore;
#ifdef DEBUG
  printf(" testingAIM/aimCleanup   instance = %d!\n", aimStore->instance);
#endif

  AIM_FREE(aimStore->tess);
  if (aimStore->mesh != NULL) {
    aim_freeMeshRef(aimStore->mesh);
    EG_free(aimStore->mesh);
  }
  EG_free(aimStore);
}


int
aimInitialize(int inst, /*@unused@*/ const char *unitSys,
              /*@unused@*/ void *aimInfo, /*@unused@*/ void **instStore,
              /*@unused@*/ int *major, /*@unused@*/ int *minor, int *nIn,
              int *nOut, int *nFields, char ***fnames, int **franks,
              int **fInOut)
{
  int        *ints, i;
  char       **strs;
  aimStorage *aimStore;

#ifdef DEBUG
/*@-nullpass@*/
  printf("\n testingAIM/aimInitialize  instance = %d  unitSys = %s!\n",
         inst, unitSys);
/*@+nullpass@*/
#endif

  /* specify the number of analysis input and out "parameters" */
  *major  = 1;
  *minor  = 0;
  *nIn    = 4;
  *nOut   = 3;
  if (inst == -1) return CAPS_SUCCESS;

  /* setup our AIM specific state */
  aimStore = (aimStorage *) EG_alloc(sizeof(aimStorage));
  if (aimStore == NULL) return EGADS_MALLOC;
  aimStore->instance = inst;
  aimStore->nBody    = 0;
  aimStore->tess     = NULL;
  aimStore->mesh     = NULL;
  
  /* only instance 0 can write */
  if (inst == 0) {
    aimStore->mesh = (aimMeshRef *) EG_alloc(sizeof(aimMeshRef));
    if (aimStore->mesh == NULL) {
      EG_free(aimStore);
      return EGADS_MALLOC;
    }
    aim_initMeshRef(aimStore->mesh);
  }

  /* specify the field variables this analysis can generate */
  *nFields = 4;
  ints     = (int *) EG_alloc((*nFields)*sizeof(int));
  if (ints == NULL) {
    aimCleanup(aimStore);
    return EGADS_MALLOC;
  }
  ints[0]  = 1;
  ints[1]  = 3;
  ints[2]  = 1;
  ints[3]  = 3;
  *franks   = ints;
  strs     = (char **) EG_alloc((*nFields)*sizeof(char *));
  if (strs == NULL) {
    EG_free(*franks);
    *franks = NULL;
    aimCleanup(aimStore);
    return EGADS_MALLOC;
  }
  strs[0]  = EG_strdup("scalar");
  strs[1]  = EG_strdup("vector");
  strs[2]  = EG_strdup("scalar");
  strs[3]  = EG_strdup("vector");
  *fnames  = strs;

  ints     = (int *) EG_alloc((*nFields)*sizeof(int));
  if (ints == NULL) {
    for (i = 0; i < *nFields; i++)
      EG_free((*fnames)[i]);
    EG_free(*fnames);
    *fnames = NULL;
    EG_free(*franks);
    *franks = NULL;
    aimCleanup(aimStore);
    return EGADS_MALLOC;
  }
  ints[0] = FieldOut;
  ints[1] = FieldOut;
  ints[2] = FieldIn;
  ints[3] = FieldIn;
  *fInOut = ints;

  *instStore = aimStore;

  return CAPS_SUCCESS;
}


void
aimFreeDiscrPtr(void *ptr)
{
  /* free up the user pointer */
  EG_free(ptr);
}


int
aimDiscr(char *tname, capsDiscr *discr)
{
  int           ibody, iface, i, vID, itri, ielem, status, found;
  int           nFace, atype, alen, tlen, nBodyDisc;
  int           ntris, nBody, state, nGlobal, global;
  int           *vid = NULL;
  double        size, box[6], params[3];
  const int     *ints, *ptype, *pindex, *tris, *nei;
  const double  *reals, *xyz, *uv;
  const char    *string, *intents;
  ego           body, *bodies, *faces = NULL, *tess = NULL;
  capsBodyDiscr *discBody;
  aimStorage    *aimStore = NULL;

#ifdef DEBUG
  printf(" testingAIM/aimDiscr: tname = %s, instance = %d!\n",
         tname, aim_getInstance(discr->aInfo));
#endif
  aimStore = (aimStorage *) discr->instStore;

  /* create the discretization structure for this capsBound */
  status = aim_getBodies(discr->aInfo, &intents, &nBody, &bodies);
  if (status != CAPS_SUCCESS) goto cleanup;

  aimStore->nBody = nBody;
  AIM_REALL(aimStore->tess, aimStore->nBody, ego, discr->aInfo, status);
  tess = aimStore->tess;

  for (ibody = 0; ibody < nBody; ibody++) {
    tess[ibody] = NULL;
    status = EG_getBoundingBox(bodies[ibody], box);
    if (status != EGADS_SUCCESS) {
      printf(" testingAIM: getBoundingBox = %d for Body %d!\n",
             status, ibody+1);
      continue;
    }
                              size = box[3]-box[0];
    if (size < box[4]-box[1]) size = box[4]-box[1];
    if (size < box[5]-box[2]) size = box[5]-box[2];
    params[0] =  0.025*size;
    params[1] =  0.001*size;
    params[2] = 15.0;
    status = EG_makeTessBody(bodies[ibody], params, &tess[ibody]);
    if (status != EGADS_SUCCESS || tess[ibody] == NULL) {
      printf(" testingAIM: makeTessBody = %d for Body %d!\n", status, ibody+1);
      continue;
    }
    /* store the tessellation back in CAPS */
    status = aim_newTess(discr->aInfo, tess[ibody]);
    if (status != EGADS_SUCCESS) {
      printf(" testingAIM: aim_setTess = %d for Body %d!\n", status, ibody+1);
      continue;
    }
    printf(" testingAIM/aimDiscr: Tessellation set!\n");
  }

  /* find any faces with our boundary marker */
  for (nBodyDisc = ibody = 0; ibody < nBody; ibody++) {
    if (tess[ibody] == NULL) continue;
    status = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nFace, &faces);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (faces == NULL) {
      status = EGADS_TOPOERR;
      goto cleanup;
    }
    found = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;
#ifdef DEBUG
      printf(" testingAIM/aimDiscr: Body %d/Face %d matches %s!\n",
             ibody+1, iface+1, tname);
#endif
      /* count the number of face with capsBound */
      found = 1;
    }
    if (found == 1) nBodyDisc++;
    AIM_FREE(faces);
  }
  if (nBodyDisc == 0) {
    printf(" testingAIM/aimDiscr: No Faces match %s!\n", tname);
    return CAPS_SUCCESS;
  }

  /* specify our single triangle element type */
  discr->nTypes = 1;
  AIM_ALLOC(discr->types, discr->nTypes, capsEleType, discr->aInfo, status);

  discr->types[0].nref  = 3;
  discr->types[0].ndata = 0;         /* data at geom reference positions
                                        (i.e. vertex centered/iso-parametric) */
  discr->types[0].ntri  = 1;
  discr->types[0].nseg  = 3;
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
  
  AIM_ALLOC(discr->types[0].segs, 2*discr->types[0].nref, int, discr->aInfo,
            status);
  discr->types[0].segs[0] = 1;
  discr->types[0].segs[1] = 2;
  discr->types[0].segs[2] = 2;
  discr->types[0].segs[3] = 3;
  discr->types[0].segs[4] = 3;
  discr->types[0].segs[5] = 1;

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
  AIM_ALLOC(discr->bodys, nBodyDisc, capsBodyDiscr, discr->aInfo, status);

  /* get the tessellation and
   make up a linear continuous triangle discretization */
  vID = nBodyDisc = 0;
  for (ibody = 0; ibody < nBody; ibody++) {
    if (tess[ibody] == NULL) continue;
    ntris  = 0;
    AIM_FREE(faces);
    status = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nFace, &faces);
    if ((status != EGADS_SUCCESS) || (faces == NULL)) {
      printf(" testingAIM: getBodyTopos (Face) = %d for Body %d!\n",
             status, ibody+1);
      status = EGADS_TOPOERR;
      goto cleanup;
    }
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
        printf(" testingAIM: EG_getTessFace %d = %d for Body %d!\n",
               iface+1, status, ibody+1);
        continue;
      }
      ntris += tlen;
      found = 1;
    }
    if (found == 0) continue;
#ifdef DEBUG
    printf(" testingAIM/aimDiscr: ntris = %d!\n", ntris);
#endif
    if (ntris == 0) {
      status = CAPS_SOURCEERR;
      goto cleanup;
    }

    discBody = &discr->bodys[nBodyDisc];
    aim_initBodyDiscr(discBody);

    discBody->nElems = ntris;
    discBody->tess   = tess[ibody];

    AIM_ALLOC(discBody->elems,      ntris, capsElement, discr->aInfo, status);
    AIM_ALLOC(discBody->gIndices, 6*ntris, int,         discr->aInfo, status);

    status = EG_statusTessBody(tess[ibody], &body, &state, &nGlobal);
    AIM_STATUS(discr->aInfo, status);

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
  discr->nBodys  = nBodyDisc;
#ifdef DEBUG
  printf(" testingAIM/aimDiscr: npts = %d!\n", vID);
#endif

  status = CAPS_SUCCESS;

cleanup:

  /* free up our stuff */
  AIM_FREE(faces);
  AIM_FREE(vid);

#ifdef DEBUG
  if (status != CAPS_SUCCESS) {
    printf(" testingAIM/aimDiscr: status = %d!\n", status);
  }
#endif

  return status;
}


int
aimLocateElement(capsDiscr *discr, double *params, double *param,
                 int *bIndex, int *eIndex, double *bary)
{
  int    i, ib, in[3], stat, ibsmall, ismall;
  double we[3], w, smallw = -1.e300;

  if (discr == NULL) return CAPS_NULLOBJ;

  ibsmall = ismall = 0;
  for (ib = 0; ib < discr->nBodys; ib++) {
    for (i = 0; i < discr->bodys[ib].nElems; i++) {
      in[0] = discr->bodys[ib].elems[i].gIndices[0] - 1;
      in[1] = discr->bodys[ib].elems[i].gIndices[2] - 1;
      in[2] = discr->bodys[ib].elems[i].gIndices[4] - 1;
      stat  = EG_inTriExact(&params[2*in[0]], &params[2*in[1]], &params[2*in[2]],
                            param, we);
      if (stat == EGADS_SUCCESS) {
        *eIndex = i+1;
        *bIndex = ib+1;
        bary[0] = we[1];
        bary[1] = we[2];
  /*    {
          double uv[2];
          uv[0] = we[0]*params[2*in[0]  ] + we[1]*params[2*in[1]  ] +
                  we[2]*params[2*in[2]  ];
          uv[1] = we[0]*params[2*in[0]+1] + we[1]*params[2*in[1]+1] +
                  we[2]*params[2*in[2]+1];
          printf(" %d: %lf %lf  %lf %lf\n", i+1, param[0], param[1], uv[0], uv[1]);
          printf("        %lf %lf\n", bary[0], bary[1]);
        }  */
        return CAPS_SUCCESS;
      }
      w = we[0];
      if (we[1] < w) w = we[1];
      if (we[2] < w) w = we[2];
      if (w > smallw) {
        ismall = i+1;
        ibsmall = ib+1;
        smallw = w;
      }
    }
  }

  /* must extrapolate! */
  if (ismall == 0) return CAPS_NOTFOUND;
  in[0] = discr->bodys[ibsmall-1].elems[ismall-1].gIndices[0] - 1;
  in[1] = discr->bodys[ibsmall-1].elems[ismall-1].gIndices[2] - 1;
  in[2] = discr->bodys[ibsmall-1].elems[ismall-1].gIndices[4] - 1;
  EG_inTriExact(&params[2*in[0]], &params[2*in[1]], &params[2*in[2]],
                param, we);
  *eIndex = ismall;
  *bIndex = ibsmall;
  bary[0] = we[1];
  bary[1] = we[2];
/*
  printf(" aimLocateElement: extropolate to %d (%lf %lf %lf)  %lf\n",
         ismall, we[0], we[1], we[2], smallw);
*/
  return CAPS_SUCCESS;
}


int
aimInputs(void *instStore, /*@unused@*/ void *aimStruc, int index, char **ainame,
          capsValue *defval)
{
  int        inst = -1;
  capsTuple  *tuple;
  aimStorage *aimStore;

  aimStore = (aimStorage *) instStore;
  if (aimStore != NULL) inst = aimStore->instance;
#ifdef DEBUG
  printf(" testingAIM/aimInputs instance = %d  index = %d!\n", inst, index);
#endif

  if (index == 1) {

    *ainame = EG_strdup("testingAIMin");
    if (*ainame == NULL) return EGADS_MALLOC;
    defval->type      = Double;
    defval->vals.real = 5.0 + inst;
    defval->units     = EG_strdup("cm");

  } else if (index == 2) {

    *ainame = EG_strdup("table");
    if (*ainame == NULL) return EGADS_MALLOC;
    tuple = (capsTuple *) EG_alloc(3*sizeof(capsTuple));
    if (tuple == NULL) return EGADS_MALLOC;
    tuple[0].name  = EG_strdup("Entry1");
    tuple[1].name  = EG_strdup("Entry2");
    tuple[2].name  = EG_strdup("Entry3");
    tuple[0].value = EG_strdup("Value1");
    tuple[1].value = EG_strdup("Value2");
    tuple[2].value = EG_strdup("Value3");

    defval->type       = Tuple;
    defval->dim        = Vector;
    defval->nrow       = 1;
    defval->ncol       = 3;
    defval->vals.tuple = tuple;

  } else if (index == 3) {

    *ainame = EG_strdup("tessIn");
    if (*ainame == NULL) return EGADS_MALLOC;
    defval->type        = Pointer;
    defval->vals.AIMptr = NULL;
    defval->nullVal     = IsNull;
    defval->units       = EG_strdup("ego");
    if (defval->units == NULL) {
      EG_free(*ainame);
      return EGADS_MALLOC;
    }

  } else {
    
    *ainame = EG_strdup("meshFile");
    if (*ainame == NULL) return EGADS_MALLOC;
    defval->type        = PointerMesh;
    defval->vals.AIMptr = NULL;
    defval->nullVal     = IsNull;
    defval->units       = EG_strdup("writer");
    if (defval->units == NULL) {
      EG_free(*ainame);
      return EGADS_MALLOC;
    }
    defval->meshWriter  = EG_strdup("testingWriter");
    if (defval->meshWriter == NULL) {
      EG_free(defval->units);
      EG_free(*ainame);
      return EGADS_MALLOC;
    }
    
  }

  return CAPS_SUCCESS;
}


int
aimPreAnalysis(void *instStore, void *aimStruc, /*@null@*/ capsValue *inputs)
{
  int              nBname, stat, state, i, npts, rank, len;
  char             **bNames, *units, *full;
  double           *data;
  capsDiscr        *discr;
  capsObject       *vobj;
  capsValue        *value;
  enum capsdMethod method;
  aimStorage       *aimStore;
  aimMeshRef       *mesh;
  ego              tess, body;
  FILE             *fp;

  aimStore = (aimStorage *) instStore;
#ifdef DEBUG
  printf(" testingAIM/aimPreAnalysis instance = %d!\n", aimStore->instance);
#endif
  
  fp = aim_fopen(aimStruc, "inputFile", "w");
  if (fp == NULL) {
    printf(" testingAIM/aimPreAnalysis fileopen = NULL!\n");
  } else {
    fprintf(fp, "Put something in the file\n");
    fclose(fp);
  }

  stat  = aim_getBounds(aimStruc, &nBname, &bNames);
  printf(" testingAIM/aimPreAnalysis aim_getBounds = %d\n", stat);
  for (i = 0; i < nBname; i++)
    printf("   Analysis in Bound = %s\n", bNames[i]);
  if (bNames != NULL) EG_free(bNames);

  stat = aim_newGeometry(aimStruc);
  printf("     aim_newGeometry = %d!\n", stat);

  if (aimStore->instance == 0) {
    /* look for parent's dependency */
    if (inputs != NULL) {
      tess = (ego) inputs[2].vals.AIMptr;
      if (tess == NULL) {
        printf("   tess is NULL!\n");
      } else {
        stat = EG_statusTessBody(tess, &body, &state, &npts);
        printf("   tess State = %d  %d   npts = %d\n", stat, state, npts);
      }
    }
    
    stat = aim_getDiscr(aimStruc, "Interface", &discr);
    printf("   getDiscr = %d\n", stat);
    if (stat == CAPS_SUCCESS) {
      stat = aim_getDataSet(discr, "scalar", &method, &npts, &rank, &data, &units);
      printf("   getDataSet = %d, rank = %d, method = %d\n", stat, rank, method);
      if (npts == 1) {
        printf("   scalar = %lf\n", data[0]);
      } else {
        printf("   %d scalars!\n", npts);
      }
    }
    
    if (aimStore->mesh->fileName == NULL) {
      printf("   meshRef = NULL\n");
    } else {
      stat = aim_deleteMeshes(aimStruc, aimStore->mesh);
      printf("   aim_deleteMeshes = %d\n", stat);
    }
  } else if (aimStore->instance == 1) {
    /* look for child's dependency */
    if (inputs != NULL) {
      tess = (ego) inputs[2].vals.AIMptr;
      if (tess == NULL) {
        printf("   tess is NULL!\n");
      } else {
        stat = EG_statusTessBody(tess, &body, &state, &npts);
        printf("   tess State = %d  %d   npts = %d\n", stat, state, npts);
      }
      
      mesh = (aimMeshRef *) inputs[3].vals.AIMptr;
      /* special internal linking -- CAPS normally sets this up */
      if (mesh == NULL)
        if (inputs[3].link != NULL) {
          vobj  = inputs[3].link;
          value = vobj->blind;
          if (value != NULL) mesh = (aimMeshRef *) value->vals.AIMptr;
        }
      
      if (mesh != NULL) {
        printf("   mesh file = %s\n", mesh->fileName);
        len  = strlen(mesh->fileName) + 5;
        full = (char *) EG_alloc(len*sizeof(char));
        if (full != NULL) {
          for (i = 0; i < len-5; i++) full[i] = mesh->fileName[i];
          full[len-5] = '.';
          full[len-4] = 't';
          full[len-3] = 'x';
          full[len-2] = 't';
          full[len-1] =  0;
          stat = aim_symLink(aimStruc, full, NULL);
          printf("   symLink  = %d\n", stat);
          EG_free(full);
        } else {
          printf(" Malloc Error on %d characters!\n", len);
        }
      }
    }

    stat = aim_getDiscr(aimStruc, "Interface", &discr);
    printf("   getDiscr = %d\n", stat);
    if (stat == CAPS_SUCCESS) {
      stat = aim_getDataSet(discr, "vector", &method, &npts, &rank, &data, &units);
      printf("   getDataSet = %d, rank = %d, method = %d\n", stat, rank, method);
      if (npts == 1) {
        printf("   vector = %lf %lf %lf\n", data[0], data[1], data[2]);
      } else if (stat != CAPS_SUCCESS) {
        printf("   vectors not setup yet!\n");
      } else {
        printf("   %d vectors!\n", npts);
      }
    }
  }

  return CAPS_SUCCESS;
}


int
aimOutputs(/*@null@*/ void *instStore, /*@unused@*/ void *aimStruc,
           int index, char **aoname, capsValue *form)
{
  int        inst = -1;
  aimStorage *aimStore;

  aimStore = (aimStorage *) instStore;
  if (aimStore != NULL) inst = aimStore->instance;
#ifdef DEBUG
  printf(" testingAIM/aimOutputs instance = %d  index = %d!\n", inst, index);
#endif

  if (index == 1) {

    *aoname = EG_strdup("testingAIMout");
    if (*aoname == NULL) return EGADS_MALLOC;
    form->type  = Double;
    form->units = EG_strdup("cm");

  } else if (index == 2) {

    *aoname = EG_strdup("tessOut");
    if (*aoname == NULL) return EGADS_MALLOC;
    form->type    = Pointer;
    form->nullVal = NotNull;
    form->units   = EG_strdup("ego");
    if (form->units == NULL) {
      EG_free(*aoname);
      return EGADS_MALLOC;
    }

  } else {
    
    *aoname = EG_strdup("meshFile");
    if (*aoname == NULL) return EGADS_MALLOC;
    form->type    = PointerMesh;
    form->nullVal = NotNull;
    form->units   = EG_strdup("writer");
    if (form->units == NULL) {
      EG_free(*aoname);
      return EGADS_MALLOC;
    }
    
  }

  return CAPS_SUCCESS;
}


int
aimExecute(void *instStore, /*@unused@*/ void *aimStruc, int *state)
{
  aimStorage *aimStore;

  aimStore = (aimStorage *) instStore;
#ifdef DEBUG
  printf(" testingAIM/aimExecute instance = %d!\n", aimStore->instance);
#endif

  *state = 0;
  return CAPS_SUCCESS;
}


int
aimPostAnalysis(void *instStore, void *aimStruc, int restart,
                /*@null@*/ capsValue *inputs)
{
  int        i, n, stat;
  capsValue  dynOut;
  aimStorage *aimStore;

  aimStore = (aimStorage *) instStore;
#ifdef DEBUG
  printf(" testingAIM/aimPostAnalysis instance = %d  restart = %d!\n",
         aimStore->instance, restart);
#endif
  
  if (restart == 1) {
    /* do pre on restarted post */
    printf(" testingAIM/aimPostAnalysis: restart -- going into pre!\n");
    stat = aimPreAnalysis(instStore, aimStruc, inputs);
    if (stat != CAPS_SUCCESS) {
      printf(" testingAIM/aimPostAnalysis: aimPreAnalysis = %d\n", stat);
      return stat;
    }
  } else {
    stat = aim_initValue(&dynOut);
    if (stat != CAPS_SUCCESS) {
      printf(" testingAIM/aimPostAnalysis: aim_initValue = %d\n", stat);
      return stat;
    }
    dynOut.vals.integer = 42;
    stat = aim_makeDynamicOutput(aimStruc, "Everything", &dynOut);
    if (stat != CAPS_SUCCESS) {
      printf(" testingAIM/aimPostAnalysis: aim_makeDynamicOutput = %d\n", stat);
      return stat;
    }
  }

  if (aimStore->instance == 0) {
    n = aim_getIndex(aimStruc, NULL, GEOMETRYIN);
    if (n < CAPS_SUCCESS)
      printf(" testingAIM/aimPostAnalysis: aim_getIndex = %d\n", n);
    for (i = 1; i <= n; i++) {
      stat = aim_getGeomInType(aimStruc, i);
      if (stat < CAPS_SUCCESS) {
        printf(" testingAIM/aimPostAnalysis: %d aim_getGeomInType = %d\n",
               i, stat);
      } else if (stat == 1) {
        printf(" testingAIM/aimPostAnalysis: %d -- %d is Config Parameter\n",
               aimStore->instance, i);
      } else if (stat == 2) {
        printf(" testingAIM/aimPostAnalysis: %d -- %d is Constant Parameter\n",
               aimStore->instance, i);
      }
    }
  }

  return CAPS_SUCCESS;
}


int
aimCalcOutput(void *instStore, void *aimStruc, int index, capsValue *val)
{
  int              i, stat, rank, npts, len;
  double           *dval;
  enum capsdMethod method;
  capsDiscr        *discr;
  char             *units, relative[PATH_MAX];
  const char       *name;
#ifdef SENSITIVITY
  int              nBody,   nFace, j, ig;
  double           *dxyz;
  ego              *bodies, *faces;
  const char       *intents;
#endif
  aimStorage       *aimStore;
  aimMesh          meshStruc;

  aimStore = (aimStorage *) instStore;
#ifdef DEBUG
  stat = aim_getName(aimStruc, index, ANALYSISOUT, &name);
  printf(" testingAIM/aimCalcOutput instance = %d  index = %d %s %d!\n",
         aimStore->instance, index, name, stat);
#endif

  if (index == 2) {
    
    val->vals.AIMptr    = NULL;
    val->nullVal        = IsNull;
    if (aimStore->tess != NULL) {
      val->nullVal      = NotNull;
      val->vals.AIMptr  = aimStore->tess[0];
    }
#ifdef DEBUG
#ifdef WIN32
    printf(" tessPtr = %llx\n", (long long) val->vals.AIMptr);
#else
    printf(" tessPtr = %lx\n", (long) val->vals.AIMptr);
#endif
#endif
    return CAPS_SUCCESS;
    
  } else if (index == 3) {
    
    if (aimStore->mesh != NULL) {
      if (aimStore->mesh->fileName == NULL) {
        stat = aim_file(aimStruc, "meshFile", relative);
        if (stat != CAPS_SUCCESS) return stat;
        len  = strlen(relative) + 1;
        aimStore->mesh->fileName = (char *) EG_alloc(len*sizeof(char));
        if (aimStore->mesh->fileName == NULL) return EGADS_MALLOC;
        for (i = 0; i < len; i++) aimStore->mesh->fileName[i] = relative[i];
/*      printf(" mesh fileName = %s\n", aimStore->mesh->fileName);  */
      }
      meshStruc.meshData = NULL;
      meshStruc.meshRef  = aimStore->mesh;
      stat = aim_writeMeshes(aimStruc, index, &meshStruc);
      if (stat == CAPS_NOTFOUND) return CAPS_SUCCESS;
      if (stat != CAPS_SUCCESS)  return stat;
      val->nullVal     = NotNull;
/*@-kepttrans@*/
      val->vals.AIMptr = aimStore->mesh;
/*@+kepttrans@*/
    }
    return CAPS_SUCCESS;
    
  }

  val->vals.real = 12.34;

  /* get a dataset */
  stat = aim_getDiscr(aimStruc, "Interface", &discr);
#ifdef DEBUG
  printf(" aim_getDiscr %d on Interface = %d\n", aimStore->instance, stat);
#endif
  if (stat == CAPS_SUCCESS) {
    stat = aim_getDataSet(discr, "scalar", &method, &npts, &rank, &dval, &units);
    if (stat == CAPS_SUCCESS) {
      printf(" aim_getDataSet %d for scalar = %d %d %d\n",
             aimStore->instance, method, npts, rank);
    } else {
      printf(" aim_getDataSet %d on Interface = %d\n",
             aimStore->instance, stat);
    }
  }

#ifdef SENSITIVITY
  /* get sensitivity on tessellation */
  ig   = 2;
  name = NULL;
  stat = aim_getName(aimStruc, ig, GEOMETRYIN, &name);
  if (stat == CAPS_SUCCESS) {
    printf(" instance = %d, %d geomIn = %s %d %d\n", inst, ig, name,
           instance[inst].gIn[ig-1].nrow, instance[inst].gIn[ig-1].ncol);
    if (inst == 0) {
      printf("\n");
      stat = aim_setSensitivity(aimStruc, name, 1, 1);
      if (stat != CAPS_SUCCESS) {
        printf(" testingAIM/aimCalcOutput: aim_setSensitivity = %d\n", stat);
      } else {
        stat = aim_getBodies(aimStruc, &intents, &nBody, &bodies);
        if (stat != CAPS_SUCCESS) {
          printf(" testingAIM/aimCalcOutput: aim_getBodies = %d!\n", stat);
        } else {
          for (i = 0; i < nBody; i++) {
            stat = EG_getBodyTopos(bodies[i], NULL, FACE, &nFace, &faces);
            if (stat != EGADS_SUCCESS) {
              printf(" testingAIM/aimCalcOutput: EG_getBodyTopos = %d!\n", stat);
              continue;
            }
            for (j = 0; j < nFace; j++) {
              stat = aim_getSensitivity(aimStruc, bodies[i+nBody], 2, j+1,
                                        &npts, &dxyz);
              if (stat != CAPS_SUCCESS) {
                printf(" testingAIM/aimCalcOutput: %d/%d aim_getSensitivity = %d!\n",
                       i+1, j+1, stat);
                continue;
              }
              printf(" sensitivity: Body %d/Face %d  nPts = %d\n",
                     i+1, j+1, npts);
              EG_free(dxyz);
            }
            EG_free(faces);
          }
        }
      }
      printf("\n");
    }
  }
#endif

  return CAPS_SUCCESS;
}


int
aimTransfer(capsDiscr *discr, /*@unused@*/ const char *name, int npts, int rank,
            double *data, /*@unused@*/ char **units)
{
  int        stat = CAPS_SUCCESS;
  int        i, j, bIndex, ptype, pindex, global;
  double     xyz[3];
  aimStorage *aimStore;

  if (discr == NULL) return CAPS_NULLOBJ;

  aimStore = (aimStorage *) discr->instStore;
#ifdef DEBUG
  printf(" testingAIM/aimTransfer name = %s  instance = %d  npts = %d/%d!\n",
         name, aimStore->instance, npts, rank);
#endif

  /* fill in with our coordinates -- for now */
  for (i = 0; i < npts; i++) {
    bIndex = discr->tessGlobal[2*i  ];
    global = discr->tessGlobal[2*i+1];
    stat   = EG_getGlobal(discr->bodys[bIndex-1].tess, global, &ptype, &pindex,
                          xyz);
    AIM_STATUS(discr->aInfo, stat);
    for (j = 0; j < rank; j++) data[rank*i+j] = xyz[j];
  }

cleanup:
  return stat;
}


int
aimInterpolation(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                 double *bary, int rank, double *data, double *result)
{
  int    in[3], i;
  double we[3];

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" testingAIM/Interpolation: %s bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_NOBODIES;
  }

  if ((eIndex <= 0) || (eIndex > discr->bodys[bIndex-1].nElems)) {
    printf(" testingAIM/Interpolation: %s eIndex = %d [1-%d]!\n",
           name, eIndex, discr->bodys[bIndex-1].nElems);
    return CAPS_BADINDEX;
  }

  we[0] = 1.0 - bary[0] - bary[1];
  we[1] = bary[0];
  we[2] = bary[1];
  in[0] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[0] - 1;
  in[1] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[2] - 1;
  in[2] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[4] - 1;
  for (i = 0; i < rank; i++)
    result[i] = data[rank*in[0]+i]*we[0] + data[rank*in[1]+i]*we[1] +
                data[rank*in[2]+i]*we[2];

  return CAPS_SUCCESS;
}


int
aimInterpolateBar(capsDiscr *discr, const char *name, int bIndex,
                  int eIndex, double *bary, int rank, double *r_bar,
                  double *d_bar)
{
  int    in[3], i;
  double we[3];

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" testingAIM/InterpolateBar: %s bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_NOBODIES;
  }

  if ((eIndex <= 0) || (eIndex > discr->bodys[bIndex-1].nElems)) {
    printf(" testingAIM/InterpolateBar: %s eIndex = %d [1-%d]!\n",
           name, eIndex, discr->bodys[bIndex-1].nElems);
    return CAPS_BADINDEX;
  }

  we[0] = 1.0 - bary[0] - bary[1];
  we[1] = bary[0];
  we[2] = bary[1];
  in[0] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[0] - 1;
  in[1] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[2] - 1;
  in[2] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[4] - 1;
  for (i = 0; i < rank; i++) {
/*  result[i] = data[rank*in[0]+i]*we[0] + data[rank*in[1]+i]*we[1] +
                data[rank*in[2]+i]*we[2];  */
    d_bar[rank*in[0]+i] += we[0]*r_bar[i];
    d_bar[rank*in[1]+i] += we[1]*r_bar[i];
    d_bar[rank*in[2]+i] += we[2]*r_bar[i];
  }

  return CAPS_SUCCESS;
}


int
aimIntegration(capsDiscr *discr, const char *name, int bIndex, int eIndex,
               int rank, double *data, double *result)
{
  int    i, in[3], stat, ptype, pindex, global[3];
  double x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area;

  if ((bIndex <= 0) || (bIndex > discr->nBodys) || (discr->bodys == NULL)) {
    printf(" testingAIM/aimIntegration: %s bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_NOBODIES;
  }

  if ((eIndex <= 0) || (eIndex > discr->bodys[bIndex-1].nElems)) {
    printf(" testingAIM/aimIntegration: %s eIndex = %d [1-%d]!\n",
           name, eIndex, discr->bodys[bIndex-1].nElems);
    return CAPS_BADINDEX;
  }

  /* Element indices */
  in[0] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[0] - 1;
  in[1] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[2] - 1;
  in[2] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++)
    global[i] = discr->tessGlobal[2*in[i]+1];

  /* get coordinates */
  stat = EG_getGlobal(discr->bodys[bIndex-1].tess, global[0], &ptype, &pindex,
                      xyz1);
  AIM_STATUS(discr->aInfo, stat);
  stat = EG_getGlobal(discr->bodys[bIndex-1].tess, global[1], &ptype, &pindex,
                      xyz2);
  AIM_STATUS(discr->aInfo, stat);
  stat = EG_getGlobal(discr->bodys[bIndex-1].tess, global[2], &ptype, &pindex,
                      xyz3);
  AIM_STATUS(discr->aInfo, stat);

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
  return stat;
}


int
aimIntegrateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                int rank, double *r_bar, double *d_bar)
{
  int    i, in[3], stat, ptype, pindex, global[3];
  double x1[3], x2[3], x3[3], xyz1[3], xyz2[3], xyz3[3], area;

  if ((bIndex <= 0) || (bIndex > discr->nBodys)) {
    printf(" testingAIM/aimIntegration: %s bIndex = %d [1-%d]!\n",
           name, bIndex, discr->nBodys);
    return CAPS_NOBODIES;
  }

  if ((eIndex <= 0) || (eIndex > discr->bodys[bIndex-1].nElems)) {
    printf(" testingAIM/aimIntegration: %s eIndex = %d [1-%d]!\n",
           name, eIndex, discr->bodys[bIndex-1].nElems);
    return CAPS_BADINDEX;
  }

  /* Element indices */
  in[0] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[0] - 1;
  in[1] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[2] - 1;
  in[2] = discr->bodys[bIndex-1].elems[eIndex-1].gIndices[4] - 1;

  /* convert to global indices */
  for (i = 0; i < 3; i++)
    global[i] = discr->tessGlobal[2*in[i]+1];

  /* get coordinates */
  stat = EG_getGlobal(discr->bodys[bIndex-1].tess, global[0], &ptype, &pindex,
                      xyz1);
  AIM_STATUS(discr->aInfo, stat);
  stat = EG_getGlobal(discr->bodys[bIndex-1].tess, global[1], &ptype, &pindex,
                      xyz2);
  AIM_STATUS(discr->aInfo, stat);
  stat = EG_getGlobal(discr->bodys[bIndex-1].tess, global[2], &ptype, &pindex,
                      xyz3);
  AIM_STATUS(discr->aInfo, stat);

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
  return stat;
}


int
aimBackdoor(void *instStore, /*@unused@*/ void *aimStruc, const char *JSONin,
            char **JSONout)
{
  aimStorage *aimStore;

  aimStore = (aimStorage *) instStore;
#ifdef DEBUG
  printf(" testingAIM/aimBackdoor instance = %d: %s!\n",
         aimStore->instance, JSONin);
#endif

  *JSONout = EG_strdup("aimBackdoor Output");

  return CAPS_SUCCESS;
}
