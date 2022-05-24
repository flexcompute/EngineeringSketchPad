/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Interference AIM
 *
 *      Copyright 2020-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"
#include "caps.h"
#include "cloud.h"

#ifdef WIN32
#define unlink _unlink
#else
#include <unistd.h>
#endif


//#define DEBUG


typedef struct {
  const char *name;
  int        bIndex;
  int        oml;          /* 0 inner body, 1 outer body, 2 "plug" */
  Cloud      cloud;
} clouds;


enum aimInputs
{
  Attr_Name = 1,               /* index is 1-based */
  OML,
  Tess_Params,
  NUMINPUT = Tess_Params       /* Total number of inputs */
};

enum aimOutputs
{
  Names = 1,                   /* index is 1-based */
  Distances,
  Volumes,
  Areas,
  CGs,
  Inertias,
  NUMOUT = Inertias            /* Total number of outputs */
};

/*!\mainpage Introduction
 * \tableofcontents
 *
 * \section overviewInter Interference AIM Overview
 * One can build and place components (Bodies) parametrically and if done
 * correctly then no single Body penetrates another. This cannot always be
 * accomplished, for example when a component is imported or the level of
 * geometric complexity makes building in these constraints very difficult.
 * Under these circumstances it is important to determine if the final
 * placements of Bodies do not intersect.
 *
 * This AIM takes a collection of Solid Bodies and returns the minimum
 * distance found between Bodies (if not intersecting) or as a negative
 * number, the penetration depth when the Bodies interfere. This is accomplished
 * by using a discrete representation of the Bodies (using the EGADS tessellator)
 * so the accuracy of the values returned is a function of how good of an
 * approximation the tessellation is to the actual BRep.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsInter
 * and \ref aimOutputsInter, respectively.
 *
 */


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys,
                  /*@unused@*/ void *aimInfo, void **inStore, int *major,
                  int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{

#ifdef DEBUG
  printf("\n interferenceAIM/aimInitialize   instance = %d!\n", inst);
#endif
  
  *major   = 1;
  *minor   = 0;

  /* specify the number of analysis input and out "parameters" */
  *nIn     = NUMINPUT;
  *nOut    = NUMOUT;
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 0;
  *fnames  = NULL;
  *franks  = NULL;
  *fInOut  = NULL;

  *inStore = NULL;    /* no internal state */

  return CAPS_SUCCESS;
}


int aimInputs(/*@unused@*/ void *instStore, void *aimInfo, int index,
              char **ainame, capsValue *defval)
{
  /*! \page aimInputsInter AIM Inputs
   * The following list outlines the Interference inputs along with their
   * default value available through the AIM interface.
   */

  int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" interferenceAIM/aimInputs  instance = %d  index = %d!\n",
         aim_getInstance(aimInfo), index);
#endif
  *ainame = NULL;

  if (index == Attr_Name) {
    *ainame             = EG_strdup("Attr_Name");
    defval->type        = String;
    defval->vals.string = EG_strdup("_name");
    AIM_NOTNULL(defval->vals.string, aimInfo, status);

    /*! \page aimInputsInter
     * - <B> Attr_Name = "_name" </B> <br>
     *  Attribute Name use to collect and label Bodies.
     */

  } else if (index == OML) {
    *ainame              = EG_strdup("OML");
    defval->type         = Boolean;
    defval->dim          = Scalar;
    defval->nrow         = 1;
    defval->ncol         = 1;
    defval->vals.integer = True;

    /*! \page aimInputsInter
     * - <B> OML = True </B> <br>
     *  Use the Body with the largest Bounding Box as a container (if True).
     *  False indicates that the Bodies are not contained.
     */

  } else if (index == Tess_Params) {
    *ainame               = EG_strdup("Tess_Params");
    defval->type          = Double;
    defval->dim           = Vector;
    defval->nrow          = 3;
    defval->ncol          = 1;
    defval->units         = NULL;
    defval->lfixed        = Fixed;
    defval->vals.reals    = (double *) EG_alloc(defval->nrow*sizeof(double));
    if (defval->vals.reals != NULL) {
      defval->vals.reals[0] = 0.025;
      defval->vals.reals[1] = 0.001;
      defval->vals.reals[2] = 15.00;
    } else return EGADS_MALLOC;

    /*! \page aimInputsInter
     * - <B> Tess_Params = [0.025, 0.001, 15.0]</B> <br>
     * Body tessellation parameters used to discretize all Bodies.
     * Tess_Params[0] and Tess_Params[1] get scaled by the bounding
     * box of the largest body. (From the EGADS manual) A set of 3 parameters
     * that drive the EDGE discretization and the FACE triangulation. The
     * first is the maximum length of an EDGE segment or triangle side
     * (in physical space). A zero is flag that allows for any length. The
     * second is a curvature-based value that looks locally at the deviation
     * between the centroid of the discrete object and the underlying
     * geometry. Any deviation larger than the input value will cause the
     * tessellation to be enhanced in those regions. The third is the maximum
     * interior dihedral angle (in degrees) between triangle facets (or Edge
     * segment tangents), note that a zero ignores this phase.
     */
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  if (status != CAPS_SUCCESS) {
    if (index == Attr_Name) AIM_FREE(defval->vals.string);
    AIM_FREE(*ainame);
  }
  return status;
}


int aimPreAnalysis(/*@unused@*/ void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
  int          i, n, numBody, status, atype, alen;
  char         file[PATH_MAX];
  const int    *ints;
  const char   *string, *intents;
  const double *reals;
  ego          *bodies;

#ifdef DEBUG
  printf(" interferenceAIM/aimPreAnalysis   instance = %d!\n",
         aim_getInstance(aimInfo));
#endif

  if (aimInputs == NULL) {
#ifdef DEBUG
    printf(" interferenceAIM/aimPreAnalysis aimInputs == NULL!\n");
#endif
    return CAPS_NULLVALUE;
  }

  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

  if ((numBody <= 1) || (bodies == NULL)) {
    AIM_ERROR(aimInfo, "AIM is given %d Bodies!", numBody);
    status = CAPS_SOURCEERR;
    goto cleanup;
  }

  /* Loop over bodies and look for our attribute */
  for (n = i = 0; i < numBody; i++) {
    status = EG_attributeRet(bodies[i], aimInputs[0].vals.string,
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS)
      if (atype == ATTRSTRING) n++;
  }
  if (n <= 1) {
    AIM_ERROR(aimInfo, "AIM is given %d Bodies with %s!",
              n, aimInputs[0].vals.string);
    status = CAPS_SOURCEERR;
    goto cleanup;
  }
  
  /* remove any old file */
  status = aim_file(aimInfo, "interference.dat", file);
  AIM_STATUS(aimInfo, status, "aim_file = %d!", status);
  unlink(file);
  status = CAPS_SUCCESS;

cleanup:
    return status;
}


int aimExecute(/*@unused@*/ void *instStore, void *aimInfo, int *state)
{
  int          i, ioml, j, m, n, numBody, status, atype, alen;
#ifdef DEBUG
  int          type;
#endif
  double       v, vol, size, params[3], box[6], xyzs[3], xyzt[3], *dist = NULL;
  double       props[14];
  const int    *ints;
  const char   *string, *intents;
  const double *reals;
  ego          *bodies;
  capsValue    *value, *nameVal;
  clouds       *cData = NULL;
  cloudPair    pair;
  FILE         *fp = NULL;
  size_t       len;

  *state = 0;
#ifdef DEBUG
  printf(" interferenceAIM/aimExecute   instance = %d!\n",
         aim_getInstance(aimInfo));
#endif

  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

  if ((numBody <= 1) || (bodies == NULL)) {
      AIM_ERROR(aimInfo, "AIM is given %d Bodies!", numBody);
      status = CAPS_SOURCEERR;
      goto cleanup;
  }
  
  status = aim_getValue(aimInfo, Attr_Name, ANALYSISIN, &nameVal);
  AIM_STATUS(aimInfo, status, "aim_getValue on Attr_Name!");

  /* Loop over bodies and look for our attribute */
  size = vol = 0.0;
  ioml = -1;
  for (n = i = 0; i < numBody; i++) {
    status = EG_attributeRet(bodies[i], nameVal->vals.string,
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS)
      if (atype == ATTRSTRING) {
        status = EG_getBoundingBox(bodies[i], box);
        AIM_STATUS(aimInfo, status, "EG_attributeRet on Body %d!", i+1);
        v = (box[3]-box[0])*(box[4]-box[1])*(box[5]-box[2]);
        if (v > vol) {
          vol  = v;
          ioml = i;
                                    size = box[3]-box[0];
          if (size < box[4]-box[1]) size = box[4]-box[1];
          if (size < box[5]-box[2]) size = box[5]-box[2];
        }
        n++;
      }
  }
  if (n <= 1) {
    AIM_ERROR(aimInfo, "AIM is given %d Bodies with %s!",
              n, nameVal->vals.string);
    status = CAPS_SOURCEERR;
    goto cleanup;
  }

#ifdef DEBUG
  printf(" Reference size = %le\n", size);
#endif
  /* set tess parameters */
  status = aim_getValue(aimInfo, Tess_Params, ANALYSISIN, &value);
  AIM_STATUS(aimInfo, status, "aim_getValue on Tess_Params!");
  params[0] = value->vals.reals[0]*size;
  params[1] = value->vals.reals[1]*size;
  params[2] = value->vals.reals[2];
  
  /* get OML state */
  status = aim_getValue(aimInfo, OML, ANALYSISIN, &value);
  AIM_STATUS(aimInfo, status, "aim_getValue on OML!");
  
  /* allocate an array to store our per Body informtion */
  cData = (clouds *) EG_alloc(n*sizeof(clouds));
  if (cData == NULL) {
    AIM_ERROR(aimInfo, " Malloc Error on %d cData!\n", n);
    status = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i < n; i++) {
    cData[i].name  = NULL;
    cData[i].oml   = 0;
  }
  if (value->vals.integer == True) {
    cData[ioml].oml = 1;
  } else {
    ioml = -1;
  }
 
  /* set up structures for interference checking */
  for (n = i = 0; i < numBody; i++) {
    status = EG_attributeRet(bodies[i], nameVal->vals.string,
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS)
      if (atype == ATTRSTRING) {
        cData[n].name   = string;
        cData[n].bIndex = i;
#ifdef DEBUG
        if (cData[n].oml == 1) {
          printf("Initialize for %s -- OML\n", cData[n].name);
        } else {
          printf("Initialize for %s\n", cData[n].name);
        }
#endif
        j = 0;
        if (cData[n].oml == 1) j = 1;
        status = initializeCloud(bodies[i], params, j, &cData[n].cloud);
        AIM_STATUS(aimInfo, status, "initializeClouds (%d) for Body %d!",
                   j, i+1);
        n++;
      }
  }
  
  /* allocte the memory for the "distance" matrix */
  dist = (double *) EG_alloc(n*n*sizeof(double));
  if (dist == NULL) {
    AIM_ERROR(aimInfo, " Malloc Error on %d^2 matrix!\n", n);
    for (j = 0; j < n; j++) freeCloud(&cData[j].cloud);
    status = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i < n*n; i++) dist[i] = 0.0;
  if (ioml != -1) dist[ioml*n+ioml] = 1.0;
  
  /* look at interference for all inner bodies */
  for (i = 0; i < n-1; i++) {
    if (cData[i].oml == 1) continue;
    for (j = i+1; j < n; j++) {
      if (cData[j].oml == 1) continue;
      
      /* 1) classify the pair of Bodies */
      status = classifyClouds(&cData[i].cloud, &cData[j].cloud, &pair);
#ifdef DEBUG
      type = 0;
      if (status == EGADS_SUCCESS) type = pair.type;
      printf(" classify = %d  type = %d  %s  %s\n",
             status, type, cData[i].name, cData[j].name);
#endif
      if (status != EGADS_SUCCESS) continue;
      
      /* 2) find the distances */
      status = minimizeClouds(&pair);
      if (status != EGADS_SUCCESS) {
#ifdef DEBUG
        printf(" ERROR: minimize = %d!\n", status);
#endif
        freeCloudPair(&pair);
        continue;
      }
      
      /* 3) report the distance and endpoints of the line segment */
      status = endpointClouds(&pair, &v, xyzs, xyzt);
      if (status != EGADS_SUCCESS) {
#ifdef DEBUG
        printf(" ERROR: endpoint = %d!\n", status);
#endif
        freeCloudPair(&pair);
        continue;
      }
      m       = j*n + i;
      dist[m] = v;
      m       = i*n + j;
      dist[m] = v;
      
      /* 4) deallocate the memory used for the pair */
      freeCloudPair(&pair);
    }
  }
  
  /* look at interference with the OML */
  for (i = 0; i < n-1; i++)
    for (j = i+1; j < n; j++) {
      if ((cData[i].oml != 1) && (cData[j].oml != 1)) continue;
      
      /* 1) classify the pair of Bodies -- OML must be first */
      if (cData[i].oml == 1) {
        if (cData[j].oml == 2) continue;
        status = classifyClouds(&cData[i].cloud, &cData[j].cloud, &pair);
      } else {
        if (cData[i].oml == 2) continue;
        status = classifyClouds(&cData[j].cloud, &cData[i].cloud, &pair);
      }
#ifdef DEBUG
      type = 0;
      if (status == EGADS_SUCCESS) type = pair.type;
      printf(" outer classify = %d  type = %d  %s  %s\n",
             status, type, cData[i].name, cData[j].name);
#endif
      if (status != EGADS_SUCCESS) continue;
      
      /* 2) find the distances */
      status = minimizeClouds(&pair);
      if (status != EGADS_SUCCESS) {
#ifdef DEBUG
        printf(" ERROR: minimize = %d!\n", status);
#endif
        freeCloudPair(&pair);
        continue;
      }
      
      /* 3) report the distance and endpoints of the line segment
            Note: xyzs is from 1st, xyzt is from second arg to classify */
      status = endpointClouds(&pair, &v, xyzs, xyzt);
      if (status != EGADS_SUCCESS) {
#ifdef DEBUG
        printf(" ERROR: endpoint = %d!\n", status);
#endif
        freeCloudPair(&pair);
        continue;
      }
      m       = j*n + i;
      dist[m] = v;
      m       = i*n + j;
      dist[m] = v;
      
      /* 4) deallocate the memory used for the pair */
      freeCloudPair(&pair);
    }

  /* output the file */
  fp = aim_fopen(aimInfo, "interference.dat", "wb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file for writing!\n");
    status = CAPS_DIRERR;
    EG_free(dist);
    for (j = 0; j < n; j++) freeCloud(&cData[j].cloud);
    goto cleanup;
  }
  
  status = CAPS_IOERR;
  len = fwrite(&n, sizeof(int), 1, fp);
  if (len != 1) goto cleanup;

  for (i = 0; i < n; i++) {
    AIM_NOTNULL(cData[i].name, aimInfo, status);
    j   = strlen(cData[i].name) + 1;
    len = fwrite(&j, sizeof(int), 1, fp);
    if (len != 1) goto cleanup;
    len = fwrite(cData[i].name, sizeof(char), j, fp);
    if (len != j) goto cleanup;
    status = EG_tessMassProps(cData[i].cloud.tess, props);
    AIM_STATUS(aimInfo, status, "EG_tessMassProps on Body %d!", i+1);
    len = fwrite(props, sizeof(double), 14, fp);
    if (len != 14) goto cleanup;
  }
  len = fwrite(dist, sizeof(double), n*n, fp);
  if (len != n*n) goto cleanup;

  fclose(fp);
  fp     = NULL;
  status = CAPS_SUCCESS;

cleanup:
  if (fp != NULL) fclose(fp);
  if (cData != NULL) {
    for (j = 0; j < n; j++) freeCloud(&cData[j].cloud);
    EG_free(cData);
  }
  if (dist != NULL) EG_free(dist);
  
  return status;
}


int aimPostAnalysis(/*@unused@*/ void *instStore, void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
#ifdef DEBUG
  printf(" interferenceAIM/aimPostAnalysis   instance = %d!\n",
         aim_getInstance(aimInfo));
#endif
  
  return aim_isFile(aimInfo, "interference.dat");
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{
  /*! \page aimOutputsInter AIM Outputs
   * The following list outlines the Interference outputs available through
   * the AIM interface.
   */

#ifdef DEBUG
  printf(" interferenceAIM/aimOutputs  instance = %d  index = %d!\n",
         aim_getInstance(aimInfo), index);
#endif

  if (index == Names) {
    *aoname = EG_strdup("Names");
    form->type        = String;
    form->dim         = Vector;
    form->nrow        = 1;
    form->ncol        = 1;
    form->units       = NULL;
    form->lfixed      = Change;
    form->nullVal     = IsNull;
    form->vals.string = NULL;

    /*! \page aimOutputsInter
     * - <B> Names </B> = A list of Attr_Name "value"s indicating the order
     * of the Bodies found in the rest of the outputs.
     */
    
  } else if (index == Distances) {
    *aoname = EG_strdup("Distances");
    form->type       = Double;
    form->dim        = Array2D;
    form->nrow       = 1;
    form->ncol       = 1;
    form->units      = NULL;
    form->lfixed     = Change;
    form->nullVal    = IsNull;
    form->vals.reals = NULL;

    /*! \page aimOutputsInter
     * - <B> Distances </B> = Distances. A symmetric NxN double array of
     * returned distances.  If it exists the OML can be found as a non-zero
     * diagonal entry. All other diagonal entries are zero.
     */
    
  } else if (index == Volumes) {
    *aoname = EG_strdup("Volumes");
    form->type       = Double;
    form->dim        = Vector;
    form->nrow       = 1;
    form->ncol       = 1;
    form->units      = NULL;
    form->lfixed     = Change;
    form->nullVal    = IsNull;
    form->vals.reals = NULL;

    /*! \page aimOutputsInter
     * - <B> Volumes </B> = Volume.
     */

  } else if (index == Areas) {
    *aoname = EG_strdup("Areas");
    form->type       = Double;
    form->dim        = Vector;
    form->nrow       = 1;
    form->ncol       = 1;
    form->units      = NULL;
    form->lfixed     = Change;
    form->nullVal    = IsNull;
    form->vals.reals = NULL;

    /*! \page aimOutputsInter
     * - <B> Areas </B> = Surface Area.
     */
    
  } else if (index == CGs) {
    *aoname = EG_strdup("CGs");
    form->type       = Double;
    form->dim        = Array2D;
    form->nrow       = 1;
    form->ncol       = 1;
    form->units      = NULL;
    form->lfixed     = Change;
    form->nullVal    = IsNull;
    form->vals.reals = NULL;

    /*! \page aimOutputsInter
     * - <B> qc/2V </B> = Center of Gravity (3 in length).
     */
    
  } else if (index == Inertias) {
    *aoname = EG_strdup("Inertias");
    form->type       = Double;
    form->dim        = Array2D;
    form->nrow       = 1;
    form->ncol       = 1;
    form->units      = NULL;
    form->lfixed     = Change;
    form->nullVal    = IsNull;
    form->vals.reals = NULL;

    /*! \page aimOutputsAInter
     * - <B> Inertias </B> = Inertial Matrix (9 in length).
     */

  } else {

    return CAPS_NOTFOUND;
    
  }

  return CAPS_SUCCESS;
}


int aimCalcOutput(/*@unused@*/ void *instStore, void *aimInfo, int index,
                  capsValue *val)
{

  int    status, i, j, k, m, n;
  double props[14];
  size_t lenx;
  FILE   *fp;

#ifdef DEBUG
  const char *name;

  status = aim_getName(aimInfo, index, ANALYSISOUT, &name);
  printf(" interferenceAIM/aimCalcOutput index %d (%s) = %d!\n",
         index, name, status);
#endif
  
  fp = aim_fopen(aimInfo, "interference.dat", "rb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file for reading!\n");
    return CAPS_DIRERR;
  }
  
  status = CAPS_IOERR;
  lenx   = fread(&n, sizeof(int), 1, fp);
  if (lenx != 1) goto cleanup;

  switch (index) {
  case Names:
      m = 1;
      for (i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        m += j;
        k = fseek(fp,  j*sizeof(char),   SEEK_CUR);
        if (k != 0) goto cleanup;
        k = fseek(fp, 14*sizeof(double), SEEK_CUR);
        if (k != 0) goto cleanup;
      }
      val->length      = n;
      val->nrow        = 1;
      val->ncol        = n;
      val->vals.string = (char *) EG_alloc(m*sizeof(char));
      if (val->vals.string == NULL) {
        AIM_ERROR(aimInfo, " Malloc Error Names (%d)!\n", m);
        status = EGADS_MALLOC;
        goto cleanup;
      }
      val->nullVal = NotNull;
      for (i = 0; i < m; i++) val->vals.string[i] = 0;
      rewind(fp);
      lenx = fread(&n, sizeof(int), 1, fp);
      if (lenx != 1) goto cleanup;
      for (m = i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        lenx = fread(&val->vals.string[m], sizeof(char), j, fp);
        if (lenx != j) goto cleanup;
        m += j;
        k = fseek(fp, 14*sizeof(double), SEEK_CUR);
        if (k != 0) goto cleanup;
      }
      break;
      
  case Distances:
      val->length  = n*n;
      val->nrow    = val->ncol = n;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
      val->nullVal = NotNull;
      status       = CAPS_IOERR;
      for (i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        k = fseek(fp,  j*sizeof(char),   SEEK_CUR);
        if (k != 0) goto cleanup;
        k = fseek(fp, 14*sizeof(double), SEEK_CUR);
        if (k != 0) goto cleanup;
      }
      lenx = fread(val->vals.reals, sizeof(double), n*n, fp);
      if (lenx != n*n) goto cleanup;
      break;
      
  case Volumes:
      val->length  = n;
      val->nrow    = 1;
      val->ncol    = n;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
      val->nullVal = NotNull;
      status       = CAPS_IOERR;
      for (i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        k = fseek(fp,  j*sizeof(char), SEEK_CUR);
        if (k != 0) goto cleanup;
        lenx = fread(props, sizeof(double), 14, fp);
        if (lenx != 14) goto cleanup;
        val->vals.reals[i] = props[0];
      }
      break;
      
  case Areas:
      val->length  = n;
      val->nrow    = 1;
      val->ncol    = n;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
      val->nullVal = NotNull;
      status       = CAPS_IOERR;
      for (i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        k = fseek(fp,  j*sizeof(char), SEEK_CUR);
        if (k != 0) goto cleanup;
        lenx = fread(props, sizeof(double), 14, fp);
        if (lenx != 14) goto cleanup;
        val->vals.reals[i] = props[1];
      }
      break;
      
  case CGs:
      val->length  = 3*n;
      val->nrow    = 3;
      val->ncol    = n;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
      val->nullVal = NotNull;
      status       = CAPS_IOERR;
      for (i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        k = fseek(fp,  j*sizeof(char), SEEK_CUR);
        if (k != 0) goto cleanup;
        lenx = fread(props, sizeof(double), 14, fp);
        if (lenx != 14) goto cleanup;
        val->vals.reals[3*i  ] = props[2];
        val->vals.reals[3*i+1] = props[3];
        val->vals.reals[3*i+2] = props[4];
      }
      break;
      
  case Inertias:
      val->length  = 9*n;
      val->nrow    = 9;
      val->ncol    = n;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
      val->nullVal = NotNull;
      status       = CAPS_IOERR;
      for (i = 0; i < n; i++) {
        lenx = fread(&j, sizeof(int), 1, fp);
        if (lenx != 1) goto cleanup;
        k = fseek(fp,  j*sizeof(char), SEEK_CUR);
        if (k != 0) goto cleanup;
        lenx = fread(props, sizeof(double), 14, fp);
        if (lenx != 14) goto cleanup;
        for (j = 0; j < 9; j++) val->vals.reals[9*i+j] = props[5+j];
      }
      break;
  }

  status = CAPS_SUCCESS;

cleanup:
  fclose(fp);
  return status;
}


void aimCleanup(/*@unused@*/ void *instStore)
{

}
