#include "caps.h"

static int setValueStr(capsObj aobj, const char *name, const char *string)
{
  int            status, nErr;
  capsObj        vobj;
  capsErrs       *errors;
  
  status = caps_childByName(aobj, VALUE, ANALYSISIN, name, &vobj, &nErr, &errors);
  if (status != CAPS_SUCCESS) return status;
  
  return caps_setValue(vobj, String, 1, 1, string, NULL, NULL, &nErr, &errors);
}


static int setValueDbl(capsObj aobj, const char *name, double real)
{
  int            status, nErr;
  capsObj        vobj;
  capsErrs       *errors;
  
  status = caps_childByName(aobj, VALUE, ANALYSISIN, name, &vobj, &nErr, &errors);
  if (status != CAPS_SUCCESS) return status;
  
  return caps_setValue(vobj, Double, 1, 1, &real, NULL, NULL, &nErr, &errors);
}


int aflr4egads(ego model, ego *tesses)
{
  char           *name, *units;
  int            i, j, oclass, mtype, status, exec, nbody, nface, nErr, *senses;
  double         size, bbox[6];
  ego            geom, tess, *bodies, *faces;
  capsObj        pobj, aobj, link, parent;
  capsErrs       *errors;
  capsOwn        current;
  enum capsoType type;
  enum capssType subtype;
  
  status = EG_getBoundingBox(model, bbox);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getBoundingBox = %d\n\n", status);
    return status;
  }
                              size = bbox[3]-bbox[0];
  if (bbox[4]-bbox[1] > size) size = bbox[4]-bbox[1];
  if (bbox[5]-bbox[2] > size) size = bbox[5]-bbox[2];
  
  /* mark all bodies */
  status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &senses);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n\n", status);
    return status;
  }
  printf(" Number of Bodies = %d\n\n", nbody);
  for (i = 0; i < nbody; i++) {
    status = EG_attributeAdd(bodies[i], "capsAIM", ATTRSTRING, 9, NULL, NULL,
                             "aflr4AIM");
    if (status != EGADS_SUCCESS)
      printf(" EG_attributeAdd capsAim %d = %d\n", i, status);
    status = EG_attributeAdd(bodies[i], "capsMeshLength", ATTRREAL, 1, NULL,
                             &size, NULL);
    if (status != EGADS_SUCCESS)
      printf(" EG_attributeAdd capsMeshLength %d = %d\n", i, status);
    status = EG_getBodyTopos(bodies[i], NULL, FACE, &nface, &faces);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos %d = %d\n\n", i, status);
      return status;
    } else {
      for (j = 0; j < nface; j++) {
        status = EG_attributeAdd(faces[j], "capsGroup", ATTRSTRING, 6, NULL,
                                 NULL, "Body");
        if (status != EGADS_SUCCESS)
          printf(" EG_attributeAdd capsGroup %d  %d = %d\n", i, j, status);
      }
      EG_free(faces);
    }
  }
  
  /* open up CAPS */
  status = caps_open("cgtTest", NULL, oMODL, model, 1, &pobj, &nErr, &errors);
  if (status != CAPS_SUCCESS) {
    printf(" caps_start = %d\n\n", status);
    return status;
  }
  
  status = caps_makeAnalysis(pobj, "aflr4AIM", "aflr4", NULL, NULL, &exec, &aobj,
                             &nErr, &errors);
  if (status != CAPS_SUCCESS) {
    printf(" caps_makeAnalysis = %d\n\n", status);
    caps_close(pobj, 0, NULL);
    return status;
  }
  
  /* set the inputs */
  status = setValueStr(aobj, "Proj_Name",          "AFLR4surfaceMeshing");
  if (status != CAPS_SUCCESS)
    printf(" setValueStr = %d for Proj_Name!\n", status);
  status = setValueStr(aobj, "Mesh_Format",        "NoOutput");
  if (status != CAPS_SUCCESS)
    printf(" setValueStr = %d for Mesh_Format!\n", status);
  status = setValueDbl(aobj, "Mesh_Length_Factor", 1.0);
  if (status != CAPS_SUCCESS)
    printf(" setValueDbl = %d for Mesh_Length_Factor!\n", status);
  status = setValueDbl(aobj, "max_scale",          0.0125);
  if (status != CAPS_SUCCESS)
    printf(" setValueDbl = %d for max_scale!\n", status);
  status = setValueDbl(aobj, "min_scale",          0.0001);
  if (status != CAPS_SUCCESS)
    printf(" setValueDbl = %d for min_scale!\n", status);
  status = setValueDbl(aobj, "curv_factor",        1.0);
  if (status != CAPS_SUCCESS)
    printf(" setValueDbl = %d for curv_factor!\n", status);
  status = setValueDbl(aobj, "erw_all",            1.0);
  if (status != CAPS_SUCCESS)
    printf(" setValueDbl = %d for erw_all!\n", status);
  
  status = caps_preAnalysis(aobj, &nErr, &errors);
  if (status != CAPS_SUCCESS) {
    printf(" caps_preAnalysis = %d\n\n", status);
    caps_close(pobj, 0, NULL);
    return status;
  }
  /* set current */
  status = caps_info(pobj, &name, &type, &subtype, &link, &parent, &current);
  if (status != CAPS_SUCCESS) {
    printf(" caps_info on Problem = %d\n", status);
    caps_close(pobj, 0, NULL);
    return status;
  }
  status = caps_postAnalysis(aobj, &nErr, &errors);
  if (status != CAPS_SUCCESS) {
    printf(" caps_postAnalysis = %d\n", status);
    caps_close(pobj, 0, NULL);
    return status;
  }
  
  /* get all of the tessellation objects */
  for (i = 1; i <= nbody; i++) {
    status = caps_bodyByIndex(aobj, -i, &tess, &units);
    if (status != CAPS_SUCCESS) {
      printf(" caps_bodyByIndex %d = %d\n", i, status);
      caps_close(pobj, 0, NULL);
      return status;
    }
    printf(" oclass = %d  mtype = %d\n", tess->oclass, tess->mtype);
    /* copy tessellation object */
    status = EG_copyObject(tess, NULL, &tesses[i-1]);
    if (status != CAPS_SUCCESS) {
      printf(" EG_copyObject %d = %d\n", i, status);
      caps_close(pobj, 0, NULL);
      return status;
    }
  }
  
  status = caps_close(pobj, 0, NULL);
  if (status != CAPS_SUCCESS)
    printf(" caps_close = %d\n", status);
  
  return EGADS_SUCCESS;
}
