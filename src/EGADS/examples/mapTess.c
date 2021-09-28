#include <math.h>
#include "egads.h"


#define SCONE
#define BOXCYL


extern int EG_sameBodyTopo(const egObject *bod1, const egObject *bod2);



int main(int argc, char *argv[])
{
  int          i, j, stat, nedge, nface, err, np1, np2, nt1, nt2;
  const int    *pt1, *pi1, *pt2, *pi2, *ts1, *tc1, *ts2, *tc2;
  double       data[7], params[3], dx[3];
  const double *x1, *x2, *p1, *p2;
  ego          context, body1, body2, tess1, tess2, *dum;
#ifdef BOXCYL
  int          oclass, mtype, *sens;
  ego          ref, mdl1, mdl2;
#else
  int          stype;
#endif

  dx[0] = 0.0;
  dx[1] = 0.0;
  dx[2] = 1.0;
  
  /* create an EGADS context */
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open return = %d\n", stat);
    return 1;
  }

#ifdef BOXCYL

  data[0] = data[1] = data[2] = -1.0;
  data[3] = data[4] = data[5] =  2.0;
  stat = EG_makeSolidBody(context, BOX, data, &tess1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody box1 return = %d\n", stat);
    return 1;
  }
  data[0] = data[2] = data[3] = data[5] = 0.0;
  data[1] = -2.0;
  data[4] =  2.0;
  data[6] =  0.5;
  stat = EG_makeSolidBody(context, CYLINDER, data, &tess2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody cyl1 return = %d\n", stat);
    return 1;
  }
  stat = EG_solidBoolean(tess1, tess2, SUBTRACTION, &mdl1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_solidBoolean 1     return = %d\n", stat);
    return 1;
  }
  printf(" EG_deleteObject tess1 = %d\n", EG_deleteObject(tess1));
  printf(" EG_deleteObject tess2 = %d\n", EG_deleteObject(tess2));
  stat = EG_getTopology(mdl1, &ref, &oclass, &mtype, data, &i, &dum, &sens);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology  1     return = %d\n", stat);
    return 1;
  }
  body1 = dum[0];

  data[0]  = data[1] = data[2] = -1.0;
  data[3]  = data[4] = data[5] =  2.0;
  data[0] += dx[0];
  data[1] += dx[1];
  data[2] += dx[2];
  stat = EG_makeSolidBody(context, BOX, data, &tess1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody box2 return = %d\n", stat);
    return 1;
  }
  data[0]  = data[2] = data[3] = data[5] = 0.0;
  data[1]  = -2.0;
  data[4]  =  2.0;
  data[6]  =  0.5;
  data[0] += dx[0];
  data[1] += dx[1];
  data[2] += dx[2];
  data[3] += dx[0];
  data[4] += dx[1];
  data[5] += dx[2];
  stat = EG_makeSolidBody(context, CYLINDER, data, &tess2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody cyl2 return = %d\n", stat);
    return 1;
  }
  stat = EG_solidBoolean(tess1, tess2, SUBTRACTION, &mdl2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_solidBoolean 2     return = %d\n", stat);
    return 1;
  }
  printf(" EG_deleteObject tess1 = %d\n", EG_deleteObject(tess1));
  printf(" EG_deleteObject tess2 = %d\n", EG_deleteObject(tess2));
  stat = EG_getTopology(mdl2, &ref, &oclass, &mtype, data, &i, &dum, &sens);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology  2     return = %d\n", stat);
    return 1;
  }
  body2 = dum[0];
  
#else
  
#ifdef SBOX
  stype   = BOX;
  data[0] = data[1] = data[2] = 0.0;
  data[3] = data[4] = data[5] = 1.0;
#endif
#ifdef SCYLINDER
  stype   = CYLINDER;
  data[0] = data[1] = data[2] = 0.0;
  data[3] = data[5] = 0.0;
  data[4] = 2.0;
  data[6] = 1.0;
#endif
#ifdef SCONE
  stype   = CONE;
  data[0] = data[1] = data[2] = 0.0;
  data[3] = data[5] = 0.0;
  data[4] = 2.0;
  data[6] = 1.0;
#endif
  stat = EG_makeSolidBody(context, stype, data, &body1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody bod1 return = %d\n", stat);
    return 1;
  }

#ifdef SBOX
  data[0] += dx[0];
  data[1] += dx[1];
  data[2] += dx[2];
#endif
#ifdef SCYLINDER
  data[0] += dx[0];
  data[1] += dx[1];
  data[2] += dx[2];
  data[3] += dx[0];
  data[4] += dx[1];
  data[5] += dx[2];
#endif
#ifdef SCONE
  data[0] += dx[0];
  data[1] += dx[1];
  data[2] += dx[2];
  data[3] += dx[0];
  data[4] += dx[1];
  data[5] += dx[2];
#endif
  stat = EG_makeSolidBody(context, stype, data, &body2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody bod2 return = %d\n", stat);
    return 1;
  }

#endif

  stat = EG_sameBodyTopo(body1, body2);
  printf(" EG_sameBodyTopo       = %d\n", stat);
  
  params[0] =  0.05;
  params[1] =  0.001;
  params[2] = 12.0;
  stat = EG_makeTessBody(body1, params, &tess1);
  printf(" EG_makeTessBody       = %d\n", stat);
  stat = EG_mapTessBody(tess1, body2, &tess2);
  printf(" EG_mapTessBody        = %d\n", stat);
  
  /* get the tessellations */
  stat = EG_getBodyTopos(body1, NULL, EDGE, &nedge, &dum);
  printf(" EG_getBodyTopos E     = %d\n", stat);
  EG_free(dum);
  
  stat = EG_getBodyTopos(body1, NULL, FACE, &nface, &dum);
  printf(" EG_getBodyTopos F     = %d\n", stat);
  EG_free(dum);
  printf("\n Number of Edges = %d   Number of Face = %d\n", nedge, nface);
 
  /* compare edges */
  err = 0;
  for (i = 1; i <= nedge; i++) {
    stat = EG_getTessEdge(tess1, i, &np1, &x1, &p1);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getTessEdge1 %d     = %d\n", i, stat);
      continue;
    }
    stat = EG_getTessEdge(tess2, i, &np2, &x2, &p2);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getTessEdge2 %d     = %d\n", i, stat);
      continue;
    }
    for (j = 0; j < np1; j++) {
      if ((fabs(x1[3*j  ]+dx[0]-x2[3*j  ]) < 1.e-14) &&
          (fabs(x1[3*j+1]+dx[1]-x2[3*j+1]) < 1.e-14) &&
          (fabs(x1[3*j+2]+dx[2]-x2[3*j+2]) < 1.e-14)) continue;
      printf(" Edge %d %d\n", i, j+1);
      err++;
    }
  }

  /* compare faces */
  if (err != 0) nface = 0;
  for (i = 1; i <= nface; i++) {
    stat = EG_getTessFace(tess1, i, &np1, &x1, &p1, &pt1, &pi1,
                                    &nt1, &ts1, &tc1);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getTessFace1 %d     = %d\n", i, stat);
      continue;
    }
    stat = EG_getTessFace(tess2, i, &np2, &x2, &p2, &pt2, &pi2,
                                    &nt2, &ts2, &tc2);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getTessFace2 %d     = %d\n", i, stat);
      continue;
    }
    printf(" Face %d: npts = %d %d,  ntris = %d %d\n", i, np1, np2, nt1, nt2);
    for (j = 0; j < np1; j++) {
      if ((fabs(x1[3*j  ]+dx[0]-x2[3*j  ]) < 1.e-14) &&
          (fabs(x1[3*j+1]+dx[1]-x2[3*j+1]) < 1.e-14) &&
          (fabs(x1[3*j+2]+dx[2]-x2[3*j+2]) < 1.e-14)) continue;
      printf(" Face %d %d  ptype = %d %d  pindex = %d %d\n",
             i, j+1, pt1[j], pt2[j], pi1[j], pi2[j]);
      err++;
    }
  }

  printf("\n");
  printf(" EG_deleteObject tess1 = %d\n", EG_deleteObject(tess1));
  printf(" EG_deleteObject tess2 = %d\n", EG_deleteObject(tess2));
#ifdef BOXCYL
  printf(" EG_deleteObject mdl1  = %d\n", EG_deleteObject(mdl1));
  printf(" EG_deleteObject mdl2  = %d\n", EG_deleteObject(mdl2));
#else
  printf(" EG_deleteObject body1 = %d\n", EG_deleteObject(body1));
  printf(" EG_deleteObject body2 = %d\n", EG_deleteObject(body2));
#endif
  printf(" EG_close the context  = %d\n", EG_close(context));

  return 0;
}
