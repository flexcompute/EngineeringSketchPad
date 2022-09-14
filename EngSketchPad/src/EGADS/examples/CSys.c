/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Coordinate System test
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include "egads.h"

#define PI               3.1415926535897931159979635


int main(int argc, char *argv[])
{
  int    stat, n, aType, aLen, nobj, oclass, mtype, *senses;
  double cosz, sinz, angle, data[9], matrix[12], uvbox[4];
  ego    context, body1, body2, oform, ref, *objs, *children;
  const int    *ints;
  const char   *str;
  const double *reals, *csys;

  angle = 30.0;
  cosz = cos(angle*PI/180.0);
  sinz = sin(angle*PI/180.0);
  
  matrix[ 0] = cosz; matrix[ 1] =-sinz; matrix[ 2] = 0.0; matrix[ 3] = 0.0;
  matrix[ 4] = sinz; matrix[ 5] = cosz; matrix[ 6] = 0.0; matrix[ 7] = 0.0;
  matrix[ 8] = 0;    matrix[ 9] = 0;    matrix[10] = 1.0; matrix[11] = 0.0;

  /* create an EGADS context */
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open return = %d\n", stat);
    return 1;
  }
  stat = EG_makeTransform(context, matrix, &oform);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTransform = %d\n", stat);
    return 1;
  }
  
  /* make a cylinder */
  data[0] = data[2] = data[3] = data[5] = 0.0;
  data[1] = -2.0;
  data[4] =  2.0;
  data[6] =  0.5;
  stat = EG_makeSolidBody(context, CYLINDER, data, &body1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody return = %d\n", stat);
    return 1;
  }
  stat = EG_getBodyTopos(body1, NULL, NODE, &nobj, &objs);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getBodyTopos return = %d\n", stat);
    return 1;
  }
  stat = EG_getTopology(objs[0], &ref, &oclass, &mtype, data, &n, &children,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology Node return = %d\n", stat);
    return 1;
  }

  /* make a CSys on body and then Node */
  data[3] = 1.0; data[4] = 0.0; data[5] = 0.0;
  data[6] = 0.0; data[7] = 1.0; data[8] = 0.0;
  stat = EG_attributeAdd(body1, "bodyCSys", ATTRCSYS, 9, NULL, data, NULL);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeAdd body return = %d\n", stat);
    return 1;
  }
  data[0] = 1.0; data[1] = 0.0; data[2] = 0.0;
  data[3] = 0.0; data[4] = 1.0; data[5] = 0.0;
  stat = EG_attributeAdd(objs[0], "nodeCSys", ATTRCSYS, 6, NULL, data, NULL);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeAdd node return = %d\n", stat);
    return 1;
  }
  EG_free(objs);
  stat = EG_attributeRet(body1, "bodyCSys", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet body return = %d\n", stat);
    return 1;
  }
  csys = &reals[aLen];
  printf(" aLen = %d  aType = %d -- bodyCSys\n", aLen, aType);
  printf(" data = %lf %lf %lf\n", reals[0], reals[1], reals[2]);
  printf("        %lf %lf %lf\n", reals[3], reals[4], reals[5]);
  printf("        %lf %lf %lf\n", reals[6], reals[7], reals[8]);
  printf(" CSys = %lf %lf %lf\n", csys[0],  csys[ 1], csys[ 2]);
  printf("        %lf %lf %lf\n", csys[3],  csys[ 4], csys[ 5]);
  printf("        %lf %lf %lf\n", csys[6],  csys[ 7], csys[ 8]);
  printf("        %lf %lf %lf\n", csys[9],  csys[10], csys[11]);

  /* get a Face */
  stat = EG_getBodyTopos(body1, NULL, FACE, &nobj, &objs);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getBodyTopos Face return = %d\n", stat);
    return 1;
  }
  stat = EG_getTopology(objs[0], &ref, &oclass, &mtype, uvbox, &n, &children,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology Node return = %d\n", stat);
    return 1;
  }
  data[0] = 0.5*(uvbox[0] + uvbox[1]);
  data[1] = 0.5*(uvbox[2] + uvbox[3]);
  data[2] = 2.0;
  stat = EG_attributeAdd(objs[0], "faceCSys", ATTRCSYS, 3, NULL, data, NULL);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeAdd face return = %d\n", stat);
    return 1;
  }
  stat = EG_attributeRet(objs[0], "faceCSys", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet body return = %d\n", stat);
    return 1;
  }
  csys = &reals[aLen];
  printf(" aLen = %d  aType = %d -- faceCSys\n", aLen, aType);
  printf(" data = %lf %lf %lf\n", reals[0], reals[1], reals[2]);
  printf(" CSys = %lf %lf %lf\n", csys[0],  csys[ 1], csys[ 2]);
  printf("        %lf %lf %lf\n", csys[3],  csys[ 4], csys[ 5]);
  printf("        %lf %lf %lf\n", csys[6],  csys[ 7], csys[ 8]);
  printf("        %lf %lf %lf\n", csys[9],  csys[10], csys[11]);
  printf("\n");
  EG_free(objs);
  
  /* transform the object */
  stat = EG_copyObject(body1, oform, &body2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_copyObject return = %d\n", stat);
    return 1;
  }
  
  stat = EG_attributeRet(body2, "bodyCSys", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet body return = %d\n", stat);
    return 1;
  }
  csys = &reals[aLen];
  printf(" aLen = %d  aType = %d -- bodyCSys\n", aLen, aType);
  printf(" data = %lf %lf %lf\n", reals[0], reals[1], reals[2]);
  printf("        %lf %lf %lf\n", reals[3], reals[4], reals[5]);
  printf("        %lf %lf %lf\n", reals[6], reals[7], reals[8]);
  printf(" CSys = %lf %lf %lf\n", csys[0],  csys[ 1], csys[ 2]);
  printf("        %lf %lf %lf\n", csys[3],  csys[ 4], csys[ 5]);
  printf("        %lf %lf %lf\n", csys[6],  csys[ 7], csys[ 8]);
  printf("        %lf %lf %lf\n", csys[9],  csys[10], csys[11]);

  stat = EG_getBodyTopos(body2, NULL, NODE, &nobj, &objs);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getBodyTopos return = %d\n", stat);
    return 1;
  }
  stat = EG_attributeRet(objs[0], "nodeCSys", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet node return = %d\n", stat);
    return 1;
  }
  csys = &reals[aLen];
  printf(" aLen = %d  aType = %d -- nodeCSys\n", aLen, aType);
  printf(" data = %lf %lf %lf\n", reals[0], reals[1], reals[2]);
  printf("        %lf %lf %lf\n", reals[3], reals[4], reals[5]);
  printf(" CSys = %lf %lf %lf\n", csys[0],  csys[ 1], csys[ 2]);
  printf("        %lf %lf %lf\n", csys[3],  csys[ 4], csys[ 5]);
  printf("        %lf %lf %lf\n", csys[6],  csys[ 7], csys[ 8]);
  printf("        %lf %lf %lf\n", csys[9],  csys[10], csys[11]);
  EG_free(objs);
  
  stat = EG_getBodyTopos(body2, NULL, FACE, &nobj, &objs);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getBodyTopos return = %d\n", stat);
    return 1;
  }
  stat = EG_attributeRet(objs[0], "faceCSys", &aType, &aLen, &ints, &reals, &str);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_attributeRet face return = %d\n", stat);
    return 1;
  }
  csys = &reals[aLen];
  printf(" aLen = %d  aType = %d -- faceCSys\n", aLen, aType);
  printf(" data = %lf %lf %lf\n", reals[0], reals[1], reals[2]);
  printf(" CSys = %lf %lf %lf\n", csys[0],  csys[ 1], csys[ 2]);
  printf("        %lf %lf %lf\n", csys[3],  csys[ 4], csys[ 5]);
  printf("        %lf %lf %lf\n", csys[6],  csys[ 7], csys[ 8]);
  printf("        %lf %lf %lf\n", csys[9],  csys[10], csys[11]);
  EG_free(objs);
  
  printf(" EG_deleteObject oform = %d\n", EG_deleteObject(oform));
  printf(" EG_deleteObject body1 = %d\n", EG_deleteObject(body1));
  printf(" EG_deleteObject body2 = %d\n", EG_deleteObject(body2));
  printf(" EG_close the context  = %d\n", EG_close(context));
  
  return 0;
}
