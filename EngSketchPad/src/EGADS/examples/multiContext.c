/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             An example of multithread/context execution
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"
#include "emp.h"


static ego context, bodies[2];


static void
otherThread(/*@unused@*/ void *struc)
{
  int    stat, oclass, mtype, len, *senses;
  double data[7];
  ego    other, sbo, bodyb, bodyc, geom, *objs;

  /* open an additional context for this thread */
  stat = EG_open(&other);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open2 return = %d\n", stat);
    EMP_ThreadExit();
  }
  
  /* do an SBO */
  data[0] = data[1] = data[2] = -1.0;
  data[3] = data[4] = data[5] =  2.0;
  stat = EG_makeSolidBody(other, BOX, data, &bodyb);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody box2 return = %d\n", stat);
    EMP_ThreadExit();
  }
  data[0] = data[2] = data[3] = data[5] = 0.0;
  data[1] = -2.0;
  data[4] =  2.0;
  data[6] =  0.5;
  stat = EG_makeSolidBody(other, CYLINDER, data, &bodyc);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody cyl2 return = %d\n", stat);
    EMP_ThreadExit();
  }
  stat = EG_solidBoolean(bodyb, bodyc, FUSION, &sbo);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_solidBoolean 2     return = %d\n", stat);
    EMP_ThreadExit();
  }
  stat = EG_getTopology(sbo, &geom, &oclass, &mtype, NULL, &len,
                        &objs, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology  2     return = %d\n", stat);
    EMP_ThreadExit();
  }
  stat = EG_attributeAdd(objs[0], "From", ATTRSTRING, 6, NULL, NULL, "Thread");
  if (stat != EGADS_SUCCESS)
    printf(" EG_attributeAdd 2     return = %d\n", stat);

  /* perform the cross context copy */
  stat = EG_copyObject(objs[0], context, &bodies[1]);
  if (stat != EGADS_SUCCESS)
    printf(" EG_copyObject 2       return = %d\n", stat);

  /* cleanup and get out */
  EG_deleteObject(sbo);
  EG_close(other);
  EMP_ThreadExit();
}


int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  int    i, stat, oclass, mtype, len, *senses;
  double data[7];
  void   *thread;
  ego    bodyb, bodyc, model, sbo, geom, *objs;
  
  bodies[0] = bodies[1] = NULL;
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open return = %d\n", stat);
    return 1;
  }
  
  /* create our other thread */
  thread = EMP_ThreadCreate(otherThread, NULL);
  if (thread == NULL) {
    printf(" Error creating thread!\n");
    return 1;
  }
  
  /* do our SBO */
  data[0] = data[1] = data[2] = -1.0;
  data[3] = data[4] = data[5] =  2.0;
  stat = EG_makeSolidBody(context, BOX, data, &bodyb);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody box1 return = %d\n", stat);
    return 1;
  }
  data[0] = data[2] = data[3] = data[5] = 0.0;
  data[1] = -2.0;
  data[4] =  2.0;
  data[6] =  0.5;
  stat = EG_makeSolidBody(context, CYLINDER, data, &bodyc);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody cyl1 return = %d\n", stat);
    return 1;
  }
  stat = EG_solidBoolean(bodyb, bodyc, SUBTRACTION, &sbo);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_solidBoolean 1     return = %d\n", stat);
    return 1;
  }
  stat = EG_getTopology(sbo, &geom, &oclass, &mtype, NULL, &len,
                        &objs, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology  1     return = %d\n", stat);
    return 1;
  }
  stat = EG_attributeAdd(objs[0], "From", ATTRSTRING, 4, NULL, NULL, "Main");
  if (stat != EGADS_SUCCESS)
    printf(" EG_attributeAdd 1     return = %d\n", stat);
  stat = EG_copyObject(objs[0], NULL, &bodies[0]);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_copyObject 1       return = %d\n", stat);
    return 1;
  }
  EG_deleteObject(sbo);
  
  /* wait for the thread to exit and destroy afterward */
  EMP_ThreadWait(thread);
  EMP_ThreadDestroy(thread);
  
  /* make the model of the 2 bodies and write it out */
  i = 1;
  if (bodies[1] != NULL) i = 2;
  printf(" Found %d Bodies!\n", i);
  stat = EG_makeTopology(context, NULL, MODEL, 0, NULL, i, bodies, NULL,
                         &model);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTopology       return = %d\n", stat);
    return 1;
  }
  stat = EG_saveModel(model, "multiContext.egads");
  if (stat != EGADS_SUCCESS)
    printf(" EG_saveModel          return = %d\n", stat);
  
  EG_close(context);

  return 0;
}
