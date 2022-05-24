/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             An example using the ruled function
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"


//#define OPENLOOP


int main(int argc, char *argv[])
{
  int    i, beg, end, status, senses[2], per, rev = 1;
  double xyz[3], data[10], range[2], xform[12];
  ego    context, nodes[3], curve, edges[2], objs[2], oform, secs[5];
  ego    others[2], body, model;
  static char *types[3] = {"Open", "Node", "Face"};

  if (argc != 3) {
    printf("\n Usage: ruled 0/1/2 0/1/2 -- 0-open, 1-node, 2-face\n\n");
    return 1;
  }
  sscanf(argv[1], "%d", &beg);
  sscanf(argv[2], "%d", &end);
  if (beg <   0) beg = 0;
  if (beg >   2) beg = 2;
  if (end <   0) end = 0;
  if (end >   2) end = 2;
  printf("\n Ruled with start %s and end %s\n\n", types[beg], types[end]);
  
  printf(" EG_open            = %d\n", EG_open(&context));

  /* make the Nodes */
  xyz[0] =  1.0;
  xyz[1] = xyz[2] = 0.0;
  printf(" EG_makeTopology N0 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[0]));
  xyz[0] = -1.0;
  printf(" EG_makeTopology N1 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[1]));
#ifdef OPENLOOP
  xyz[0] = 0.0;
  xyz[1] = 1.0;
  printf(" EG_makeTopology N2 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[2]));
#endif

  /* make the Curve */
  data[0] = data[1] = data[2] = data[4] = data[5] = data[6] = data[8] = 0.0;
  data[3] = data[7] = data[9] = 1.0;
  printf(" EG_makeGeometry C0 = %d\n", EG_makeGeometry(context, CURVE, CIRCLE,
                                                       NULL, NULL, data,
                                                       &curve));
  status = EG_getRange(curve, range, &per);
  printf(" EG_getRange     C0 = %d -  %lf %lf\n", status, range[0], range[1]);

  /* construct the Edges */
#ifdef OPENLOOP
  objs[0] = nodes[0];
  objs[1] = nodes[2];
  data[0] = range[0];
  data[1] = range[0] + 0.25*(range[1]-range[0]);
  printf(" EG_makeTopology E0 = %d\n", EG_makeTopology(context, curve, EDGE,
                                                       TWONODE, data, 2,
                                                       objs, NULL, &edges[0]));
  objs[0] = nodes[2];
  objs[1] = nodes[1];
  data[0] = data[1];
  data[1] = range[0] + 0.5*(range[1]-range[0]);
  printf(" EG_makeTopology E1 = %d\n", EG_makeTopology(context, curve, EDGE,
                                                       TWONODE, data, 2,
                                                       objs, NULL, &edges[1]));
  /* make the Loop */
  senses[0] = rev;
  senses[1] = rev;
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, LOOP,
                                                       OPEN, NULL, 2, edges,
                                                       senses, &secs[0]));
#else
  objs[0] = nodes[0];
  objs[1] = nodes[1];
  data[0] = range[0];
  data[1] = range[0] + 0.5*(range[1]-range[0]);
  printf(" EG_makeTopology E0 = %d\n", EG_makeTopology(context, curve, EDGE,
                                                       TWONODE, data, 2,
                                                       objs, NULL, &edges[0]));
  objs[0] = nodes[1];
  objs[1] = nodes[0];
  data[0] = data[1];
  data[1] = range[1];
  printf(" EG_makeTopology E1 = %d\n", EG_makeTopology(context, curve, EDGE,
                                                       TWONODE, data, 2,
                                                       objs, NULL, &edges[1]));
  /* make the Loop */
  senses[0] = rev;
  senses[1] = rev;
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, LOOP,
                                                       CLOSED, NULL, 2, edges,
                                                       senses, &secs[0]));
#endif

  /* make a transform */
  for (i = 1; i < 11; i++) xform[i] = 0.0;
  xform[0]  = xform[5] = xform[10]  = 1.1;
  xform[11] = 1.0;
  printf(" EG_makeTransform   = %d\n", EG_makeTransform(context, xform,
                                                        &oform));

  /* make the sections */
  for (i = 1; i < 5; i++)
    printf(" EG_copyObject %d    = %d\n", i, EG_copyObject( secs[i-1], oform,
                                                           &secs[i]));
  EG_deleteObject(oform);
  
  /* deal with the ends */
  others[0] = others[1] = NULL;
  if (beg == 1) {
    EG_deleteObject(secs[0]);
    xyz[0] = xyz[1] = xyz[2] = 0.0;
    printf(" EG_makeTopology Nb = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                         xyz, 0, NULL, NULL,
                                                         &secs[0]));
  } else if (beg == 2) {
    others[0] = secs[0];
    printf(" EG_makeFace beg    = %d\n", EG_makeFace(others[0], SREVERSE*rev,
                                                     NULL, &secs[0]));
  }
  if (end == 1) {
    EG_deleteObject(secs[4]);
    xyz[0] = xyz[1] = 0.0;
    xyz[2] = 5.0;
    printf(" EG_makeTopology Ne = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                         xyz, 0, NULL, NULL,
                                                         &secs[4]));
  } else if (end == 2) {
    others[1] = secs[4];
    printf(" EG_makeFace end    = %d\n", EG_makeFace(others[1], SFORWARD*rev,
                                                     NULL, &secs[4]));
  }
  
  /* blend & save */
  printf(" EG_ruled           = %d\n", EG_ruled(5, secs, &body));
  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                       NULL, 1, &body,
                                                       NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "ruled.egads"));
  printf("\n");

  /* cleanup */
  EG_deleteObject(model);
  for (i = 4; i >= 0; i--) EG_deleteObject(secs[i]);
  if (others[0] != NULL) EG_deleteObject(others[0]);
  if (others[1] != NULL) EG_deleteObject(others[1]);
  EG_deleteObject(edges[1]);
  EG_deleteObject(edges[0]);
  EG_deleteObject(curve);
  EG_deleteObject(nodes[1]);
  EG_deleteObject(nodes[0]);
#ifdef OPENLOOP
  EG_deleteObject(nodes[2]);
#endif

  EG_setOutLevel(context, 2);
  printf(" EG_close           = %d\n", EG_close(context));

  return 0;
}
