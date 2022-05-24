/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             An example of "bottom-up" construction
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"


int main(/*@unused@*/int argc, /*@unused@*/ char *argv[])
{
  int    senses[3];
  double xyz[3], data[6], dum[3], range[2];
  ego    context, nodes[3], lines[3], edges[3], objs[2], loop, face, body, other;
  ego    model = NULL, copyModel = NULL;

  printf(" EG_open            = %d\n", EG_open(&context));

  /* make the Nodes */
  xyz[0] = xyz[1] = xyz[2] = 0.0;
  printf(" EG_makeTopology N0 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[0]));
  xyz[0] = 1.0;
  printf(" EG_makeTopology N1 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[1]));
  xyz[0] = 0.0;
  xyz[1] = 2.0;
  printf(" EG_makeTopology N2 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[2]));

  /* make the Curves */
  data[0] = data[1] = data[2] = data[4] = data[5] = 0.0;
  data[3] = 1.0;
  printf(" EG_makeGeometry L0 = %d\n", EG_makeGeometry(context, CURVE, LINE,
                                                       NULL, NULL, data,
                                                       &lines[0]));
  data[3] = 0.0;
  data[4] = 2.0;
  printf(" EG_makeGeometry L1 = %d\n", EG_makeGeometry(context, CURVE, LINE,
                                                       NULL, NULL, data,
                                                       &lines[1]));
  data[0] =  1.0;
  data[3] = -1.0;
  printf(" EG_makeGeometry L2 = %d\n", EG_makeGeometry(context, CURVE, LINE,
                                                       NULL, NULL, data,
                                                       &lines[2]));

  /* construct the Edges */
  xyz[0] = xyz[1] = xyz[2] = 0.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[0], xyz, &range[0],
                                                     dum));
  xyz[0] = 1.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[0], xyz, &range[1],
                                                     dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[0];
  objs[1] = nodes[1];
  printf(" EG_makeTopology E0 = %d\n", EG_makeTopology(context, lines[0], EDGE,
                                                       TWONODE, range, 2, 
                                                       objs, NULL, &edges[0]));

  xyz[0] = 0.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[1], xyz, &range[0],
                                                     dum));
  xyz[1] = 2.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[1], xyz, &range[1],
                                                     dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[0];
  objs[1] = nodes[2];
  printf(" EG_makeTopology E1 = %d\n", EG_makeTopology(context, lines[1], EDGE,
                                                       TWONODE, range, 2, 
                                                       objs, NULL, &edges[1]));

  xyz[0] = 1.0;
  xyz[1] = 0.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[2], xyz, &range[0],
                                                     dum));
  xyz[0] = 0.0;
  xyz[1] = 2.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[2], xyz, &range[1],
                                                     dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[1];
  objs[1] = nodes[2];
  printf(" EG_makeTopology E2 = %d\n", EG_makeTopology(context, lines[2], EDGE,
                                                       TWONODE, range, 2, 
                                                       objs, NULL, &edges[2]));

  /* make the Loop and then the Face */
  senses[0] = -1;
  senses[1] =  1;
  senses[2] = -1;
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, LOOP,
                                                       CLOSED, NULL, 3, edges,
                                                       senses, &loop));
  
#ifdef WIREBOD
  printf(" EG_makeTopology B  = %d\n", EG_makeTopology(context, NULL, BODY,
                                                       WIREBODY, NULL, 1, &loop,
                                                       NULL, &body));
  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                       NULL, 1, &body,
                                                       NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "wirebody.egads"));
#else
  printf(" EG_makeFace        = %d\n", EG_makeFace(loop, SREVERSE, NULL,
                                                   &face));

  /* revolve */
  data[0] = data[1] = data[2] = data[3] = data[5] = 0.0;
  data[4] = 1.0;
  printf(" EG_rotate          = %d\n", EG_rotate(face, 180.0, data, &body));
  
  /* make a copy */
  printf(" EG_copyObject      = %d\n", EG_copyObject(body, NULL, &model));
  printf(" EG_isEquivalent    = %d\n", EG_isEquivalent(body, model));
  printf(" EG_deleteObject    = %d\n", EG_deleteObject(model));
  printf(" EG_attributeAdd    = %d\n", EG_attributeAdd(body, "myAttribute",
                                                       ATTRSTRING, 7, NULL, NULL,
                                                       "myValue"));
  
  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                       NULL, 1, &body,
                                                       NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "mkCone.egads"));
  printf("\n");
#endif
  
  printf(" EG_open  Other     = %d\n", EG_open(&other));
  printf(" EG_copyObject      = %d\n", EG_copyObject(model, other, &copyModel));
  
  EG_setOutLevel(context, 2);
  printf(" EG_close           = %d\n", EG_close(context));
  
  printf(" EG_saveModel Other = %d\n", EG_saveModel(copyModel, "mkCono.egads"));
  if (copyModel != NULL) EG_deleteObject(copyModel);
  printf(" EG_close Other     = %d\n", EG_close(other));

  return 0;
}
