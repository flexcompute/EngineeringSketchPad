/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Test the non-planar "makeFace" options
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"

#define NEDGE 4


/* entry point not in egads.h */
extern int EG_isPlanar(ego object);


int main(int argc, char *argv[])
{
  int    i, senses[4], style;
  double xyz[3], data[10], dum[3], range[2], t[20];
  ego    context, nodes[4], lines[4], edges[8], objs[2], loop, surface;
  ego    face, body, model, *bedges = NULL;

  style = 0.0;
  if (argc > 2) {
    printf("\n Usage: makeFace3D [-1/0/1]\n\n");
    return 1;
  }
  if (argc == 2) {
    sscanf(argv[1], "%d", &style);
    printf(" style = %d\n\n", style);
  }
  
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
  xyz[0] = 1.0;
  printf(" EG_makeTopology N3 = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                       xyz, 0, NULL, NULL,
                                                       &nodes[3]));
  
  /* make the Curves */
  
  data[0] = 0.5;
  data[1] = data[2] = data[4] = data[5] = data[6] = data[7] = 0.0;
  data[3] = data[8] = -1.0;
  data[9] = 0.5;
  printf(" EG_makeGeometry C0 = %d\n", EG_makeGeometry(context, CURVE, CIRCLE,
                                                       NULL, NULL, data,
                                                       &lines[0]));
  data[0] = data[1] = data[2] = data[4] = data[5] = 0.0;
  data[3] = 1.0;
/*
  printf(" EG_makeGeometry L0 = %d\n", EG_makeGeometry(context, CURVE, LINE,
                                                       NULL, NULL, data,
                                                       &lines[0]));
*/
  data[3] = 0.0;
  data[4] = 2.0;
  printf(" EG_makeGeometry L1 = %d\n", EG_makeGeometry(context, CURVE, LINE,
                                                       NULL, NULL, data,
                                                       &lines[1]));
  data[0] =  1.0;
#if NEDGE == 3
  data[3] = -1.0;
#endif
  printf(" EG_makeGeometry L2 = %d\n", EG_makeGeometry(context, CURVE, LINE,
                                                       NULL, NULL, data,
                                                       &lines[2]));
  data[0] = 0.5;
  data[2] = data[4] = data[5] = data[6] = data[7] = 0.0;
  data[1] = 2.0;
  data[3] = data[8] = -1.0;
  data[9] = 0.5;
  printf(" EG_makeGeometry C3 = %d\n", EG_makeGeometry(context, CURVE, CIRCLE,
                                                       NULL, NULL, data,
                                                       &lines[3]));
  
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

#if NEDGE == 3
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
  
#else

  xyz[0] = 1.0;
  xyz[1] = 2.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[3], xyz, &range[0],
                                                      dum));
  xyz[0] = 0.0;
  xyz[1] = 2.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[3], xyz, &range[1],
                                                      dum));
  if (range[1] < range[0]) range[1] = 2.0*range[0];
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[3];
  objs[1] = nodes[2];
  printf(" EG_makeTopology E2 = %d\n", EG_makeTopology(context, lines[3], EDGE,
                                                       TWONODE, range, 2,
                                                       objs, NULL, &edges[2]));
  xyz[0] = 1.0;
  xyz[1] = 0.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[2], xyz, &range[0],
                                                      dum));
  xyz[0] = 1.0;
  xyz[1] = 2.0;
  printf(" EG_invEvaluate     = %d\n", EG_invEvaluate(lines[2], xyz, &range[1],
                                                      dum));
  printf("                      range = %lf %lf\n", range[0], range[1]);
  objs[0] = nodes[1];
  objs[1] = nodes[3];
  printf(" EG_makeTopology E3 = %d\n", EG_makeTopology(context, lines[2], EDGE,
                                                       TWONODE, range, 2,
                                                       objs, NULL, &edges[3]));
#endif
  
  /* make the loop -- without a surface */
  senses[0] = -1;
  senses[1] =  1;
  senses[2] = -1;
  senses[3] = -1;
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, LOOP,
                                                       CLOSED, NULL, NEDGE,
                                                       edges, senses, &loop));

  if (EG_isPlanar(loop) == 0) {

    printf(" EG_makeFace        = %d\n", EG_makeFace(loop, SFORWARD, NULL,
                                                     &face));
  } else {
    
    /*
     * makeFace like procedure:
     *   get the surface that fits the loop
     */
    printf(" EG_isoCline        = %d\n", EG_isoCline(loop, style, 0.0,
                                                     &surface));
  
    /*   update the loop to include the surface & pcurves   */
    printf(" EG_deleteObject L  = %d\n", EG_deleteObject(loop));
    for (i = 0; i < NEDGE; i++)
      printf(" EG_otherCurve  PC%d = %d\n",
                                      i, EG_otherCurve(surface, edges[i], 0.0,
                                                       &edges[i+NEDGE]));
    
    printf(" EG_makeTopology Ls = %d\n", EG_makeTopology(context, surface, LOOP,
                                                         CLOSED, NULL, NEDGE,
                                                         edges, senses, &loop));
    /*   make the face the "hard way" */
    printf(" EG_makeTopology F  = %d\n", EG_makeTopology(context, surface, FACE,
                                                         SFORWARD, NULL, 1,
                                                         &loop, senses, &face));
  }
  
  /* make a FACEBODY and save the geometry */
  printf(" EG_makeTopology B  = %d\n", EG_makeTopology(context, NULL, BODY,
                                                       FACEBODY, NULL, 1, &face,
                                                       NULL, &body));

  printf(" EG_getBodyTopos E  = %d\n", EG_getBodyTopos(body, NULL, EDGE, &i,
                                                       &bedges));
  for (i = 1; i <= 20; i++) t[i-1] = i/10.5;
  EG_attributeAdd(bedges[0], ".tPos", ATTRREAL, 20, NULL, t, NULL);
  EG_attributeAdd(bedges[2], ".tPos", ATTRREAL, 20, NULL, t, NULL);
  EG_free(bedges);

  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                       NULL, 1, &body,
                                                       NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "makeFace3D.egads"));

  /* cleanup */
  EG_setOutLevel(context, 0);
  printf(" EG_deleteObject  C = %d\n", EG_deleteObject(context));
  printf(" EG_deleteObject  M = %d\n", EG_deleteObject(model));
  EG_setOutLevel(context, 2);
  printf(" EG_close           = %d\n", EG_close(context));

  return 0;
}
