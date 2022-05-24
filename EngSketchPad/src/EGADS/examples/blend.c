/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             An example using the blend function
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <assert.h>
#include <math.h>
#include "egads.h"

//#define PERIODIC
//#define OPENLOOP
//#define SAVEMODEL


int main(int argc, char *argv[])
{
  int    i, j, beg, end, brnd, ernd, status, senses[2], per, nface, rev = 1;
  int    btgt, etgt, atype, alen, nsec;
  double xyz[3], data[10], range[2], xform[12], rc[8], *rc1, *rcN;
  double area, btan[4], etan[4];
  ego    context, nodes[3], curve, edges[2], objs[2], oform, secs[7];
  ego    others[2], body, *faces;
#ifdef SAVEMODEL
  ego    model;
#endif
  const  char   *str;
  const  int    *ints;
  const  double *reals;
  static char   *types[4] = {"Open", "Node", "Face", "Tip"};

  rc1   = rcN = NULL;
  rc[0] = 0.05;
  rc[4] = 0.4;
  rc[1] = rc[6] = 1.0;
  rc[2] = rc[3] = rc[5] = rc[7] = 0.0;
  brnd  = ernd = 0;
  btgt  = etgt = 0;
  if (argc != 3) {
    printf("\n Usage: blend 0/1/2/3 0/1/2/3 -- 0-open, 1-node, 2-face, 3-tip\n\n");
    return 1;
  }
  sscanf(argv[1], "%d", &beg);
  sscanf(argv[2], "%d", &end);
  if (beg == -1) beg = brnd = 1;
  if (end == -1) end = ernd = 1;
  if (beg == -2) beg = btgt = 2;
  if (end == -2) end = etgt = 2;
  if (beg <   0) beg = 0;
  if (beg >   3) beg = 3;
  if (end <   0) end = 0;
  if (end >   3) end = 3;
  printf("\n Blend with start %s and end %s\n\n", types[beg], types[end]);
  if (btgt == 2) {
    btan[0] = 1.0;
    printf(" Enter tangent at beginning: ");
    scanf(" %lf %lf %lf", &btan[1], &btan[2], &btan[3]);
    printf("\n");
  }
  if (etgt == 2) {
    etan[0] = 1.0;
    printf(" Enter tangent at end: ");
    scanf(" %lf %lf %lf", &etan[1], &etan[2], &etan[3]);
    printf("\n");
  }
  if (beg == 3) {
    btan[0] = 0.0;
    btan[1] = 1.0;
    rc1     = btan;
  }
  if (end == 3) {
    etan[0] = 0.0;
    etan[1] = 1.0;
    rcN     = etan;
  }
  
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
#ifdef PERIODIC
  objs[0] = nodes[0];
  objs[1] = nodes[1];
  data[0] = range[0];
  data[1] = range[1];
  printf(" EG_makeTopology E0 = %d\n", EG_makeTopology(context, curve, EDGE,
                                                       ONENODE, data, 1,
                                                       objs, NULL, &edges[0]));
  /* make the Loop */
  senses[0] = rev;
  printf(" EG_makeTopology L  = %d\n", EG_makeTopology(context, NULL, LOOP,
                                                       CLOSED, NULL, 1, edges,
                                                       senses, &secs[0]));
#else
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
    if (brnd == 1) rc1 = rc;
    EG_deleteObject(secs[0]);
    xyz[0] = xyz[1] = xyz[2] = 0.0;
    printf(" EG_makeTopology Nb = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                         xyz, 0, NULL, NULL,
                                                         &secs[0]));
  } else if (beg >= 2) {
    if (btgt == 2) rc1 = btan;
    others[0] = secs[0];
    per       = SFORWARD;
    area      = 0.0;
    printf(" EG_getArea  beg    = %d  area = %lf\n", EG_getArea(others[0], NULL,
                                                                &area), area);
    if (area < 0.0) per = SREVERSE;
    printf(" EG_makeFace beg    = %d\n", EG_makeFace(others[0], per, NULL,
                                                     &secs[0]));
  }
  if (end == 1) {
    if (ernd == 1) rcN = rc;
    EG_deleteObject(secs[4]);
    xyz[0] = xyz[1] = 0.0;
    xyz[2] = 5.0;
    printf(" EG_makeTopology Ne = %d\n", EG_makeTopology(context, NULL, NODE, 0,
                                                         xyz, 0, NULL, NULL,
                                                         &secs[4]));
  } else if (end >= 2) {
    if (etgt == 2) rcN = etan;
    others[1] = secs[4];
    per       = SFORWARD;
    area      = 0.0;
    printf(" EG_getArea  end    = %d  area = %lf\n", EG_getArea(others[1], NULL,
                                                                &area), area);
    if (area < 0.0) per = SREVERSE;
    printf(" EG_makeFace end    = %d\n", EG_makeFace(others[1], per, NULL,
                                                     &secs[4]));
  }
  nsec = 5;
  
/* multiplicity of sections
  for (i = nsec-1; i > 2; i--) secs[i+1] = secs[i];
  secs[3] = secs[2];
  nsec++;
*/
/*
  for (i = nsec-1; i > 2; i--) secs[i+2] = secs[i];
  secs[3] = secs[4] = secs[2];
  nsec += 2;
*/

  /* blend & save */
  printf(" EG_blend           = %d\n", EG_blend(nsec, secs, rc1, rcN, &body));
  printf(" EG_getBodyTopos    = %d\n", EG_getBodyTopos(body, NULL, FACE,
                                                       &nface, &faces));

  for (i = 0; i < nface; i++) {
    status = EG_attributeRet(faces[i], ".blendSamples", &atype, &alen,
                             &ints, &reals, &str);
    if (status == EGADS_SUCCESS) {
      printf("   Face %d/%d: blendSamples =", i+1, nface);
      if (atype == ATTRREAL) {
        for (j = 0; j < alen; j++) printf(" %lf", reals[j]);
      } else {
        printf(" atype = %d, alen = %d", atype, alen);
      }
      printf("\n");
    }
    status = EG_attributeRet(faces[i], ".blendSenses", &atype, &alen,
                             &ints, &reals, &str);
    if (status == EGADS_SUCCESS) {
      printf("   Face %d/%d: blendSenses =", i+1, nface);
      if (atype == ATTRINT) {
        for (j = 0; j < alen; j++) printf(" %d", ints[j]);
      } else {
        printf(" atype = %d, alen = %d", atype, alen);
      }
      printf("\n");
    }
  }
#ifdef SAVEMODEL
  printf(" EG_makeTopology M  = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                       NULL, 1, &body,
                                                       NULL, &model));
  printf(" EG_saveModel       = %d\n", EG_saveModel(model, "blend.egads"));
  printf("\n");
  EG_free(faces);

  /* cleanup */
  EG_deleteObject(model);
#else
  printf(" EG_saveModel       = %d\n", EG_saveModel(body, "blend.egads"));
  printf("\n");
  EG_free(faces);
  
  /* cleanup */
  EG_deleteObject(body);
#endif
  for (i = nsec-1; i >= 0; i--) EG_deleteObject(secs[i]);
  if (others[0] != NULL) EG_deleteObject(others[0]);
  if (others[1] != NULL) EG_deleteObject(others[1]);
#ifndef PERIODIC
  EG_deleteObject(edges[1]);
#endif
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
