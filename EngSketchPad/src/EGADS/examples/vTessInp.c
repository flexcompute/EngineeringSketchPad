/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Tessellation using wv testing Tessellation Input
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

/* #define RELOAD         */    // open up existing tessellation


#include "egads.h"
#include "wsserver.h"

#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif
#include <winsock2.h>
#endif


/* structure to hold on to the EGADS triangulation per Body */
typedef struct {
  ego *faces;
  ego *edges;
  ego body;
  ego tess;
  int mtype;
  int nfaces;
  int nedges;
} bodyData;


/* globals used in these functions */
  static wvContext *cntxt;
  static bodyData  *bodydata;



/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *uPtr, /*@unused@*/ void *wsi,
                    char *text, /*@unused@*/ int lena)
{

  printf(" Recv'ed: %s\n", text);

}


int main(int argc, char *argv[])
{
  int          i, j, k, ibody, stat, oclass, mtype, len, ntri, sum, state, npts;
  int          nseg, *segs, *senses;
  int          nbody;
  const int    *tris, *tric, *ptype, *pindex;
  float        arg, focus[4], color[3];
  double       box[6], size, params[3];
  const double *xyzs, *uvs, *ts;
  char         gpname[34], *startapp;
  const char   *OCCrev;
  ego          context, model, geom, tess, body, *bodies, *dum;
  wvData       items[5];
  int          sides[3][2] = {{1,2}, {2,0}, {0,1}};
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vTessInp filename [angle maxlen sag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d %s\n\n", i, j, OCCrev);

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  printf(" EG_loadModel      = %d\n", EG_loadModel(context, 0, argv[1], &model));
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");

                            size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

  focus[0] = 0.5*(box[0]+box[3]);
  focus[1] = 0.5*(box[1]+box[4]);
  focus[2] = 0.5*(box[2]+box[5]);
  focus[3] = size;

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyData *) malloc(nbody*sizeof(bodyData));
  if (bodydata == NULL) {
    printf(" MALLOC Error on Body storage!\n");
    return 1;
  }

  params[0] =  0.025*size;
  params[1] =  0.001*size;
  params[2] = 15.0;
  if (argc == 5) {
    sscanf(argv[2], "%f", &arg);
    params[2] = arg;
    sscanf(argv[3], "%f", &arg);
    params[0] = arg;
    sscanf(argv[4], "%f", &arg);
    params[1] = arg;
    printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
           params[2], params[0], params[1]);
    params[0] *= size;
    params[1] *= size;
  }

  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
    mtype = 0;
    EG_getTopology(bodies[ibody], &geom, &oclass,
                   &mtype, NULL, &j, &dum, &senses);
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", ibody+1);
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", ibody+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", ibody+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", ibody+1);
    }

    stat = EG_getBodyTopos(bodies[ibody], NULL, FACE,
                           &bodydata[ibody].nfaces, &bodydata[ibody].faces);
    i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE,
                           &bodydata[ibody].nedges, &bodydata[ibody].edges);
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      return 1;
    }

#ifdef RELOAD
    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
    stat = EG_openTessBody(bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_openTessBody %d = %d\n", ibody, stat);
      continue;
    }
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len, &xyzs, &uvs,
                            &ptype, &pindex, &ntri, &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        printf("       %d EG_getTessFace %d = %d!\n", ibody, i+1, stat);
        continue;
      }
      stat = EG_setTessFace(bodydata[ibody].tess, i+1, len, xyzs, uvs,
                            ntri, tris);
      if (stat != EGADS_SUCCESS)
        printf("       %d EG_setTessFace %d = %d!\n", ibody, i+1, stat);
    }
    stat = EG_statusTessBody(bodydata[ibody].tess, &body, &state, &npts);
    printf(" statusTessBody %d = %d\n", ibody, stat);
#else
    stat = EG_makeTessBody(bodies[ibody], params, &tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
    stat = EG_initTessBody(bodies[ibody], &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_initTessBody %d = %d\n", ibody, stat);
      continue;
    }
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      stat = EG_getTessEdge(tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) {
        printf("       %d EG_getTessEdge %d = %d!\n", ibody, i+1, stat);
        continue;
      }
      stat = EG_setTessEdge(bodydata[ibody].tess, i+1, len, xyzs, ts);
      if (stat != EGADS_SUCCESS)
        printf("       %d EG_setTessEdge %d = %d!\n", ibody, i+1, stat);
    }
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(tess, i+1, &len, &xyzs, &uvs, &ptype, &pindex,
                            &ntri, &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        printf("       %d EG_getTessFace %d = %d!\n", ibody, i+1, stat);
        continue;
      }
      stat = EG_setTessFace(bodydata[ibody].tess, i+1, len, xyzs, uvs,
                            ntri, tris);
      if (stat != EGADS_SUCCESS)
        printf("       %d EG_setTessFace %d = %d!\n", ibody, i+1, stat);
    }
    EG_deleteObject(tess);
    stat = EG_statusTessBody(bodydata[ibody].tess, &body, &state, &npts);
    printf(" statusTessBody %d = %d\n", ibody, stat);
#endif
  }
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
    printf(" failed to create wvContext!\n");
    for (ibody = 0; ibody < nbody; ibody++) {
      EG_deleteObject(bodydata[ibody].tess);
      EG_free(bodydata[ibody].edges);
      EG_free(bodydata[ibody].faces);
    }
    free(bodydata);

    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
    printf(" EG_close          = %d\n", EG_close(context));
    return 1;
  }

  /* make the scene */
  for (sum = stat = ibody = 0; ibody < nbody; ibody++) {

    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      snprintf(gpname, 34, "Body %d Face %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      color[0]  = 1.0;
      color[1]  = ibody;
      color[1] /= nbody;
      color[2]  = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      for (nseg = j = 0; j < ntri; j++)
        for (k = 0; k < 3; k++)
          if (tric[3*j+k] < j+1) nseg++;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      for (nseg = j = 0; j < ntri; j++)
        for (k = 0; k < 3; k++)
          if (tric[3*j+k] < j+1) {
            segs[2*nseg  ] = tris[3*j+sides[k][0]];
            segs[2*nseg+1] = tris[3*j+sides[k][1]];
            nseg++;
          }
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", i, gpname);
      free(segs);
/*    color[0] = color[1] = color[2] = 0.8;  */
      color[0] = color[1] = color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", i, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
                         WV_ON|WV_ORIENTATION, 5, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      sum += ntri;
    }

    /* get edges */
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      stat  = EG_getTessEdge(bodydata[ibody].tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) continue;
      nseg = len-1;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Edge %d\n",
               nseg, ibody, i+1);
        continue;
      }
      for (j = 0; j < len-1; j++) {
        segs[2*j  ] = j + 1;
        segs[2*j+1] = j + 2;
      }

      snprintf(gpname, 34, "Body %d Edge %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs, WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS,   &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_INDICES,  &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(segs);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
      if (stat < 0) {
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      } else {
        if (cntxt != NULL)
          if (cntxt->gPrims != NULL) {
            cntxt->gPrims[stat].lWidth = 1.5;
            if (wv_addArrowHeads(cntxt, stat, 0.05, 1, &nseg) != 0)
              printf(" wv_addArrowHeads Error\n");
          }
      }
    }
  }
  printf(" ** %d gPrims with %d triangles **\n", stat+1, sum);

  /* start the server code */

  stat = 0;
  wv_setCallBack(cntxt, browserMessage);
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    /* we have a single valid server -- stay alive a long as we have a client */
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
    }
  }
  wv_cleanupServers();

  /* finish up */
  for (ibody = 0; ibody < nbody; ibody++) {
    EG_deleteObject(bodydata[ibody].tess);
    EG_free(bodydata[ibody].edges);
    EG_free(bodydata[ibody].faces);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
