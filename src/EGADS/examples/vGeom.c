/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Display the EGADS Geometry using wv (the WebViewer)
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#endif

#include "egads.h"
#include "wsserver.h"


#define NCOLOR 256
//#define CONVERT


/* structure to hold on to the EGADS discretization per Body */
typedef struct {
  ego *faces;
  ego *faceTess;
  ego *edges;
  ego *edgeTess;
  ego body;
  int mtype;
  int nfaces;
  int nedges;
} bodyData;


  static wvContext *cntxt;

  
int main(int argc, char *argv[])
{
  int        i, j, k, m, ibody, stat, oclass, mtype, len, ntri, sum, nbody;
  int        nloops, nseg, *segs, *senses, sizes[2], *tris, per, nnodes, ndiv;
  float      color[3], focus[4];
  double     box[6], tlimits[4], data[14], size, *xyzs;
  char       gpname[33], *startapp;
  const char *OCCrev;
  ego        context, model, geom, *bodies, *dum, *loops, *nodes;
  wvData     items[5];
  bodyData   *bodydata;
  float      eye[3]    = {0.0, 0.0, 7.0};
  float      center[3] = {0.0, 0.0, 0.0};
  float      up[3]     = {0.0, 1.0, 0.0};

  /* get our starting application line 
   * 
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");
  
  if ((argc != 2) && (argc != 3)) {
    printf("\n Usage: vGeom filename [nDiv]\n\n");
    return 1;
  }
  ndiv = 37;
  if (argc == 3) {
    ndiv = atoi(argv[2]);
    printf("\n nDiv = %d\n\n", ndiv);
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
  
  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
  
    mtype = 0;
    EG_getTopology(bodies[ibody], &geom, &oclass, 
                   &mtype, NULL, &j, &dum, &senses);
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", ibody+1);
      stat = EG_getMassProperties(bodydata[ibody].body, data);
      if (stat == EGADS_SUCCESS)
        printf("                 CoG = %lf %lf %lf\n", data[2], data[3], data[4]);
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
      free(bodydata);
      return 1;
    }
    
    printf(" EG_getBodyTopos:    %d nFaces  = %d\n", 
           ibody+1, bodydata[ibody].nfaces);
    printf(" EG_getBodyTopos:    %d nEdges  = %d\n", 
           ibody+1, bodydata[ibody].nedges);
    bodydata[ibody].faceTess = (ego *) 
                        malloc(bodydata[ibody].nfaces*sizeof(ego));
    bodydata[ibody].edgeTess = (ego *) 
                        malloc(bodydata[ibody].nedges*sizeof(ego));
    if ((bodydata[ibody].faceTess == NULL) || 
        (bodydata[ibody].edgeTess == NULL)) {
      printf(" CAN NOT Allocate Face Tessellation Objects!\n");
      free(bodydata);
      return 1;
    }
    
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      bodydata[ibody].faceTess[i] = NULL;
      stat = EG_getTopology(bodydata[ibody].faces[i], &geom, &oclass, 
                            &mtype, NULL, &nloops, &loops, &senses);
      if (stat != EGADS_SUCCESS) continue;
      printf(" EG_getTopology:     %d Face %d -- nLoops = %d\n",
             ibody+1, i+1, nloops);
      stat = EG_getMassProperties(bodydata[ibody].faces[i], data);
      if (stat == EGADS_SUCCESS) 
      printf("                 CoG = %lf %lf %lf   Area = %le\n",
             data[2], data[3], data[4], data[1]);
      sizes[0] = ndiv;
      sizes[1] = ndiv;
      if (mtype == SREVERSE) sizes[0] = -ndiv;
#ifdef CONVERT
      if (geom->mtype == BSPLINE) {
        stat = EG_getRange(bodydata[ibody].faces[i], tlimits, &per);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_getRange Face return = %d!\n", stat);
          free(bodydata);
          return 1;
        }
      } else {
        stat = EG_convertToBSpline(bodydata[ibody].faces[i], &geom);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_convertToBSpline Face return = %d!\n", stat);
          free(bodydata);
          return 1;
        }
        stat = EG_getRange(geom, tlimits, &per);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_getRange Face return = %d!\n", stat);
          free(bodydata);
          return 1;
        }
      }
#else
      stat = EG_getRange(bodydata[ibody].faces[i], tlimits, &per);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getRange Face return = %d!\n", stat);
        free(bodydata);
        return 1;
      }
#endif
      stat = EG_makeTessGeom(geom, tlimits, sizes,
                             &bodydata[ibody].faceTess[i]);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_makeTessGeom Face return = %d!\n", stat);
        free(bodydata);
        return 1;
      }
    }
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      double range[2];
      
      bodydata[ibody].edgeTess[i] = NULL;
      stat = EG_getTopology(bodydata[ibody].edges[i], &geom, &oclass, 
                            &mtype, range, &nnodes, &nodes, &senses);
      if (stat != EGADS_SUCCESS) continue;
/*    {
        int nfac, ii;
        ego *facs;
        stat = EG_getBodyTopos(bodies[ibody], bodydata[ibody].edges[i], FACE,
                               &nfac, &facs);
        if (stat != EGADS_SUCCESS) {
          printf("  EG_getBodyTopo = %d\n", stat);
        } else {
          for (ii = 0; ii < nfac; ii++) {
            stat = EG_getEdgeUV(facs[ii], nodes[0], 0, 0.0, data);
            printf(" EG_getEdgeUV Node = %d  %lf %lf\n", stat, data[0], data[1]);
            stat = EG_getEdgeUV(facs[ii], bodydata[ibody].edges[i], 0, range[0],
                                data);
            printf(" EG_getEdgeUV Edge = %d  %lf %lf\n", stat, data[0], data[1]);
          }
          EG_free(facs);
        }
      }  */
      if (mtype == DEGENERATE) continue;
      printf(" EG_getTopology:     %d Edge %d -- nNodes = %d\n",
             ibody+1, i+1, nnodes);
#ifdef CONVERT
      if (geom->mtype == BSPLINE) {
        stat = EG_getRange(bodydata[ibody].edges[i], tlimits, &per);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_getRange Edge return = %d!\n", stat);
          free(bodydata);
          return 1;
        }
      } else {
        stat = EG_convertToBSpline(bodydata[ibody].edges[i], &geom);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_convertToBSpline Edge return = %d!\n", stat);
          free(bodydata);
          return 1;
        }
        stat = EG_getRange(geom, tlimits, &per);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_getRange Edge return = %d!\n", stat);
          free(bodydata);
          return 1;
        }
      }
#else
      stat = EG_getRange(bodydata[ibody].edges[i], tlimits, &per);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getRange Edge return = %d!\n", stat);
        free(bodydata);
        return 1;
      }
#endif
      sizes[0] = ndiv;
      sizes[1] = 0;               /* not required */
      stat = EG_makeTessGeom(geom, tlimits, sizes, 
                             &bodydata[ibody].edgeTess[i]);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_makeTessGeom Edge return = %d!\n", stat);
        free(bodydata);
        return 1;
      }
    }
  }
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) printf(" failed to create wvContext!\n");
        
  /* make the scene */
  for (sum = stat = ibody = 0; ibody < nbody; ibody++) {
  
    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessGeom(bodydata[ibody].faceTess[i], sizes, &xyzs);
      if (stat != EGADS_SUCCESS) continue;
      len  = sizes[0]*sizes[1];
      ntri = 2*(sizes[0]-1)*(sizes[1]-1);
      tris = (int *) malloc(3*ntri*sizeof(int));
      if (tris == NULL) {
        printf(" Can not allocate %d triangles for Body %d Face %d\n",
               ntri, ibody, i+1);
        continue;
      }
      nseg = (sizes[0]-1)*sizes[1] + sizes[0]*(sizes[1]-1);
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Face %d\n",
               nseg, ibody, i+1);
        free(tris);
        continue;
      }
      for (m = k = 0; k < sizes[1]-1; k++)
        for (j = 0; j < sizes[0]-1; j++) {
          tris[3*m  ] = j + sizes[0]*k     + 1;
          tris[3*m+1] = j + sizes[0]*k     + 2;
          tris[3*m+2] = j + sizes[0]*(k+1) + 2;
          m++;
          tris[3*m  ] = j + sizes[0]*(k+1) + 2;
          tris[3*m+1] = j + sizes[0]*(k+1) + 1;
          tris[3*m+2] = j + sizes[0]*k     + 1;
          m++;          
        }
      for (m = k = 0; k < sizes[1]; k++)
        for (j = 0; j < sizes[0]-1; j++) {
          segs[2*m  ] = j + sizes[0]*k + 1;
          segs[2*m+1] = j + sizes[0]*k + 2;
          m++;          
        }
      for (k = 0; k < sizes[1]-1; k++)
        for (j = 0; j < sizes[0]; j++) {
          segs[2*m  ] = j + sizes[0]*k     + 1;
          segs[2*m+1] = j + sizes[0]*(k+1) + 1;
          m++;          
        }

      sprintf(gpname, "Body %d Face %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      free(tris);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      color[0]  = 1.0;
      color[1]  = ibody;
      color[1] /= nbody;
      color[2]  = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
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
      if (bodydata[ibody].edgeTess[i] == NULL) continue;
      stat  = EG_getTessGeom(bodydata[ibody].edgeTess[i], sizes, &xyzs);
      if (stat != EGADS_SUCCESS) continue;
      len  = sizes[0];
      nseg = sizes[0]-1;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Edge %d\n",
               nseg, ibody, i+1);
        continue;
      }
      for (j = 0; j < sizes[0]-1; j++) {
        segs[2*j  ] = j + 1;
        segs[2*j+1] = j + 2;          
      }

      sprintf(gpname, "Body %d Edge %d", ibody+1, i+1);
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
          if (cntxt->gPrims != NULL) cntxt->gPrims[stat].lWidth = 1.5;
      }
    }
  }
  printf(" ** %d gPrims with %d triangles **\n", stat+1, sum);

  /* start the server code */

  stat = 0;
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    /* we have a single valid server -- stay alive a long as we have a client */
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
/*    wv_broadcastText("I'm here!");  */
    }
  }
  wv_cleanupServers();

  /* finish up */
  for (ibody = 0; ibody < nbody; ibody++) {
    for (i = 0; i < bodydata[ibody].nedges; i++) 
      EG_deleteObject(bodydata[ibody].edgeTess[i]);
    free(bodydata[ibody].edgeTess);
    for (i = 0; i < bodydata[ibody].nfaces; i++) 
      EG_deleteObject(bodydata[ibody].faceTess[i]);
    free(bodydata[ibody].faceTess);
    EG_free(bodydata[ibody].faces);
    EG_free(bodydata[ibody].edges);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}


/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, char *text, /*@unused@*/ int lena)
{
  int          stat;
  static float cols[NCOLOR*3] =
      { 0.0000,    0.0000,    1.0000,     0.0000,    0.0157,    1.0000,
        0.0000,    0.0314,    1.0000,     0.0000,    0.0471,    1.0000,
        0.0000,    0.0627,    1.0000,     0.0000,    0.0784,    1.0000,
        0.0000,    0.0941,    1.0000,     0.0000,    0.1098,    1.0000,
        0.0000,    0.1255,    1.0000,     0.0000,    0.1412,    1.0000,
        0.0000,    0.1569,    1.0000,     0.0000,    0.1725,    1.0000,
        0.0000,    0.1882,    1.0000,     0.0000,    0.2039,    1.0000,
        0.0000,    0.2196,    1.0000,     0.0000,    0.2353,    1.0000,
        0.0000,    0.2510,    1.0000,     0.0000,    0.2667,    1.0000,
        0.0000,    0.2824,    1.0000,     0.0000,    0.2980,    1.0000,
        0.0000,    0.3137,    1.0000,     0.0000,    0.3294,    1.0000,
        0.0000,    0.3451,    1.0000,     0.0000,    0.3608,    1.0000,
        0.0000,    0.3765,    1.0000,     0.0000,    0.3922,    1.0000,
        0.0000,    0.4078,    1.0000,     0.0000,    0.4235,    1.0000,
        0.0000,    0.4392,    1.0000,     0.0000,    0.4549,    1.0000,
        0.0000,    0.4706,    1.0000,     0.0000,    0.4863,    1.0000,
        0.0000,    0.5020,    1.0000,     0.0000,    0.5176,    1.0000,
        0.0000,    0.5333,    1.0000,     0.0000,    0.5490,    1.0000,
        0.0000,    0.5647,    1.0000,     0.0000,    0.5804,    1.0000,
        0.0000,    0.5961,    1.0000,     0.0000,    0.6118,    1.0000,
        0.0000,    0.6275,    1.0000,     0.0000,    0.6431,    1.0000,
        0.0000,    0.6588,    1.0000,     0.0000,    0.6745,    1.0000,
        0.0000,    0.6902,    1.0000,     0.0000,    0.7059,    1.0000,
        0.0000,    0.7216,    1.0000,     0.0000,    0.7373,    1.0000,
        0.0000,    0.7529,    1.0000,     0.0000,    0.7686,    1.0000,
        0.0000,    0.7843,    1.0000,     0.0000,    0.8000,    1.0000,
        0.0000,    0.8157,    1.0000,     0.0000,    0.8314,    1.0000,
        0.0000,    0.8471,    1.0000,     0.0000,    0.8627,    1.0000,
        0.0000,    0.8784,    1.0000,     0.0000,    0.8941,    1.0000,
        0.0000,    0.9098,    1.0000,     0.0000,    0.9255,    1.0000,
        0.0000,    0.9412,    1.0000,     0.0000,    0.9569,    1.0000,
        0.0000,    0.9725,    1.0000,     0.0000,    0.9882,    1.0000,
        0.0000,    1.0000,    0.9961,     0.0000,    1.0000,    0.9804,
        0.0000,    1.0000,    0.9647,     0.0000,    1.0000,    0.9490,
        0.0000,    1.0000,    0.9333,     0.0000,    1.0000,    0.9176,
        0.0000,    1.0000,    0.9020,     0.0000,    1.0000,    0.8863,
        0.0000,    1.0000,    0.8706,     0.0000,    1.0000,    0.8549,
        0.0000,    1.0000,    0.8392,     0.0000,    1.0000,    0.8235,
        0.0000,    1.0000,    0.8078,     0.0000,    1.0000,    0.7922,
        0.0000,    1.0000,    0.7765,     0.0000,    1.0000,    0.7608,
        0.0000,    1.0000,    0.7451,     0.0000,    1.0000,    0.7294,
        0.0000,    1.0000,    0.7137,     0.0000,    1.0000,    0.6980,
        0.0000,    1.0000,    0.6824,     0.0000,    1.0000,    0.6667,
        0.0000,    1.0000,    0.6510,     0.0000,    1.0000,    0.6353,
        0.0000,    1.0000,    0.6196,     0.0000,    1.0000,    0.6039,
        0.0000,    1.0000,    0.5882,     0.0000,    1.0000,    0.5725,
        0.0000,    1.0000,    0.5569,     0.0000,    1.0000,    0.5412,
        0.0000,    1.0000,    0.5255,     0.0000,    1.0000,    0.5098,
        0.0000,    1.0000,    0.4941,     0.0000,    1.0000,    0.4784,
        0.0000,    1.0000,    0.4627,     0.0000,    1.0000,    0.4471,
        0.0000,    1.0000,    0.4314,     0.0000,    1.0000,    0.4157,
        0.0000,    1.0000,    0.4000,     0.0000,    1.0000,    0.3843,
        0.0000,    1.0000,    0.3686,     0.0000,    1.0000,    0.3529,
        0.0000,    1.0000,    0.3373,     0.0000,    1.0000,    0.3216,
        0.0000,    1.0000,    0.3059,     0.0000,    1.0000,    0.2902,
        0.0000,    1.0000,    0.2745,     0.0000,    1.0000,    0.2588,
        0.0000,    1.0000,    0.2431,     0.0000,    1.0000,    0.2275,
        0.0000,    1.0000,    0.2118,     0.0000,    1.0000,    0.1961,
        0.0000,    1.0000,    0.1804,     0.0000,    1.0000,    0.1647,
        0.0000,    1.0000,    0.1490,     0.0000,    1.0000,    0.1333,
        0.0000,    1.0000,    0.1176,     0.0000,    1.0000,    0.1020,
        0.0000,    1.0000,    0.0863,     0.0000,    1.0000,    0.0706,
        0.0000,    1.0000,    0.0549,     0.0000,    1.0000,    0.0392,
        0.0000,    1.0000,    0.0235,     0.0000,    1.0000,    0.0078,
        0.0078,    1.0000,    0.0000,     0.0235,    1.0000,    0.0000,
        0.0392,    1.0000,    0.0000,     0.0549,    1.0000,    0.0000,
        0.0706,    1.0000,    0.0000,     0.0863,    1.0000,    0.0000,
        0.1020,    1.0000,    0.0000,     0.1176,    1.0000,    0.0000,
        0.1333,    1.0000,    0.0000,     0.1490,    1.0000,    0.0000,
        0.1647,    1.0000,    0.0000,     0.1804,    1.0000,    0.0000,
        0.1961,    1.0000,    0.0000,     0.2118,    1.0000,    0.0000,
        0.2275,    1.0000,    0.0000,     0.2431,    1.0000,    0.0000,
        0.2588,    1.0000,    0.0000,     0.2745,    1.0000,    0.0000,
        0.2902,    1.0000,    0.0000,     0.3059,    1.0000,    0.0000,
        0.3216,    1.0000,    0.0000,     0.3373,    1.0000,    0.0000,
        0.3529,    1.0000,    0.0000,     0.3686,    1.0000,    0.0000,
        0.3843,    1.0000,    0.0000,     0.4000,    1.0000,    0.0000,
        0.4157,    1.0000,    0.0000,     0.4314,    1.0000,    0.0000,
        0.4471,    1.0000,    0.0000,     0.4627,    1.0000,    0.0000,
        0.4784,    1.0000,    0.0000,     0.4941,    1.0000,    0.0000,
        0.5098,    1.0000,    0.0000,     0.5255,    1.0000,    0.0000,
        0.5412,    1.0000,    0.0000,     0.5569,    1.0000,    0.0000,
        0.5725,    1.0000,    0.0000,     0.5882,    1.0000,    0.0000,
        0.6039,    1.0000,    0.0000,     0.6196,    1.0000,    0.0000,
        0.6353,    1.0000,    0.0000,     0.6510,    1.0000,    0.0000,
        0.6667,    1.0000,    0.0000,     0.6824,    1.0000,    0.0000,
        0.6980,    1.0000,    0.0000,     0.7137,    1.0000,    0.0000,
        0.7294,    1.0000,    0.0000,     0.7451,    1.0000,    0.0000,
        0.7608,    1.0000,    0.0000,     0.7765,    1.0000,    0.0000,
        0.7922,    1.0000,    0.0000,     0.8078,    1.0000,    0.0000,
        0.8235,    1.0000,    0.0000,     0.8392,    1.0000,    0.0000,
        0.8549,    1.0000,    0.0000,     0.8706,    1.0000,    0.0000,
        0.8863,    1.0000,    0.0000,     0.9020,    1.0000,    0.0000,
        0.9176,    1.0000,    0.0000,     0.9333,    1.0000,    0.0000,
        0.9490,    1.0000,    0.0000,     0.9647,    1.0000,    0.0000,
        0.9804,    1.0000,    0.0000,     0.9961,    1.0000,    0.0000,
        1.0000,    0.9882,    0.0000,     1.0000,    0.9725,    0.0000,
        1.0000,    0.9569,    0.0000,     1.0000,    0.9412,    0.0000,
        1.0000,    0.9255,    0.0000,     1.0000,    0.9098,    0.0000,
        1.0000,    0.8941,    0.0000,     1.0000,    0.8784,    0.0000,
        1.0000,    0.8627,    0.0000,     1.0000,    0.8471,    0.0000,
        1.0000,    0.8314,    0.0000,     1.0000,    0.8157,    0.0000,
        1.0000,    0.8000,    0.0000,     1.0000,    0.7843,    0.0000,
        1.0000,    0.7686,    0.0000,     1.0000,    0.7529,    0.0000,
        1.0000,    0.7373,    0.0000,     1.0000,    0.7216,    0.0000,
        1.0000,    0.7059,    0.0000,     1.0000,    0.6902,    0.0000,
        1.0000,    0.6745,    0.0000,     1.0000,    0.6588,    0.0000,
        1.0000,    0.6431,    0.0000,     1.0000,    0.6275,    0.0000,
        1.0000,    0.6118,    0.0000,     1.0000,    0.5961,    0.0000,
        1.0000,    0.5804,    0.0000,     1.0000,    0.5647,    0.0000,
        1.0000,    0.5490,    0.0000,     1.0000,    0.5333,    0.0000,
        1.0000,    0.5176,    0.0000,     1.0000,    0.5020,    0.0000,
        1.0000,    0.4863,    0.0000,     1.0000,    0.4706,    0.0000,
        1.0000,    0.4549,    0.0000,     1.0000,    0.4392,    0.0000,
        1.0000,    0.4235,    0.0000,     1.0000,    0.4078,    0.0000,
        1.0000,    0.3922,    0.0000,     1.0000,    0.3765,    0.0000,
        1.0000,    0.3608,    0.0000,     1.0000,    0.3451,    0.0000,
        1.0000,    0.3294,    0.0000,     1.0000,    0.3137,    0.0000,
        1.0000,    0.2980,    0.0000,     1.0000,    0.2824,    0.0000,
        1.0000,    0.2667,    0.0000,     1.0000,    0.2510,    0.0000,
        1.0000,    0.2353,    0.0000,     1.0000,    0.2196,    0.0000,
        1.0000,    0.2039,    0.0000,     1.0000,    0.1882,    0.0000,
        1.0000,    0.1725,    0.0000,     1.0000,    0.1569,    0.0000,
        1.0000,    0.1412,    0.0000,     1.0000,    0.1255,    0.0000,
        1.0000,    0.1098,    0.0000,     1.0000,    0.0941,    0.0000,
        1.0000,    0.0784,    0.0000,     1.0000,    0.0627,    0.0000,
        1.0000,    0.0471,    0.0000,     1.0000,    0.0314,    0.0000,
        1.0000,    0.0157,    0.0000,     1.0000,    0.0000,    0.0000 };

  /* ping it back
  wv_sendText(wsi, text); */
  if (strcmp(text,"coarser") == 0) {
    stat = wv_setKey(cntxt, 0,      NULL, 0.0, 0.0, NULL);
  } else {
    stat = wv_setKey(cntxt, NCOLOR, cols, 0.0, 1.0, text);
  }
  printf(" setKey %s = %d\n", text, stat);
  
}
