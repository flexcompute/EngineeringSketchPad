/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Display the EGADS Tessellation using wv (the WebViewer)
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

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

//#define INNER


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



/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *uPtr, /*@unused@*/ void *wsi,
                    char *text, /*@unused@*/ int lena)
{

  printf(" RX: %s\n", text);
  /* ping it back
  wv_sendText(wsi, text); */

}


int main(int argc, char *argv[])
{
  int          i, j, k, m, n, ibody, stat, oclass, mtype, len, ntri, sum, nbody;
  int          nloops, nseg, nledges, *segs, *lsenses, *senses, *esenses;
  int          npatch, nq, ip, n1, n2, nh, *trix, *heads;
  const int    *tris, *tric, *ptype, *pindex, *pvi, *pvb;
  float        arg, color[3], focus[4], *lsegs;
  double       box[6], size, params[3], qparam[3] = {0.0, 0.0, 0.0};
#ifdef INNER
  double       uvbox[4];
#endif
  const double *xyzs, *uvs, *ts;
  char         gpname[46], *startapp;
  const char   *OCCrev;
  ego          context, model, geom, *bodies, *dum, *ledges, *loops;
  wvContext    *cntxt;
  wvData       items[5];
  bodyData     *bodydata;
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vQuad filename [angle maxlen sag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

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
      free(bodydata);
      return 1;
    }

#ifdef INNER
    /* set the innerLoop attribute */
    for (m = 0; m < bodydata[ibody].nfaces; m++) {
      stat = EG_getTopology(bodydata[ibody].faces[m], &geom, &oclass, &mtype,
                            uvbox, &len, &loops, &lsenses);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 1) continue;
      stat = EG_attributeAdd(bodydata[ibody].faces[m], ".innerLoops", ATTRSTRING,
                             3, NULL, NULL, "Off");
      if (stat != EGADS_SUCCESS)
        printf(" EG_attributeAdd = %d\n", stat);
    }
#endif

    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
    /* make the quads */
    printf("\n");
    EG_setOutLevel(context, 0);
    for (i = 0 ; i < bodydata[ibody].nfaces; i++)
      EG_makeQuads(bodydata[ibody].tess, qparam, i+1);
    EG_setOutLevel(context, 1);
  }
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) printf(" failed to create wvContext!\n");

  /* make the scene */
  for (sum = stat = ibody = 0; ibody < nbody; ibody++) {

    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getQuads(bodydata[ibody].tess, i+1, &len, &xyzs, &uvs, &ptype,
                         &pindex, &npatch);
      if (stat != EGADS_SUCCESS) continue;

      /* got quads */
      if (npatch != 0) {
        for (nq = ip = 0; ip < npatch; ip++) {
          stat = EG_getPatch(bodydata[ibody].tess, i+1, ip+1, &n1, &n2,
                             &pvi, &pvb);
          if (stat != EGADS_SUCCESS) {
            printf(" Face %d: patch %d (%d) status = %d\n",
                   i+1, ip+1, npatch, stat);
            continue;
          }
          nq += (n1-1)*(n2-1);
        }
        printf(" face %d: npatch = %d, nquad = %d\n", i+1, npatch, nq);
        if (nq == 0) continue;

        snprintf(gpname, 34, "Body %d Face %d", ibody+1, i+1);
        stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
        wv_adjustVerts(&items[0], focus);
        trix = (int *) malloc(3*2*nq*sizeof(int));
        if (trix == NULL) continue;
        for (ntri = ip = 0; ip < npatch; ip++) {
          stat = EG_getPatch(bodydata[ibody].tess, i+1, ip+1, &n1, &n2,
                             &pvi, &pvb);
          if (stat != EGADS_SUCCESS) continue;
          for (k = 0; k < n2-1; k++)
            for (j = 0; j < n1-1; j++) {
              trix[3*ntri  ] = pvi[n1*k     + j  ];
              trix[3*ntri+1] = pvi[n1*k     + j+1];
              trix[3*ntri+2] = pvi[n1*(k+1) + j+1];
              ntri++;
              trix[3*ntri  ] = pvi[n1*k     + j  ];
              trix[3*ntri+1] = pvi[n1*(k+1) + j+1];
              trix[3*ntri+2] = pvi[n1*(k+1) + j  ];
              ntri++;
            }
        }
        stat = wv_setData(WV_INT32, 3*ntri, (void *) trix, WV_INDICES, &items[1]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
        free(trix);
        color[0]  = 1.0;
        color[1]  = ibody;
        color[1] /= nbody;
        color[2]  = 0.0;
        stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
        segs = (int *) malloc(2*4*nq*sizeof(int));
        if (segs == NULL) {
          printf(" Can not allocate %d Sides!\n", 2*4*nq);
          continue;
        }
        for (nseg = ip = 0; ip < npatch; ip++) {
          stat = EG_getPatch(bodydata[ibody].tess, i+1, ip+1, &n1, &n2,
                             &pvi, &pvb);
          if (stat != EGADS_SUCCESS) continue;
          for (k = 0; k < n2-1; k++)
            for (j = 0; j < n1-1; j++) {
              segs[2*nseg  ] = pvi[n1*k     + j  ];
              segs[2*nseg+1] = pvi[n1*k     + j+1];
              nseg++;
              segs[2*nseg  ] = pvi[n1*k     + j+1];
              segs[2*nseg+1] = pvi[n1*(k+1) + j+1];
              nseg++;
              segs[2*nseg  ] = pvi[n1*(k+1) + j+1];
              segs[2*nseg+1] = pvi[n1*(k+1) + j  ];
              nseg++;
              segs[2*nseg  ] = pvi[n1*(k+1) + j  ];
              segs[2*nseg+1] = pvi[n1*k     + j  ];
              nseg++;
            }
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
        continue;
      }
      /* plot triangles */
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

    /* get loops */
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTopology(bodydata[ibody].faces[i], &geom, &oclass,
                            &mtype, NULL, &nloops, &loops, &lsenses);
      if (stat != EGADS_SUCCESS) continue;
      for (nh = j = 0; j < nloops; j++) {
        stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, NULL,
                              &nledges, &ledges, &esenses);
        if (stat != EGADS_SUCCESS) continue;

        /* count */
        for (ntri = nseg = k = 0; k < nledges; k++) {
          m = 0;
          while (ledges[k] != bodydata[ibody].edges[m]) {
            m++;
            if (m == bodydata[ibody].nedges) break;
          }
          /* assume that the edge is degenerate and removed */
          if (m == bodydata[ibody].nedges) continue;
          stat = EG_getTessEdge(bodydata[ibody].tess, m+1, &len,
                                &xyzs, &ts);
          if (stat != EGADS_SUCCESS) {
            printf(" EG_getTessEdge %d = %d!\n", m+1, stat);
            nseg = 0;
            break;
          }
          if (len == 2)
            if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
                (xyzs[2] == xyzs[5])) continue;
          nh++;
#ifdef NONINDEXED
          nseg += 2*(len-1);
#else
          nseg += len;
          ntri += 2*(len-1);
#endif
        }
        if (nseg == 0) continue;
        lsegs = (float *) malloc(3*nseg*sizeof(float));
        if (lsegs == NULL) {
          printf(" Can not allocate %d Segments!\n", nseg);
          continue;
        }
#ifndef NONINDEXED
        segs = (int *) malloc(ntri*sizeof(int));
        if (segs == NULL) {
          printf(" Can not allocate %d Line Segments!\n", ntri);
          free(lsegs);
          continue;
        }
#endif
        heads = (int *) malloc(nh*sizeof(int));
        if (heads == NULL) {
          printf(" Can not allocate %d Heads!\n", nh);
#ifndef NONINDEXED
          free(segs);
#endif
          free(lsegs);
          continue;
        }

        /* fill */
        for (nh = ntri = nseg = k = 0; k < nledges; k++) {
          m = 0;
          while (ledges[k] != bodydata[ibody].edges[m]) {
            m++;
            if (m == bodydata[ibody].nedges) break;
          }
          /* assume that the edge is degenerate and removed */
          if (m == bodydata[ibody].nedges) continue;
          EG_getTessEdge(bodydata[ibody].tess, m+1, &len, &xyzs, &ts);
          if (len == 2)
            if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
                (xyzs[2] == xyzs[5])) continue;
#ifdef NONINDEXED
          if (esenses[k] == -1) heads[nh] = -nseg/2 - 1;
          for (n = 0; n < len-1; n++) {
            lsegs[3*nseg  ] = xyzs[3*n  ];
            lsegs[3*nseg+1] = xyzs[3*n+1];
            lsegs[3*nseg+2] = xyzs[3*n+2];
            nseg++;
            lsegs[3*nseg  ] = xyzs[3*n+3];
            lsegs[3*nseg+1] = xyzs[3*n+4];
            lsegs[3*nseg+2] = xyzs[3*n+5];
            nseg++;
          }
          if (esenses[k] ==  1) heads[nh] = nseg/2;
#else
          if (esenses[k] == -1) heads[nh] = -ntri/2 - 1;
          for (n = 0; n < len-1; n++) {
            segs[ntri] = n+nseg+1;
            ntri++;
            segs[ntri] = n+nseg+2;
            ntri++;
          }
          if (esenses[k] ==  1) heads[nh] = ntri/2;
          for (n = 0; n < len; n++) {
            lsegs[3*nseg  ] = xyzs[3*n  ];
            lsegs[3*nseg+1] = xyzs[3*n+1];
            lsegs[3*nseg+2] = xyzs[3*n+2];
            nseg++;
          }
#endif
          nh++;
        }
        snprintf(gpname, 46, "Body %d Loop %d/%d", ibody+1, i+1, j+1);
        stat = wv_setData(WV_REAL32, nseg, (void *) lsegs,  WV_VERTICES, &items[0]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
        wv_adjustVerts(&items[0], focus);
        free(lsegs);
        stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[1]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
#ifdef NONINDEXED
        stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 2, items);
#else
        stat = wv_setData(WV_INT32, ntri, (void *) segs, WV_INDICES, &items[2]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
        free(segs);
        stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
#endif
        if (stat < 0) {
          printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
        } else {
          if (cntxt != NULL)
            if (cntxt->gPrims != NULL) {
              cntxt->gPrims[stat].lWidth = 1.0;
              n = wv_addArrowHeads(cntxt, stat, 0.05, nh, heads);
              if (n != 0) printf(" wv_addArrowHeads = %d\n", n);
            }
        }
        free(heads);
/*      wv_printGPrim(cntxt, stat);  */
      }
    }
  }
  printf("\n ** %d gPrims with %d triangles **\n", stat+1, sum);

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
/*    wv_broadcastText("I'm here!");  */
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
