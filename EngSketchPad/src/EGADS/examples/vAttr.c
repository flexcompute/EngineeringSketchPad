/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Tessellation using wv w/ Attribute reporting
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

//#define DISJOINTQUADS
//#define REGULAR


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
  int          i, j, iBody, ient, stat, nattr, atype, alen, oclass, mtype;
  int          *senses;
  const int    *pints;
  char         tag[5];
  const char   *name, *pstr;
  const double *preals;
  ego          obj, geom, *objs;

  if (strncmp(text,"Picked: ", 8) == 0) {
    iBody = 0;
    sscanf(&text[12], "%d %s %d", &iBody, tag, &ient);
    if (iBody != 0) {
      printf(" Picked: iBody = %d, type = %s, index = %d\n", iBody, tag, ient);
      if (strcmp(tag,"Face") == 0) {
        obj = bodydata[iBody-1].faces[ient-1];
      } else {
        obj = bodydata[iBody-1].edges[ient-1];
      }
      stat = EG_getTopology(obj, &geom, &oclass, &mtype, NULL, &i, &objs,
                            &senses);
      if ((stat == EGADS_SUCCESS) && (geom != NULL))
        printf("         Geom type = %d\n", geom->mtype);
      nattr = 0;
      stat  = EG_attributeNum(obj, &nattr);
      if ((stat == EGADS_SUCCESS) && (nattr != 0)) {
        for (i = 1; i <= nattr; i++) {
          stat = EG_attributeGet(obj, i, &name, &atype, &alen,
                                 &pints, &preals, &pstr);
          if (stat != EGADS_SUCCESS) continue;
          printf("   %s: ", name);
          if ((atype == ATTRREAL) || (atype == ATTRCSYS)) {
            for (j = 0; j < alen; j++) printf("%lf ", preals[j]);
          } else if (atype == ATTRSTRING) {
            printf("%s", pstr);
          } else {
            for (j = 0; j < alen; j++) printf("%d ", pints[j]);
          }
          printf("\n");
        }
      }
    }
  }

}


int main(int argc, char *argv[])
{
  int          i, j, k, ibody, stat, oclass, mtype, mbody, len, ntri, sum;
  int          nattr, nseg, *segs, *senses;
  int          nbody, ngp, atype, alen, quad;
  const int    *tris, *tric, *ptype, *pindex, *ints;
  float        arg, focus[4], color[3];
  double       box[6], size, tol, params[3];
  const double *xyzs, *uvs, *ts, *reals;
  char         gpname[34], *startapp;
  const char   *OCCrev, *string, *name, *str;
  ego          context, model, geom, *bodies, *dum;
#ifdef DISJOINTQUADS
  ego          tess;
#endif
  wvData       items[5];
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = {  0,     3,     4,     2};

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vAttr filename [angle maxlen sag]\n\n");
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
  stat = EG_getTopology(model, &geom, &oclass, &mbody, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d %d\n", nbody, mbody);
  if (nbody < mbody) nbody = mbody;
  bodydata = (bodyData *) malloc(nbody*sizeof(bodyData));
  if (bodydata == NULL) {
    printf(" MALLOC Error on Body storage!\n");
    return 1;
  }

  params[0] =  0.025*size;
  params[1] =  0.001*size;
  params[2] = 15.0;
#ifdef REGULAR
  params[0] =  0.050*size;
  params[1] =  0.002*size;
  params[2] = 20.0;
#endif
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
  printf(" Reference size = %le\n", size);

  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
    bodydata[ibody].body  = NULL;
    bodydata[ibody].mtype = 0;
    bodydata[ibody].tess  = NULL;
    bodydata[ibody].faces = NULL;
    bodydata[ibody].edges = NULL;
    stat = EG_getTopology(bodies[ibody], &geom, &oclass, &mtype,
                          NULL, &j, &dum, &senses);
    if (stat != EGADS_SUCCESS) continue;
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    if (mtype == WIREBODY) {
      printf(" Body %2d:  Type = WireBody ", ibody+1);
    } else if (mtype == FACEBODY) {
      printf(" Body %2d:  Type = FaceBody ", ibody+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %2d:  Type = SheetBody", ibody+1);
    } else {
      printf(" Body %2d:  Type = SolidBody", ibody+1);
    }

    if (oclass == EBODY) {
      stat = EG_getBodyTopos(bodies[ibody], NULL, EFACE,
                             &bodydata[ibody].nfaces, &bodydata[ibody].faces);
      i    = EG_getBodyTopos(bodies[ibody], NULL, EEDGE,
                             &bodydata[ibody].nedges, &bodydata[ibody].edges);
      printf("  Effective Topology\n");
    } else {
      stat = EG_getBodyTopos(bodies[ibody], NULL, FACE,
                             &bodydata[ibody].nfaces, &bodydata[ibody].faces);
      i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE,
                             &bodydata[ibody].nedges, &bodydata[ibody].edges);
      stat = EG_tolerance(bodies[ibody], &tol);
      printf("  tol = %le\n", tol);
    }
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      return 1;
    }
    printf("           nFaces = %d   nEdges = %d\n", bodydata[ibody].nfaces,
           bodydata[ibody].nedges);

    /* do we have a tessellation? */
    if (mbody > 0) {
      for (i = 0; i < nbody; i++) {
        if (i == ibody) continue;
        if (bodies[i]->oclass != TESSELLATION) continue;
        stat = EG_statusTessBody(bodies[i], &geom, &j, &len);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_statusTessBody %d = %d\n", i+1, stat);
          continue;
        }
        if (geom != bodies[ibody]) continue;
        printf("           Found Tessellation %d for Body %d\n", i+1, ibody+1);
        stat = EG_copyObject(bodies[i], NULL, &bodydata[ibody].tess);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_copyObject %d = %d\n", i+1, stat);
          continue;
        } else {
          nattr = 0;
          stat  = EG_attributeNum(bodydata[ibody].tess, &nattr);
          if ((stat == EGADS_SUCCESS) && (nattr != 0)) {
            for (k = 1; k <= nattr; k++) {
              stat = EG_attributeGet(bodydata[ibody].tess, k, &name, &atype,
                                     &alen, &ints, &reals, &str);
              if (stat != EGADS_SUCCESS) continue;
              printf("           %s: ", name);
              if ((atype == ATTRREAL) || (atype == ATTRCSYS)) {
                for (j = 0; j < alen; j++) printf("%lf ", reals[j]);
              } else if (atype == ATTRSTRING) {
                printf("%s", str);
              } else {
                for (j = 0; j < alen; j++) printf("%d ", ints[j]);
              }
              printf("\n");
            }
          }
        }
        break;
      }
    }
    if (bodydata[ibody].tess == NULL) {
      if (mbody > 0)
        printf("           Tessellating Body %d\n", ibody+1);
      stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_makeTessBody %d = %d\n", ibody, stat);
        continue;
      }
    }
#ifdef DISJOINTQUADS
    tess = bodydata[ibody].tess;
#ifndef REGULAR
    /* disable regularization in EGADS */
    stat = EG_attributeAdd(tess, ".qRegular", ATTRSTRING, 3, NULL, NULL, "Off");
    if (stat != EGADS_SUCCESS)
      printf(" EG_attributeAdd qRegular %d = %d\n", ibody, stat);
#endif
    stat = EG_quadTess(tess, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_quadTess %d = %d  -- reverting...\n", ibody, stat);
      bodydata[ibody].tess = tess;
      continue;
    }
    EG_deleteObject(tess);
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
  for (ngp = sum = stat = ibody = 0; ibody < nbody; ibody++) {
    if (bodydata[ibody].tess == NULL) continue;
    
    quad = 0;
    stat = EG_attributeRet(bodydata[ibody].tess, ".tessType", &atype,
                           &alen, &ints, &reals, &string);
    if (stat == EGADS_SUCCESS)
      if (atype == ATTRSTRING)
        if (strcmp(string, "Quad") == 0) quad = 1;

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
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) nseg++;
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) nseg++;
      }
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) {
              segs[2*nseg  ] = tris[3*j+sides[k][0]];
              segs[2*nseg+1] = tris[3*j+sides[k][1]];
              nseg++;
            }
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) {
              segs[2*nseg  ] = tris[6*j+sideq[k][0]];
              segs[2*nseg+1] = tris[6*j+sideq[k][1]];
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
      if (stat > 0) ngp = stat+1;
      sum += ntri;
    }

    /* get edges */
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      stat  = EG_getTessEdge(bodydata[ibody].tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
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
        ngp = stat+1;
      }
    }
  }
  printf(" ** %d gPrims with %d triangles **\n", ngp, sum);

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
    if (bodydata[ibody].tess == NULL) continue;
    EG_deleteObject(bodydata[ibody].tess);
    EG_free(bodydata[ibody].edges);
    EG_free(bodydata[ibody].faces);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
