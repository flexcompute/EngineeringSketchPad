/**
 * ESP client of libxddm. For XDDM documentation, see
 * $CART3D/doc/xddm/xddm.html. The library uses XML Path Language (XPath) to
 * navigate the elements of XDDM documents. For XPath tutorials, see the web,
 * e.g. www.developer.com/net/net/article.php/3383961/NET-and-XML-XPath-Queries.htm
 *
 * Dependency: libxml2, www.xmlsoft.org. This library is usually present on
 * most systems, check existence of 'xml2-config' script
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "writeTrix.h"

#include "bodyTess.h"

/* OpenCSM Defines & Includes */
#include "common.h"
#include "OpenCSM.h"

/* external functions not declared in OpenCSM.h */
extern int velocityOfNode(modl_T *modl, int ibody, int inode, double dxyz[]);
extern int velocityOfEdge(modl_T *modl, int ibody, int iedge, int npnt,
                          /*@null@*/double t[], double dxyz[]);
extern int velocityOfFace(modl_T *modl, int ibody, int iface, int npnt,
                          /*@null@*/double uv[], double dxyz[]);



static int
parsePmtr(char *name, char *pname, int *irow, int *icol)
{
  int i, open, len;
  
  *irow = 1;
  *icol = 1;
  
  open = -1;
  len  = strlen(name);
  for (i = 0; i <= len; i++) {
    if ((name[i] == '~') || (name[i] == '[')) {
      open = i;
      break;
    }
    pname[i] = name[i];
  }
  if (open == -1) return 0;
  
  pname[open] = 0;
  sscanf(&name[open+1], "%d,%d", irow, icol);
  return 1;
}


int
main(int argc, char *argv[])
{
  int          i, j, k, m, ipmtr, stat, major, minor, type, nrow, ncol, opts;
  int          irow, icol, buildTo, builtTo, ibody, nvert, ntri, noI;
  int          nbody, nface, nedge, plen, tlen, kb, nb, *tris;
  char         *filename, name[MAX_NAME_LEN], pname[129];
  double       value, dot, size, lower, upper;
  double       global[3], tparam[3], box[6], **dvar, *xyzs, *psens;
  const char   *occ_rev;
  const int    *tri, *tric, *ptype, *pindex;
  const double *points, *uv, *pcsens;
  ego          context, body;
  void         *modl;
  modl_T       *MODL;
  verTags      *vtags;
  p_tsXddm     p_xddm = NULL;
  p_tsXmParent p_p;
  
  if (argc != 3) {
    printf("Uasge: ESPxddm <xddm_filename> <xpath_expression>\n\n");
    return 1;
  }
  
  ocsmVersion(&major, &minor);
  printf("\n Using OpenCSM %d.%02d\n", major, minor);
  EG_revision(&major, &minor, &occ_rev);
  printf(" Using EGADS   %d.%02d %s\n\n", major, minor, occ_rev);
  
  opts = 1;  /* be verbose */
  
  p_xddm = xddm_readFile(argv[1], argv[2], opts);
  if (NULL == p_xddm) {
    printf("xddm_readFile failed to parse\n");
    return 1;
  }
  xddm_echo(p_xddm);
  
  /* get filename -- attached to ID */
  filename = NULL;
  if (p_xddm->p_parent) {
    for (i = 0; i < p_xddm->p_parent->nAttr; i++) {
      p_p = p_xddm->p_parent;
      if (strcmp(p_p->p_attr[i].p_name, "ID") == 0) {
        filename = p_p->p_attr[i].p_value;
        break;
      }
    }
  }
  if (filename == NULL) {
    printf("ID not found!\n");
    xddm_free(p_xddm);
    return 1;
  }
  
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS failed to Open: %d\n", stat);
    xddm_free(p_xddm);
    return 1;
  }

  /* do an OpenCSM load */
  stat = ocsmLoad(filename, &modl);
  if (stat < SUCCESS) {
    printf(" ocsmLoad failed: %d\n", stat);
    EG_close(context);
    xddm_free(p_xddm);
    return 1;
  }
  MODL = (modl_T *) modl;
  MODL->context   = context;
  MODL->tessAtEnd = 0;

  /* check that Branches are properly ordered */
  stat = ocsmCheck(modl);
  if (stat < SUCCESS) {
    printf(" ocsmCheck failed: %d\n", stat);
    ocsmFree(modl);
    ocsmFree(NULL);
    EG_close(context);
    xddm_free(p_xddm);
    return 1;
  }

  printf("\n");
  for (i = 0; i < p_xddm->nv; i++) {
    if (p_xddm->a_v[i].p_id == NULL)  continue;
    if (p_xddm->a_v[i].val  == UNSET) continue;
    noI = parsePmtr(p_xddm->a_v[i].p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      if (type != OCSM_DESPMTR) continue;
      if (strcmp(pname, name) == 0) {
        ipmtr = j+1;
        break;
      }
    }
    if (ipmtr != 0) {
      if ((noI == 0) && ((ncol > 1) || (nrow > 1))) {
        printf(" Variable %s not indexed!\n", pname);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      if ((irow < 1) || (irow > nrow) || (icol < 1) || (icol > ncol)) {
        printf(" Variable %s not in range [%d,%d]!\n",
               p_xddm->a_v[i].p_id, nrow, ncol);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      stat = ocsmGetValu( modl, ipmtr, irow, icol, &value, &dot);
      stat = ocsmSetValuD(modl, ipmtr, irow, icol, p_xddm->a_v[i].val);
      printf(" Setting Variable %s from %lf to %lf  stat = %d\n",
             p_xddm->a_v[i].p_id, value, p_xddm->a_v[i].val, stat);
      stat = ocsmGetBnds(modl, ipmtr, irow, icol, &lower, &upper);
      if (stat == SUCCESS) {
        if (lower != -HUGEQ) p_xddm->a_v[i].minVal = lower;
        if (upper !=  HUGEQ) p_xddm->a_v[i].maxVal = upper;
      }
    } else {
      printf(" Variable %s not found!\n", pname);
      ocsmFree(modl);
      ocsmFree(NULL);
      EG_close(context);
      xddm_free(p_xddm);
      return 1;
    }
  }
  
  printf("\n");
  for (i = 0; i < p_xddm->nc; i++) {
    if (p_xddm->a_c[i].p_id == NULL)  continue;
    if (p_xddm->a_c[i].val  == UNSET) continue;
    noI = parsePmtr(p_xddm->a_c[i].p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      if (type != OCSM_DESPMTR) continue;
      if (strcmp(pname, name) == 0) {
        ipmtr = j+1;
        break;
      }
    }
    if (ipmtr != 0) {
      if ((noI == 0) && ((ncol > 1) || (nrow > 1))) {
        printf(" Constant %s not indexed!\n", pname);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      if ((irow < 1) || (irow > nrow) || (icol < 1) || (icol > ncol)) {
        printf(" Constant %s not in range [%d,%d]!\n",
               p_xddm->a_c[i].p_id, nrow, ncol);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      stat = ocsmGetValu( modl, ipmtr, irow, icol, &value, &dot);
      stat = ocsmSetValuD(modl, ipmtr, irow, icol, p_xddm->a_c[i].val);
      printf(" Setting Constant %s from %lf to %lf  stat = %d\n",
             p_xddm->a_c[i].p_id, value, p_xddm->a_c[i].val, stat);
    } else {
      printf(" Constant %s not found!\n", pname);
      ocsmFree(modl);
      ocsmFree(NULL);
      EG_close(context);
      xddm_free(p_xddm);
      return 1;
    }
  }
  printf("\n");

  buildTo = 0;                          /* all */
  nbody   = 0;
  stat    = ocsmBuild(modl, buildTo, &builtTo, &nbody, NULL);
  EG_deleteObject(context);             /* clean up after build */
  if (stat != SUCCESS) {
    printf(" ocsmBuild failed: %d\n", stat);
    ocsmFree(modl);
    ocsmFree(NULL);
    EG_close(context);
    xddm_free(p_xddm);
    return 1;
  }
  nbody = 0;
  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    nbody++;
  }
  printf("\n nBody = %d\n\n", nbody);
  if (nbody <= 0) {
    ocsmFree(modl);
    ocsmFree(NULL);
    EG_close(context);
    xddm_free(p_xddm);
    return 1;
  }
  /* NOTE: currently will only work for a single body! */
  if (nbody != 1) {
    printf(" ESPxddm: currently functions for a single body!\n");
    ocsmFree(modl);
    ocsmFree(NULL);
    EG_close(context);
    xddm_free(p_xddm);
    return 1;
  }
  
  /* set analysis values */
  
  for (i = 0; i < p_xddm->na; i++) {
    if (p_xddm->a_ap[i].p_id == NULL)  continue;
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      if (name[0] != '@') continue;
      if (strcmp(p_xddm->a_ap[i].p_id, &name[1]) == 0) {
        ipmtr = j+1;
        break;
      }
    }
    if (ipmtr != 0) {
      stat = ocsmGetValu(modl, ipmtr, 1, 1, &value, &dot);
      p_xddm->a_ap[i].val = value;
      printf(" Setting Analysis Parameter %s to %lf\n", p_xddm->a_ap[i].p_id,
             value);
    } else {
      printf(" Analysis Parameter %s not found!\n", p_xddm->a_ap[i].p_id);
      ocsmFree(modl);
      ocsmFree(NULL);
      EG_close(context);
      xddm_free(p_xddm);
      return 1;
    }
  }
  printf("\n");
  
  /* get global tesselate parameters */
  global[0] = global[1] = global[2] = 0.0;
  for (i = 0; i < p_xddm->ntess; i++) {
    if (p_xddm->a_tess[i].p_id != NULL) continue;
    for (j = 0; j < p_xddm->a_tess[i].nAttr; j++)
      if (strcmp(p_xddm->a_tess[i].p_attr[j].p_name, "MaxEdge") == 0) {
        sscanf(p_xddm->a_tess[i].p_attr[j].p_value, "%lf", &global[0]);
      } else if (strcmp(p_xddm->a_tess[i].p_attr[j].p_name, "Sag") == 0) {
        sscanf(p_xddm->a_tess[i].p_attr[j].p_value, "%lf", &global[1]);
      } else if (strcmp(p_xddm->a_tess[i].p_attr[j].p_name, "Angle") == 0) {
        sscanf(p_xddm->a_tess[i].p_attr[j].p_value, "%lf", &global[2]);
      } else {
        printf(" Tessellation (global) Attribute %s not Understood!\n",
               p_xddm->a_tess[i].p_attr[j].p_name);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
  }
  if ((global[0] == 0.0) && (global[1] == 0.0) && (global[2] == 0.0)) {
    global[0] =  0.025;
    global[1] =  0.001;
    global[2] = 12.0;
  }

  /* tessellate */
  kb = 0;
  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    kb++;
    body      = MODL->body[ibody].ebody;
    tparam[0] = global[0];
    tparam[1] = global[1];
    tparam[2] = global[2];
    for (i = 0; i < p_xddm->ntess; i++) {
      if (p_xddm->a_tess[i].p_id == NULL) continue;
      sscanf(p_xddm->a_tess[i].p_id, "%d", &stat);
      if (stat != kb) continue;
      for (j = 0; j < p_xddm->a_tess[i].nAttr; j++)
        if (strcmp(p_xddm->a_tess[i].p_attr[j].p_name, "MaxEdge") == 0) {
          sscanf(p_xddm->a_tess[i].p_attr[j].p_value, "%lf", &tparam[0]);
        } else if (strcmp(p_xddm->a_tess[i].p_attr[j].p_name, "Sag") == 0) {
          sscanf(p_xddm->a_tess[i].p_attr[j].p_value, "%lf", &tparam[1]);
        } else if (strcmp(p_xddm->a_tess[i].p_attr[j].p_name, "Angle") == 0) {
          sscanf(p_xddm->a_tess[i].p_attr[j].p_value, "%lf", &tparam[2]);
        } else {
          printf(" Tessellation (ID=%d) Attribute %s not Understood!\n",
                 kb, p_xddm->a_tess[i].p_attr[j].p_name);
          ocsmFree(modl);
          ocsmFree(NULL);
          EG_close(context);
          xddm_free(p_xddm);
          return 1;
        }
    }
    printf(" Tessellating %d with  MaxEdge = %lf   Sag = %lf   Angle = %lf\n",
           kb, tparam[0], tparam[1], tparam[2]);
    stat = EG_getBoundingBox(body, box);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getBoundingBox failed: %d!\n", stat);
      ocsmFree(modl);
      ocsmFree(NULL);
      EG_close(context);
      xddm_free(p_xddm);
      return 1;
    }
    size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                (box[4]-box[1])*(box[4]-box[1]) +
                (box[5]-box[2])*(box[5]-box[2]));
    tparam[0] *= size;
    tparam[1] *= size;
    stat = EG_makeTessBody(body, tparam, &MODL->body[ibody].etess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody failed: %d!\n", stat);
      ocsmFree(modl);
      ocsmFree(NULL);
      EG_close(context);
      xddm_free(p_xddm);
      return 1;
    }
    stat = bodyTess(MODL->body[ibody].etess, &nface, &nedge, &nvert, &xyzs,
                    &vtags, &ntri, &tris);
    if (stat != EGADS_SUCCESS) {
      printf(" bodyTess failed: %d!\n", stat);
      ocsmFree(modl);
      ocsmFree(NULL);
      EG_close(context);
      xddm_free(p_xddm);
      return 1;
    }
    printf("    nVerts = %d   nTris = %d\n", nvert, ntri);
    
    dvar = NULL;
    if (p_xddm->nv > 0) {
      dvar = (double **) EG_alloc(p_xddm->nv*sizeof(double *));
      if (dvar == NULL) {
        printf(" sensitivity allocation failed!\n");
        EG_free(tris);
        EG_free(vtags);
        EG_free(xyzs);
        ocsmFree(modl);
        ocsmFree(NULL);
        EG_close(context);
        xddm_free(p_xddm);
        return 1;
      }
      for (i = 0; i < p_xddm->nv; i++) {
        dvar[i] = EG_alloc(3*nvert*sizeof(double));
        if (dvar[i] == NULL) {
          printf(" sensitivity %d allocation failed!\n", i+1);
          for (j = 0; j < i; j++) EG_free(dvar[j]);
          EG_free(dvar);
          EG_free(tris);
          EG_free(vtags);
          EG_free(xyzs);
          ocsmFree(modl);
          ocsmFree(NULL);
          EG_close(context);
          xddm_free(p_xddm);
          return 1;
        }
        for (j = 0; j < 3*nvert; j++) dvar[i][j] = 0.0;
      }
    }
    if (dvar == NULL) goto write;
    
    /* compute the sensitivities */
    printf("\n");
    for (i = 0; i < p_xddm->nv; i++) {
      if (p_xddm->a_v[i].p_id == NULL)  continue;
      if (p_xddm->a_v[i].val  == UNSET) continue;
      noI = parsePmtr(p_xddm->a_v[i].p_id, pname, &irow, &icol);
      for (ipmtr = j = 0; j < MODL->npmtr; j++) {
        stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
        if (strcmp(pname, name) == 0) {
          ipmtr = j+1;
          break;
        }
      }
      /* clear all then set */
      ocsmSetDtime(modl, 0);
      ocsmSetVelD(modl, 0,     0,    0,    0.0);
      ocsmSetVelD(modl, ipmtr, irow, icol, 1.0);
      buildTo = 0;
      nb      = 0;
      stat    = ocsmBuild(modl, buildTo, &builtTo, &nb, NULL);
      if (noI == 0) {
        printf("\n*** compute parameter %d (%s) sensitivity status = %d (%d)***\n",
               ipmtr, p_xddm->a_v[i].p_id, stat, nb);
      } else {
        printf("\n*** compute parameter %d [%d,%d] (%s) sensitivity = %d (%d)***\n",
               ipmtr, irow, icol, p_xddm->a_v[i].p_id, stat, nb);
      }
      if (MODL->perturb == NULL) {
        if (p_xddm->a_v[i].p_comment != NULL) {
          if (strcmp(p_xddm->a_v[i].p_comment, "FD") == 0) {
            stat = ocsmSetDtime(MODL, 0.001);
            printf("\n*** forced finite differencing for %s (%d) ***\n",
                   pname, stat);
          }
          if (strcmp(p_xddm->a_v[i].p_comment, "oFD") == 0) {
            stat = ocsmSetDtime(MODL, 0.001);
            printf("\n*** forcing OpenCSM finite differencing for %s (%d) ***\n",
                   pname, stat);
          }
        }
      }
      
      /* do we need to do a full vertex-based Finite Difference? */
      if (MODL->perturb != NULL)
        if (p_xddm->a_v[i].p_comment != NULL)
          if (strcmp(p_xddm->a_v[i].p_comment, "FD") == 0) {
            /* fill in the Nodes */
            for (j = 0; j < nvert; j++) {
              if (vtags[j].ptype != 0) continue;
              stat = velocityOfNode(modl, ibody, vtags[j].pindex, &dvar[i][3*j]);
              if (stat != EGADS_SUCCESS) {
                printf(" velocityOfNode Parameter %d vert %d failed: %d  (Node = %d)!\n",
                       i+1, j+1, stat, vtags[j].pindex);
                for (k = 0; k < p_xddm->nv; k++) EG_free(dvar[k]);
                EG_free(dvar);
                EG_free(tris);
                EG_free(vtags);
                EG_free(xyzs);
                ocsmFree(modl);
                ocsmFree(NULL);
                EG_close(context);
                xddm_free(p_xddm);
                return 1;
              }
/*            printf(" Node %d: %lf %lf %lf\n", vtags[j].pindex, dvar[i][3*j  ],
                     dvar[i][3*j+1], dvar[i][3*j+2]);  */
            }
            /* do the edges */
            for (j = 1; j <= nedge; j++) {
              EG_getTessEdge(MODL->body[ibody].etess, j, &plen, &points, &uv);
              if (plen <= 2) continue;
              psens = (double *) EG_alloc(3*plen*sizeof(double));
              if (psens == NULL) {
                printf(" Can not get sens vec for Edge %d (%d)!\n", j, plen);
                continue;
              }
              stat = velocityOfEdge(modl, ibody, j, plen, NULL, psens);
              if (stat != EGADS_SUCCESS) {
                printf(" velocityOfEdge Parameter %d Edge %d failed: %d!\n",
                       i+1, j, stat);
                EG_free(psens);
                for (k = 0; k < p_xddm->nv; k++) EG_free(dvar[k]);
                EG_free(dvar);
                EG_free(tris);
                EG_free(vtags);
                EG_free(xyzs);
                ocsmFree(modl);
                ocsmFree(NULL);
                EG_close(context);
                xddm_free(p_xddm);
                return 1;
              }
              for (k = 0; k < nvert; k++)
                if ((vtags[k].ptype > 0) && (vtags[k].pindex == j)) {
                  m              = vtags[k].ptype - 1;
                  dvar[i][3*k  ] = psens[3*m  ];
                  dvar[i][3*k+1] = psens[3*m+1];
                  dvar[i][3*k+2] = psens[3*m+2];
/*                printf(" Edge %d: %d %lf %lf %lf\n", j, m,
                         psens[3*m  ], psens[3*m+1], psens[3*m+2]);  */
                }
              EG_free(psens);
            }
            /* do the faces */
            for (j = 1; j <= nface; j++) {
              EG_getTessFace(MODL->body[ibody].etess, j, &plen, &points, &uv,
                             &ptype, &pindex, &tlen, &tri, &tric);
              if (plen <= 2) continue;
              psens = (double *) EG_alloc(3*plen*sizeof(double));
              if (psens == NULL) {
                printf(" Can not get sens vec for Face %d (%d)!\n", j, plen);
                continue;
              }
              stat = velocityOfFace(modl, ibody, j, plen, NULL, psens);
              if (stat != EGADS_SUCCESS) {
                printf(" velocityOfFace Parameter %d Face %d failed: %d!\n",
                       i+1, j, stat);
                EG_free(psens);
                for (k = 0; k < p_xddm->nv; k++) EG_free(dvar[k]);
                EG_free(dvar);
                EG_free(tris);
                EG_free(vtags);
                EG_free(xyzs);
                ocsmFree(modl);
                ocsmFree(NULL);
                EG_close(context);
                xddm_free(p_xddm);
                return 1;
              }
              for (k = 0; k < nvert; k++)
                if ((vtags[k].ptype < 0) && (vtags[k].pindex == j)) {
                  m              = -vtags[k].ptype - 1;
                  dvar[i][3*k  ] = psens[3*m  ];
                  dvar[i][3*k+1] = psens[3*m+1];
                  dvar[i][3*k+2] = psens[3*m+2];
                }
              EG_free(psens);
            }
            printf("\n");
            continue;
          }
  
      /* fill in the Nodes */
      for (j = 0; j < nvert; j++) {
        if (vtags[j].ptype != 0) continue;
        stat = ocsmGetTessVel(modl, ibody, OCSM_NODE, vtags[j].pindex, &pcsens);
        if (stat != EGADS_SUCCESS) {
          printf(" ocsmGetTessVel Parameter %d vert %d failed: %d  (Node = %d)!\n",
                 i+1, j+1, stat, vtags[j].pindex);
          for (k = 0; k < p_xddm->nv; k++) EG_free(dvar[k]);
          EG_free(dvar);
          EG_free(tris);
          EG_free(vtags);
          EG_free(xyzs);
          ocsmFree(modl);
          ocsmFree(NULL);
          EG_close(context);
          xddm_free(p_xddm);
          return 1;
        }
/*      printf(" Node %d: %lf %lf %lf\n", vtags[j].pindex, pcsens[0], pcsens[1],
               pcsens[2]);  */
        dvar[i][3*j  ] = pcsens[0];
        dvar[i][3*j+1] = pcsens[1];
        dvar[i][3*j+2] = pcsens[2];
      }

      /* next do all of the edges */
      for (j = 1; j <= nedge; j++) {
        stat = ocsmGetTessVel(modl, ibody, OCSM_EDGE, j, &pcsens);
        if (stat != EGADS_SUCCESS) {
          printf(" ocsmGetTessVel Parameter %d Edge %d failed: %d!\n",
                 i+1, j, stat);
          for (k = 0; k < p_xddm->nv; k++) EG_free(dvar[k]);
          EG_free(dvar);
          EG_free(tris);
          EG_free(vtags);
          EG_free(xyzs);
          ocsmFree(modl);
          ocsmFree(NULL);
          EG_close(context);
          xddm_free(p_xddm);
          return 1;
        }
        for (k = 0; k < nvert; k++)
          if ((vtags[k].ptype > 0) && (vtags[k].pindex == j)) {
            m              = vtags[k].ptype - 1;
            dvar[i][3*k  ] = pcsens[3*m  ];
            dvar[i][3*k+1] = pcsens[3*m+1];
            dvar[i][3*k+2] = pcsens[3*m+2];
/*          printf(" Edge %d: %d %lf %lf %lf\n", j, m,
                   pcsens[3*m  ], pcsens[3*m+1], pcsens[3*m+2]);  */
          }
      }
      
      /* do all of the faces */
      for (j = 1; j <= nface; j++) {
        stat = ocsmGetTessVel(modl, ibody, OCSM_FACE, j, &pcsens);
        if (stat != EGADS_SUCCESS) {
          printf(" ocsmGetTessVel Parameter %d Face %d failed: %d!\n",
                 i+1, j, stat);
          for (k = 0; k < p_xddm->nv; k++) EG_free(dvar[k]);
          EG_free(dvar);
          EG_free(tris);
          EG_free(vtags);
          EG_free(xyzs);
          ocsmFree(modl);
          ocsmFree(NULL);
          EG_close(context);
          xddm_free(p_xddm);
          return 1;
        }
        for (k = 0; k < nvert; k++)
          if ((vtags[k].ptype < 0) && (vtags[k].pindex == j)) {
            m              = -vtags[k].ptype - 1;
            dvar[i][3*k  ] = pcsens[3*m  ];
            dvar[i][3*k+1] = pcsens[3*m+1];
            dvar[i][3*k+2] = pcsens[3*m+2];
          }
      }

      printf("\n");
    }
    
write:
    /* write trix files -- Note: move out of body loop for multiple bodies */
    stat = writeTrix("Components.i.tri", p_xddm, kb, nvert, xyzs, p_xddm->nv,
                     dvar, ntri, tris);
    if (dvar != NULL) {
      for (j = 0; j < p_xddm->nv; j++) EG_free(dvar[j]);
      EG_free(dvar);
    }
    EG_free(tris);
    EG_free(vtags);
    EG_free(xyzs);
  }
  printf("\n");

  xddm_updateAnalysisParams(argv[1], p_xddm, opts);
  
  ocsmFree(modl);
  ocsmFree(NULL);
  EG_close(context);
  xddm_free(p_xddm);

  return 0;
}
