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

#include "xddmInternals.h"


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
  int          stat = EGADS_SUCCESS;
  int          i, j, k, ipmtr, major, minor, type, nrow, ncol;
  unsigned     opts;
  int          irow, icol, buildTo, builtTo, ibody, nvert, noI, ntri;
  int          nbody=0, nface, kb, nb, iglobal, iface, outLevel;
  char         *filename, name[MAX_NAME_LEN], pname[129];
  double       value, dot, size, lower, upper, step;
  double       global[3], tparam[3], box[6], ***dvar=NULL, *psens=NULL;
  const char   *occ_rev;
  const double *pcsens, *xyz, *uv;
  const int    *ptype, *pindx, *tris, *tric;
  ego          context=NULL, body, *tess=NULL;
  void         *modl = NULL;
  modl_T       *MODL = NULL;
  verTags      *vtags = NULL;
  p_tsXddm     p_xddm = NULL;
  
  if (argc != 3) {
    printf("Usage: ESPxddm <xddm_filename> <xpath_expression>\n\n");
    return 1;
  }
  
  ocsmVersion(&major, &minor);
  printf("\n Using OpenCSM %d.%02d\n", major, minor);
  EG_revision(&major, &minor, &occ_rev);
  printf(" Using EGADS   %d.%02d %s\n\n", major, minor, occ_rev);
  
  opts = 1;  /* be verbose */
  
  p_xddm = xddm_readFile(argv[1], argv[2], &opts);
  if (NULL == p_xddm) {
    printf("xddm_readFile failed to parse\n");
    goto cleanup;
  }
  xddm_echo(p_xddm, 2); // 2 is indent
  
  /* get filename -- attached to ID */
  filename = p_xddm->p_e->p_id;
  if (filename == NULL) {
    printf("ID not found!\n");
    stat = 1;
    goto cleanup;
  }
  
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS failed to Open: %d\n", stat);
    goto cleanup;
  }

  /* do an OpenCSM load */
  stat = ocsmLoad(filename, &modl);
  if (stat < SUCCESS || modl == NULL) {
    printf(" ocsmLoad failed: %d\n", stat);
    goto cleanup;
  }
  MODL = (modl_T *) modl;
  MODL->context   = context;
  MODL->tessAtEnd = 0;

  /* check that Branches are properly ordered */
  stat = ocsmCheck(modl);
  if (stat < SUCCESS) {
    printf(" ocsmCheck failed: %d\n", stat);
    goto cleanup;
  }

  printf("\n");
  for (i = 0; i < p_xddm->nv; i++) {
    if (p_xddm->a_v[i].p_e->p_id == NULL)  continue;
    if (p_xddm->a_v[i].val  == UNSET) continue;
    noI = parsePmtr(p_xddm->a_v[i].p_e->p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        goto cleanup;
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
        goto cleanup;
      }
      if ((irow < 1) || (irow > nrow) || (icol < 1) || (icol > ncol)) {
        printf(" Variable %s not in range [%d,%d]!\n",
               p_xddm->a_v[i].p_e->p_id, nrow, ncol);
        goto cleanup;
      }
      stat = ocsmGetValu( modl, ipmtr, irow, icol, &value, &dot);
      stat = ocsmSetValuD(modl, ipmtr, irow, icol, p_xddm->a_v[i].val);
      printf(" Setting Variable %s from %lf to %lf  stat = %d\n",
             p_xddm->a_v[i].p_e->p_id, value, p_xddm->a_v[i].val, stat);
      stat = ocsmGetBnds(modl, ipmtr, irow, icol, &lower, &upper);
      if (stat == SUCCESS) {
        if (lower != -HUGEQ) p_xddm->a_v[i].minVal = lower;
        if (upper !=  HUGEQ) p_xddm->a_v[i].maxVal = upper;
      }
    } else {
      printf(" Variable %s not found!\n", pname);
      stat = EGADS_NOTFOUND;
      goto cleanup;
    }
  }
  
  printf("\n");
  for (i = 0; i < p_xddm->nc; i++) {
    if (p_xddm->a_c[i].p_e->p_id == NULL)  continue;
    if (p_xddm->a_c[i].val  == UNSET) continue;
    noI = parsePmtr(p_xddm->a_c[i].p_e->p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        goto cleanup;
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
        goto cleanup;
      }
      if ((irow < 1) || (irow > nrow) || (icol < 1) || (icol > ncol)) {
        printf(" Constant %s not in range [%d,%d]!\n",
               p_xddm->a_c[i].p_e->p_id, nrow, ncol);
        goto cleanup;
      }
      stat = ocsmGetValu( modl, ipmtr, irow, icol, &value, &dot);
      stat = ocsmSetValuD(modl, ipmtr, irow, icol, p_xddm->a_c[i].val);
      printf(" Setting Constant %s from %lf to %lf  stat = %d\n",
             p_xddm->a_c[i].p_e->p_id, value, p_xddm->a_c[i].val, stat);
    } else {
      printf(" Constant %s not found!\n", pname);
      stat = EGADS_NOTFOUND;
      goto cleanup;
    }
  }
  
  printf("\n");
  for (i = 0; i < p_xddm->nf; i++) {
    if (p_xddm->a_f[i].p_e->p_id == NULL)  continue;
    if (p_xddm->a_f[i].val  == UNSET) continue;
    noI = parsePmtr(p_xddm->a_f[i].p_e->p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        goto cleanup;
      }
      if (type != OCSM_DESPMTR) continue;
      if (strcmp(pname, name) == 0) {
        ipmtr = j+1;
        break;
      }
    }
    if (ipmtr != 0) {
      if ((noI == 0) && ((ncol > 1) || (nrow > 1))) {
        printf(" Function %s not indexed!\n", pname);
        goto cleanup;
      }
      if ((irow < 1) || (irow > nrow) || (icol < 1) || (icol > ncol)) {
        printf(" Function %s not in range [%d,%d]!\n",
               p_xddm->a_f[i].p_e->p_id, nrow, ncol);
        goto cleanup;
      }
      stat = ocsmGetValu( modl, ipmtr, irow, icol, &value, &dot);
      stat = ocsmSetValuD(modl, ipmtr, irow, icol, p_xddm->a_f[i].val);
      printf(" Setting Function %s from %lf to %lf  stat = %d\n",
             p_xddm->a_f[i].p_e->p_id, value, p_xddm->a_f[i].val, stat);
    } else {
      printf(" Function %s not found!\n", pname);
      stat = EGADS_NOTFOUND;
      goto cleanup;
    }
  }
  printf("\n");

  buildTo = 0;                          /* all */
  nbody   = 0;
  stat    = ocsmBuild(modl, buildTo, &builtTo, &nbody, NULL);
  if (stat != SUCCESS) goto cleanup;

  nbody = 0;
  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    nbody++;
  }
  printf("\n nBody = %d\n\n", nbody);
  if (nbody <= 0) {
    stat = 1;
    goto cleanup;
  }
  tess = (ego*)EG_alloc(nbody*sizeof(ego));
  if (tess == NULL) {
    printf(" Memory allocation failed!\n");
    stat = 1;
    goto cleanup;
  }
  for (i = 0; i < nbody; i++) tess[i] = NULL;
  
  /* set analysis values */
  
  for (i = 0; i < p_xddm->na; i++) {
    if (p_xddm->a_ap[i].p_e->p_id == NULL)  continue;
    noI = parsePmtr(p_xddm->a_ap[i].p_e->p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (stat != SUCCESS) {
        printf(" ocsmGetPmtr %d failed: %d\n", ipmtr, stat);
        goto cleanup;
      }
      if (type != OCSM_OUTPMTR) continue;
//      if (name[0] != '@') continue;
//      if (strcmp(p_xddm->a_ap[i].p_e->p_id, &name[1]) == 0) {
      if (name[0] == '@') {
        char* substr = name + 1;
        memmove(name, substr, strlen(substr) + 1);
      }
      if (strcmp(pname, name) == 0) {
        ipmtr = j+1;
        break;
      }
    }
    if (ipmtr != 0) {
      if ((noI == 0) && ((ncol > 1) || (nrow > 1))) {
        printf(" Analysis Parameter %s not indexed!\n", pname);
        goto cleanup;
      }
      if ((irow < 1) || (irow > nrow) || (icol < 1) || (icol > ncol)) {
        printf(" Analysis Parameter %s not in range [%d,%d]!\n",
               p_xddm->a_v[i].p_e->p_id, nrow, ncol);
        goto cleanup;
      }
      stat = ocsmGetValu(modl, ipmtr, irow, icol, &value, &dot);
      p_xddm->a_ap[i].val = value;
      printf(" Setting Analysis Parameter %s to %lf\n", p_xddm->a_ap[i].p_e->p_id,
             value);

      // allocate analysis parameter arrays for variables and sensitivities
      p_xddm->a_ap[i].ndvs = p_xddm->nv;
      p_xddm->a_ap[i].pa_dvs = malloc( p_xddm->nv * sizeof(*(p_xddm->a_ap[i].pa_dvs)) );
      p_xddm->a_ap[i].a_lin = malloc( p_xddm->nv * sizeof(*p_xddm->a_ap[i].a_lin) );
      for (j = 0; j < p_xddm->nv; j++) {
        if (p_xddm->a_v[j].p_e->p_id == NULL)  continue;
        if (p_xddm->a_v[j].val  == UNSET) continue;
        p_xddm->a_ap[i].pa_dvs[j] = fillXmlString( (xmlChar *) p_xddm->a_v[j].p_e->p_id);
        noI = parsePmtr(p_xddm->a_v[j].p_e->p_id, pname, &irow, &icol);
        int kpmtr = 0;
        for (k = 0; k < MODL->npmtr; k++) {
          stat = ocsmGetPmtr(modl, k+1, &type, &nrow, &ncol, name);
          if (strcmp(pname, name) == 0) {
            kpmtr = k+1;
            break;
          }
        }
        if (kpmtr != 0) {
          /* clear all then set */
          ocsmSetVelD(modl, 0, 0, 0, 0.0);
          ocsmSetVelD(modl, kpmtr, irow, icol, 1.0);
          buildTo = 0;
          nb      = 0;
          printf(" CAPS Info: Building sensitivity information for: %s[%d,%d]\n",
                name, irow, icol);
          stat = ocsmBuild(modl, buildTo, &builtTo, &nb, NULL);
          stat = ocsmGetValu(modl, ipmtr, 1, 1, &value, &dot);
          p_xddm->a_ap[i].a_lin[j] = dot;
          if (noI == 0) {
            printf("\n*** compute parameter %d (%s) sensitivity status = %d (%d)***\n",
                  ipmtr, p_xddm->a_v[j].p_e->p_id, stat, nb);
          } else {
            printf("\n*** compute parameter %d [%d,%d] (%s) sensitivity = %d (%d)***\n",
                  ipmtr, irow, icol, p_xddm->a_v[j].p_e->p_id, stat, nb);
          }
        }
      }

    } else {
      printf(" Analysis Parameter %s not found!\n", p_xddm->a_ap[i].p_e->p_id);
      stat = EGADS_NOTFOUND;
      goto cleanup;
    }
  }
  printf("\n");
  
  /* get global tesselate parameters */
  global[0] = global[1] = global[2] = 0.0;
  for (i = 0; i < p_xddm->nt; i++) {
    if (p_xddm->a_t[i].p_id != NULL) continue;
    for (j = 0; j < p_xddm->a_t[i].nAttr; j++)
      if (strcmp(p_xddm->a_t[i].p_attr[j].p_name, "MaxEdge") == 0) {
        sscanf(p_xddm->a_t[i].p_attr[j].p_value, "%lf", &global[0]);
      } else if (strcmp(p_xddm->a_t[i].p_attr[j].p_name, "Sag") == 0) {
        sscanf(p_xddm->a_t[i].p_attr[j].p_value, "%lf", &global[1]);
      } else if (strcmp(p_xddm->a_t[i].p_attr[j].p_name, "Angle") == 0) {
        sscanf(p_xddm->a_t[i].p_attr[j].p_value, "%lf", &global[2]);
      } else {
        printf(" Tessellation (global) Attribute %s not Understood!\n",
               p_xddm->a_t[i].p_attr[j].p_name);
        goto cleanup;
      }
  }
  if ((global[0] == 0.0) && (global[1] == 0.0) && (global[2] == 0.0)) {
    global[0] =  0.025;
    global[1] =  0.001;
    global[2] = 12.0;
  }

  /* allocate sensitvity arrays */
  if (p_xddm->nv > 0) {
    dvar = (double ***) EG_alloc(nbody*sizeof(double **));
    if (dvar == NULL) {
      printf(" sensitivity allocation failed!\n");
      goto cleanup;
    }
    for (i = 0; i < nbody; i++) dvar[i] = NULL;

    for (i = 0; i < nbody; i++) {
      dvar[i] = (double**)EG_alloc(p_xddm->nv*sizeof(double*));
      if (dvar[i] == NULL) {
        printf(" sensitivity %d allocation failed!\n", i+1);
        goto cleanup;
      }
      for (j = 0; j < p_xddm->nv; j++) dvar[i][j] = NULL;
    }
  }


  /* tessellate all bodies */
  kb = 0;
  for (ibody = 1; ibody <= MODL->nbody; ibody++) {
    if (MODL->body[ibody].onstack != 1) continue;
    if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
    body      = MODL->body[ibody].ebody;
    tparam[0] = global[0];
    tparam[1] = global[1];
    tparam[2] = global[2];
    for (i = 0; i < p_xddm->nt; i++) {
      if (p_xddm->a_t[i].p_id == NULL) continue;
      sscanf(p_xddm->a_t[i].p_id, "%d", &stat);
      if (stat != kb+1) continue;
      for (j = 0; j < p_xddm->a_t[i].nAttr; j++)
        if (strcmp(p_xddm->a_t[i].p_attr[j].p_name, "MaxEdge") == 0) {
          sscanf(p_xddm->a_t[i].p_attr[j].p_value, "%lf", &tparam[0]);
        } else if (strcmp(p_xddm->a_t[i].p_attr[j].p_name, "Sag") == 0) {
          sscanf(p_xddm->a_t[i].p_attr[j].p_value, "%lf", &tparam[1]);
        } else if (strcmp(p_xddm->a_t[i].p_attr[j].p_name, "Angle") == 0) {
          sscanf(p_xddm->a_t[i].p_attr[j].p_value, "%lf", &tparam[2]);
        } else {
          printf(" Tessellation (ID=%d) Attribute %s not Understood!\n",
                 kb+1, p_xddm->a_t[i].p_attr[j].p_name);
          goto cleanup;
        }
    }
    printf(" Tessellating %d with  MaxEdge = %lf   Sag = %lf   Angle = %lf\n",
           kb+1, tparam[0], tparam[1], tparam[2]);
    stat = EG_getBoundingBox(body, box);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getBoundingBox failed: %d!\n", stat);
      goto cleanup;
    }
    size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                (box[4]-box[1])*(box[4]-box[1]) +
                (box[5]-box[2])*(box[5]-box[2]));
    tparam[0] *= size;
    tparam[1] *= size;
    stat = EG_makeTessBody(body, tparam, &MODL->body[ibody].etess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody failed: %d!\n", stat);
      goto cleanup;
    }
    tess[kb] = MODL->body[ibody].etess;
    
    stat = EG_statusTessBody(tess[kb], &body, &i, &nvert);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (p_xddm->nv > 0) {
      if (dvar     == NULL) { stat = EGADS_MALLOC; goto cleanup; }
      if (dvar[kb] == NULL) { stat = EGADS_MALLOC; goto cleanup; }
      for (i = 0; i < p_xddm->nv; i++) {
        dvar[kb][i] = (double*) EG_alloc(3*nvert*sizeof(double));
        if (dvar[kb][i] == NULL) {
          printf(" sensitivity %d allocation failed!\n", i+1);
          goto cleanup;
        }
        for (j = 0; j < 3*nvert; j++) dvar[kb][i][j] = 0.0;
      }
    }
    kb++;
  }
    
  /* compute the sensitivities */
  printf("\n");
  for (i = 0; i < p_xddm->nv; i++) {
    if (p_xddm->a_v[i].p_e->p_id == NULL)  continue;
    if (p_xddm->a_v[i].val  == UNSET) continue;
    noI = parsePmtr(p_xddm->a_v[i].p_e->p_id, pname, &irow, &icol);
    for (ipmtr = j = 0; j < MODL->npmtr; j++) {
      stat = ocsmGetPmtr(modl, j+1, &type, &nrow, &ncol, name);
      if (strcmp(pname, name) == 0) {
        ipmtr = j+1;
        break;
      }
    }
    /* clear all then set */
    step = 0.0;
    ocsmSetVelD(modl, 0,     0,    0,    0.0);
    ocsmSetVelD(modl, ipmtr, irow, icol, 1.0);
    if (p_xddm->a_v[i].p_e->p_comment != NULL) {
      if (strcmp(p_xddm->a_v[i].p_e->p_comment, "FD") == 0) {
        step = 0.001;
        printf("\n*** forced finite differencing for %s (%d) ***\n",
               pname, stat);
      }
      if (strcmp(p_xddm->a_v[i].p_e->p_comment, "oFD") == 0) {
        step = 0.001;
        printf("\n*** forcing OpenCSM finite differencing for %s (%d) ***\n",
               pname, stat);
      }
    }
    outLevel = ocsmSetOutLevel(0);
    stat = ocsmSetDtime(modl, step);
    if (step == 0.0) {
      buildTo = 0;
      nb      = 0;
      printf(" CAPS Info: Building sensitivity information for: %s[%d,%d]\n",
             name, irow, icol);
      stat    = ocsmBuild(modl, buildTo, &builtTo, &nb, NULL);
    }
    ocsmSetOutLevel(outLevel);
    if (noI == 0) {
      printf("\n*** compute parameter %d (%s) sensitivity status = %d (%d)***\n",
             ipmtr, p_xddm->a_v[i].p_e->p_id, stat, nb);
    } else {
      printf("\n*** compute parameter %d [%d,%d] (%s) sensitivity = %d (%d)***\n",
             ipmtr, irow, icol, p_xddm->a_v[i].p_e->p_id, stat, nb);
    }

    kb = 0;
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
      if (MODL->body[ibody].onstack != 1) continue;
      if (MODL->body[ibody].botype  == OCSM_NULL_BODY) continue;
      body = MODL->body[ibody].ebody;

      stat = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* do all of the faces */
      for (iface = 1; iface <= nface; iface++) {
        outLevel = ocsmSetOutLevel(0);
        stat = ocsmGetTessVel(MODL, ibody, OCSM_FACE, iface, &pcsens);
        ocsmSetOutLevel(outLevel);
        if (stat != EGADS_SUCCESS) {
          printf(" ocsmGetTessVel Parameter %d Face %d failed: %d!\n",
                 i+1, j, stat);
          goto cleanup;
        }
        stat = EG_getTessFace(tess[kb], iface, &nvert, &xyz, &uv,
                              &ptype, &pindx, &ntri, &tris, &tric);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_getTessFace Parameter %d Face %d failed: %d!\n",
                 i+1, j, stat);
          goto cleanup;
        }
        if (dvar        == NULL) { stat = EGADS_MALLOC; goto cleanup; }
        if (dvar[kb]    == NULL) { stat = EGADS_MALLOC; goto cleanup; }
        if (dvar[kb][i] == NULL) { stat = EGADS_MALLOC; goto cleanup; }
        for (k = 1; k <= nvert; k++) {
          stat = EG_localToGlobal(tess[kb], iface, k, &iglobal);
          if (stat != EGADS_SUCCESS) goto cleanup;

          dvar[kb][i][3*iglobal-3] = pcsens[3*k-3];
          dvar[kb][i][3*iglobal-2] = pcsens[3*k-2];
          dvar[kb][i][3*iglobal-1] = pcsens[3*k-1];
        }
      }
      kb++;
    }
  }
  printf("\n");

  /* write tri file */
  stat = writeTrix("Components.i.tri", nbody, tess,
                   p_xddm, p_xddm->nv, dvar);

  printf("\n");

  xddm_updateAnalysisParams(argv[1], p_xddm, opts);
  
cleanup:

  if (dvar != NULL && p_xddm != NULL) {
    for (i = 0; i < nbody; i++) {
      if (dvar[i] != NULL) {
        for (j = 0; j < p_xddm->nv; j++) EG_free(dvar[i][j]);
        EG_free(dvar[i]);
      }
    }
  }
  EG_free(dvar);
  EG_free(tess);
  EG_free(vtags);
  EG_free(psens);
  ocsmFree(modl);
  ocsmFree(NULL);
  if (context != NULL) EG_close(context);
  xddm_free(p_xddm, opts);

  if (stat == EGADS_SUCCESS)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
