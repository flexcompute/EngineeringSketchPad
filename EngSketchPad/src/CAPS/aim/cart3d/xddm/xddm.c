/*
 * $Id: xddm.c,v 1.4 2014/05/12 22:05:21 mnemec Exp $
 */

/**
 * Reader and writer of XML files that parses the
 * Extensible-Design-Description-Markup (XDDM) elements.
 *
 * Dependency: libxml2, www.xmlsoft.org. This library is usually present on
 * most systems, check existence of 'xml2-config' script
 *
 * Marian Nemec, STC, June 2013, latest update May 2014
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h> /* for fsync */

#include "xddm.h"
#include "xddmInternals.h"

static int
xddm_readElement(const xmlChar *str, xmlXPathContextPtr xpathCtx, p_tsXddmElem *p_e)
{
  int rc, i, ne;
  xmlXPathObjectPtr  xddmObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  xddmObj = xmlXPathEvalExpression(str, xpathCtx);
  if (xddmObj) {
    ne = countNodes(xddmObj); /* number of elements */

    if ( ne > 1 ) {
      ERR("more than one element while parsing %s\n", str);
      return -1;
    }

    if ( ne > 0 ) {
      /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr || pp_xmlAttr[0] == NULL) {
        ERR("problem in parsing %s\n", str);
        rc = -1;
        goto cleanup;
      }

      (*p_e) = xddm_allocElement(1);
      if ((*p_e) == NULL) {
        ERR("malloc problem\n");
        rc = -1;
        goto cleanup;
      }
      const p_tsXddmXmlAttr p_xA = pp_xmlAttr[0];

      /* first look for the comment */
      for (i=0; i<p_xA->n; i++) {
        if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
          (*p_e)->p_comment = fillXmlString(p_xA->a_values[i]);
        } else {
          rc = addXmlAttribute(p_xA->a_names[i], p_xA->a_values[i], &(*p_e)->nAttr, &(*p_e)->p_attr);
          if (rc != 0) {
            ERR("addXmlAttribute problem\n");
            goto cleanup;
          }
        }
      } /* end for attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if variables */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  } /* end if xddmObj */

  rc = 0;
cleanup:
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  return rc;
}

static int
xddm_readVariable(const xmlChar *str, xmlXPathContextPtr xpathCtx, int *nv, p_tsXddmVar *a_v)
{
  int rc = 0, i, inode=0;
  xmlXPathObjectPtr  xddmObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  xddmObj = xmlXPathEvalExpression(str, xpathCtx);
  if (xddmObj) {
    (*nv) = countNodes(xddmObj); /* number of variables */

    if ( (*nv) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing %s\n", str);
        rc = 1;
        goto cleanup;
      }

      (*a_v) = xddm_allocVariable((*nv));
      if ( NULL == (*a_v) ) {
        ERR("allocVariable failure\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*nv)) { rc = -1; goto cleanup; }
        const p_tsXddmVar p_v = &(*a_v)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_v->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_v->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
            p_v->val = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "min" ) ) {
            p_v->minVal = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "max" ) ) {
            p_v->maxVal = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "typicalsize" ) ) {
            p_v->typicalSize = fillDouble(p_xA->a_values[i]);
          }
        } /* end for attributes */

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if variables */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  } /* end if xddmObj */

  rc = 0;
cleanup:
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  return rc;
}


static int
xddm_readSensitivity(xmlXPathContextPtr xpathCtx, int *ndvs, double **a_lin, char ***pa_dvs)
{
  int rc = 0, i, j;
  xmlXPathObjectPtr linObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL, *pp_sensit  = NULL;

  /* check if the node has a sensitivity array */
  linObj = xmlXPathEvalExpression(BAD_CAST "./SensitivityArray/Sensitivity", xpathCtx);
  if (linObj) {
    if ( countNodes(linObj) > 0 ) {
      pp_sensit = xddmParseXpathObj(linObj);

      (*ndvs)=0;
      while (NULL != pp_sensit[(*ndvs)]) {
        (*ndvs)++;
      }

      (*pa_dvs) = malloc( (*ndvs) * sizeof(*(*pa_dvs)));
      if (NULL == (*pa_dvs)) {
        ERR("malloc failed\n");
        rc = 1;
        goto cleanup;
      }

      (*a_lin) = malloc( (*ndvs) * sizeof(*(*a_lin)));
      if (NULL == (*a_lin)) {
        ERR("malloc failed\n");
        rc = 1;
        goto cleanup;
      }

      i=0;
      while (NULL != pp_sensit[i]) {
        const p_tsXddmXmlAttr p_xA = pp_sensit[i];

        for (j=0; j<p_xA->n; j++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[j], BAD_CAST "p" ) ) {
            (*pa_dvs)[i] = fillXmlString(p_xA->a_values[j]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[j], BAD_CAST "value" ) ) {
            (*a_lin)[i] = fillDouble(p_xA->a_values[j]);
          }
        }
        i++;
      } /* end while sensitivities */

      freeXddmAttr(pp_sensit);
      pp_sensit = NULL;
    }

    xmlXPathFreeObject(linObj);
    linObj = NULL;
  }

  rc = 0;
cleanup:
  if (linObj    ) xmlXPathFreeObject(linObj);
  if (pp_sensit ) freeXddmAttr(pp_sensit);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  return rc;
}

static int
xddm_readFunctional(const xmlChar *str, xmlXPathContextPtr xpathCtx, int *nf, p_tsXddmFun *a_f)
{
  int rc = 0, i, inode=0;
  xmlXPathObjectPtr  xddmObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  xmlNodePtr p_nodeCur=NULL;

  (*nf) = 0;
  (*a_f) = NULL;

  /* save the current node so it can be reset
   */
  p_nodeCur = xpathCtx->node;

  xddmObj = xmlXPathEvalExpression(str, xpathCtx);
  if (xddmObj) {
    (*nf) = countNodes(xddmObj); /* number of variables */

    if ( (*nf) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing %s\n", str);
        rc = 1;
        goto cleanup;
      }

      (*a_f) = xddm_allocFunctinoal((*nf));
      if ( NULL == (*a_f) ) {
        ERR("allocVariable failure\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*nf)) { rc = -1; goto cleanup; }
        const p_tsXddmFun p_f = &(*a_f)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_f->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_f->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "expr" ) ) {
            p_f->p_expr = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
            p_f->val = fillDouble(p_xA->a_values[i]);
          }
          else {
            rc = addXmlAttribute(p_xA->a_names[i], p_xA->a_values[i], &p_f->nAttr, &p_f->p_attr);
            if (rc != 0) {
              ERR("addXmlAttribute problem\n");
              goto cleanup;
            }
          }

          xpathCtx->node = p_xA->p_node;

          /* check if the node has a sensitivity array */
          rc = xddm_readSensitivity(xpathCtx, &p_f->ndvs, &p_f->a_lin, &p_f->pa_dvs);
          if (rc != 0) {
            ERR("addXmlAttribute problem\n");
            goto cleanup;
          }

        } /* end for attributes */

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if variables */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  } /* end if xddmObj */

  rc = 0;
cleanup:
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur;

  return rc;
}

static int
xddm_readAeroFun(xmlXPathContextPtr xpathCtx, int *na, p_tsXddmAFun *a_afun)
{
  int rc = 0, i, inode=0;
  xmlXPathObjectPtr  xddmObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  xmlNodePtr p_nodeCur=NULL;

  /* save the current node so it can be reset
   */
  p_nodeCur = xpathCtx->node;

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./AeroFun", xpathCtx);
  if (xddmObj) {
    (*na) = countNodes(xddmObj); /* number of analysis parameters */

    if ( (*na) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in analysis parameters parsing\n");
        rc = 1;
        goto cleanup;
      }

      (*a_afun) = xddm_allocAeroFun(*na);
      if ( NULL == (*a_afun) ) {
        ERR("allocAnalysis failed\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*na)) { rc = -1; goto cleanup; }
        const p_tsXddmAFun p_a = &(*a_afun)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_a->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_a->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "options" ) ) {
            p_a->p_options = fillXmlString(p_xA->a_values[i]);
          }
        } /* end for attributes */

        xmlChar *p_xmlstr;
        p_xmlstr = xmlNodeGetContent(p_xA->p_node);
        if (NULL == p_xmlstr) {
          ERR("xmlNodeGetContent failed\n");
          rc = 1;
          goto cleanup;
        }
        char* p_str = (char *) p_xmlstr;
        char* pch = strtok(p_str, "\r\n");

        while (pch != NULL)
        {
          if (p_a->nt == 0) {
            p_a->pa_text = malloc(sizeof(*p_a->pa_text));
          } else {
            p_a->pa_text = realloc(p_a->pa_text, (p_a->nt+1)*sizeof(*p_a->pa_text));
          }
          p_a->pa_text[p_a->nt] = fillString(pch);
          p_a->nt++;
          pch = strtok(NULL, "\r\n");
        }

        free(p_xmlstr);

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if analysis parameters */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  rc = 0;
cleanup:
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur;

  return rc;
}

static int
xddm_readAnalysis(xmlXPathContextPtr xpathCtx, int *na, p_tsXddmAPar *a_ap)
{
  int rc = 0, i, inode=0, naero=0;
  xmlXPathObjectPtr  xddmObj = NULL, linObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  xmlNodePtr p_nodeCur=NULL;

  /* save the current node so it can be reset
   */
  p_nodeCur = xpathCtx->node;

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Analysis", xpathCtx);
  if (xddmObj) {
    (*na) = countNodes(xddmObj); /* number of analysis parameters */

    if ( (*na) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in analysis parameters parsing\n");
        rc = 1;
        goto cleanup;
      }

      (*a_ap) = xddm_allocAnalysis(*na);
      if ( NULL == (*a_ap) ) {
        ERR("allocAnalysis failed\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*na)) { rc = -1; goto cleanup; }
        const p_tsXddmAPar p_a = &(*a_ap)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_a->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_a->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
            p_a->val = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "discretizationerror" ) ) {
            p_a->discrErr = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "sensitivity" ) ) {
            if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "required" ) ) {
              p_a->lin = 1;
            } else if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "none" ) ) {
              p_a->lin = 0;
            }
          }
        } /* end for attributes */

        /* move the node into the analysis element */
        xpathCtx->node = pp_xmlAttr[inode]->p_node;

        /* look for an AeroFun element */
        rc = xddm_readAeroFun(xpathCtx, &naero, &p_a->p_afun);
        if (rc != 0) goto cleanup;
        if ( naero > 1 ) {
          ERR("Analysis may only contain one AeroFun!\n");
          rc = 1;
          goto cleanup;
        }

        /* check if the node has a sensitivity array */
        rc = xddm_readSensitivity(xpathCtx, &p_a->ndvs, &p_a->a_lin, &p_a->pa_dvs);
        if (rc != 0) {
          ERR("addXmlAttribute problem\n");
          goto cleanup;
        }

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if analysis parameters */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  rc = 0;
cleanup:
  if (linObj    ) xmlXPathFreeObject(linObj );
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur;

  return rc;
}


static int
xddm_readDesignPoint(xmlXPathContextPtr xpathCtx, int *nd, p_tsXddmDesP *a_dp)
{
  int rc = 0, i, inode=0;
  xmlXPathObjectPtr  xddmObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  xmlNodePtr p_nodeCur=NULL;

  /* save the current node so it can be reset
   */
  p_nodeCur = xpathCtx->node;

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./DesignPoint", xpathCtx);
  if (xddmObj) {
    (*nd) = countNodes(xddmObj); /* number of design points */

    if ( (*nd) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in design point parsing\n");
        rc = 1;
        goto cleanup;
      }

      (*a_dp) = xddm_allocDesignPoint(*nd);
      if ( NULL == (*a_dp) ) {
        ERR("allocDesignPoint failure\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*nd)) { rc = -1; goto cleanup; }
        const p_tsXddmDesP p_dp = &(*a_dp)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_dp->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_dp->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "geometry" ) ) {
            p_dp->p_geometry = fillXmlString(p_xA->a_values[i]);
          } else {
            rc = addXmlAttribute((*pp_xmlAttr)->a_names[i], (*pp_xmlAttr)->a_values[i], &p_dp->nAttr, &p_dp->p_attr);
            if (rc != 0) {
              ERR("addXmlAttribute problem\n");
              goto cleanup;
            }
          }
        } /* end for attributes */

        /* move the node into the design point element */
        xpathCtx->node = pp_xmlAttr[inode]->p_node;

        /* ---------------
         * parse Variables
         * ---------------
         */

        rc = xddm_readVariable(BAD_CAST "./Variable", xpathCtx, &p_dp->nv, &p_dp->a_v);
        if (rc != 0) {
          ERR("xddm_readVariable Variable failed\n");
          goto cleanup;
        }

        /* ---------------
         * parse Constants
         * ---------------
         */

        rc = xddm_readVariable(BAD_CAST "./Constant", xpathCtx, &p_dp->nc, &p_dp->a_c);
        if (rc != 0) {
          ERR("xddm_readVariable Constant failed\n");
          goto cleanup;
        }

        /* -------------------------
         * parse Analysis parameters
         * -------------------------
         */

        rc = xddm_readAnalysis(xpathCtx, &p_dp->na, &p_dp->a_ap);
        if (rc != 0) {
          ERR("xddm_readAnalysis failed\n");
          goto cleanup;
        }

        /* -------------------------
         * parse Objective function
         * -------------------------
         */

        int nobj = 0;
        rc = xddm_readFunctional(BAD_CAST "./Objective", xpathCtx, &nobj, &p_dp->p_obj);
        if (rc != 0 || nobj > 1) {
          ERR("xddm_readFunctional failed\n");
          freeFunctinoal(nobj, p_dp->p_obj);
          rc = -1;
          goto cleanup;
        }

        /* -------------------------
         * parse Constraint functions
         * -------------------------
         */

        rc = xddm_readFunctional(BAD_CAST "./Constraint", xpathCtx, &p_dp->ncr, &p_dp->a_cr);
        if (rc != 0 || nobj > 1) {
          ERR("xddm_readFunctional failed\n");
          rc = -1;
          goto cleanup;
        }

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if analysis parameters */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  rc = 0;
cleanup:
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur;

  return rc;
}

static int
xddm_readComponent(xmlXPathContextPtr xpathCtx, int *ncmp, p_tsXddmComp *a_cmp)
{
  int rc = 0, i, inode=0;
  xmlXPathObjectPtr  xddmObj = NULL, dataObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL, *pp_data = NULL;
  xmlNodePtr p_nodeCur=NULL;

  /* save the current node so it can be reset
   */
  p_nodeCur = xpathCtx->node;

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Component", xpathCtx);
  if (xddmObj) {
    (*ncmp) = countNodes(xddmObj); /* number of variables */

    if ( (*ncmp) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing Component\n");
        rc = 1;
        goto cleanup;
      }

      (*a_cmp) = xddm_allocComponent((*ncmp));
      if ( NULL == (*a_cmp) ) {
        ERR("allocVariable failure\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*ncmp)) { rc = -1; goto cleanup; }
        const p_tsXddmComp p_cmp = &(*a_cmp)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "name" ) ) {
            p_cmp->p_name = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_cmp->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "parent" ) ) {
            p_cmp->p_parent = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "type" ) ) {
            p_cmp->p_type = fillXmlString(p_xA->a_values[i]);
          }
          else {
            rc = addXmlAttribute(p_xA->a_names[i], p_xA->a_values[i], &p_cmp->nAttr, &p_cmp->p_attr);
            if (rc != 0) {
              ERR("addXmlAttribute problem\n");
              goto cleanup;
            }
          }

          /* move the node into the comonent element */
          xpathCtx->node = p_xA->p_node;

          /* look for an AeroFun element */
          dataObj = xmlXPathEvalExpression(BAD_CAST "./Data", xpathCtx);
          if (dataObj) {
            if ( countNodes(dataObj) > 1 ) {
              ERR("Component may only contain one Data!\n");
              rc = 1;
              goto cleanup;
            }
            if ( countNodes(dataObj) == 1 ) {
              pp_data = xddmParseXpathObj(dataObj);

              const p_tsXddmXmlAttr p_xA = pp_data[0];

              xmlChar *p_xmlstr;
              p_xmlstr = xmlNodeGetContent(p_xA->p_node);
              if (NULL == p_xmlstr) {
                ERR("xmlNodeGetContent failed\n");
                rc = 1;
                goto cleanup;
              }
              p_cmp->p_data = fillXmlString(p_xmlstr);
              free(p_xmlstr);

              freeXddmAttr(pp_data);
              pp_data = NULL;
            }

            xmlXPathFreeObject(dataObj);
            dataObj = NULL;
          }

        } /* end for attributes */

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if variables */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  } /* end if xddmObj */

  rc = 0;
cleanup:
  if (dataObj   ) xmlXPathFreeObject(dataObj);
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);
  if (pp_data   ) freeXddmAttr(pp_data);

  xpathCtx->node = p_nodeCur;

  return rc;
}


static int
xddm_readTessellate(xmlXPathContextPtr xpathCtx, int *ntess, p_tsXmTess *a_tess)
{
  int rc = 0, i, inode=0;
  xmlXPathObjectPtr  xddmObj = NULL;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  xmlNodePtr p_nodeCur=NULL;

  /* save the current node so it can be reset
   */
  p_nodeCur = xpathCtx->node;

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Tessellate", xpathCtx);
  if (xddmObj) {
    (*ntess) = countNodes(xddmObj); /* number of variables */

    if ( (*ntess) > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing Component\n");
        rc = 1;
        goto cleanup;
      }

      (*a_tess) = xddm_allocTessellate((*ntess));
      if ( NULL == (*a_tess) ) {
        ERR("allocVariable failure\n");
        rc = 1;
        goto cleanup;
      }

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        if (inode >= (*ntess)) { rc = -1; goto cleanup; }
        p_tsXmTess p_t = &(*a_tess)[inode];
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_t->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_t->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "sensitivity" ) ) {
            if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "required" ) ) {
              p_t->lin = 1;
            } else if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "none" ) ) {
              p_t->lin = 0;
            }

          }
          else {
            if ( 0 == p_t->nAttr ) {
              p_t->nAttr  = 1;
              p_t->p_attr = malloc(sizeof(*p_t->p_attr));
              if (NULL == p_t->p_attr) {
                ERR("malloc failed\n");
                exit(1);
              }
            }
            else {
              p_tsXddmAttr p_tmp = NULL;
              p_t->nAttr  += 1;
              p_tmp = realloc(p_t->p_attr, p_t->nAttr*sizeof(tsXddmAttr));
              if (NULL == p_tmp) {
                ERR("realloc failed\n");
                exit(1);
              }
              else {
                p_t->p_attr = p_tmp;
              }
            }
            p_t->p_attr[p_t->nAttr-1].p_name  = fillXmlString(p_xA->a_names[i]);
            p_t->p_attr[p_t->nAttr-1].p_value = fillXmlString(p_xA->a_values[i]);
          }
        } /* end for attributes */

        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
      pp_xmlAttr = NULL;
    } /* end if variables */

    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  } /* end if xddmObj */

  rc = 0;
cleanup:
  if (xddmObj   ) xmlXPathFreeObject(xddmObj);
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur;

  return rc;
}

/**
 * xddm_readFile(): parse an XPath expression from an XDDM file and return
 * structure containing all elements in the path
 */
p_tsXddm
xddm_readFile(const char *const p_fileName, const char *const xpathExpr, const int options)
{
  xmlDocPtr          doc = NULL;
  xmlParserCtxtPtr   ctxt = NULL;
  xmlXPathObjectPtr  xddmObj = NULL;
  xmlXPathContextPtr xpathCtx = NULL;

  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  p_tsXddm p_xddm = NULL; /* return object */

  p_tsXmParent p_p = NULL;

  int parserOpts=0;
  int rc, i;
                                     /* initialize and check libxml2 version */
  { LIBXML_TEST_VERSION };

  xmlXPathInit();  /* initialize xpath */

  if ( XDDM_VERBOSE & options ) {
    printf(" o Parsing file \"%s\" with libxml2\n", p_fileName);
  }

  ctxt = xmlNewParserCtxt(); /* create a document parser context */
  if (ctxt == NULL) {
    ERR(" xddm_readFile failed to allocate parser context\n");
    xmlCleanupParser();
    return NULL;
  }

  doc = xmlCtxtReadFile(ctxt, p_fileName, NULL, parserOpts);
  if (doc == NULL) {
    if ( XDDM_VERBOSE & options ) {
      ERR("%s is not valid XML\n", p_fileName);
    }
    xmlCleanupParser();
    return NULL;
  }

  xmlFreeParserCtxt(ctxt); /* done with document parser context */

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("xmlXPathNewContext failed to create context\n");
    rc = 1;
    goto cleanup;
  }

  /* looks like we have a valid xml document that we are ready to parse
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if (NULL == xddmObj) {
    ERR("no elements found for expression %s\n", xpathExpr);
    rc = 1;
    goto cleanup;
  }

  if ( 0 == countNodes(xddmObj) ) {
    ERR("no elements found for expression \'%s\'\n", xpathExpr);
    rc = 1;
    goto cleanup;
  }

  p_xddm = xddm_alloc(); /* return object */

  p_xddm->fileName  = fillString(p_fileName);
  p_xddm->xpathExpr = fillString(xpathExpr);

  pp_xmlAttr = xddmParseXpathObj(xddmObj);
  if ( NULL == pp_xmlAttr || pp_xmlAttr[0] == NULL ) {
    ERR("problem in parsing variables\n");
    rc = 1;
    goto cleanup;
  }

  p_p = p_xddm->p_parent;
  p_p->name   = fillXmlString((*pp_xmlAttr)->p_node->name);
  p_p->nAttr  = (*pp_xmlAttr)->n;
  p_p->p_attr = xddm_allocAttribute(p_p->nAttr);
  if (NULL == p_p->p_attr) {
    ERR("malloc failed\n");
    rc = 1;
    goto cleanup;
  }

  for (i=0; i<p_p->nAttr; i++) {
    p_p->p_attr[i].p_name  = fillXmlString((*pp_xmlAttr)->a_names[i]);
    p_p->p_attr[i].p_value = fillXmlString((*pp_xmlAttr)->a_values[i]);
  }

  freeXddmAttr(pp_xmlAttr);
  pp_xmlAttr = NULL;

  xpathCtx->node = xddmObj->nodesetval->nodeTab[0];
  xmlXPathFreeObject(xddmObj);
  xddmObj = NULL;

  /* ---------------
   * parse Configure
   * ---------------
   */

  rc = xddm_readElement(BAD_CAST "./Configure", xpathCtx, &p_xddm->p_config);
  if (rc != 0) {
    ERR("xddm_readElement Configure failed\n");
    goto cleanup;
  }

  /* ---------------
   * parse Intersect
   * ---------------
   */

  rc = xddm_readElement(BAD_CAST "./Intersect", xpathCtx, &p_xddm->p_inter);
  if (rc != 0) {
    ERR("xddm_readElement Intersect failed\n");
    goto cleanup;
  }

  /* ---------------
   * parse Variables
   * ---------------
   */

  rc = xddm_readVariable(BAD_CAST "./Variable", xpathCtx, &p_xddm->nv, &p_xddm->a_v);
  if (rc != 0) {
    ERR("xddm_readVariable Variable failed\n");
    goto cleanup;
  }

  /* ---------------
   * parse Constants
   * ---------------
   */

  rc = xddm_readVariable(BAD_CAST "./Constant", xpathCtx, &p_xddm->nc, &p_xddm->a_c);
  if (rc != 0) {
    ERR("xddm_readVariable Constant failed\n");
    goto cleanup;
  }

  /* -------------------------
   * parse Analysis parameters
   * -------------------------
   */

  rc = xddm_readAnalysis(xpathCtx, &p_xddm->na, &p_xddm->a_ap);
  if (rc != 0) {
    ERR("xddm_readAnalysis failed\n");
    goto cleanup;
  }

  /* -------------------------
   * parse Design Point
   * -------------------------
   */

  rc = xddm_readDesignPoint(xpathCtx, &p_xddm->nd, &p_xddm->a_dp);
  if (rc != 0) {
    ERR("xddm_readDesignPoint failed\n");
    goto cleanup;
  }

  /* -------------------------
   * parse Component
   * -------------------------
   */

  rc = xddm_readComponent(xpathCtx, &p_xddm->ncmp, &p_xddm->a_cmp);
  if (rc != 0) {
    ERR("xddm_readComponent failed\n");
    goto cleanup;
  }

  /* -------------------------
   * parse AeroFun
   * -------------------------
   */

  rc = xddm_readAeroFun(xpathCtx, &p_xddm->nafun, &p_xddm->a_afun);
  if (rc != 0) {
    ERR("xddm_readAeroFun failed\n");
    goto cleanup;
  }

  /* -------------------------
   * parse Tessellate
   * -------------------------
   */

  rc = xddm_readTessellate(xpathCtx, &p_xddm->ntess, &p_xddm->a_tess);
  if (rc != 0) {
    ERR("xddm_readTessellate failed\n");
    goto cleanup;
  }

cleanup:
  /* all done, clean up and return
   */
  if (rc != 0) {
    xddm_free(p_xddm);
    p_xddm = NULL;
  }
  if (pp_xmlAttr) freeXddmAttr(pp_xmlAttr);
  if (xddmObj)    xmlXPathFreeObject(xddmObj);
  if (xpathCtx)   xmlXPathFreeContext(xpathCtx);
  if (doc)        shutdown_libxml(doc);

  return p_xddm;
}

int
countNodes(const xmlXPathObjectPtr xpathObj)
{
  int n=0;
  int nNodes, i;
  xmlNodeSetPtr nodes;

  assert(xpathObj);

  if(xmlXPathNodeSetIsEmpty(xpathObj->nodesetval)){
    nNodes = 0;
  }
  else {
    nodes = xpathObj->nodesetval;
    nNodes= (nodes) ? nodes->nodeNr : 0;
  }
                                            /* find number of element nodes */
  for (i=0; i < nNodes; i++) {
    if ( nodes->nodeTab[i]->type != XML_ELEMENT_NODE ) continue;
    n ++;
  }

  return n;
}

p_tsXddmXmlAttr *
xddmParseXpathObj(const xmlXPathObjectPtr xpathObj)
{
  int i;

  const xmlNodeSetPtr nodes = xpathObj->nodesetval;
  const int nNodes= (nodes) ? nodes->nodeNr : 0;

  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  pp_xmlAttr = malloc( (nNodes+1) * sizeof(*pp_xmlAttr));
  if (NULL == pp_xmlAttr) {
    ERR("malloc failed\n");
    exit(1);
  }
  for (i=0; i < nNodes; i++) pp_xmlAttr[i] = NULL;

  for (i=0; i < nNodes; i++) {
    if ( nodes->nodeTab[i]->type != XML_ELEMENT_NODE ) continue;

    pp_xmlAttr[i] = xddmParseNode(nodes->nodeTab[i]);
    if (NULL == pp_xmlAttr[i]) {
      ERR("node has no attributes\n");
      exit(1);
    }
  }

  pp_xmlAttr[nNodes] = NULL; /* sentinel */

  return pp_xmlAttr;
}

p_tsXddmXmlAttr
xddmParseNode(const xmlNodePtr p_node)
{
  int             i;
  xmlAttr        *p_attribute;
  p_tsXddmXmlAttr p_xmlAttr;

  assert(p_node);

  p_xmlAttr = malloc(sizeof(*p_xmlAttr));
  if (NULL == p_xmlAttr) {
    ERR("malloc failed\n");
    exit(1);
  }

  p_xmlAttr->n = 0;

  p_attribute = p_node->properties;
  while (p_attribute && p_attribute->name && p_attribute->children) {
    p_attribute = p_attribute->next;
    p_xmlAttr->n++;
  }

  /* printf("\nElement %s has %d attributes\n",(char *) p_node->name, p_xmlAttr->n); */

  p_xmlAttr->a_names = malloc(p_xmlAttr->n * sizeof(*p_xmlAttr->a_names));
  if (NULL == p_xmlAttr->a_names) {
    ERR("malloc failed\n");
    exit(1);
  }

  p_xmlAttr->a_values = malloc(p_xmlAttr->n * sizeof(*p_xmlAttr->a_values));
  if (NULL == p_xmlAttr->a_values) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i = 0; i < p_xmlAttr->n; i++) {
    p_xmlAttr->a_names[i]  = NULL;
    p_xmlAttr->a_values[i] = NULL;
  }

  i = 0;
  p_attribute = p_node->properties;
  while (p_attribute && p_attribute->name && p_attribute->children) {
    p_xmlAttr->a_names[i]  = p_attribute->name;
    p_xmlAttr->a_values[i] = xmlGetProp(p_node, p_attribute->name);
    if (NULL == p_xmlAttr->a_values[i]) {
      ERR("xddmParseNode xmlGetProp failed for %s\n",(char *)p_attribute->name);
      freeAttributes(p_xmlAttr);
      fflush(stdout);
      fflush(stderr);
      return NULL;
    }
    /* printf(" Attribute %s Value %s\n",(char *)p_xmlAttr->a_names[i], */
    /*        (char *)p_xmlAttr->a_values[i]); */
    p_attribute = p_attribute->next;
    i++;
  }

  p_xmlAttr->p_node = p_node;

  return p_xmlAttr;
}

char *
fillXmlString(const xmlChar * p_xml)
{
  const int len = xmlStrlen( p_xml ) + 1;
  char *p_c = NULL;

  if ( len > MAX_STR_LEN ) {
    ERR("huge xml string\n");
    exit(1);
  }

  p_c = malloc( len * sizeof(*p_c));
  if ( NULL == p_c ) {
    ERR("malloc failed\n");
    exit(1);
  }

  memset((void*) p_c,'\0',len);

  (void) strncpy(p_c, (char *) p_xml, len);

  p_c[len - 1] = '\0'; /* force end string */

  return p_c;
}

char *
fillString(const char* p_string)
{
  const int len = strlen( p_string ) + 1;
  char *p_c = NULL;

  if ( len > MAX_STR_LEN ) {
    ERR("huge string\n");
    exit(1);
  }

  p_c = malloc( len * sizeof(*p_c));
  if ( NULL == p_c ) {
    ERR("malloc failed\n");
    exit(1);
  }

  memset((void*) p_c,'\0',len);

  (void) strncpy(p_c, (char *) p_string, len);

  p_c[len - 1] = '\0'; /* force end string */

  return p_c;
}

double
fillDouble(xmlChar * p_xml)
{
  double val = xmlXPathCastStringToNumber(p_xml);
                                                       /* basic error checks */
  if ( 0 != xmlXPathIsInf(val) ) {
    ERR("Value is INFINITE\n");
    exit(1);
  }
  if ( 1 == xmlXPathIsNaN(val) ) {
    ERR("Value is NaN\n");
    exit(1);
  }

  if ( UNSET == val ) {
    ERR("Value is UNSET -- conflict with internal defaults\n");
    exit(1);
  }

  return val;
}

p_tsXddmAttr
xddm_allocAttribute(const int n)
{
  int i;
  p_tsXddmAttr p_a=NULL;

  p_a = malloc(n * sizeof(*p_a));
  if (NULL == p_a) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_a[i].p_name  = NULL;
    p_a[i].p_value = NULL;
  }

  return p_a;
}

int
xddm_addAttribute(const char *const name, const char *const value, int *n, p_tsXddmAttr *p_a)
{
  if (0 == (*n)) {
    (*p_a) = malloc(sizeof(**p_a));
  } else {
    (*p_a) = realloc((*p_a), ((*n)+1) * sizeof(**p_a));
  }
  if (NULL == p_a) {
    ERR("malloc failed\n");
    return -1;
  }

  (*p_a)[(*n)].p_name  = fillString(name);
  (*p_a)[(*n)].p_value = fillString(value);
  if (NULL == (*p_a)[(*n)].p_name ||
      NULL == (*p_a)[(*n)].p_value) {
    ERR("malloc failed\n");
    return -1;
  }
  (*n)++;

  return 0;
}

int
addXmlAttribute(const xmlChar *name, const xmlChar *value, int *n, p_tsXddmAttr *p_a)
{
  if (0 == (*n)) {
    (*p_a) = malloc(sizeof(**p_a));
  } else {
    (*p_a) = realloc((*p_a), ((*n)+1) * sizeof(**p_a));
  }
  if (NULL == p_a) {
    ERR("malloc failed\n");
    return -1;
  }

  (*p_a)[(*n)].p_name  = fillXmlString(name);
  (*p_a)[(*n)].p_value = fillXmlString(value);
  if (NULL == (*p_a)[(*n)].p_name ||
      NULL == (*p_a)[(*n)].p_value) {
    ERR("malloc failed\n");
    return -1;
  }
  (*n)++;

  return 0;
}

void
freeAttribute(const int n, p_tsXddmAttr p_a)
{
  int i;
  if (!p_a) return;

  for (i=0; i<n; i++) {
    if (p_a[i].p_name ) free(p_a[i].p_name);
    if (p_a[i].p_value) free(p_a[i].p_value);
  }

  free(p_a);
  return;
}

p_tsXmParent
xddm_allocParent()
{
  p_tsXmParent p_p=NULL;

  p_p = malloc(sizeof(*p_p));
  if (NULL == p_p) {
    ERR("malloc failed\n");
    exit(1);
  }

  p_p->name   = NULL;
  p_p->nAttr  = 0;
  p_p->p_attr = NULL;

  return p_p;
}

void
freeParent(p_tsXmParent p_p)
{
  if (!p_p) return;

  if (p_p->name) free(p_p->name);
  freeAttribute(p_p->nAttr, p_p->p_attr);

  free(p_p);
  return;
}

p_tsXddmElem
xddm_allocElement(const int n)
{
  int i;
  p_tsXddmElem p_e=NULL;

  p_e = malloc(n*sizeof(*p_e));
  if (NULL == p_e) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_e[i].p_comment = NULL;
    p_e[i].nAttr     = 0;
    p_e[i].p_attr    = NULL;
  }

  return p_e;
}

void
freeElement(p_tsXddmElem p_p)
{
  if (!p_p) return;

  if (p_p->p_comment) free(p_p->p_comment);
  freeAttribute(p_p->nAttr, p_p->p_attr);

  free(p_p);
  return;
}

p_tsXddmVar
xddm_allocVariable(const int n)
{
  int i;
  p_tsXddmVar p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_id        = NULL;
    p_v[i].p_comment   = NULL;
    /* p_v[i].p_parent    = NULL; */
    p_v[i].val         = UNSET;
    p_v[i].minVal      = UNSET;
    p_v[i].maxVal      = UNSET;
    p_v[i].typicalSize = UNSET;
  }

  return p_v;
}

void
freeVariable(const int n, p_tsXddmVar p_v)
{
  int i;
  if (!p_v) return;

  for (i=0; i<n; i++) {
    if (p_v[i].p_id     ) free(p_v[i].p_id);
    if (p_v[i].p_comment) free(p_v[i].p_comment);
  }
  free(p_v);
  return;
}

p_tsXddmFun
xddm_allocFunctinoal(const int n)
{
  int i;
  p_tsXddmFun p_f=NULL;

  p_f = malloc(n * sizeof(*p_f));
  if (NULL == p_f) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_f[i].p_id        = NULL;
    p_f[i].p_comment   = NULL;
    p_f[i].p_expr      = NULL;
    p_f[i].val         = UNSET;
    p_f[i].lin         = -1;
    p_f[i].ndvs        = 0;
    p_f[i].a_lin       = NULL;
    p_f[i].pa_dvs      = NULL;
    p_f[i].nAttr     = 0;
    p_f[i].p_attr    = NULL;
  }

  return p_f;
}

void
freeFunctinoal(const int n, p_tsXddmFun p_f)
{
  int i, j;
  if (!p_f) return;

  for (i=0; i<n; i++) {
    if (p_f[i].p_id     ) free(p_f[i].p_id);
    if (p_f[i].p_comment) free(p_f[i].p_comment);
    if (p_f[i].p_expr   ) free(p_f[i].p_expr);
    for (j=0; j<p_f[i].ndvs; j++)
      if (p_f[i].pa_dvs[j]) free(p_f[i].pa_dvs[j]);
    if (p_f[i].a_lin) free(p_f[i].a_lin);
  }
  free(p_f);
  return;
}

p_tsXddmAFun
xddm_allocAeroFun(const int n)
{
  int i;
  p_tsXddmAFun p_afun=NULL;

  p_afun = malloc(n * sizeof(*p_afun));
  if (NULL == p_afun) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_afun[i].p_id      = NULL;
    p_afun[i].p_comment = NULL;
    p_afun[i].p_options = NULL;
    p_afun[i].nt        = 0;
    p_afun[i].pa_text   = NULL;
    p_afun[i].nAttr     = 0;
    p_afun[i].p_attr    = NULL;
  }

  return p_afun;
}

void
freeAeroFun(const int n, p_tsXddmAFun p_afun)
{
  int i, j;
  if (!p_afun) return;

  for (i=0; i<n; i++) {
    if (p_afun[i].p_id     ) free(p_afun[i].p_id);
    if (p_afun[i].p_comment) free(p_afun[i].p_comment);
    if (p_afun[i].p_options) free(p_afun[i].p_options);
    if (p_afun[i].pa_text) {
      for (j=0; j<p_afun[i].nt; j++)
        if (p_afun[i].pa_text[j]) free(p_afun[i].pa_text[j]);
      free(p_afun[i].pa_text);
    }
    freeAttribute(p_afun[i].nAttr, p_afun[i].p_attr);
  }

  free(p_afun);
  return;
}

p_tsXddmAPar
xddm_allocAnalysis(const int n)
{
  int i;
  p_tsXddmAPar p_a=NULL;

  p_a = malloc(n * sizeof(*p_a));
  if (NULL == p_a) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_a[i].p_id      = NULL;
    p_a[i].p_comment = NULL;
    /* p_v[i].p_parent  = NULL; */
    p_a[i].val       = UNSET;
    p_a[i].discrErr  = UNSET;
    p_a[i].lin       = -1;
    p_a[i].ndvs      = 0;
    p_a[i].a_lin     = NULL;
    p_a[i].pa_dvs    = NULL;

    /* AeroFun Element */
    p_a[i].p_afun    = NULL;
  }

  return p_a;
}

void
freeAnalysis(const int n, p_tsXddmAPar p_a)
{
  int i, j;
  if (!p_a) return;

  for (i=0; i<n; i++) {
    if (p_a[i].p_id     ) free(p_a[i].p_id);
    if (p_a[i].p_comment) free(p_a[i].p_comment);

    if (p_a[i].a_lin) free(p_a[i].a_lin);
    if (p_a[i].ndvs) {
      for (j=0; j<p_a[i].ndvs; j++)
        if (p_a[i].pa_dvs[j]) free(p_a[i].pa_dvs[j]);
      free(p_a[i].pa_dvs);
    }

    freeAeroFun(1, p_a[i].p_afun);
  }

  free(p_a);
  return;
}

p_tsXddmDesP
xddm_allocDesignPoint(const int n)
{
  int i;
  p_tsXddmDesP p_dp=NULL;

  p_dp = malloc(n * sizeof(*p_dp));
  if (NULL == p_dp) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_dp[i].p_id       = NULL;
    p_dp[i].p_comment  = NULL;
    p_dp[i].p_geometry = NULL;
    p_dp[i].nv         = 0;
    p_dp[i].a_v        = NULL;
    p_dp[i].nc         = 0;
    p_dp[i].a_c        = NULL;
    p_dp[i].na         = 0;
    p_dp[i].a_ap       = NULL;
    p_dp[i].p_obj      = NULL;
    p_dp[i].ncr        = 0;
    p_dp[i].a_cr       = NULL;
    p_dp[i].nAttr      = 0;
    p_dp[i].p_attr     = NULL;
  }

  return p_dp;
}

void
freeDesignPoint(const int n, p_tsXddmDesP p_dp)
{
  int i;
  if (!p_dp) return;

  for (i=0; i<n; i++) {
    if (p_dp[i].p_id      ) free(p_dp[i].p_id);
    if (p_dp[i].p_comment ) free(p_dp[i].p_comment);
    if (p_dp[i].p_geometry) free(p_dp[i].p_geometry);
    freeVariable(p_dp[i].nv, p_dp[i].a_v);
    freeVariable(p_dp[i].nc, p_dp[i].a_c);
    freeAnalysis(p_dp[i].na, p_dp[i].a_ap);
    freeFunctinoal(1          , p_dp[i].p_obj);
    freeFunctinoal(p_dp[i].ncr, p_dp[i].a_cr);
    freeAttribute(p_dp[i].nAttr, p_dp[i].p_attr);
  }

  free(p_dp);
  return;
}

p_tsXddmComp
xddm_allocComponent(const int n)
{
  int i;
  p_tsXddmComp p_cmp=NULL;

  p_cmp = malloc(n * sizeof(*p_cmp));
  if (NULL == p_cmp) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_cmp[i].p_name     = NULL;
    p_cmp[i].p_comment  = NULL;
    p_cmp[i].p_parent   = NULL;
    p_cmp[i].p_type     = NULL;
    p_cmp[i].p_data     = NULL;
    p_cmp[i].nAttr      = 0;
    p_cmp[i].p_attr     = NULL;
  }

  return p_cmp;
}

void
freeComponent(const int n, p_tsXddmComp p_cmp)
{
  int i;
  if (!p_cmp) return;

  for (i=0; i<n; i++) {
    if (p_cmp[i].p_name    ) free(p_cmp[i].p_name);
    if (p_cmp[i].p_comment ) free(p_cmp[i].p_comment);
    if (p_cmp[i].p_parent  ) free(p_cmp[i].p_parent);
    if (p_cmp[i].p_type    ) free(p_cmp[i].p_type);
    if (p_cmp[i].p_data    ) free(p_cmp[i].p_data);
    freeAttribute(p_cmp[i].nAttr, p_cmp[i].p_attr);
  }

  free(p_cmp);
  return;
}

p_tsXmTess
xddm_allocTessellate(const int n)
{
  int i;
  p_tsXmTess p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_id      = NULL;
    p_v[i].p_comment = NULL;
    p_v[i].lin       = -1;
    p_v[i].nAttr     = 0;
    p_v[i].p_attr    = NULL;
  }

  return p_v;
}

void
freeTessellate(const int n, p_tsXmTess p_v)
{
  int i;
  if (!p_v) return;

  for (i=0; i<n; i++) {
    if (p_v[i].p_id     ) free(p_v[i].p_id);
    if (p_v[i].p_comment) free(p_v[i].p_comment);
    freeAttribute(p_v[i].nAttr, p_v[i].p_attr);
  }

  free(p_v);
  return;
}

void
freeAttributes(p_tsXddmXmlAttr p_xmlAttr)
{
  int j;
  for (j=0; j<p_xmlAttr->n; j++) {
    xmlFree(p_xmlAttr->a_values[j]);
  }
  free(p_xmlAttr->a_values);
  free(p_xmlAttr->a_names);
  free(p_xmlAttr);
  return;
}

void
freeXddmAttr(p_tsXddmXmlAttr *pp_xmlAttr)
{
  int i=0;

  if (pp_xmlAttr) {
    while (NULL != pp_xmlAttr[i]) {
      if (pp_xmlAttr[i]) freeAttributes(pp_xmlAttr[i]);
      i++;
    }
    free(pp_xmlAttr);
  }

  return;
}

static void
shutdown_libxml(xmlDocPtr doc)
{
  xmlFreeDoc(doc);
  xmlCleanupParser(); /* shutdown libxml */
  xmlMemoryDump(); /* this is to debug memory for regression tests */
  return;
}

static void
echoAttributes(const int nAttr, const p_tsXddmAttr p_attr)
{
  int i;
  if (nAttr == 0) return;
  printf("   Attributes:\n");
  for (i=0; i<nAttr; i++) {
    printf("      %s = %s\n", p_attr[i].p_name, p_attr[i].p_value);
  }
}

static void
echoElement(const p_tsXddmElem p_e)
{
  if (!p_e) return;

  echoAttributes(p_e->nAttr, p_e->p_attr);
  if (NULL != p_e->p_comment) printf("        %s\n",p_e->p_comment);
}

static void
echoVariable(const p_tsXddmVar p_v)
{
  if (!p_v) return;

  if (NULL != p_v->p_id) printf(" ID=%s",p_v->p_id);
  if (UNSET != p_v->val) printf("  Value=%g",p_v->val);
  printf("\n");
  if (UNSET != p_v->minVal) printf("        Min=%g",p_v->minVal);
  if (UNSET != p_v->maxVal) printf(" Max=%g",p_v->maxVal);
  if (UNSET != p_v->typicalSize) printf(" TypicalSize=%g",p_v->typicalSize);
  if (UNSET != p_v->minVal || UNSET != p_v->maxVal ||
      UNSET != p_v->typicalSize ) printf("\n");
  if (NULL != p_v->p_comment) printf("        %s\n",p_v->p_comment);
}

static void
echoFunctional(const p_tsXddmFun p_f)
{
  int j;
  if (!p_f) return;

  if (NULL  != p_f->p_id) printf(" ID=%s",p_f->p_id);
  if (UNSET != p_f->val) printf("  Value=%g",p_f->val);
  if (NULL  != p_f->p_expr) printf("  Expr=%s",p_f->p_expr);
  printf("\n");

  if ( p_f->ndvs > 0 ) {
    printf("        Sensitivity array\n");
    for (j=0; j<p_f->ndvs; j++) {
      printf("           DV=%s Value=%g\n",p_f->pa_dvs[j], p_f->a_lin[j]);
    }
  }
  echoAttributes(p_f->nAttr, p_f->p_attr);
  if (NULL != p_f->p_comment) printf("        %s\n",p_f->p_comment);
}

static void
echoAeroFun(const p_tsXddmAFun p_afun)
{
  int j;
  if (!p_afun) return;

  if (NULL != p_afun->p_id     ) printf(" ID=%s",p_afun->p_id);
  if (NULL != p_afun->p_options) printf(" Options=%s",p_afun->p_options);
  printf("\n");
  if (NULL != p_afun->p_comment) printf("        %s\n",p_afun->p_comment);
  echoAttributes(p_afun->nAttr, p_afun->p_attr);

  if (p_afun->nt > 0) {
    for (j=0; j<p_afun->nt; j++) {
      printf("\n%s",p_afun->pa_text[j]);
    }
    printf("\n");
  }
}

static void
echoAnalysis(const p_tsXddmAPar p_a)
{
  int j;
  if (!p_a) return;

  if (NULL != p_a->p_id) printf(" ID=%s",p_a->p_id);
  if (UNSET != p_a->val) printf("  Value=%g",p_a->val);
  printf("\n");
  if (p_a->lin == 1) printf("        Linearization Required\n");
  if (p_a->lin == 0) printf("        Linearization Dissabled\n");
  if (NULL != p_a->p_comment) printf("        %s\n",p_a->p_comment);

  if (p_a->p_afun) {
    printf("        AeroFun:");
    echoAeroFun(p_a->p_afun);
    printf("\n");
  }

  if ( p_a->ndvs > 0 ) {
    printf("        Sensitivity array\n");
    for (j=0; j<p_a->ndvs; j++) {
      printf("           DV=%s Value=%g\n",p_a->pa_dvs[j], p_a->a_lin[j]);
    }
  }
}

static void
echoDesignPoint(const p_tsXddmDesP p_dp)
{
  int i;
  if (!p_dp) return;

  if (NULL != p_dp->p_id) printf(" ID=%s",p_dp->p_id);
  if (NULL != p_dp->p_geometry) printf("  Geometry=%s",p_dp->p_geometry);
  printf("\n");
  if (NULL != p_dp->p_comment) printf("        %s\n",p_dp->p_comment);

  printf(" o Number of variables = %d\n",p_dp->nv);
  for (i=0; i<p_dp->nv; i++) {
    printf("\n   %4d",i);
    echoVariable(&p_dp->a_v[i]);
  }

  printf("\n");
  printf(" o Number of constants = %d\n",p_dp->nc);
  for (i=0; i<p_dp->nc; i++) {
    printf("\n   %4d",i);
    echoVariable(&p_dp->a_c[i]);
  }

  printf("\n");
  printf(" o Number of analysis parameters = %d\n",p_dp->na);
  for (i=0; i<p_dp->na; i++) {
    printf("\n   %4d",i);
    echoAnalysis(&p_dp->a_ap[i]);
  }

  if (p_dp->p_obj) {
    printf("\n");
    printf(" o Objective functional:\n");
    echoFunctional(p_dp->p_obj);
  }

  if (p_dp->ncr > 0) {
    printf("\n");
    printf(" o Number of Constraint functional = %d\n", p_dp->ncr);
    for (i=0; i<p_dp->ncr; i++) {
      printf("\n   %4d",i);
      echoFunctional(&p_dp->a_cr[i]);
    }
  }
}

static void
echoComponent(const p_tsXddmComp p_cmp)
{
  if (!p_cmp) return;

  if (NULL  != p_cmp->p_name) printf(" Name=%s",p_cmp->p_name);
  if (NULL  != p_cmp->p_parent) printf("  Parent=%s",p_cmp->p_parent);
  if (NULL  != p_cmp->p_type) printf("  Parent=%s",p_cmp->p_type);
  printf("\n");
  if (NULL  != p_cmp->p_data) printf("  Data=%s",p_cmp->p_data);
  echoAttributes(p_cmp->nAttr, p_cmp->p_attr);
  if (NULL != p_cmp->p_comment) printf("        %s\n",p_cmp->p_comment);
}

void
xddm_echo(const p_tsXddm p_xddm)
{
  int i, j;

  printf("\n");

  if (p_xddm->xpathExpr && p_xddm->fileName) {
    printf(" o Evaluated \'%s\' for file \'%s\'\n", p_xddm->xpathExpr,
           p_xddm->fileName);
  }

  if (p_xddm->p_parent) {
    if (p_xddm->p_parent->name) printf(" o Name of parent element = %s\n",
                                       p_xddm->p_parent->name);
    echoAttributes(p_xddm->p_parent->nAttr, p_xddm->p_parent->p_attr);
    printf("\n");
  }

  if (p_xddm->p_config) {
    printf(" o Configure\n");
    echoElement(p_xddm->p_config);
    printf("\n");
  }

  if (p_xddm->p_inter) {
    printf(" o Intersect\n");
    echoElement(p_xddm->p_inter);
    printf("\n");
  }

  printf(" o Number of variables = %d\n",p_xddm->nv);
  for (i=0; i<p_xddm->nv; i++) {
    printf("\n   %4d",i);
    echoVariable(&p_xddm->a_v[i]);
  }

  printf("\n");
  printf(" o Number of constants = %d\n",p_xddm->nc);
  for (i=0; i<p_xddm->nc; i++) {
    printf("\n   %4d",i);
    echoVariable(&p_xddm->a_c[i]);
  }

  printf("\n");
  printf(" o Number of analysis parameters = %d\n",p_xddm->na);
  for (i=0; i<p_xddm->na; i++) {
    printf("\n   %4d",i);
    echoAnalysis(&p_xddm->a_ap[i]);
  }

  printf("\n");
  printf(" o Number of design points = %d\n",p_xddm->nd);
  for (i=0; i<p_xddm->nd; i++) {
    printf("\n   %4d",i);
    echoDesignPoint(&p_xddm->a_dp[i]);
  }

  printf("\n");
  printf(" o Number of components = %d\n",p_xddm->ncmp);
  for (i=0; i<p_xddm->ncmp; i++) {
    printf("\n   %4d",i);
    echoComponent(&p_xddm->a_cmp[i]);
  }

  printf("\n");
  printf(" o Number of tessellate elements = %d\n",p_xddm->ntess);
  for (i=0; i<p_xddm->ntess; i++) {
    printf("\n   %4d",i);
    if (NULL != p_xddm->a_tess[i].p_id) printf(" ID=%s",p_xddm->a_tess[i].p_id);
    printf("\n");
    if (p_xddm->a_tess[i].lin == 1) printf("        Linearization Required\n");
    if (NULL != p_xddm->a_tess[i].p_comment) printf("        %s\n",p_xddm->a_tess[i].p_comment);
    if ( 0 < p_xddm->a_tess[i].nAttr ) {
      printf("        Attributes: %d\n",p_xddm->a_tess[i].nAttr);
      for (j=0; j<p_xddm->a_tess[i].nAttr; j++) {
        printf("          %s %s\n", p_xddm->a_tess[i].p_attr[j].p_name, p_xddm->a_tess[i].p_attr[j].p_value);
      }
    }
  }
  printf("\n");

  return;
}

int
xddm_addAeroFunForce(p_tsXddmAFun p_afun,
                     const char *const name,
                     const int force,
                     const int frame,
                     const int J,
                     const int N,
                     const double target,
                     const double weight,
                     const int bnd,
                     const char *const comp)
{
  char tmp[512];
  const char *cmp = comp;

  if (p_afun->pa_text == NULL) {
    /* add the header information */
    p_afun->pa_text = malloc(3*sizeof(*p_afun->pa_text));
    if (p_afun->pa_text == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
    p_afun->pa_text[0] = fillString("#         Name    Force   Frame    J      N    Target   Weight  Bound  GMP Comp");
    p_afun->pa_text[1] = fillString("#        (String) (0,1,2) (0,1) (0,1,2) (int)  (dble)   (dble)   (0)");
    p_afun->pa_text[2] = fillString("#------------------------------------------------------------------------------");
    p_afun->nt = 3;
    if (p_afun->pa_text[0] == NULL ||
        p_afun->pa_text[1] == NULL ||
        p_afun->pa_text[2] == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
  }

  if (comp == NULL) {
    cmp = "entire";
  }

  p_afun->pa_text = realloc(p_afun->pa_text, (p_afun->nt+1)*sizeof(*p_afun->pa_text));
  if (p_afun->pa_text == NULL) {
    ERR("malloc failed\n");
    return -1;
  }

  snprintf(tmp, 512, "optForce   %7s   %d      %d      %d      %d      %6g  %6g   %d    %s",
                     name, force, frame, J, N, target, weight, bnd, cmp);

  p_afun->pa_text[p_afun->nt] = fillString(tmp);
  p_afun->nt++;

  return 0;
}

int
xddm_addAeroFunMoment_Point(p_tsXddmAFun p_afun,
                            const char *const name,
                            const int index,
                            const int moment,
                            const int frame,
                            const int J,
                            const int N,
                            const double target,
                            const double weight,
                            const int bnd,
                            const char *const comp)
{
  char tmp[512];
  const char *cmp = comp;

  if (p_afun->pa_text == NULL) {
    /* add the header information */
    p_afun->pa_text = malloc(3*sizeof(*p_afun->pa_text));
    if (p_afun->pa_text == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
    p_afun->pa_text[0] = fillString("#                  Name   Index  Moment  Frame   J     N   Target  Weight  Bound  GMP_Comp");
    p_afun->pa_text[1] = fillString("#                (String) (int) (0,1,2)  (0,1) (0,1) (int) (dble)  (dble)  (0)");
    p_afun->pa_text[2] = fillString("#---------------------------------------------------------------------------------------");
    p_afun->nt = 3;
    if (p_afun->pa_text[0] == NULL ||
        p_afun->pa_text[1] == NULL ||
        p_afun->pa_text[2] == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
  }

  if (comp == NULL) {
    cmp = "entire";
  }

  p_afun->pa_text = realloc(p_afun->pa_text, (p_afun->nt+1)*sizeof(*p_afun->pa_text));
  if (p_afun->pa_text == NULL) {
    ERR("malloc failed\n");
    return -1;
  }

  snprintf(tmp, 512, "optMoment_Point %7s    %d      %d        %d      %d    %d  %6g %6g     %d    %s",
                     name, index, moment, frame, J, N, target, weight, bnd, cmp);

  p_afun->pa_text[p_afun->nt] = fillString(tmp);
  p_afun->nt++;

  return 0;
}

int
xddm_addAeroFunLoD(p_tsXddmAFun p_afun,
                   const char *const name,
                   const int frame,
                   const int J,
                   const int N,
                   const double A,
                   const double bias,
                   const double target,
                   const double weight,
                   const int bnd,
                   const char *const comp)
{
  char tmp[512];
  const char *cmp = comp;

  if (p_afun->pa_text == NULL) {
    /* add the header information */
    p_afun->pa_text = malloc(6*sizeof(*p_afun->pa_text));
    if (p_afun->pa_text == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
    p_afun->pa_text[0] = fillString("# L/D -> SIGN(CL)*ABS(CL)^A/(CD+Bias) in Aero Frame");
    p_afun->pa_text[1] = fillString("#     -> SIGN(CN)*ABS(CN)^A/(CA+Bias) in Body Frame");
    p_afun->pa_text[2] = fillString("# Format:");
    p_afun->pa_text[3] = fillString("#      Name   Frame   J     N     A     Bias  Target  Weight  Bound  GMP_Comp");
    p_afun->pa_text[4] = fillString("#    (String) (0,1) (0,1) (int) (dble) (dble) (dble)  (dble)   (0)");
    p_afun->pa_text[5] = fillString("#----------------------------------------------------------------------------");
    p_afun->nt = 6;
    if (p_afun->pa_text[0] == NULL ||
        p_afun->pa_text[1] == NULL ||
        p_afun->pa_text[2] == NULL ||
        p_afun->pa_text[3] == NULL ||
        p_afun->pa_text[4] == NULL ||
        p_afun->pa_text[5] == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
  }

  if (comp == NULL) {
    cmp = "entire";
  }

  p_afun->pa_text = realloc(p_afun->pa_text, (p_afun->nt+1)*sizeof(*p_afun->pa_text));
  if (p_afun->pa_text == NULL) {
    ERR("malloc failed\n");
    return -1;
  }

  snprintf(tmp, 512, "optLD  %7s   %d      %d   %d   %6g  %6g   %6g    %6g   %d    %s",
                     name, frame, J, N, A, bias, target, weight, bnd, cmp);

  p_afun->pa_text[p_afun->nt] = fillString(tmp);
  p_afun->nt++;

  return 0;
}

#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

#define MY_ENCODING "ISO-8859-1"

int
xddm_writeAttr(xmlTextWriterPtr writer, p_tsXmParent p_node)
{
  int j, rc;

  rc = xmlTextWriterStartElement(writer, BAD_CAST p_node->name);
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  for (j=0; j<p_node->nAttr; j++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_node->p_attr[j].p_name,
                                             BAD_CAST p_node->p_attr[j].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  /* do not end element here in case we need to nest another one */

  return 0;
}

int
xddm_writeElement(xmlTextWriterPtr writer, const char *name, const p_tsXddmElem p_e)
{
  int i, rc;

  rc = xmlTextWriterStartElement(writer, BAD_CAST name);
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  if (p_e->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_e->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  for (i=0; i<p_e->nAttr; i++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_e->p_attr[i].p_name,
                                             BAD_CAST p_e->p_attr[i].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  rc = xmlTextWriterEndElement(writer);         /* end configure */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  return 0;
}

int
xddm_writeDouble(xmlTextWriterPtr writer, const char *p_name, const double val)
{
  int rc;

  if (UNSET != val) {
    if ( 0 != xmlXPathIsInf(val) ) {
      ERR("Value is INFINITE\n");
      return -1;
    }
    if ( 1 == xmlXPathIsNaN(val) ) {
      ERR("Value is NaN\n");
      return -1;
    }

    if ( 0 == strcasecmp(p_name, "value") ) {
      rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST p_name,"%.17e", val);
      if (rc < 0) {
        ERR("xmlTextWriterWriteFormatAttribute failed\n");
        return -1;
      }
    }
    else {
      rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST p_name,"%g", val);
      if (rc < 0) {
        ERR("xmlTextWriterWriteFormatAttribute failed\n");
        return -1;
      }
    }
  }

  return 0;
}

static int
xddm_writeSensitivity(xmlTextWriterPtr writer, const int ndvs, const double *a_lin, char **pa_dvs)
{
  int i, rc;

  if ( ndvs > 0 ) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "SensitivityArray");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    for (i=0; i<ndvs; i++) {
      rc = xmlTextWriterStartElement(writer, BAD_CAST "Sensitivity");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return -1;
      }

      rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "P", BAD_CAST pa_dvs[i]);
      if (rc < 0) {
        ERR("xmlTextWriterWriteAttribute failed\n");
        return -1;
      }

      rc = xddm_writeDouble(writer, "Value", a_lin[i]);
      if (rc) {
        ERR("xddm_writeDouble failed\n");
        return -1;
      }

      rc = xmlTextWriterEndElement(writer);         /* end sensitivity */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
    }

    rc = xmlTextWriterEndElement(writer);         /* end sensitivity array */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  return 0;
}

int
xddm_writeVariable(xmlTextWriterPtr writer, const p_tsXddmVar p_v)
{
  int rc;

  if (p_v->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_v->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_v->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  rc = xddm_writeDouble(writer, "Min", p_v->minVal);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  rc = xddm_writeDouble(writer, "Max", p_v->maxVal);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  rc = xddm_writeDouble(writer, "TypicalSize", p_v->typicalSize);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  if (p_v->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_v->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  return 0;
}


static int
xddm_writeFunctional(xmlTextWriterPtr writer, const p_tsXddmFun p_f)
{
  int rc, j;

  if (p_f->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_f->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_f->p_expr) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Expr", BAD_CAST p_f->p_expr);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_f->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  for (j=0; j<p_f->nAttr; j++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_f->p_attr[j].p_name,
                                             BAD_CAST p_f->p_attr[j].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_f->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_f->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  rc = xddm_writeSensitivity(writer, p_f->ndvs, p_f->a_lin, p_f->pa_dvs);
  if (rc < 0) {
    ERR("xddm_writeSensitivity failed\n");
    return -1;
  }

  return 0;
}

int
xddm_writeAeroFun(xmlTextWriterPtr writer, const p_tsXddmAFun p_afun)
{
  int i, rc;

  if (p_afun->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_afun->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_afun->p_options) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Options", BAD_CAST p_afun->p_options);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_afun->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_afun->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  for (i=0; i<p_afun->nAttr; i++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_afun->p_attr[i].p_name,
                                             BAD_CAST p_afun->p_attr[i].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  for (i=0; i<p_afun->nt; i++) {
    rc = xmlTextWriterWriteFormatString(writer, "\n%s", p_afun->pa_text[i]);
    if (rc < 0) {
      ERR("xmlTextWriterWriteFormatString failed\n");
      return -1;
    }
  }

  return 0;
}

int
xddm_writeAnalysis(xmlTextWriterPtr writer, const p_tsXddmAPar p_a)
{
  int rc;

  if (p_a->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_a->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_a->val);
  if (rc < 0) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  rc = xddm_writeDouble(writer, "DiscretizationError", p_a->discrErr);
  if (rc < 0) {
    ERR("xddm_writeDouble failed\n");
    return -1;
  }

  if (p_a->lin == 1) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity", BAD_CAST "Required");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }
  else if (p_a->lin == 0) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity", BAD_CAST "None");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_a->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_a->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  /* AeroFun Element */
  if (p_a->p_afun) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "AeroFun");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeAeroFun(writer, p_a->p_afun);
    if (rc < 0) {
      ERR("xddm_writeAeroFun failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end AeroFun */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  /* Sensitivities */
  if (p_a->ndvs > 0) {
    rc = xddm_writeSensitivity(writer, p_a->ndvs, p_a->a_lin, p_a->pa_dvs);
    if (rc < 0) {
      ERR("xddm_writeSensitivity failed\n");
      return -1;
    }
  }

  return 0;
}

static int
xddm_writeDesignPoint(xmlTextWriterPtr writer, const p_tsXddmDesP p_dp)
{
  int i, rc;

  if (p_dp->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_dp->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_dp->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_dp->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_dp->p_geometry) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Geometry", BAD_CAST p_dp->p_geometry);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  for (i=0; i<p_dp->nAttr; i++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_dp->p_attr[i].p_name,
                                             BAD_CAST p_dp->p_attr[i].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if ( p_dp->nv > 0 ) {
    for (i=0; i<p_dp->nv; i++) {
      rc = xmlTextWriterStartElement(writer, BAD_CAST "Variable");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return -1;
      }

      rc = xddm_writeVariable(writer, &p_dp->a_v[i]);
      if (rc < 0) {
        ERR("xddm_writeVariable failed\n");
        return -1;
      }
      
      rc = xmlTextWriterEndElement(writer);         /* end variable */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
    }
  }

  if ( p_dp->nc > 0 ) {
    for (i=0; i<p_dp->nc; i++) {
      rc = xmlTextWriterStartElement(writer, BAD_CAST "Constant");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return -1;
      }

      rc = xddm_writeVariable(writer, &p_dp->a_c[i]);
      if (rc < 0) {
        ERR("xddm_writeVariable failed\n");
        return -1;
      }
      
      rc = xmlTextWriterEndElement(writer);         /* end constant */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
    }
  }

  if ( p_dp->na > 0 ) {
    for (i=0; i<p_dp->na; i++) {

      rc = xmlTextWriterStartElement(writer, BAD_CAST "Analysis");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return -1;
      }

      rc = xddm_writeAnalysis(writer, &p_dp->a_ap[i]);
      if (rc < 0) {
        ERR("xddm_writeAnalysis failed\n");
        return -1;
      }
      
      rc = xmlTextWriterEndElement(writer);         /* end analysis */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
    }
  }

  if ( p_dp->p_obj ) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Objective");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeFunctional(writer, p_dp->p_obj);
    if (rc < 0) {
      ERR("xddm_writeFunctional failed\n");
      return -1;
    }
    
    rc = xmlTextWriterEndElement(writer);         /* end objective */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  if ( p_dp->ncr > 0 ) {
    for (i=0; i<p_dp->ncr; i++) {

      rc = xmlTextWriterStartElement(writer, BAD_CAST "Constraint");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return -1;
      }

      rc = xddm_writeFunctional(writer, &p_dp->a_cr[i]);
      if (rc < 0) {
        ERR("xddm_writeFunctional failed\n");
        return -1;
      }
      
      rc = xmlTextWriterEndElement(writer);         /* end constraint */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
    }
  }

  return 0;
}

static int
xddm_writeComponent(xmlTextWriterPtr writer, const p_tsXddmComp p_cmp)
{
  int i, rc;

  if (p_cmp->p_name) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name", BAD_CAST p_cmp->p_name);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_cmp->p_parent) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Parent", BAD_CAST p_cmp->p_parent);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_cmp->p_type) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type", BAD_CAST p_cmp->p_type);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_cmp->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_cmp->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  for (i=0; i<p_cmp->nAttr; i++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_cmp->p_attr[i].p_name,
                                             BAD_CAST p_cmp->p_attr[i].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_cmp->p_data) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Data");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xmlTextWriterWriteFormatString(writer, "%s", p_cmp->p_data);
    if (rc < 0) {
      ERR("xmlTextWriterWriteFormatString failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end data */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  return 0;
}

static int
xddm_writeTessellate(xmlTextWriterPtr writer, const p_tsXmTess p_tess)
{
  int i, rc;

  if (p_tess->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_tess->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_tess->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_tess->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  if (p_tess->lin == 1) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity", BAD_CAST "Required");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }
  else if (p_tess->lin == 0) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity", BAD_CAST "None");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  for (i=0; i<p_tess->nAttr; i++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_tess->p_attr[i].p_name,
                                             BAD_CAST p_tess->p_attr[i].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  return 0;
}

int
xddm_writeFile(const char *const p_fileName, const p_tsXddm p_xddm, const int options)
{
  xmlTextWriterPtr writer;
  xmlBufferPtr buf;

  int i, fd, rc;
  FILE *p_strm;

  { LIBXML_TEST_VERSION };

  /* create a new XML buffer, to which the XML document will be written */
  buf = xmlBufferCreate();
  if (buf == NULL) {
    ERR("xmlBufferCreate failed to create xml buffer\n");
    return -1;
  }

  /* create a new XmlWriter for memory, with no compression */
  writer = xmlNewTextWriterMemory(buf, 0);
  if (writer == NULL) {
    ERR("xmlNewTextWriterMemory failed to create xml writer\n");
    return -1;
  }

  /* start the document with the xml default and encoding ISO 8859-1 */
  rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
  if (rc < 0) {
    ERR("xmlTextWriterStartDocument failed\n");
    return -1;
  }

  rc = xmlTextWriterSetIndent(writer, 1); /* pretty indent */
  if (rc < 0) {
    ERR("xmlTextWriterSetIndent failed\n");
    return -1;
  }

  /* write parent element */
  rc = xddm_writeAttr(writer, p_xddm->p_parent);
  if (rc < 0) {
    ERR("xddm_writeAttr p_xddm->p_parent failed\n");
    return -1;
  }

  if (p_xddm->p_config) {
    rc = xddm_writeElement(writer, "Configure", p_xddm->p_config);
    if (rc < 0) {
      ERR("xddm_writeElement failed\n");
      return -1;
    }
  }

  if (p_xddm->p_inter) {
    rc = xddm_writeElement(writer, "Intersect", p_xddm->p_inter);
    if (rc < 0) {
      ERR("xddm_writeElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->nv; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "Variable");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeVariable(writer, &p_xddm->a_v[i]);
    if (rc < 0) {
      ERR("xddm_writeVariable failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end variable */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->nc; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "Constant");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeVariable(writer, &p_xddm->a_c[i]);
    if (rc < 0) {
      ERR("xddm_writeVariable failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end constant */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->na; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "Analysis");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeAnalysis(writer, &p_xddm->a_ap[i]);
    if (rc < 0) {
      ERR("xddm_writeAnalysis failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end analysis */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->nd; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "DesignPoint");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeDesignPoint(writer, &p_xddm->a_dp[i]);
    if (rc < 0) {
      ERR("xddm_writeDesignPoint failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end design point */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->ncmp; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "Component");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeComponent(writer, &p_xddm->a_cmp[i]);
    if (rc < 0) {
      ERR("xddm_writeComponent failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end component */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->nafun; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "AeroFun");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeAeroFun(writer, &p_xddm->a_afun[i]);
    if (rc < 0) {
      ERR("xddm_writeAeroFun failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end AeroFun */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  for (i=0; i<p_xddm->ntess; i++) {

    rc = xmlTextWriterStartElement(writer, BAD_CAST "Tessellate");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xddm_writeTessellate(writer, &p_xddm->a_tess[i]);
    if (rc < 0) {
      ERR("xddm_writeComponent failed\n");
      return -1;
    }

    rc = xmlTextWriterEndElement(writer);         /* end tessellate */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  }

  rc = xmlTextWriterEndElement(writer);         /* end model element */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  xmlFreeTextWriter(writer);
                                                     /* write buffer to disk */
  p_strm = fopen(p_fileName, "w");
  if (p_strm == NULL) {
    printf("fopen failed\n");
    return -1;
  }
  fprintf(p_strm, "%s", (const char *) buf->content);
  fflush(p_strm);
                                                            /* force to disk */
  fd = fileno(p_strm);
  rc = fsync(fd);
  if ( 0 != rc ) {
    WARN("fsync on mesh io failed\n");
  }
  fclose(p_strm);

  xmlBufferFree(buf);
  xmlCleanupParser(); /* shutdown libxml2 */
  xmlMemoryDump(); /* this is to debug memory for regression tests */

  if ( XDDM_VERBOSE & options ) {
    printf(" o Wrote \'%s\'\n",p_fileName);
  }

  return 0;
}

/**
 * xddm_updateAnalysisParams(): update values and sensitivities of analysis parameters
 */
int
xddm_updateAnalysisParams(const char *const p_fileName, const p_tsXddm p_xddm,
                          const int options)
{
  int parserOpts=0, i, ia;
  size_t len=0;
  char *xpathExpr;

  xmlDocPtr          doc;
  xmlParserCtxtPtr   ctxt;
  xmlXPathContextPtr xpathCtx = NULL;
  xmlXPathObjectPtr  xddmObj  = NULL;
  xmlNodePtr         node, sarray_node, new_node;
  xmlChar           *xmlstr;
  xmlAttrPtr         new_attr;

  FILE *p_f;

  /* initialize and check libxml2 version */
  { LIBXML_TEST_VERSION };

  xmlXPathInit();  /* initialize xpath */

  if ( XDDM_VERBOSE & options ) {
    printf("    o  Parsing file \"%s\" with libxml2\n", p_fileName);
  }

  ctxt = xmlNewParserCtxt(); /* create a document parser context */
  if (ctxt == NULL) {
    ERR(" xddm_updateAnalysisParams failed to allocate parser context\n");
    xmlCleanupParser();
    return -1;
  }

  doc = xmlCtxtReadFile(ctxt, p_xddm->fileName, NULL, parserOpts);
  if (doc == NULL) {
    if ( XDDM_VERBOSE & options ) {
      ERR("%s is not valid XML\n", p_fileName);
    }
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    return -1;
  }

  xmlFreeParserCtxt(ctxt); /* done with document parser context */

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("xmlXPathNewContext failed to create context\n");
    shutdown_libxml(doc);
    return -1;
  }

  if (NULL == p_xddm->xpathExpr) {
    ERR("missing xpath expression\n");
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return -1;
  }

  xddmObj = xmlXPathEvalExpression(BAD_CAST p_xddm->xpathExpr, xpathCtx);
  if (NULL == xddmObj) {
    WARN("no elements found for expression %s\n", p_xddm->xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return -1;
  }

  if ( 0 == countNodes(xddmObj) ) {
    WARN("no elements found for expression \'%s\'\n", p_xddm->xpathExpr);
    xmlXPathFreeObject(xddmObj);
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return -1;
  }

  xpathCtx->node = xddmObj->nodesetval->nodeTab[0];
  xmlXPathFreeObject(xddmObj);

  for (ia=0; ia<p_xddm->na; ia++) {
    const p_tsXddmAPar p_a = &p_xddm->a_ap[ia];

    len  = strlen("./Analysis[@ID=\"\"]") + 1;
    len += strlen(p_a->p_id);

    xpathExpr = malloc( len * sizeof(*xpathExpr) );
    if ( NULL == xpathExpr ) {
      ERR("malloc failed");
      return -1;
    }
    memset((void*) xpathExpr,'\0',len);

    snprintf(xpathExpr,len,"./Analysis[@ID=\"%s\"]", p_a->p_id);

    xddmObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);

    if (NULL == xddmObj) {
      WARN("no elements found for expression %s\n", xpathExpr);
      xmlXPathFreeContext(xpathCtx);
      shutdown_libxml(doc);
      return -1;
    }

    if ( 1 != countNodes(xddmObj) ) {
      ERR("analysis parameter not unique\n");
      xmlXPathFreeContext(xpathCtx);
      xmlXPathFreeObject(xddmObj);
      shutdown_libxml(doc);
      return -1;
    }

    if (xpathExpr) free(xpathExpr);
    xpathExpr = NULL;

    node = xddmObj->nodesetval->nodeTab[0];

    xmlstr = xmlXPathCastNumberToString(p_a->val);
    xmlNewProp(node, BAD_CAST "Value", xmlstr);
    xmlFree(xmlstr);

                                                        /* add sensitivities */
    if (p_a->ndvs > 0) {
      sarray_node = xmlNewChild(node, NULL, BAD_CAST "SensitivityArray", NULL);
      for (i=0; i<p_a->ndvs; i++) {
        new_node = xmlNewChild(sarray_node, NULL, BAD_CAST "Sensitivity", NULL);
        if (NULL==new_node) {
          ERR("xmlNewChild failed\n");
          exit(1);
        }
        new_attr = xmlNewProp(new_node, BAD_CAST "P", BAD_CAST p_a->pa_dvs[i]);
        if (NULL==new_attr) {
          ERR("xmlNewProp failed\n");
          exit(1);
        }
        xmlstr   = xmlXPathCastNumberToString(p_a->a_lin[i]);
        if (NULL==xmlstr) {
          ERR("xmlXPathCastNumberToString failed\n");
          exit(1);
        }
        new_attr = xmlNewProp(new_node, BAD_CAST "Value", xmlstr);
        if (NULL==new_attr) {
          ERR("xmlNewProp failed\n");
          exit(1);
        }
        xmlFree(xmlstr);
      }
    }

    xmlXPathFreeObject(xddmObj);
  }

  xmlXPathFreeContext(xpathCtx);

  p_f = fopen("model.output.xml", "w");
  if (p_f == NULL) {
    printf("fopen failed\n");
    return(-1);
  }

  xmlKeepBlanksDefault(0);
  xmlDocFormatDump(p_f, doc, 1);
  xmlSaveFormatFile("fileout.xml", doc, 1);
  fclose(p_f);

  /* all done, clean up and return
   */

  shutdown_libxml(doc);

  return 0;
}

#undef MY_ENCODING

#else

int
xddm_writeFile(const char *const p_fileName, const int options, p_tsXddm p_xddm);
{
  ERR("LIBXML2 writer not compiled\n");
  return;
}

#endif

/**
 * xddm_alloc(): constructor of XDDM datastructure with basic initialization
 */
p_tsXddm
xddm_alloc(void)
{
  p_tsXddm p_xddm = NULL; /* return object */

  p_xddm = malloc(sizeof(*p_xddm));
  if (NULL == p_xddm) {
    ERR("malloc failed\n");
    exit(1);
  }

  p_xddm->fileName  = NULL;
  p_xddm->xpathExpr = NULL;

  p_xddm->p_parent = xddm_allocParent();

  p_xddm->p_config = NULL;
  p_xddm->p_inter  = NULL;

  p_xddm->nv       = 0;
  p_xddm->nc       = 0;
  p_xddm->na       = 0;
  p_xddm->nd       = 0;
  p_xddm->ncmp     = 0;
  p_xddm->nafun    = 0;
  p_xddm->ntess    = 0;
  p_xddm->a_v      = NULL;
  p_xddm->a_c      = NULL;
  p_xddm->a_ap     = NULL;
  p_xddm->a_dp     = NULL;
  p_xddm->a_cmp    = NULL;
  p_xddm->a_afun   = NULL;
  p_xddm->a_tess   = NULL;

  assert(p_xddm);

  return p_xddm;
}

/**
 * xddm_free(): destructor of XDDM datastructure
 */
void
xddm_free(p_tsXddm p_xddm)
{
  if (p_xddm) {
    if (p_xddm->fileName ) free(p_xddm->fileName);
    if (p_xddm->xpathExpr) free(p_xddm->xpathExpr);

    freeParent(p_xddm->p_parent);
    freeElement(p_xddm->p_config);
    freeElement(p_xddm->p_inter);
    freeVariable(p_xddm->nv, p_xddm->a_v);
    freeVariable(p_xddm->nc, p_xddm->a_c);
    freeAnalysis(p_xddm->na, p_xddm->a_ap);
    freeDesignPoint(p_xddm->nd, p_xddm->a_dp);
    freeComponent(p_xddm->ncmp, p_xddm->a_cmp);
    freeAeroFun(p_xddm->nafun, p_xddm->a_afun);
    freeTessellate(p_xddm->ntess, p_xddm->a_tess);

    free(p_xddm);
  }

  return;
}
