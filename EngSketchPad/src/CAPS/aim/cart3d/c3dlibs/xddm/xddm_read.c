
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: xddm_read.c,v 1.3 2022/11/07 23:01:39 mnemec Exp $
 */

/* open source */

/**
 * XML file reader that parses the
 * Extensible-Design-Description-Markup (XDDM) elements.
 *
 * Dependency: libxml2
 * https://gitlab.gnome.org/GNOME/libxml2/-/wikis/homewww.xmlsoft.org
 * This library is usually present on most systems, check via
 * existence of 'xml2-config' script.
 */

#include "xddm.h"
#include "xddmInternals.h"

void
xddm_freeElement(p_tsXddmElem p_e);

/**
 * Local prototypes
 */
static int xddm_readVariable(xmlXPathObjectPtr xddmObj, size_t *p_nv,
                             p_tsXddmVar *a_v, const unsigned options);
static int xddm_readConfigure(xmlXPathObjectPtr xddmObj, unsigned *p_config);
static int xddm_readElement(xmlXPathObjectPtr xddmObj, size_t *p_ne,
                            p_tsXddmElem *a_e, const unsigned options);
static int xddm_parseBasicTypes(xmlXPathContextPtr xpathCtx, p_tsXddm p_xddm,
                                const unsigned options);
static int xddm_readAnalysis(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
                             size_t *p_na, p_tsXddmAPar *a_ap,
                             const unsigned options);
static int xddm_readFunction(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
                             size_t *p_nf, p_tsXddmFun *a_f,
                             const unsigned options);
static int xddm_readAeroFun(xmlXPathObjectPtr aero, p_tsXddmAFun *pp_afun,
                            size_t *p_naf, const unsigned options);
static int xddm_readIntersect(xmlXPathObjectPtr xddmObj, p_tsXddmGeom *pp_geom,
                              size_t *p_ng, const unsigned options);
static int xddm_readSum(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
                        size_t *p_ns, p_tsXddmSum *a_s, const unsigned options);
static int xddm_fillElement(p_tsXddmElem p_e, const xmlChar *p_name,
                            const xmlChar *p_value);
static int xddm_addXmlAttribute(const xmlChar * p_name, const xmlChar * p_value,
                                size_t *p_nAttr, p_tsXddmAttr *pp_attr);
static int xddm_readXDDM(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
                         size_t *p_ndp, p_tsXddm *a_dp, const unsigned options);

/**
 * xddm_readFile(): parse an XPath expression from an XDDM file and return
 * structure containing all elements in the path
 */
p_tsXddm
xddm_readFile(const char *const p_fileName, const char *const xpathExpr,
              unsigned *p_options)
{
  xmlDocPtr          doc      = NULL;
  xmlParserCtxtPtr   ctxt     = NULL;
  xmlXPathObjectPtr  xddmObj  = NULL;
  xmlXPathContextPtr xpathCtx = NULL;

  p_tsXddm p_xddm  = NULL; /* return object */

  unsigned libxmlOpts=0;

  size_t np=0;
  int rc=0;
                                     /* initialize and check libxml2 version */
  { LIBXML_TEST_VERSION };

  xmlXPathInit();                                        /* initialize xpath */

  if ( XDDM_VERBOSE & *p_options ) {
    printf(" o Parsing file \"%s\" with libxml2\n", p_fileName);
  }

  ctxt = xmlNewParserCtxt();             /* create a document parser context */
  if (ctxt == NULL) {
    ERR(" xddm_readFile failed to allocate parser context\n");
    xmlCleanupParser();
    return NULL;
  }

  doc = xmlCtxtReadFile(ctxt, p_fileName, NULL, libxmlOpts);
  if (doc == NULL) {
    if ( XDDM_VERBOSE & *p_options ) {
      ERR("%s is not valid XML\n", p_fileName);
    }
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    return NULL;
  }

  xmlFreeParserCtxt(ctxt);              /* done with document parser context */

  xpathCtx = xmlXPathNewContext(doc);     /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("xmlXPathNewContext failed to create context\n");
    rc = 1;
    goto cleanup;
  }
       /* looks like we have a valid xml document that we are ready to parse */
  xddmObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if (NULL == xddmObj) {
    ERR("Unable to parse expression \'%s\'\n", xpathExpr);
    rc = 1;
    goto cleanup;
  }

  if ( xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
    printf(" o No elements found for expression \'%s\'\n", xpathExpr);
    rc = 0;
    goto cleanup;
  }
                                   /* xpath returned valid nodes, parse xddm */
  p_xddm = xddm_new(1);

  if (p_xddm->p_e) {
    xddm_freeElement(p_xddm->p_e);
    free(p_xddm->p_e);
    p_xddm->p_e = NULL;
  }

  rc = xddm_readElement(xddmObj, &np, &p_xddm->p_e, *p_options);
  if ( rc ) {
    ERR("xddm_readElement failed for root node\n");
    goto cleanup;
  }

  if ( 1 != np ) {
    WARN("Expecting one root element for '%s'\n",xpathExpr);
    rc = 1;
    goto cleanup;
  }
                                      /* set xpath root to xddm root element */
  xpathCtx->node = xddmObj->nodesetval->nodeTab[0];

  xmlXPathFreeObject(xddmObj);
  xddmObj = NULL;

  /* ---------------
   * parse Configure
   * ---------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "//Configure", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readConfigure(xddmObj, p_options);
      if ( rc ) {
        ERR("xddm_readConfigure failed\n");
        goto cleanup;
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* --------------------------
   * parse all XDDM basic types
   * --------------------------
   */

  rc = xddm_parseBasicTypes(xpathCtx, p_xddm, *p_options);
  if (rc) {
    ERR("xddm_parseBasicTypes failed\n");
    goto cleanup;
  }

  /* --------------------------
   * parse DesignPoint elements
   * --------------------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./DesignPoint", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readXDDM(xddmObj, xpathCtx, &p_xddm->nk, &p_xddm->a_kids,
                         *p_options);
      if ( rc ) {
        ERR("xddm_readXDDM failed for DesignPoints\n");
        goto cleanup;
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ------------------------------
   * Optional: read everything else
   * ------------------------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./*[not(self::Configure or self::Variable or self::Constant or self::Analysis or self::Function or self::Sum or self::Objective or self::Constraint or self::Intersect or self::Tessellate or self::DesignPoint )]", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      int i;
      xmlNodeSetPtr nodes = xddmObj->nodesetval;
      int nNodes= (nodes) ? nodes->nodeNr : 0;
      for (i=0; i < nNodes; i++) {
        if ( nodes->nodeTab[i]->type != XML_ELEMENT_NODE ) continue;
        printf("node name %s\n",nodes->nodeTab[i]->name);
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

 cleanup:
                                            /* all done, clean up and return */
  if (0 != rc) {
    xddm_free(p_xddm, XDDM_VERBOSE);
    p_xddm = NULL;
  }

  if (xddmObj)  xmlXPathFreeObject(xddmObj);
  if (xpathCtx) xmlXPathFreeContext(xpathCtx);
  if (doc)      xmlFreeDoc(doc);
  shutdown_libxml();

  return p_xddm;
}

/**
 * xddm_readConfigure(): parse Configure node
 */
static int
xddm_readConfigure(xmlXPathObjectPtr xddmObj, unsigned *p_config)
{
  size_t inode=0, n=0;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  n = countNodes(xddmObj);

  if ( n != 1 ) {
    ERR("problem in number of Configure elements\n");
    return 1;
  }
                                                     /* parse node attributes */
  pp_xmlAttr = xddmParseXpathObj(xddmObj);
  if ( NULL == pp_xmlAttr ) {
    ERR("Problem parsing Configure\n");
    return 1;
  }

  inode = 0;                     /* loop over each node, sentinel is NULL */
  while (NULL != pp_xmlAttr[inode]) {
    if (inode >= n) {
      ERR("Overflow in number of Configure elements\n");
      freeXddmAttr(pp_xmlAttr);
      return 1;
    }
    else {
      int i;
      const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
                               /* loop over attributes and save relevant ones */
      for (i=0; i<p_xA->n; i++) {
        if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "sensitivity" ) ) {
          if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "required" ) ) {
            *p_config |= (XDDM_LIN);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "none" ) ) {
            *p_config |= (XDDM_NOLIN);
          }
        }
      } /* end for attributes */
    }
    inode++;
  } /* end while attributes */
  freeXddmAttr(pp_xmlAttr);
  pp_xmlAttr = NULL;

  return 0;
}

/**
 * xddm_parseBasicTypes(): parse basic XDDM elements
 */
static int
xddm_parseBasicTypes(xmlXPathContextPtr xpathCtx, p_tsXddm p_xddm,
                     const unsigned options)
{
  int rc = 0;
  xmlXPathObjectPtr xddmObj  = NULL;

  /* ---------------
   * parse Variables
   * ---------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Variable", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readVariable(xddmObj, &p_xddm->nv, &p_xddm->a_v, options);
      if ( rc ) {
        ERR("xddm_readVariable Variable failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ---------------
   * parse Constants
   * ---------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Constant", xpathCtx);
  if (xddmObj) {
   if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
     rc = xddm_readVariable(xddmObj, &p_xddm->nc, &p_xddm->a_c, options);
     if ( rc ) {
       ERR("xddm_readVariable Constant failed\n");
     }
   }
   xmlXPathFreeObject(xddmObj);
   xddmObj = NULL;
  }

  /* -------------------------
   * parse Analysis parameters
   * -------------------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Analysis", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readAnalysis(xddmObj, xpathCtx, &p_xddm->na, &p_xddm->a_ap,
                             options);
      if ( rc ) {
        ERR("xddm_readAnalysis failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ---------------
   * parse Functions
   * ---------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Function", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readFunction(xddmObj, xpathCtx, &p_xddm->nf, &p_xddm->a_f,
                             options);
      if ( rc ) {
        ERR("xddm_readFunction failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ----------
   * parse Sums
   * ----------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Sum", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readSum(xddmObj, xpathCtx, &p_xddm->ns, &p_xddm->a_s, options);
      if ( rc ) {
        ERR("xddm_readFunction failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ----------------
   * parse Objectives
   * ----------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Objective", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readFunction(xddmObj, xpathCtx, &p_xddm->nj, &p_xddm->a_j,
                             options);
      if ( rc ) {
        ERR("xddm_readFunction failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* -----------------
   * parse Constraints
   * -----------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Constraint", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readFunction(xddmObj, xpathCtx, &p_xddm->ncon, &p_xddm->a_con,
                             options);
      if ( rc ) {
        ERR("xddm_readFunction failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ------------------------
   * parse Intersect elements
   * ------------------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Intersect", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readIntersect(xddmObj, &p_xddm->a_geo, &p_xddm->ng, options);
      if ( rc ) {
        ERR("xddm_readIntersect failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ----------------------
   * parse AeroFun elements
   * ----------------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./AeroFun", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readAeroFun(xddmObj, &p_xddm->a_afun, &p_xddm->naf, options);
      if ( rc ) {
        ERR("xddm_readAeroFun failed\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  /* ----------------
   * parse Tessellate
   * ----------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Tessellate", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      rc = xddm_readElement(xddmObj, &p_xddm->nt, &p_xddm->a_t, options);
      if ( rc ) {
        ERR("xddm_readElement failed for tessellate node\n");
      }
    }
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  return rc;
}

/**
 * xddm_readXDDM(): read any XDDM nodes
 */
static int
xddm_readXDDM(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
              size_t *p_ndp, p_tsXddm *a_dp, const unsigned options)
{
  xmlNodePtr p_nodeCur = xpathCtx->node; /* save current node to return */
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  size_t inode;
  int rc=0;

  assert(xddmObj);

  *p_ndp = countNodes(xddmObj);
  if ( 0 == *p_ndp ) {
    ERR("Problem parsing DesignPoint parameters\n");
    return 1;
  }
                                    /* parse node attributes of all elements */
  pp_xmlAttr = xddmParseXpathObj(xddmObj);

  (*a_dp) = xddm_new(*p_ndp);

  inode = 0;                        /* loop over each node, sentinel is NULL */
  while (NULL != pp_xmlAttr[inode]) {
    if (inode >= (*p_ndp)) {
      ERR("overflow in number of design point params\n");
      freeXddmAttr(pp_xmlAttr);
      return 1;
    }
    else {
      int i;
      p_tsXddm               p_d = (*a_dp) + inode;
      const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

      p_d->p_e->p_xn = xddm_getName(p_xA->p_node, options);
      p_d->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
      for (i=0; i<p_xA->n; i++) {
        rc = xddm_fillElement(p_d->p_e, p_xA->a_names[i], p_xA->a_values[i]);
        if (rc) {
          ERR("Failed read element attributes\n");
          freeXddmAttr(pp_xmlAttr);
          return 1;
        }
      } /* end for attributes */

                             /* move xpath node into the design point element */
      xpathCtx->node = p_xA->p_node;

      //      p_d->p_xd = xddm_new(1);

      rc = xddm_parseBasicTypes(xpathCtx, p_d, options);
      if (rc) {
        ERR("xddm_parseBasicTypes Failed to read design point\n");
        freeXddmAttr(pp_xmlAttr);
        xddm_free(p_d, XDDM_VERBOSE);
        p_d = NULL;
        return rc;
      }
    } /* end if */

    inode++;
  } /* end while attributes */

  freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur; /* reset root node */

  return rc;
}

/**
 * xddm_readVariable(): read XDDM variable
 */
static int
xddm_readVariable(xmlXPathObjectPtr xddmObj, size_t *p_nv, p_tsXddmVar *a_v,
                  const unsigned options)
{
  size_t inode=0;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  assert(xddmObj);

  *p_nv = countNodes(xddmObj); /* number of variables */

  if ( *p_nv > 0 ) {
                                    /* parse node attributes of all elements */
    pp_xmlAttr = xddmParseXpathObj(xddmObj);
    if ( NULL == pp_xmlAttr ) {
      ERR("Parsing problem\n");
      return 1;
    }

    (*a_v) = xddm_newVariable(*p_nv);

    inode = 0;                      /* loop over each node, sentinel is NULL */
    while (NULL != pp_xmlAttr[inode]) {
      if (inode >= (*p_nv)) {
        ERR("Overflow in number of Variables or Constants\n");
        freeXddmAttr(pp_xmlAttr);
        return 1;
      }
      else {
        size_t i;
        int rc;
        const p_tsXddmVar p_v      = (*a_v) + inode;
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        p_v->p_e->p_xn = xddm_getName(p_xA->p_node, options);
        p_v->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
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
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "fdstep" ) ) {
            p_v->fdstep = fillDouble(p_xA->a_values[i]);
          }
          else {
            rc = xddm_fillElement(p_v->p_e, p_xA->a_names[i], p_xA->a_values[i]);
            if (rc) {
              ERR("Failed read element attributes\n");
              freeXddmAttr(pp_xmlAttr);
              return 1;
            }
          }
        } /* end for attributes */
      }
      inode++;
    } /* end while attributes */

    freeXddmAttr(pp_xmlAttr);
    pp_xmlAttr = NULL;
  } /* end if variables */

  return 0;
}

static int
xddm_readSensitivity(xmlXPathObjectPtr linObj, size_t *p_ndvs, double **pa_lin,
                     char ***ppa_dvs)
{
  size_t i, j;
  p_tsXddmXmlAttr *pp_sensit  = NULL;
  char **pa_dvs = NULL; /* ragged array of dvs names */
  double *a_lin = NULL; /* derivative values */

  if ( 0 == countNodes(linObj) ) {
    ERR("Syntax error likely in SensitivityArray\n");
    return 1;
  }

  pp_sensit = xddmParseXpathObj(linObj);

  *p_ndvs = 0;
  while (NULL != pp_sensit[*p_ndvs]) {
    (*p_ndvs)++;
  }

  if ( 0 == (*p_ndvs) ) {
    WARN("Syntax error likely in SensitivityArray\n");
    freeXddmAttr(pp_sensit);
    return 1;
  }

  pa_dvs = malloc( (*p_ndvs) * sizeof(*pa_dvs) );
  if (NULL == pa_dvs) {
    ERR("malloc failed\n");
    exit(1);
  }

  a_lin = malloc( (*p_ndvs) * sizeof(*a_lin) );
  if (NULL == a_lin) {
    ERR("malloc failed\n");
    exit(1);
  }

  i=0;
  while (NULL != pp_sensit[i]) {
    const p_tsXddmXmlAttr p_s = pp_sensit[i];

    for (j=0; j<p_s->n; j++) {
      if ( 0 == xmlStrcasecmp( p_s->a_names[j], BAD_CAST "p" ) ) {
        pa_dvs[i] = fillXmlString(p_s->a_values[j]);
      }
      else if ( 0 == xmlStrcasecmp( p_s->a_names[j], BAD_CAST "value" ) ) {
        a_lin[i] = fillDouble(p_s->a_values[j]);
      }
    }
    i++;
  } /* end while sensitivities */

  *ppa_dvs = pa_dvs; /* hook into xddm data struct */
  *pa_lin  = a_lin;

  freeXddmAttr(pp_sensit);

  return 0;
}

static int
xddm_readAeroFun(xmlXPathObjectPtr aero, p_tsXddmAFun *pp_afun, size_t *p_naf,
                 const unsigned options)
{
  size_t inode=0;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  xmlChar *p_xmlstr;
  p_tsXddmAFun p_afun = NULL;

  *p_naf = countNodes(aero);

  if ( *p_naf < 1 ) {
    WARN("Syntax error likely in AeroFun: %zu\n", *p_naf);
    return 1;
  }

  p_afun = xddm_newAeroFun( *p_naf );
  pp_xmlAttr = xddmParseXpathObj(aero);

  inode=0;
  while (NULL != pp_xmlAttr[inode]) {
    size_t i;
    int rc;
    const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
    p_tsXddmAFun           p_f = p_afun + inode;

    p_f->p_e->p_xn = xddm_getName(p_xA->p_node, options);
    p_f->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
    for (i=0; i<p_xA->n; i++) {
      rc = xddm_fillElement(p_f->p_e, p_xA->a_names[i], p_xA->a_values[i]);
      if (rc) {
        ERR("Failed read element attributes\n");
        freeXddmAttr(pp_xmlAttr);
        return 1;
      }
    } /* end for attributes */

    p_xmlstr = xmlNodeGetContent(p_xA->p_node);
    if (NULL == p_xmlstr) {
      ERR("xmlNodeGetContent failed in xddm_readAeroFun\n");
      freeXddmAttr(pp_xmlAttr);
      return 1;
    }

    if (! p_xA->p_node->ns ) {
      p_f->p_text = fillXmlString(p_xmlstr);
    }
    else {
      ERR("Unable to read TEXT node\n");
    }

    xmlFree(p_xmlstr);
    inode++;
  }

  freeXddmAttr(pp_xmlAttr);

  *pp_afun = p_afun; /* hook into xddm data struct */

  return 0;
}

static int
xddm_readAnalysis(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
                  size_t *p_na, p_tsXddmAPar *a_ap, const unsigned options)
{
  xmlNodePtr p_nodeCur = xpathCtx->node; /* save current node to return */
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  size_t inode;
  int rc=0;

  assert(xddmObj);

  *p_na = countNodes(xddmObj); /* number of analysis parameters */
  if ( 0 == *p_na ) {
    ERR("Problem parsing Analysis parameters\n");
    return 1;
  }
                                    /* parse node attributes of all elements */
  pp_xmlAttr = xddmParseXpathObj(xddmObj);

  (*a_ap) = xddm_newAnalysis(*p_na);

  inode = 0;                        /* loop over each node, sentinel is NULL */
  while (NULL != pp_xmlAttr[inode]) {
    if (inode >= (*p_na)) {
      ERR("overflow in number of analysis params\n");
      freeXddmAttr(pp_xmlAttr);
      return 1;
    }
    else {
      size_t i;
      const p_tsXddmAPar    p_a  = (*a_ap) + inode;
      const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
      xmlXPathObjectPtr     xObj = NULL;

      p_a->p_e->p_xn = xddm_getName(p_xA->p_node, options);
      p_a->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
      for (i=0; i<p_xA->n; i++) {
        if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
          p_a->val = fillDouble(p_xA->a_values[i]);
        }
        else if ( 0 == xmlStrcasecmp( p_xA->a_names[i],
                                      BAD_CAST "discretizationerror" ) ) {
          p_a->derr = fillDouble(p_xA->a_values[i]);
        }
        else {
          rc = xddm_fillElement(p_a->p_e, p_xA->a_names[i], p_xA->a_values[i]);
          if (rc) {
            ERR("xddm_fillElement failed\n");
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        }
      } /* end for attributes */
                                                   /* look for sensitivities */
                                /* move xpath node into the analysis element */
      xpathCtx->node = pp_xmlAttr[inode]->p_node;
                                /* check if the node has a sensitivity array */
      xObj = xmlXPathEvalExpression(BAD_CAST "./SensitivityArray/Sensitivity",
                                        xpathCtx);
      if (xObj) {
        if ( ! xmlXPathNodeSetIsEmpty(xObj->nodesetval) ) {
          rc = xddm_readSensitivity(xObj, &p_a->ndvs, &p_a->a_lin, &p_a->pa_dvs);
          if ( rc ) {
            ERR("xddm_readSensitivity failed for %s\n",p_a->p_e->p_id);
            xmlXPathFreeObject(xObj);
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        }
        xmlXPathFreeObject(xObj);
        xObj = NULL;
      }
      xpathCtx->node = p_nodeCur; /* reset root node */

                                                /* look for AeroFun elements */
      xpathCtx->node = pp_xmlAttr[inode]->p_node;
      xObj = xmlXPathEvalExpression(BAD_CAST "./AeroFun", xpathCtx);
      if (xObj) {
        if ( ! xmlXPathNodeSetIsEmpty(xObj->nodesetval) ) {
          size_t naf = 0;
          rc = xddm_readAeroFun(xObj, &p_a->p_afun, &naf, options);
          if ( rc ) {
            ERR("xddm_readAeroFun failed for %s\n",p_a->p_e->p_id);
            xmlXPathFreeObject(xObj);
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
          if ( naf > 1 ) {
            ERR("expecting only one AeroFun node\n");
            return 1;
          }
        }
        xmlXPathFreeObject(xObj);
        xObj = NULL;
      }

      xpathCtx->node = p_nodeCur; /* reset root node */
    } /* end if inode */

    inode++;
  } /* end while attributes */

  freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur; /* reset root node */

  return 0;
}

static int
xddm_readFunction(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
                  size_t *p_nf, p_tsXddmFun *a_f, const unsigned options)
{
  xmlNodePtr p_nodeCur = xpathCtx->node; /* save current node to return */
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  size_t inode;
  int rc=0;

  assert(xddmObj);

  *p_nf = countNodes(xddmObj); /* number of analysis parameters */
  if ( 0 == *p_nf ) {
    ERR("xddm_readFunction countNodes failed\n");
    return 1;
  }
                                    /* parse node attributes of all elements */
  pp_xmlAttr = xddmParseXpathObj(xddmObj);

  (*a_f) = xddm_newFunction(*p_nf);

  inode = 0;                        /* loop over each node, sentinel is NULL */
  while (NULL != pp_xmlAttr[inode]) {
    if (inode >= (*p_nf)) {
      ERR("overflow in number of function params\n");
      freeXddmAttr(pp_xmlAttr);
      return 1;
    }
    else {
      size_t i;
      const p_tsXddmFun     p_f  = (*a_f) + inode;
      const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
      xmlXPathObjectPtr     xObj = NULL;

      p_f->p_e->p_xn = xddm_getName(p_xA->p_node, options);
      p_f->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
      for (i=0; i<p_xA->n; i++) {
        if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
          p_f->val = fillDouble(p_xA->a_values[i]);
        }
        else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "expr" ) ) {
          p_f->p_expr = fillXmlString(p_xA->a_values[i]);
        }
        else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "min" ) ) {
          p_f->min = fillDouble(p_xA->a_values[i]);
        }
        else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "max" ) ) {
          p_f->max = fillDouble(p_xA->a_values[i]);
        }
        else {
          rc = xddm_fillElement(p_f->p_e, p_xA->a_names[i], p_xA->a_values[i]);
          if (rc) {
            ERR("xddm_fillElement failed\n");
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        }
      } /* end for attributes */
                                                   /* look for sensitivities */
                                /* move xpath node into the analysis element */
      xpathCtx->node = pp_xmlAttr[inode]->p_node;
                                /* check if the node has a sensitivity array */
      xObj = xmlXPathEvalExpression(BAD_CAST "./SensitivityArray/Sensitivity",
                                        xpathCtx);
      if (xObj) {
        if ( ! xmlXPathNodeSetIsEmpty(xObj->nodesetval) ) {
          rc = xddm_readSensitivity(xObj, &p_f->ndvs, &p_f->a_lin, &p_f->pa_dvs);
          if ( rc ) {
            ERR("xddm_readSensitivity failed for %s\n",p_f->p_e->p_id);
            xmlXPathFreeObject(xObj);
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        }
        xmlXPathFreeObject(xObj);
        xObj = NULL;
      }
      xpathCtx->node = p_nodeCur; /* reset root node */
    } /* end if inode */

    inode++;
  } /* end while attributes */

  freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur; /* reset root node */

  return 0;
}

static int
xddm_readSum(xmlXPathObjectPtr xddmObj, xmlXPathContextPtr xpathCtx,
             size_t *p_ns, p_tsXddmSum *a_s, const unsigned options)
{
  xmlNodePtr p_nodeCur = xpathCtx->node; /* save current node to return */
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  size_t inode;
  int rc=0;

  if (!xddmObj) return 1;

  *p_ns = countNodes(xddmObj); /* number of analysis parameters */
  if ( 0 == *p_ns ) {
    ERR("xddm_readSum countNodes failed\n");
    return 1;
  }
                                    /* parse node attributes of all elements */
  pp_xmlAttr = xddmParseXpathObj(xddmObj);

  (*a_s) = xddm_newSum(*p_ns);

  inode = 0;                        /* loop over each node, sentinel is NULL */
  while (NULL != pp_xmlAttr[inode]) {
    if (inode >= (*p_ns)) {
      ERR("overflow in number of sum nodes\n");
      freeXddmAttr(pp_xmlAttr);
      return 1;
    }
    else {
      size_t i;
      const p_tsXddmSum     p_s  = (*a_s) + inode;
      const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
      xmlXPathObjectPtr     xObj = NULL;

      p_s->p_e->p_xn = xddm_getName(p_xA->p_node, options);
      p_s->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
      for (i=0; i<p_xA->n; i++) {
        if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
          p_s->val = fillDouble(p_xA->a_values[i]);
        }
        else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "expr" ) ) {
          p_s->p_expr = fillXmlString(p_xA->a_values[i]);
        }
        // else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "min" ) ) {
        //   p_f->min = fillDouble(p_xA->a_values[i]);
        // }
        // else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "max" ) ) {
        //   p_f->max = fillDouble(p_xA->a_values[i]);
        // }
        else {
          rc = xddm_fillElement(p_s->p_e, p_xA->a_names[i], p_xA->a_values[i]);
          if (rc) {
            ERR("xddm_fillElement failed\n");
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        }
      } /* end for attributes */
                                                   /* look for sensitivities */
                                /* move xpath node into the analysis element */
      xpathCtx->node = pp_xmlAttr[inode]->p_node;
                                /* check if the node has a sensitivity array */
      xObj = xmlXPathEvalExpression(BAD_CAST "./SensitivityArray/Sensitivity",
                                        xpathCtx);
      if (xObj) {
        if ( ! xmlXPathNodeSetIsEmpty(xObj->nodesetval) ) {
          rc = xddm_readSensitivity(xObj, &p_s->ndvs, &p_s->a_lin, &p_s->pa_dvs);
          if ( rc ) {
            ERR("xddm_readSensitivity failed for %s\n",p_s->p_e->p_id);
            xmlXPathFreeObject(xObj);
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        }
        xmlXPathFreeObject(xObj);
        xObj = NULL;
      }
      xpathCtx->node = p_nodeCur; /* reset root node */
    } /* end if inode */

    inode++;
  } /* end while attributes */

  freeXddmAttr(pp_xmlAttr);

  xpathCtx->node = p_nodeCur; /* reset root node */

  return 0;
}

static int
xddm_readIntersect(xmlXPathObjectPtr xddmObj, p_tsXddmGeom *pp_geom,
                   size_t *p_ng, const unsigned options)
{
  size_t inode=0;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;
  p_tsXddmGeom p_geom = NULL;

  *p_ng = countNodes(xddmObj);

  if ( *p_ng == 0 ) {
    WARN("Syntax error likely in Intersect: %zu\n", *p_ng);
    return 1;
  }

  p_geom = xddm_newIntersect( *p_ng );
  pp_xmlAttr = xddmParseXpathObj(xddmObj);

  inode=0;
  while (NULL != pp_xmlAttr[inode]) {
    size_t i;
    int rc;
    const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
    p_tsXddmGeom           p_g = p_geom + inode;

    p_g->p_e->p_xn = xddm_getName(p_xA->p_node, options);
    p_g->p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
    for (i=0; i<p_xA->n; i++) {
      if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "parts" ) ) {
        p_g->p_parts = fillXmlString(p_xA->a_values[i]);
      }
      else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comp2tri" ) ) {
        p_g->p_comp2tri = fillXmlString(p_xA->a_values[i]);
      }
      else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "cutout" ) ) {
        p_g->p_cutout = fillXmlString(p_xA->a_values[i]);
      }
      else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "overlap" ) ) {
        p_g->p_overlap = fillXmlString(p_xA->a_values[i]);
      }
      else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "ps" ) ) {
        p_g->p_ps = fillXmlString(p_xA->a_values[i]);
      }
      else {
        rc = xddm_fillElement(p_g->p_e, p_xA->a_names[i], p_xA->a_values[i]);
        if (rc) {
          ERR("Failed read element attributes\n");
          freeXddmAttr(pp_xmlAttr);
          return 1;
        }
      }
    } /* end for attributes */
    inode++;
  }

  freeXddmAttr(pp_xmlAttr);

  *pp_geom = p_geom; /* hook into xddm data struct */

  return 0;
}

static int
xddm_addXmlAttribute(const xmlChar * p_name, const xmlChar * p_value,
                     size_t *p_nAttr, p_tsXddmAttr *pp_attr)
{
  if ( 0 == *p_nAttr ) {
    *p_nAttr = 1;
    *pp_attr = xddm_newAttribute(*p_nAttr);
  }
  else {
    p_tsXddmAttr p_tmp = NULL;
    (*p_nAttr) += 1;
    p_tmp = realloc(*pp_attr, (*p_nAttr)*sizeof(*(*pp_attr)));
    if (NULL == p_tmp) {
      ERR("realloc failed\n");
      exit(1);
    }
    else {
      *pp_attr = p_tmp;
    }
  }
  (*pp_attr)[(*p_nAttr)-1].p_name  = fillXmlString(p_name);
  (*pp_attr)[(*p_nAttr)-1].p_value = fillXmlString(p_value);

  // printf("%zu %s %s\n",*p_nAttr,(*pp_attr)[(*p_nAttr)-1].p_name,(*pp_attr)[(*p_nAttr)-1].p_value);

  return 0;
}

static int
xddm_fillElement(p_tsXddmElem p_e, const xmlChar *p_name,
                 const xmlChar *p_value)
{
  int rc = 0;

  if ( 0 == xmlStrcasecmp( p_name, BAD_CAST "id" ) ) {
    p_e->p_id = fillXmlString(p_value);
  }
  else if ( 0 == xmlStrcasecmp( p_name, BAD_CAST "comment" ) ) {
    p_e->p_comment = fillXmlString(p_value);
  }
  else if ( 0 == xmlStrcasecmp( p_name, BAD_CAST "sensitivity" ) ) {
    if ( 0 == xmlStrcasecmp( p_value, BAD_CAST "required" ) ) {
      p_e->flags |= XDDM_LIN;
    }
    else if ( 0 == xmlStrcasecmp( p_value, BAD_CAST "none" ) ) {
      p_e->flags |= XDDM_NOLIN;
    }
  }
  else if ( 0 == xmlStrcasecmp( p_name, BAD_CAST "bound" ) ) {
    if ( 0 == xmlStrcasecmp( p_value, BAD_CAST "upper" ) ) {
      p_e->flags |= XDDM_BOUND_UPPER;
    }
    else if ( 0 == xmlStrcasecmp( p_value, BAD_CAST "lower" ) ) {
      p_e->flags |= XDDM_BOUND_LOWER;
    }
  }
  else {                               /* extra attributes outside xddm spec */
    rc = xddm_addXmlAttribute(p_name, p_value, &p_e->nAttr, &p_e->p_attr);
    if (rc) {
      ERR("Failed read of extra attributes\n");
    }
  }

  return rc;
}

static int
xddm_readElement(xmlXPathObjectPtr xddmObj, size_t *p_ne, p_tsXddmElem *a_e,
                 const unsigned options)
{
  size_t inode;
  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  assert(xddmObj);

  *p_ne = countNodes(xddmObj);

  if ( *p_ne > 0 ) {
                                    /* parse node attributes of all elements */
    pp_xmlAttr = xddmParseXpathObj(xddmObj);
    if ( NULL == pp_xmlAttr ) {
      ERR("Parsing problem in tessellate\n");
      return 1;
    }

    (*a_e) = xddm_newElement(*p_ne);

    inode = 0;                      /* loop over each node, sentinel is NULL */
    while (NULL != pp_xmlAttr[inode]) {
      if (inode >= (*p_ne)) {
        ERR("Overflow in number of nodes\n");
        freeXddmAttr(pp_xmlAttr);
        return 1;
      }
      else {
        size_t i;
        int rc;
        p_tsXddmElem           p_e = (*a_e) + inode;
        const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];

        p_e->p_xn = xddm_getName(p_xA->p_node, options);
        p_e->p_nn = fillXmlString(p_xA->p_node->name);

           /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<p_xA->n; i++) {
          rc = xddm_fillElement(p_e, p_xA->a_names[i], p_xA->a_values[i]);
          if (rc) {
            ERR("Failed read element attributes\n");
            freeXddmAttr(pp_xmlAttr);
            return 1;
          }
        } /* end for attributes */
      } /* end if inode */
      inode++;
    } /* end while attributes */

    freeXddmAttr(pp_xmlAttr);
  } /* end if tessellate nodes */

  return 0;
}

/**
 * countNodes(): count XML nodes in xpath object
 */
size_t
countNodes(const xmlXPathObjectPtr xpathObj)
{
  size_t n=0;
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
    n++;
  }

  return n;
}

/**
 * xddmParseXpathObj(): filter xpath nodes to store only element nodes
 * in local data struct to simplify parsing
 */
p_tsXddmXmlAttr *
xddmParseXpathObj(const xmlXPathObjectPtr xpathObj)
{
  int i,j;

  const xmlNodeSetPtr nodes = xpathObj->nodesetval;
  int nNodes= (nodes) ? nodes->nodeNr : 0;

  p_tsXddmXmlAttr *pp_xmlAttr = NULL;

  j = 0;
  for (i=0; i < nNodes; i++) {
    if ( nodes->nodeTab[i]->type != XML_ELEMENT_NODE ) continue;
    j++;
  }

  nNodes = j;
  j = 0;

  if ( nNodes > 0 ) {
    pp_xmlAttr = malloc( (nNodes+1) * sizeof(*pp_xmlAttr));
    if (NULL == pp_xmlAttr) {
      ERR("malloc failed\n");
      exit(1);
    }

    for (i=0; i < nNodes; i++) {
      if ( nodes->nodeTab[i]->type != XML_ELEMENT_NODE ) continue;

      pp_xmlAttr[j] = xddmParseNode(nodes->nodeTab[i]);
      if (NULL == pp_xmlAttr[j]) {
        ERR("node has no attributes\n");
        exit(1);
      }
      j++;
    }
    pp_xmlAttr[nNodes] = NULL; /* sentinel */
  }
  return pp_xmlAttr;
}

/**
 * xddmParseNode(): save all element attributes
 */
p_tsXddmXmlAttr
xddmParseNode(const xmlNodePtr p_node)
{
  size_t          i;
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

  i = 0;
  p_attribute = p_node->properties;
  while (p_attribute && p_attribute->name && p_attribute->children) {
    p_xmlAttr->a_names[i]  = p_attribute->name;
    p_xmlAttr->a_values[i] = xmlGetProp(p_node, p_attribute->name);
    if (NULL == p_xmlAttr->a_values[i]) {
      ERR("xddmParseNode xmlGetProp failed for %s\n",(char *)p_attribute->name);
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
