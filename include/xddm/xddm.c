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

/**
 * xddm_readFile(): parse an XPath expression from an XDDM file and return
 * structure containing all elements in the path
 */
p_tsXddm
xddm_readFile(const char *const p_fileName, const char *const xpathExpr, const int options)
{
  xmlDocPtr          doc;
  xmlParserCtxtPtr   ctxt;
  xmlXPathObjectPtr  xddmObj = NULL, linObj = NULL;
  xmlXPathContextPtr xpathCtx = NULL;
  
  p_tsXddmXmlAttr *pp_xmlAttr = NULL, *pp_sensit  = NULL;

  p_tsXddm p_xddm = NULL; /* return object */
  
  p_tsXmParent p_p;
  
  int parserOpts=0;
  int i, j, inode=0;
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
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    return NULL;
  }
 
  xmlFreeParserCtxt(ctxt); /* done with document parser context */

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("xmlXPathNewContext failed to create context\n");
    shutdown_libxml(doc);
    return NULL;
  }

  /* looks like we have a valid xml document that we are ready to parse
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if (NULL == xddmObj) {
    WARN("no elements found for expression %s\n", xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return NULL;
  }

  if ( 0 == countNodes(xddmObj) ) {
    WARN("no elements found for expression \'%s\'\n", xpathExpr);
    xmlXPathFreeObject(xddmObj);
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return NULL;
  }

  p_xddm = xddm_alloc(); /* return object */
  
  p_xddm->fileName  = fillString(p_fileName); 
  p_xddm->xpathExpr = fillString(xpathExpr);
  
  pp_xmlAttr = xddmParseXpathObj(xddmObj);
  if ( NULL == pp_xmlAttr ) {
    ERR("problem in parsing variables\n");
    exit(1);
  }

  p_p = p_xddm->p_parent;
  p_p->name   = fillXmlString((*pp_xmlAttr)->p_node->name);
  p_p->nAttr  = (*pp_xmlAttr)->n;
  p_p->p_attr = malloc(p_p->nAttr * sizeof(*p_p->p_attr));
  if (NULL == p_p->p_attr) {
    ERR("malloc failed\n");
    exit(1);
  }
  
  for (i=0; i<p_p->nAttr; i++) {
    p_p->p_attr[i].p_name  = fillXmlString((*pp_xmlAttr)->a_names[i]);
    p_p->p_attr[i].p_value = fillXmlString((*pp_xmlAttr)->a_values[i]);
  }

  freeXddmAttr(pp_xmlAttr);
  
  xpathCtx->node = xddmObj->nodesetval->nodeTab[0];
  xmlXPathFreeObject(xddmObj);
  
  /* ---------------
   * parse Variables
   * ---------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Variable", xpathCtx);
  if (xddmObj) {
    p_xddm->nv = countNodes(xddmObj); /* number of variables */

    if ( p_xddm->nv > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing variables\n");
        exit(1);
      }
    
      p_xddm->a_v = allocVariable(p_xddm->nv);

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        const p_tsXddmVar p_v = &p_xddm->a_v[inode];
        
        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<pp_xmlAttr[inode]->n; i++) {
          const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
        
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
    } /* end if variables */
  
    xmlXPathFreeObject(xddmObj);
  } /* end if xddmObj */

  /* ---------------
   * parse Constants
   * ---------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Constant", xpathCtx);
  if (xddmObj) {
    p_xddm->nc = countNodes(xddmObj); /* number of constants */

    if ( p_xddm->nc > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing constants\n");
        exit(1);
      }
    
      p_xddm->a_c = allocVariable(p_xddm->nc);

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        const p_tsXddmVar p_c = &p_xddm->a_c[inode];

        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<pp_xmlAttr[inode]->n; i++) {
          const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
        
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_c->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_c->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
            p_c->val = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "min" ) ) {
            p_c->minVal = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "max" ) ) {
            p_c->maxVal = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "typicalsize" ) ) {
            p_c->typicalSize = fillDouble(p_xA->a_values[i]);
          }
        } /* end for attributes */
      
        inode++;
      } /* end while attributes */

      freeXddmAttr(pp_xmlAttr);
    } /* end if constants */
  
    xmlXPathFreeObject(xddmObj);
  } /* end if xddmObj */


  /* ----------------
   * parse Tessellate
   * ----------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Tessellate", xpathCtx);
  if (xddmObj) {
    p_xddm->ntess = countNodes(xddmObj); /* number of variables */
    if ( p_xddm->ntess > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in parsing tessellate\n");
        exit(1);
      }
    
      p_xddm->a_tess = allocTessellate(p_xddm->ntess);

      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        p_tsXmTess p_t = &p_xddm->a_tess[inode];
        
        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<pp_xmlAttr[inode]->n; i++) {
          const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
        
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_t->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_t->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "sensitivity" ) ) {
            if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "required" ) ) {
              p_t->lin = 1;
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
    } /* end if variables */
  
    xmlXPathFreeObject(xddmObj);
  } /* end if xddmObj */

  /* -------------------------
   * parse Analysis parameters
   * -------------------------
   */

  xddmObj = xmlXPathEvalExpression(BAD_CAST "./Analysis", xpathCtx);
  if (xddmObj) {
    p_xddm->na = countNodes(xddmObj); /* number of analysis parameters */

    if ( p_xddm->na > 0 ) {
                                    /* parse node attributes of all elements */
      pp_xmlAttr = xddmParseXpathObj(xddmObj);
      if ( NULL == pp_xmlAttr ) {
        ERR("problem in analysis parameters parsing\n");
        exit(1);
      }
    
      p_xddm->a_ap = allocAnalysis(p_xddm->na);
      
      inode = 0;                    /* loop over each node, sentinel is NULL */
      while (NULL != pp_xmlAttr[inode]) {
        const p_tsXddmAPar p_a = &p_xddm->a_ap[inode];
        
        /* for each node, loop over its attributes and save relevant ones */
        for (i=0; i<pp_xmlAttr[inode]->n; i++) {
          const p_tsXddmXmlAttr p_xA = pp_xmlAttr[inode];
        
          if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "id" ) ) {
            p_a->p_id = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "comment" ) ) {
            p_a->p_comment = fillXmlString(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "value" ) ) {
            p_a->val = fillDouble(p_xA->a_values[i]);
          }
          else if ( 0 == xmlStrcasecmp( p_xA->a_names[i], BAD_CAST "sensitivity" ) ) {
            if ( 0 == xmlStrcasecmp( p_xA->a_values[i], BAD_CAST "required" ) ) {
              p_a->lin = 1;
            }
          }
        } /* end for attributes */

                                /* check if the node has a sensitivity array */
        xpathCtx->node = pp_xmlAttr[inode]->p_node;
        linObj = xmlXPathEvalExpression(BAD_CAST "./SensitivityArray/Sensitivity", xpathCtx);
        if (linObj) {
          if ( countNodes(linObj) > 0 ) {
            pp_sensit = xddmParseXpathObj(linObj);
      
            p_a->ndvs=0;
            while (NULL != pp_sensit[p_a->ndvs]) {
              p_a->ndvs++;
            }
            
            p_a->pa_dvs = malloc( p_a->ndvs * sizeof(*p_a->pa_dvs));
            if (NULL == p_a->pa_dvs) {
              ERR("malloc failed\n");
              exit(1);
            }
      
            p_a->a_lin = malloc( p_a->ndvs * sizeof(*p_a->a_lin));
            if (NULL == p_a->a_lin) {
              ERR("malloc failed\n");
              exit(1);
            }
      
            i=0;
            while (NULL != pp_sensit[i]) {
              const p_tsXddmXmlAttr p_s = pp_sensit[i];
              
              for (j=0; j<p_s->n; j++) {
                if ( 0 == xmlStrcasecmp( p_s->a_names[j], BAD_CAST "p" ) ) {
                  p_a->pa_dvs[i] = fillXmlString(p_s->a_values[j]);
                }
                else if ( 0 == xmlStrcasecmp( p_s->a_names[j], BAD_CAST "value" ) ) {
                  p_a->a_lin[i] = fillDouble(p_s->a_values[j]);
                }
              }
              i++;
            } /* end while sensitivities */
        
            freeXddmAttr(pp_sensit);
          }
      
          xmlXPathFreeObject(linObj);
        }
        inode++;
      } /* end while attributes */
      
      freeXddmAttr(pp_xmlAttr);
    } /* end if analysis parameters */

    xmlXPathFreeObject(xddmObj);
  }
  
  /* DO NOT PROCESS NODES HERE BECAUSE xpathCtx->node HAS BEEN ASSIGNED TO
     Sensitivity ARRAY AND HAS NOT BEEN RESET */

  xmlXPathFreeContext(xpathCtx);

  /* all done, clean up and return
   */

  shutdown_libxml(doc);

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

p_tsXddmVar
allocVariable(const int n)
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

p_tsXddmAPar
allocAnalysis(const int n)
{
  int i;
  p_tsXddmAPar p_v=NULL;
  
  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_id      = NULL;
    p_v[i].p_comment = NULL;
    /* p_v[i].p_parent  = NULL; */
    p_v[i].val       = UNSET;
    p_v[i].lin       = 0;
    p_v[i].ndvs      = 0;
    p_v[i].a_lin     = NULL;
    p_v[i].pa_dvs    = NULL;
  }
  
  return p_v;
} 

p_tsXmTess
allocTessellate(const int n)
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
    p_v[i].lin       = 0;
    p_v[i].nAttr     = 0;
    p_v[i].p_attr    = NULL;
  }
  
  return p_v;
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

void
shutdown_libxml(xmlDocPtr doc)
{
  xmlFreeDoc(doc);
  xmlCleanupParser(); /* shutdown libxml */
  xmlMemoryDump(); /* this is to debug memory for regression tests */ 
  return;
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
    printf("   Attributes:\n");
    for (i=0; i<p_xddm->p_parent->nAttr; i++) {
      p_tsXmParent p_p = p_xddm->p_parent;
      printf("      %s %s\n", p_p->p_attr[i].p_name, p_p->p_attr[i].p_value);
    }
    printf("\n");
  }
  
  printf(" o Number of variables = %d\n",p_xddm->nv);
  for (i=0; i<p_xddm->nv; i++) {
    printf("\n   %4d",i);
    if (NULL != p_xddm->a_v[i].p_id) printf(" ID=%s",p_xddm->a_v[i].p_id);
    if (UNSET != p_xddm->a_v[i].val) printf("  Value=%g",p_xddm->a_v[i].val);
    printf("\n");
    if (UNSET != p_xddm->a_v[i].minVal) printf("        Min=%g",p_xddm->a_v[i].minVal);
    if (UNSET != p_xddm->a_v[i].maxVal) printf(" Max=%g",p_xddm->a_v[i].maxVal);
    if (UNSET != p_xddm->a_v[i].typicalSize) printf(" TypicalSize=%g",p_xddm->a_v[i].typicalSize);
    if (UNSET != p_xddm->a_v[i].minVal || UNSET != p_xddm->a_v[i].maxVal ||
        UNSET != p_xddm->a_v[i].typicalSize ) printf("\n");
    if (NULL != p_xddm->a_v[i].p_comment) printf("        %s\n",p_xddm->a_v[i].p_comment);    
  }

  printf("\n");
  printf(" o Number of constants = %d\n",p_xddm->nc);
  for (i=0; i<p_xddm->nc; i++) {
    printf("\n   %4d",i);
    if (NULL != p_xddm->a_c[i].p_id) printf(" ID=%s",p_xddm->a_c[i].p_id);
    if (UNSET != p_xddm->a_c[i].val) printf("  Value=%g",p_xddm->a_c[i].val);
    printf("\n");
    if (UNSET != p_xddm->a_c[i].minVal) printf("        Min=%g",p_xddm->a_c[i].minVal);
    if (UNSET != p_xddm->a_c[i].maxVal) printf(" Max=%g",p_xddm->a_c[i].maxVal);
    if (UNSET != p_xddm->a_c[i].typicalSize) printf(" TypicalSize=%g",p_xddm->a_c[i].typicalSize);
    if (UNSET != p_xddm->a_c[i].minVal || UNSET != p_xddm->a_c[i].maxVal ||
        UNSET != p_xddm->a_c[i].typicalSize ) printf("\n");
    if (NULL != p_xddm->a_c[i].p_comment) printf("        %s\n",p_xddm->a_c[i].p_comment);    
  }

  printf("\n");
  printf(" o Number of analysis parameters = %d\n",p_xddm->na);
  for (i=0; i<p_xddm->na; i++) {
    printf("\n   %4d",i);
    if (NULL != p_xddm->a_ap[i].p_id) printf(" ID=%s",p_xddm->a_ap[i].p_id);
    if (UNSET != p_xddm->a_ap[i].val) printf("  Value=%g",p_xddm->a_ap[i].val);
    printf("\n");
    if (p_xddm->a_ap[i].lin) printf("        Linearization Required\n");
    if (NULL != p_xddm->a_ap[i].p_comment) printf("        %s\n",p_xddm->a_ap[i].p_comment);
    if ( p_xddm->a_ap[i].ndvs > 0 ) {
      printf("        Sensitivity array\n");
      for (j=0; j<p_xddm->a_ap[i].ndvs; j++) {
        printf("           DV=%s Value=%g\n",p_xddm->a_ap[i].pa_dvs[j], p_xddm->a_ap[i].a_lin[j]);
      }
    }
  }

  printf("\n");
  printf(" o Number of tessellate elements = %d\n",p_xddm->ntess);
  for (i=0; i<p_xddm->ntess; i++) {
    printf("\n   %4d",i);
    if (NULL != p_xddm->a_tess[i].p_id) printf(" ID=%s",p_xddm->a_tess[i].p_id);
    printf("\n");
    if (p_xddm->a_tess[i].lin) printf("        Linearization Required\n");
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

#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

#define MY_ENCODING "ISO-8859-1"

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
    return 1;
  }

  /* create a new XmlWriter for memory, with no compression */
  writer = xmlNewTextWriterMemory(buf, 0);
  if (writer == NULL) {
    ERR("xmlNewTextWriterMemory failed to create xml writer\n");
    return 1;
  }

  /* start the document with the xml default and encoding ISO 8859-1 */
  rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
  if (rc < 0) {
    ERR("xmlTextWriterStartDocument failed\n");
    return 1;
  }

  rc = xmlTextWriterSetIndent(writer, 1); /* pretty indent */
  if (rc < 0) {
    ERR("xmlTextWriterSetIndent failed\n");
    return 1;
  }

  /* write parent element */
  xddm_writeAttr(writer, p_xddm->p_parent);

  for (i=0; i<p_xddm->nv; i++) {
  
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Variable");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    xddm_writeVariable(writer, &p_xddm->a_v[i]);
  
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

    xddm_writeVariable(writer, &p_xddm->a_c[i]);
  
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

    xddm_writeAnalysis(writer, &p_xddm->a_ap[i]);
  
    rc = xmlTextWriterEndElement(writer);         /* end analysis */
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
    return 1;
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
      return 1;
    }
  }

  /* do not end element here in case we need to nest another one */
  
  return 0;
}

int
xddm_writeDouble(xmlTextWriterPtr writer, const char *p_name, const double val)
{
  int rc;
  
  if (UNSET != val) {
    if ( 0 != xmlXPathIsInf(val) ) {
      ERR("Value is INFINITE\n");
      exit(1);
    }
    if ( 1 == xmlXPathIsNaN(val) ) {
      ERR("Value is NaN\n");
      exit(1);
    }

    if ( 0 == strcasecmp(p_name, "value") ) {
      rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST p_name,"%.17e", val);
      if (rc < 0) {
        ERR("xmlTextWriterWriteFormatAttribute failed\n");
        return 1;
      }
    }
    else {
      rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST p_name,"%g", val);
      if (rc < 0) {
        ERR("xmlTextWriterWriteFormatAttribute failed\n");
        return 1;
      }
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
      return 1;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_v->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  rc = xddm_writeDouble(writer, "Min", p_v->minVal);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  rc = xddm_writeDouble(writer, "Max", p_v->maxVal);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  rc = xddm_writeDouble(writer, "TypicalSize", p_v->typicalSize);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }
  
  if (p_v->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_v->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return 1;
    }
  }
  
  return 0;
}

int
xddm_writeAnalysis(xmlTextWriterPtr writer, const p_tsXddmAPar p_a)
{
  int i, rc;

  if (p_a->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_a->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return 1;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_a->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  if (p_a->lin) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity", BAD_CAST "Required");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return 1;
    }
  }
  else {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity", BAD_CAST "None");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return 1;
    }
  }

  if (p_a->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment", BAD_CAST p_a->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return 1;
    }
  }
  
  if ( p_a->ndvs > 0 ) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "SensitivityArray");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return 1;
    }
    
    for (i=0; i<p_a->ndvs; i++) {
      rc = xmlTextWriterStartElement(writer, BAD_CAST "Sensitivity");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return 1;
      }
      
      rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "P", BAD_CAST p_a->pa_dvs[i]);
      if (rc < 0) {
        ERR("xmlTextWriterWriteAttribute failed\n");
        return 1;
      }

      rc = xddm_writeDouble(writer, "Value", p_a->a_lin[i]);
      if (rc) {
        ERR("xddm_writeDouble failed\n");
        return 1;
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
    ERR(" xddm_readFile failed to allocate parser context\n");
    xmlCleanupParser();
    return 1;
  }

  doc = xmlCtxtReadFile(ctxt, p_xddm->fileName, NULL, parserOpts);
  if (doc == NULL) {
    if ( XDDM_VERBOSE & options ) {
      ERR("%s is not valid XML\n", p_fileName);
    }
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    return 1;
  }
 
  xmlFreeParserCtxt(ctxt); /* done with document parser context */

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("xmlXPathNewContext failed to create context\n");
    shutdown_libxml(doc);
    return 1;
  }

  if (NULL == p_xddm->xpathExpr) {
    ERR("missing xpath expression\n");
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return 1;
  }
  
  xddmObj = xmlXPathEvalExpression(BAD_CAST p_xddm->xpathExpr, xpathCtx);
  if (NULL == xddmObj) {
    WARN("no elements found for expression %s\n", p_xddm->xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return 1;
  }

  if ( 0 == countNodes(xddmObj) ) {
    WARN("no elements found for expression \'%s\'\n", p_xddm->xpathExpr);
    xmlXPathFreeObject(xddmObj);
    xmlXPathFreeContext(xpathCtx);
    shutdown_libxml(doc);
    return 1;
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
      exit(1);
    }
    memset((void*) xpathExpr,'\0',len);
           
    strncpy(xpathExpr,"./Analysis[@ID=\"",strlen("./Analysis[@ID=\""));
    strncat(xpathExpr,p_a->p_id,strlen(p_a->p_id));
    strncat(xpathExpr,"\"]",strlen("\"]"));
    
    xddmObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx); 

    if (NULL == xddmObj) { 
      WARN("no elements found for expression %s\n", xpathExpr); 
      xmlXPathFreeContext(xpathCtx); 
      shutdown_libxml(doc); 
      return 1; 
    } 

    if ( 1 != countNodes(xddmObj) ) { 
      ERR("analysis parameter not unique\n"); 
      xmlXPathFreeContext(xpathCtx); 
      xmlXPathFreeObject(xddmObj); 
      shutdown_libxml(doc); 
      return 1; 
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
        new_attr = xmlNewProp(new_node, BAD_CAST "P", BAD_CAST p_a->pa_dvs[i]);
        xmlstr   = xmlXPathCastNumberToString(p_a->a_lin[i]);
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
  
  p_xddm->p_parent = malloc(sizeof(*p_xddm->p_parent));
  if (NULL == p_xddm->p_parent) {
    ERR("malloc failed\n");
    exit(1);
  }
  p_xddm->p_parent->name   = NULL;
  p_xddm->p_parent->nAttr  = 0;
  p_xddm->p_parent->p_attr = NULL;
  
  p_xddm->nv       = 0;
  p_xddm->nc       = 0;
  p_xddm->na       = 0;
  p_xddm->ntess    = 0;
  p_xddm->a_v      = NULL;
  p_xddm->a_c      = NULL;
  p_xddm->a_ap     = NULL;
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
  int i,j;
  
  if (p_xddm) {
    if (p_xddm->fileName) free(p_xddm->fileName);
    if (p_xddm->xpathExpr) free(p_xddm->xpathExpr);
    
    if (p_xddm->p_parent) {
      if (p_xddm->p_parent->name) free(p_xddm->p_parent->name);
      for (i=0; i<p_xddm->p_parent->nAttr; i++) {
        if (p_xddm->p_parent->p_attr[i].p_name) free(p_xddm->p_parent->p_attr[i].p_name);
        if (p_xddm->p_parent->p_attr[i].p_value) free(p_xddm->p_parent->p_attr[i].p_value);
      }
      if (p_xddm->p_parent->p_attr) free(p_xddm->p_parent->p_attr);
    }
    
    for (i=0; i<p_xddm->nv; i++) {
      if (p_xddm->a_v[i].p_id) free(p_xddm->a_v[i].p_id);
      if (p_xddm->a_v[i].p_comment) free(p_xddm->a_v[i].p_comment);
    }
    for (i=0; i<p_xddm->nc; i++) {
      if (p_xddm->a_c[i].p_id) free(p_xddm->a_c[i].p_id);
      if (p_xddm->a_c[i].p_comment) free(p_xddm->a_c[i].p_comment);
    }

    for (i=0; i<p_xddm->na; i++) {
      p_tsXddmAPar p_a = &p_xddm->a_ap[i];
      
      if (p_a->p_id)      free(p_a->p_id);
      if (p_a->p_comment) free(p_a->p_comment);
      if (p_a->a_lin)     free(p_a->a_lin);
      for (j=0; j<p_a->ndvs; j++) {
        if (p_a->pa_dvs[j]) free(p_a->pa_dvs[j]);
      }
      if (p_a->pa_dvs) free(p_a->pa_dvs);
    }

    for (i=0; i<p_xddm->ntess; i++) {
      p_tsXmTess p_t = &p_xddm->a_tess[i];
      
      if (p_t->p_id)      free(p_t->p_id);
      if (p_t->p_comment) free(p_t->p_comment);
      
      for (j=0; j<p_t->nAttr; j++) {
        if (p_t->p_attr[j].p_name) free(p_t->p_attr[j].p_name);
        if (p_t->p_attr[j].p_value) free(p_t->p_attr[j].p_value);
      }
      if (p_t->p_attr) free(p_t->p_attr);
    }
    
    if (p_xddm->p_parent) free(p_xddm->p_parent);
    if (p_xddm->a_v)      free(p_xddm->a_v);
    if (p_xddm->a_c)      free(p_xddm->a_c);
    if (p_xddm->a_ap)     free(p_xddm->a_ap);
    if (p_xddm->a_tess)   free(p_xddm->a_tess);
 
    free(p_xddm);
    p_xddm = NULL;
  }
  
  return;
}
