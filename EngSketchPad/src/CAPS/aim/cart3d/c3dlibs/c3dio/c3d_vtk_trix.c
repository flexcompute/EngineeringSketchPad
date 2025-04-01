
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
 * $Id: c3d_vtk_trix.c,v 1.2 2022/11/07 18:45:42 mnemec Exp $
 */

/* open source */

/*
 * Cart3D IO lib for TRIX files
 */

/* ----- include files ------ */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h> /* for fsync */

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include "c3dio_lib.h"

/* ----- local macros ------- */

#define MY_ENCODING "ISO-8859-1"

/*                       XML_PARSE_HUGE is only defined libxml2 v2.7 and newer
 *                                   set to (524288 = 1<<19) and leave defined
 */
#if !defined( XML_PARSE_HUGE )
#define XML_PARSE_HUGE 524288
#endif

/* -------------private prototypes -------------- */
int trix_readCompInfo(xmlDocPtr doc, p_tsTriangulation *pp_config,
                      int *const p_nComps, const char *p_compName,
                      const int firstComp, const int options);
int trix_setTriInfo(p_tsTriangulation p_config, const int nComps,
                    const char *p_compName, const xmlNodePtr cur_node,
                    const int firstComp);
int trix_readTriDataInfo(xmlDocPtr doc, p_tsTriangulation p_config,
                         const int nComps, const char *p_dataName,
                         const int firstComp);
int trix_setVertInfo(p_tsTriangulation p_config, int nComps,
                     const char *p_compName, const xmlNodePtr cur_node,
                     const int firstComp);
int trix_readVertDataInfo(xmlDocPtr doc, p_tsTriangulation p_config,
                          const int nComps, const char *p_dataName,
                          const int firstComp);
int trix_readVertCoord(p_tsTriangulation p_config, const int icomp,
                       const xmlNodePtr cur_node);
int trix_readConnectivity(p_tsTriangulation p_config, const int icomp,
                          const xmlNodePtr cur_node);
int trix_readTriData(xmlXPathContextPtr xpathCtx, p_tsTriangulation p_comp);
int trix_readVertData(xmlXPathContextPtr xpathCtx,
                      const p_tsTriangulation p_comp);
int trix_fillContent(xmlDocPtr doc, const p_tsTriangulation p_config,
                     const int nComps, const int firstComp);
VTKtype   getVTKtype(const char *type, const int len);
TRIXtype getTRIXtype(const char *type, const int len);

/* ------------- private datatypes ------------------------------------ */
char *vtkTypes[11] =
  {
    "Int8",
    "Uint8",
    "Int16",
    "UInt16",
    "Int32",
    "UInt32",
    "Int64",
    "UInt64",
    "Float32",
    "Float64",
    "UNSET"
  };

char *trixTypes[5] =
  {
    "UNSET",
    "FLOW_VARIABLE",
    "SHAPE_LINEARIZATION",
    "COMPONENT_TAG",
    "OTHER",
  };


/* ==================================================================== */


/**
 * -- io_readSurfTrix(): read VTK unstructured grid ASCII file.  For format
 *         definition see www.vtk.org/VTK/img/file-formats.pdf.  Libxml2 is
 *         used to parse the data (xmlsoft.org). Note that tri-connectivity
 *         uses  zero-offset.
 *
 *       Function returns the number of components, nComps, and a
 *         triangulation for each component, p_config. nComps must be
 *         initialized to zero by the caller to read in the first component
 *         set. If p_compName is specified then only that component is
 *         retrieved. For each component, you may also request specific
 *         vert and tri data. Initializing p_compName, p_vertDataNames,
 *         p_triDataNames to ALL retrieves everything from the
 *         file. Setting 'options'=1 denotes verbose mode and this variable
 *         can be extended to include other options by encoding powers of
 *         two and using bitwise comparisons.
 *
 *       Caller must free p_config when done using c3d_freeTriangulation().
 *
 *       NOTE: If we are compiling against an old libxml2 library that does
 *         not have the XML_PARSE_HUGE parser option (versions before 2.7.0)
 *         then we set it manually. This allows us to bind to a newer
 *         dynamic library (if required) so that we can still read large
 *         (over 10MB) tri files
 */
int io_readSurfTrix(const char *const p_fileName,
                    p_tsTriangulation *pp_config,
                    int *const p_nComps,
                    const char *const p_compName,
                    const char *const p_vertDataNames,
                    const char *const p_triDataNames,
                    const int options)
{
  xmlDocPtr doc;
  xmlNode  *root = NULL;
  xmlParserCtxtPtr ctxt; /* the parser context */

  p_tsTriangulation p_comp = NULL;
  const int firstComp = *p_nComps;

  int rc, i, parserOpts;
  /* initializes and check libxml2 version */
  { LIBXML_TEST_VERSION };

  xmlXPathInit();  /* initialize xpath */

  if ( TRIX_VERBOSE & options ) {
    printf("    o  Parsing file \"%s\" with libxml2\n", p_fileName);
  }


  ctxt = xmlNewParserCtxt();  /* create a parser context */
  if (ctxt == NULL) {
    ERR(" io_readSurfTrix failed to allocate parser context\n");
    xmlCleanupParser();
    return -1;
  }

  parserOpts  = 0; /* set parser to be silent */
  parserOpts |= (XML_PARSE_NOERROR);
  parserOpts |= (XML_PARSE_NOWARNING);
  parserOpts |= (XML_PARSE_HUGE); /* relax any hardcoded limit from the parser */

                                                        /* load XML document */
  doc = xmlCtxtReadFile(ctxt, p_fileName, NULL, parserOpts);
  if (doc == NULL) {
    if ( TRIX_VERBOSE & options ) {
      CONT("%s does not appear to be valid XML (io_readSurfTrix)\n", p_fileName);
    }
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    return -1;
  }

  xmlFreeParserCtxt(ctxt); /* done with parser */

  root = xmlDocGetRootElement(doc); /* check root element node */
  if ( !root ||
       !root->name ||
       xmlStrcmp(root->name, (xmlChar *) "VTKFile") ) {
    WARN(" io_readSurfTrix root element not VTKFile\n");
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return -1;
  }
                               /* read number of Verts and Tris for compName */
  rc = trix_readCompInfo(doc, pp_config, p_nComps, p_compName, firstComp,
                         options);
  if (0 != rc) {
    ERR("trix_readCompInfo: unable to parse file \"%s\"\n", p_fileName);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return -1;
  }
                                   /* read number of tri scalars and vectors */
  rc = trix_readTriDataInfo(doc, *pp_config, *p_nComps, p_triDataNames,
                            firstComp);
  if (0 != rc) {
    ERR("readTriDataInfo: unable to parse file content \"%s\"\n", p_fileName);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return -1;
  }
                                            /* read vert scalars and vectors */
  rc = trix_readVertDataInfo(doc, *pp_config, *p_nComps, p_vertDataNames,
                             firstComp);
  if (0 != rc) {
    ERR("readVertDataInfo: unable to parse file content \"%s\"\n", p_fileName);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return -1;
  }
                                /* alloc memory for component triangulations */
  for (i=firstComp; i<*p_nComps; ++i) {
    p_comp = (*pp_config) + i;
    rc = c3d_allocTriangulation(&p_comp);
    if (rc != 0) {
      ERR(" io_readSurfTrix c3d_allocTriangulation failed\n");
      xmlFreeDoc(doc);
      xmlCleanupParser();
      return -1;
    }
  } /* end for comps */

                                        /* read verts, connectivity and data */
  rc = trix_fillContent(doc, *pp_config, *p_nComps, firstComp);
  if (rc != 0) {
    ERR("trix_fillContent failed\n");
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return -1;
  }

  xmlFreeDoc(doc);
  xmlCleanupParser(); /* shutdown libxml */
  xmlMemoryDump();    /* this is to debug memory for regression tests */

  return 0;
}

/**
 * -- trix_initLibXML2(): mid-level-wrapper for LibXML2 parser
 *                        initialization only needed by apps that
 *                        handle parser init and shutdown externally
 *                        (e.g. overgrid)
 */
void trix_initLibXML2(void)
{
  xmlInitParser(); /* init libxml */
  xmlXPathInit();  /* initialize xpath */
  LIBXML_TEST_VERSION;
}

/**
 * -- trix_stopLibXML2(): mid-level-wrapper for xmlCleanupParser()
 *                        needed by apps directly calling
 *                        trix_initLibXML2()
 */
void trix_stopLibXML2(void)
{
  xmlCleanupParser(); /* shutdown libxml */
}

/**
 * -- trix_readSurf(): Mid-Level function for reading VTK sturface
 *                     tri-files, when application needs to manage
 *                     init and shutdown of the LibXML2 parser and
 *                     XPath Lib. This is basically an exact copy of
 *                     io_readSurfTrix() but with the initialization
 *                     and shutdown being removed. If an application
 *                     calls this, it needs to handle both init and
 *                     shutdown externally.  Typically the calling
 *                     would be: trix_initLibXML2(); trix_readSurf();
 *                     trix_stopLibXML2()
 */
int trix_readSurf(const char *const p_fileName,
                  p_tsTriangulation *pp_config,
                  int *const p_nComps,
                  const char *const p_compName,
                  const char *const p_vertDataNames,
                  const char *const p_triDataNames,
                  const int options)
{
  xmlDocPtr doc;
  xmlNode  *root = NULL;
  xmlParserCtxtPtr ctxt; /* the parser context */

  p_tsTriangulation p_comp = NULL;
  const int firstComp = *p_nComps;

  int rc, i, parserOpts;

  if ( TRIX_VERBOSE & options )
    printf("    o  Parsing file \"%s\" with libxml2 version %d\n", p_fileName,
           LIBXML_VERSION);

  ctxt = xmlNewParserCtxt();  /* create a parser context */
  if (ctxt == NULL) {
    ERR(" trix_readSurfTrix failed to allocate parser context\n");
    /* xmlCleanupParser(); */
    return -1;
  }

  parserOpts  = 0; /* set parser to be silent */
  parserOpts |= (XML_PARSE_NOERROR);
  parserOpts |= (XML_PARSE_NOWARNING);
#ifdef DARWIN
  parserOpts |= (XML_PARSE_HUGE); /* needed for laptops */
  /* parserOpts |= (XML_PARSE_PEDANTIC); */
#endif

                                                        /* load XML document */
  doc = xmlCtxtReadFile(ctxt, p_fileName, NULL, parserOpts);
  if (doc == NULL) {
    if ( TRIX_VERBOSE & options ) {
      CONT("%s is not valid XML (trix_readSurfTrix)\n", p_fileName);
    }
    xmlFreeParserCtxt(ctxt);
    return -1;
  }

  if ( IS_TRIX & options ) {
    xmlFreeParserCtxt(ctxt);
    return VALID_TRIX_FILE;
  }

  xmlFreeParserCtxt(ctxt); /* done with parser */

  root = xmlDocGetRootElement(doc); /* check root element node */
  if ( !root ||
       !root->name ||
       xmlStrcmp(root->name, (xmlChar *) "VTKFile") ) {
    WARN(" trix_readSurfTrix root element not VTKFile\n");
    xmlFreeDoc(doc);
    return -1;
  }
                               /* read number of Verts and Tris for compName */
  rc = trix_readCompInfo(doc, pp_config, p_nComps, p_compName,firstComp,
                         options);
  if (0 != rc) {
    ERR("trix_readCompInfo: unable to parse file \"%s\"\n", p_fileName);
    xmlFreeDoc(doc);
    return -1;
  }
                                   /* read number of tri scalars and vectors */
  rc = trix_readTriDataInfo(doc, *pp_config, *p_nComps, p_triDataNames,
                            firstComp);
  if (0 != rc) {
    ERR("readTriDataInfo: unable to parse file content \"%s\"\n", p_fileName);
    xmlFreeDoc(doc);
    return -1;
  }
                                            /* read vert scalars and vectors */
  rc = trix_readVertDataInfo(doc, *pp_config, *p_nComps, p_vertDataNames,
                             firstComp);
  if (0 != rc) {
    ERR("readVertDataInfo: unable to parse file content \"%s\"\n", p_fileName);
    xmlFreeDoc(doc);
    return -1;
  }
                                /* alloc memory for component triangulations */
  for (i=firstComp; i<*p_nComps; ++i) {
    p_comp = (*pp_config) + i;
    rc = c3d_allocTriangulation(&p_comp);
    if (rc != 0) {
      ERR(" trix_readSurfTrix c3d_allocTriangulation failed\n");
      xmlFreeDoc(doc);
      return -1;
    }
  } /* end for comps */

                                        /* read verts, connectivity and data */
  rc = trix_fillContent(doc, *pp_config, *p_nComps, firstComp);
  if (rc != 0) {
    ERR("trix_fillContent failed\n");
    xmlFreeDoc(doc);
    return -1;
  }

  xmlFreeDoc(doc);

  return 0;
}

/**
 * -- io_readTrixIntersectDims(): Basic function that returns
 *                                dimensions of triangulation in
 *                                p_fileName. Passing pointers and
 *                                returning void so it is callable
 *                                from fortran, i.e. Intersect.  If
 *                                any returned value is -1, this
 *                                indicates an error condition.
 *                                Should pass in name of component;
 *                                future work.
 */
void io_readTrixIntersectDims(const char *const p_fileName,
                              int *p_nVerts, int *p_nTris,
                              int *p_nVertScalars, int *p_nTriScalars,
                              int *p_nVertData, int *p_nTriData,
                              const int options)
{
  xmlDocPtr doc;
  xmlNode  *root = NULL;
  xmlParserCtxtPtr ctxt; /* the parser context */

  xmlXPathContextPtr xpathCtx;
  xmlXPathObjectPtr  xpathObj;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         cur_node;
  xmlChar           *p_xmlstr;

  char               xpathExpr[STRING_LEN];
  int                rc, nNodes, parserOpts;

  *p_nVerts = *p_nTris = *p_nVertScalars = *p_nTriScalars = *p_nVertData =
    *p_nTriData = -1;

  { LIBXML_TEST_VERSION };

  xmlXPathInit();  /* initialize xpath */

  if ( TRIX_VERBOSE & options )
    printf("    o  Parsing file \"%s\" with libxml2 version %d\n", p_fileName,
           LIBXML_VERSION);

  /* #ifndef DARWIN   */
  /*   if (xmlHasFeature((xmlFeature)31)) { /\* check if there is ZLIB *\/ */
  /*     NOTE("Have ZLIB support\n"); */
  /*   } */
  /*   else { */
  /*     NOTE("No ZLIB support\n"); */
  /*   } */
  /* #endif */

  ctxt = xmlNewParserCtxt();  /* create a parser context */
  if (ctxt == NULL) {
    ERR(" io_readSurfTrix failed to allocate parser context\n");
    xmlCleanupParser();
    return;
  }

  parserOpts = 0; /* set parser to be silent */
  parserOpts |= (XML_PARSE_NOERROR);
  parserOpts |= (XML_PARSE_NOWARNING);
  parserOpts |= (XML_PARSE_HUGE); /* relax any hardcoded limit from the parser */
  /* load XML document */
  doc = xmlCtxtReadFile(ctxt, p_fileName, NULL, parserOpts);
  if (doc == NULL) {
    if ( TRIX_VERBOSE & options ) {
      CONT("%s does not appear to be valid XML (io_readTrixIntersectDims)\n", p_fileName);
    }
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    return;
  }

  xmlFreeParserCtxt(ctxt); /* done with parser */

  root = xmlDocGetRootElement(doc); /* check root element node */
  if ( !root ||
       !root->name ||
       xmlStrcmp(root->name, (xmlChar *) "VTKFile") ) {
    WARN(" io_readSurfTrix root element not VTKFile\n");
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("trix_readCompInfo: xmlXPathNewContext failed to create new XPath\n");
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
                                                   /* extract all components */
  strcpy(xpathExpr,"//UnstructuredGrid/Piece");

  xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if(xpathObj == NULL) {
    ERR("trix_readCompInfo: failed to evaluate \"%s\"\n", xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }

  nodes  = xpathObj->nodesetval;
  nNodes = (nodes) ? nodes->nodeNr : 0;
  if (0 == nNodes) {
    ERR("io_readTrixDims: Failed to find any %s elements\n", xpathExpr);
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }

  if (1 != nNodes) {
    ERR("io_readTrixDims: Too many %s elements\n", xpathExpr);
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
                                                    /* scan nverts and ntris */
  ASSERT(nodes->nodeTab[0]);
  ASSERT(nodes->nodeTab[0]->type == XML_ELEMENT_NODE);

  cur_node = nodes->nodeTab[0];

  p_xmlstr = xmlGetProp(cur_node, BAD_CAST "NumberOfPoints");
  if (NULL == p_xmlstr) {
    ERR("trix_readCompInfo: PIECE element missing NumberOfPoints\n");
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
  rc = sscanf( (const char *) p_xmlstr, "%d", p_nVerts);
  if (0 == rc) {
    ERR("Unable to convert NumberOfPoints %s\n", (char *) p_xmlstr);
    xmlFree(p_xmlstr);
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
  xmlFree(p_xmlstr);

  p_xmlstr = xmlGetProp(cur_node, BAD_CAST "NumberOfCells");
  if (NULL == p_xmlstr) {
    ERR("PIECE element missing NumberOfCells\n");
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
  rc = sscanf( (const char *) p_xmlstr,"%d", p_nTris);
  if (0 == rc) {
    ERR("Unable to convert NumberOfCells %s\n", (char *) p_xmlstr);
    xmlFree(p_xmlstr);
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
  xmlFree(p_xmlstr);

  xmlXPathFreeObject(xpathObj);
  nodes    = NULL;
  cur_node = NULL;
  /* take only COMPONENT_TAGS from triData */
  snprintf(xpathExpr,STRING_LEN*sizeof(char),
           "//UnstructuredGrid/Piece/CellData/DataArray[@TRIXtype='%s']",
           "COMPONENT_TAG");

  xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if(xpathObj == NULL) {
    ERR("Unable to evaluate expression \"%s\"\n", (char *) xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
  nodes  = xpathObj->nodesetval;
  nNodes = (nodes) ? nodes->nodeNr : 0;
  /* subtract one for IntersectComponents */
  (*p_nTriScalars) = (*p_nTriData) = nNodes - 1;

  xmlXPathFreeObject(xpathObj);
  nodes = NULL;
                              /* take only SHAPE_LINEARIZATION from vertData */
  snprintf(xpathExpr,STRING_LEN*sizeof(char),
           "//UnstructuredGrid/Piece/PointData/DataArray[@TRIXtype='%s']",
           "SHAPE_LINEARIZATION");

  xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if(xpathObj == NULL) {
    ERR("Unable to evaluate expression \"%s\"\n", (char *) xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return;
  }
  nodes  = xpathObj->nodesetval;
  nNodes = (nodes) ? nodes->nodeNr : 0;

  *p_nVertData    = nNodes;
  *p_nVertScalars = nNodes*3;

  xmlXPathFreeObject(xpathObj);
  nodes = NULL;

  xmlXPathFreeContext(xpathCtx);

  xmlFreeDoc(doc);
  xmlCleanupParser(); /* shutdown libxml */
  xmlMemoryDump();    /* this is to debug memory for regression tests */

  return;
}

/**
 * -- trix_readCompInfo(): for a requested component, verify component
 *                         name and parse out nVerts and nTris
 */
int trix_readCompInfo(xmlDocPtr doc, p_tsTriangulation *pp_config,
                      int *const p_nComps, const char *p_compName,
                      const int firstComp, const int options)
{
  xmlXPathContextPtr xpathCtx;
  xmlXPathObjectPtr  xpathObj;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         cur_node;
  xmlChar           *p_xmlstr;

  char               xpathExpr[STRING_LEN];
  int                i, rc, nNodes;

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("trix_readCompInfo: xmlXPathNewContext failed to create new XPath\n");
    return(-1);
  }

  if ( 0==strcmp("ALL",p_compName) ) { /* extract all components */
    strcpy(xpathExpr,"//UnstructuredGrid/Piece");
  }
  else {
    snprintf(xpathExpr,STRING_LEN*sizeof(char),
             "//UnstructuredGrid/Piece[@Name='%s']",
             p_compName);
  }

  xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if(xpathObj == NULL) {
    ERR("trix_readCompInfo: failed to evaluate \"%s\"\n", xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    return(-1);
  }
  nodes  = xpathObj->nodesetval;
  nNodes = (nodes) ? nodes->nodeNr : 0;
  if (0 == nNodes) {
    ERR("trix_readCompInfo: Failed to find any %s elements\n", xpathExpr);
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    return(-1);
  }
                        /* allocate and initialize tri structs for each comp */
  rc = c3d_newTriangulation(pp_config, firstComp, (*p_nComps) + nNodes);
  if (0 != rc) {
    ERR("trix_readCompInfo: c3d_newTriangulation failed\n");
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    return(-1);
  }
                                   /* scan component names, nverts and ntris */
  for (i = 0; i < nNodes; ++i) {
    ASSERT(nodes->nodeTab[i]);
    ASSERT(nodes->nodeTab[i]->type == XML_ELEMENT_NODE);

    cur_node = nodes->nodeTab[i];

    p_xmlstr = xmlGetProp(cur_node, BAD_CAST "Name");
    if (NULL == p_xmlstr) {
      ERR("trix_readCompInfo: PIECE element NOT named\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    strcpy( (*pp_config)[*p_nComps].geomName, (const char *) p_xmlstr);
    xmlFree(p_xmlstr);

    p_xmlstr = xmlGetProp(cur_node, BAD_CAST "NumberOfPoints");
    if (NULL == p_xmlstr) {
      ERR("trix_readCompInfo: PIECE element missing NumberOfPoints\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    rc = sscanf( (const char *) p_xmlstr, "%d",
                 &(*pp_config)[*p_nComps].nVerts);
    if (0 == rc) {
      ERR("Unable to convert NumberOfPoints %s\n", (char *) p_xmlstr);
      xmlFree(p_xmlstr);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    xmlFree(p_xmlstr);

    p_xmlstr = xmlGetProp(cur_node, BAD_CAST "NumberOfCells");
    if (NULL == p_xmlstr) {
      ERR("PIECE element missing NumberOfCells\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    rc = sscanf( (const char *) p_xmlstr,"%d",
                 &(*pp_config)[*p_nComps].nTris);
    if (0 == rc) {
      ERR("Unable to convert NumberOfCells %s\n", (char *) p_xmlstr);
      xmlFree(p_xmlstr);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    xmlFree(p_xmlstr);

    if ( TRIX_DP_VERTS & options ) {
      /* set infoCode to alloc double precision verts */
      (*pp_config + *p_nComps)->infoCode = DP_VERTS_CODE;
    }

    (*p_nComps)++; /* increase nComps */
  }

  xmlXPathFreeObject(xpathObj);
  nodes    = NULL;
  cur_node = NULL;
  xmlXPathFreeContext(xpathCtx);

  return 0;
}

/**
 * -- trix_setTriInfo()
 */
int trix_setTriInfo(p_tsTriangulation p_config, const int nComps,
                    const char *p_compName, const xmlNodePtr cur_node,
                    const int firstComp)
{
  int               i, rc;
  xmlChar          *p_tmpstr;
  char              name[STRING_LEN];
  p_tsTriangulation p_comp = NULL;

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "Name"); /* get name of DataArray */
  if (NULL == p_tmpstr) {
    ERR("setTriInfo: xmlGetProp failed to get NAME\n");
    return -1;
  }
  strcpy(name, (char *) p_tmpstr);
  xmlFree(p_tmpstr);

  for (i=firstComp; i<nComps; ++i) { /* find comp in p_config */
    if ( 0 == strcmp(p_compName, p_config[i].geomName) ) {
      p_comp = p_config + i;
      break;
    }
  }

  if ( NULL == p_comp ) {
    WARN("setTriInfo: Component \"%s\" cannot be set - skipping DataArray\n",
         p_compName);
    return 0;
  }

  rc = c3d_allocTriData(&p_comp, 1); /* alloc one triData struct */
  if (0 != rc) {
    ERR("c3d_allocTriData failed\n");
  }

  strcpy(p_comp->p_triData[p_comp->nTriData-1].name, name);

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "NumberOfComponents");
  if (NULL == p_tmpstr) {
    ERR("setTriInfo: xmlGetProp failed to get NumberOfComponents\n");
    return -1;
  }
  rc = sscanf( (const char *) p_tmpstr,"%d",
               &p_comp->p_triData[p_comp->nTriData-1].dim);
  xmlFree(p_tmpstr);
  if (0 == rc) {
    ERR("setTriInfo: Unable to convert NumberOfComponents %s\n",
        (char *) p_tmpstr);
    return -1;
  }

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "type");
  if (NULL == p_tmpstr) {
    ERR("setTriInfo: xmlGetProp failed to get TYPE\n");
    return -1;
  }
  p_comp->p_triData[p_comp->nTriData-1].type =
    getVTKtype( (char *) p_tmpstr, xmlStrlen(p_tmpstr));
  xmlFree(p_tmpstr);

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "TRIXtype");
  if (NULL == p_tmpstr) {
    WARN("trix_setTriInfo: xmlGetProp failed to get TRIXtype, using unset\n");
    p_comp->p_triData[p_comp->nTriData-1].info = TRIX_unset;
  }
  else {
    p_comp->p_triData[p_comp->nTriData-1].info =
      getTRIXtype( (char *) p_tmpstr, xmlStrlen(p_tmpstr));
    xmlFree(p_tmpstr);
  }

                              /* compute offset into to a_scalar_t tri array */
  p_comp->p_triData[p_comp->nTriData-1].offset = 0;
  for (i=0; i<p_comp->nTriData-1; ++i) {
    p_comp->p_triData[p_comp->nTriData-1].offset += p_comp->p_triData[i].dim;
  }

  return 0;
}

/**
 * trix_readTriDataInfo()
 */
int trix_readTriDataInfo(xmlDocPtr doc, p_tsTriangulation p_config,
                         const int nComps, const char *p_dataName,
                         const int firstComp)
{
  xmlXPathContextPtr xpathCtx;
  xmlXPathObjectPtr  xpathObj;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         parent;
  xmlChar           *p_tmpstr;
  char               xpathExpr[STRING_LEN];
  int                i, nNodes;

  if (NULL == p_dataName) {
    WARN("readTriDataInfo: p_dataName is NULL\n");
    return 0;
  }

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("Failed to create new XPath context\n");
    return(-1);
  }

  if ( 0==strcmp("ALL",p_dataName) ) { /* extract all */
    strcpy(xpathExpr,"//UnstructuredGrid/Piece/CellData/DataArray");
  }
  else {
    snprintf(xpathExpr,STRING_LEN*sizeof(char),
             "//UnstructuredGrid/Piece/CellData/DataArray[@Name='%s']",
             p_dataName);
  }

  xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if(xpathObj == NULL) {
    ERR("Unable to evaluate expression \"%s\"\n", (char *) xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    return(-1);
  }
  nodes  = xpathObj->nodesetval;
  nNodes = (nodes) ? nodes->nodeNr : 0;

              /* loop over nodes, get nTriData and save it with correct comp */
  for (i=0; i<nNodes; ++i) {
    /* this is the 'Piece' element, which holds comp name */
    parent = (nodes->nodeTab[i]->parent)->parent;

    if ( 0 != strcmp((char *)parent->name,"Piece") ) {    /* must be 'Piece' */
      ERR("readTriDataInfo: Failed to find PIECE\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    p_tmpstr = xmlGetProp(parent, BAD_CAST "Name");
    if (NULL == p_tmpstr) {
      ERR("readTriDataInfo: xmlGetProp failed to find NAME of Piece\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
                                                          /* set tri scalars */
    trix_setTriInfo(p_config, nComps, (const char *) p_tmpstr,
                    nodes->nodeTab[i], firstComp);
    xmlFree(p_tmpstr);
  } /* end for nodes */

  xmlXPathFreeObject(xpathObj);
  nodes  = NULL;
  xmlXPathFreeContext(xpathCtx);
  return 0;
}

/**
 * trix_setVertInfo()
 */
int trix_setVertInfo(p_tsTriangulation p_config, int nComps,
                     const char *p_compName, const xmlNodePtr cur_node,
                     const int firstComp)
{
  int               i, rc;
  xmlChar          *p_tmpstr;
  char              name[STRING_LEN];
  p_tsTriangulation p_comp = NULL;

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "Name");
  if (NULL == p_tmpstr) {
    ERR("setVertInfo: xmlGetProp failed to get NAME\n");
    return -1;
  }
  strcpy(name, (char *) p_tmpstr);
  xmlFree(p_tmpstr);

  for (i=firstComp; i<nComps; ++i) { /* find component in p_config */
    if ( 0 == strcmp(p_compName, p_config[i].geomName) ) {
      p_comp = p_config + i;
      break;
    }
  }

  if ( NULL == p_comp ) {
    WARN("setVertInfo: Component \"%s\" cannot be set - skipping DataArray\n",
         p_compName);
    return 0;
  }

  rc = c3d_allocVertData(&p_comp, 1);
  if (0 != rc) {
    ERR("c3d_allocVertData failed\n");
  }

  strcpy(p_comp->p_vertData[p_comp->nVertData-1].name, name);

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "NumberOfComponents");
  if (NULL == p_tmpstr) {
    ERR("setVertInfo: xmlGetProp failed to get NumberOfComponents\n");
    return -1;
  }
  rc = sscanf( (const char *) p_tmpstr,"%d",
               &p_comp->p_vertData[p_comp->nVertData-1].dim);
  xmlFree(p_tmpstr);
  if (0 == rc) {
    ERR("setVertInfo: sscanf failed to convert NumberOfComponents %s\n",
        (char *)p_tmpstr);
    return -1;
  }

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "type");
  if (NULL == p_tmpstr) {
    ERR("setVertInfo: xmlGetProp failed to get TYPE\n");
    return -1;
  }
  p_comp->p_vertData[p_comp->nVertData-1].type =
    getVTKtype( (char *) p_tmpstr, xmlStrlen(p_tmpstr));
  xmlFree(p_tmpstr);

  p_tmpstr = xmlGetProp(cur_node, BAD_CAST "TRIXtype");
  if (NULL == p_tmpstr) {
    WARN("trix_setVertInfo: xmlGetProp failed to get TRIXtype, using unset\n");
    p_comp->p_vertData[p_comp->nVertData-1].info = TRIX_unset;
  }
  else {
    p_comp->p_vertData[p_comp->nVertData-1].info =
      getTRIXtype( (char *) p_tmpstr, xmlStrlen(p_tmpstr));
    xmlFree(p_tmpstr);
  }

                              /* compute offset into to a_scalar_t tri array */
  p_comp->p_vertData[p_comp->nVertData-1].offset = 0;
  for (i=0; i<p_comp->nVertData-1; ++i) {
    p_comp->p_vertData[p_comp->nVertData-1].offset +=p_comp->p_vertData[i].dim;
  }

  return 0;
}

/**
 * trix_readVertDataInfo()
 */
int trix_readVertDataInfo(xmlDocPtr doc, p_tsTriangulation p_config,
                          const int nComps, const char *p_dataName,
                          const int firstComp)
{
  xmlXPathContextPtr xpathCtx;
  xmlXPathObjectPtr  xpathObj;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         parent;
  xmlChar           *p_tmpstr;
  char               xpathExpr[STRING_LEN];
  int                i, nNodes;

  if (NULL == p_dataName) {
    WARN("readVertDataInfo: p_dataName is NULL\n");
    return 0;
  }

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("Failed to create new XPath context\n");
    return(-1);
  }

  if ( 0==strcmp("ALL",p_dataName) ) {
    strcpy(xpathExpr,"//UnstructuredGrid/Piece/PointData/DataArray");
  }
  else {
    snprintf(xpathExpr, STRING_LEN*sizeof(char),
             "//UnstructuredGrid/Piece/PointData/DataArray[@Name='%s']",
             p_dataName);
  }

  xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if(xpathObj == NULL) {
    ERR("Failed to evaluate expression \"%s\"\n", (char *) xpathExpr);
    xmlXPathFreeContext(xpathCtx);
    return(-1);
  }
  nodes  = xpathObj->nodesetval;
  nNodes = (nodes) ? nodes->nodeNr : 0;
  /* loop over nodes, get nTriData and save it with correct comp */
  for (i=0; i<nNodes; ++i) {
    /* this is the 'Piece' element, which holds comp name */
    parent = (nodes->nodeTab[i]->parent)->parent;

    if ( 0 != strcmp((char *)parent->name,"Piece") ) {    /* must be 'Piece' */
      ERR("readVertDataInfo: node name is not Piece\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    p_tmpstr = xmlGetProp(parent, BAD_CAST "Name");
    if (NULL == p_tmpstr) {
      ERR("readVertDataInfo: xmlGetProp failed to find NAME\n");
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    /* set tri scalars */
    trix_setVertInfo(p_config, nComps, (const char *) p_tmpstr,
                     nodes->nodeTab[i], firstComp);
    xmlFree(p_tmpstr);
  } /* end for nodes */

  xmlXPathFreeObject(xpathObj);
  nodes  = NULL;
  xmlXPathFreeContext(xpathCtx);
  return 0;
}

/**
 * -- trix_readVertCoord(): read in all vert coordinates
 */
int trix_readVertCoord(p_tsTriangulation p_config, const int icomp,
                       const xmlNodePtr cur_node)
{
  int      i, j, nTokens, rc;
  char    *p_end = NULL;
  xmlChar *p_xmlstr;
  VTKtype  type;

  const p_tsTriangulation p_comp = p_config + icomp;

  ASSERT(XML_ELEMENT_NODE == cur_node->type);

  p_xmlstr = xmlGetProp(cur_node, BAD_CAST "NumberOfComponents");
  rc = sscanf( (const char *) p_xmlstr,"%d", &nTokens);
  if (0 == rc) {
    ERR("Failed to convert nverts %s\n", (char *) p_xmlstr);
    xmlFree(p_xmlstr);
    return(-1);
  }
  xmlFree(p_xmlstr);

  ASSERT(3 == nTokens);

  p_xmlstr = xmlGetProp(cur_node, BAD_CAST "type");
  if (NULL == p_xmlstr) {
    ERR("DataArray element missing type\n");
    xmlFree(p_xmlstr);
    return(-1);
  }
  type = getVTKtype( (char *) p_xmlstr, xmlStrlen(p_xmlstr));
  xmlFree(p_xmlstr);

  p_xmlstr = xmlNodeGetContent(cur_node); /* get all verts */
  if (NULL == p_xmlstr) {
    ERR("xmlNodeGetContent failed in readVertCoord\n");
    return -1;
  }

  p_end = (char *) p_xmlstr;
  if (VTK_Float64 == type && p_comp->a_dpVerts) {
    NOTE("c3dio: double precision read of vert coordinates\n");
    for (i=0; i<p_comp->nVerts; ++i) {
      for (j=0; j<nTokens; ++j) {
        p_comp->a_dpVerts[i].x[j] = strtod(p_end, &p_end);
        p_comp->a_Verts[i].x[j]   = (float) p_comp->a_dpVerts[i].x[j];
      }
    }
  }
  else if (VTK_Float32 == type && p_comp->a_dpVerts) {
    /* If the actual data type is single precision, then do not populate the
     * dpVerts array. Instead, free it and override the user request. This
     * could be detected in readCompInfo but would require more parsing there.
     */
    WARN("c3dio: double precision read requested, but file vert coordinates\n");
    CONT("are in single precision, continuing with single precision\n");
    free(p_comp->a_dpVerts);
    p_comp->a_dpVerts = NULL;
    p_comp->infoCode = p_comp->infoCode & (~ DP_VERTS_CODE); /* reset DP flag */
    for (i=0; i<p_comp->nVerts; ++i) {
      for (j=0; j<nTokens; ++j) {
        p_comp->a_Verts[i].x[j] = (float) strtod(p_end, &p_end);
      }
    }
  }
  else {
    for (i=0; i<p_comp->nVerts; ++i) {
      for (j=0; j<nTokens; ++j) {
        p_comp->a_Verts[i].x[j] = (float) strtod(p_end, &p_end);
      }
    }
  }

  xmlFree(p_xmlstr);
  return 0;
}

/**
 * -- trix_readConnectivity(): read connectivity. We assume triangles
 *                             for element type, which have
 *                             connectivity offset of 3. This way we
 *                             don't have to read the 'offsets' and
 *                             'type' arrays. VTK format uses zero
 *                             offset.
 */
int trix_readConnectivity(p_tsTriangulation p_config, const int icomp,
                          const xmlNodePtr cur_node)
{
  int      i, j;
  char    *p_end = NULL;
  xmlChar *p_xmlstr;

  const int nTokens = 3;
  const p_tsTriangulation p_comp = p_config + icomp;

  ASSERT(XML_ELEMENT_NODE == cur_node->type);

  p_xmlstr = xmlNodeGetContent(cur_node);
  if (NULL == p_xmlstr) {
    ERR("xmlNodeGetContent failed in readConnectivity\n");
    return -1;
  }

  p_end = (char *) p_xmlstr;
  for (i=0; i<p_comp->nTris; ++i) {
    for (j=0; j<nTokens; ++j) {
      p_comp->a_Tris[i].vtx[j] = (int) strtol(p_end, &p_end, 10);
    }
    p_comp->a_Tris[i].Comp = BAD_INDEX; /* initialized to BAD_INDEX */
  }

  xmlFree(p_xmlstr);
  return 0;
}

/**
 * -- trix_readTriData()
 */
int trix_readTriData(xmlXPathContextPtr xpathCtx, p_tsTriangulation p_comp)
{
  xmlXPathObjectPtr xpathObj;
  xmlNodeSetPtr     nodes;
  xmlChar          *p_xmlstr;

  int  i, j, k, nNodes=0, indx=0;
  char xpathExpr[STRING_LEN], *p_end = NULL;

  for (j=0; j<p_comp->nTriData; ++j) {
    snprintf(xpathExpr,STRING_LEN*sizeof(char),
             ".//CellData/DataArray[@Name='%s']",p_comp->p_triData[j].name);

    xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
    if(xpathObj == NULL) {
      ERR("readTriData: Failed to evaluate expression \"%s\"\n",
          (char *) xpathExpr);
      return(-1);
    }
    nodes  = xpathObj->nodesetval;
    nNodes = (nodes) ? nodes->nodeNr : 0;
    if (1 != nNodes) {
      ERR("readTriData: Too many nodes %d\n", nNodes);
      xmlXPathFreeObject(xpathObj);
      return(-1);
    }

    p_xmlstr = xmlNodeGetContent(nodes->nodeTab[0]);
    if (NULL == p_xmlstr) {
      ERR("readTriData: Failed to get content\n");
      xmlXPathFreeObject(xpathObj);
      return -1;
    }

    p_end = (char *) p_xmlstr;
    for (i=0; i<p_comp->nTris; ++i) {
      for (k=0; k<p_comp->p_triData[j].dim; ++k) {
        indx = k + i*(p_comp->p_triData[j].dim) +
          p_comp->p_triData[j].offset*p_comp->nTris;
        p_comp->a_scalar0_t[indx] = strtod(p_end, &p_end);
      }
    }

    xmlFree(p_xmlstr); /* cleanup */
    xmlXPathFreeObject(xpathObj);
    nodes = NULL;
  } /* done for TriData */

  return 0;
}

/**
 * -- trix_readVertData()
 */
int trix_readVertData(xmlXPathContextPtr xpathCtx,
                      const p_tsTriangulation p_comp)
{
  xmlXPathObjectPtr xpathObj;
  xmlNodeSetPtr     nodes;
  xmlChar          *p_xmlstr;

  int  i, j, k, nNodes=0, indx=0;
  char xpathExpr[STRING_LEN], *p_end=NULL;

  for (j=0; j<p_comp->nVertData; ++j) {
    snprintf(xpathExpr, STRING_LEN*sizeof(char),
             ".//PointData/DataArray[@Name='%s']",p_comp->p_vertData[j].name);
    xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
    if(xpathObj == NULL) {
      ERR("Unable to evaluate expression \"%s\"\n", xpathExpr);
      return(-1);
    }
    nodes  = xpathObj->nodesetval;
    nNodes = (nodes) ? nodes->nodeNr : 0;
    if (1 != nNodes) {
      ERR("Too many nodes %d\n", nNodes);
      xmlXPathFreeObject(xpathObj);
      return(-1);
    }

    p_xmlstr = xmlNodeGetContent(nodes->nodeTab[0]); /* this is the data */
    if (NULL == p_xmlstr) {
      ERR("Unable to get content\n");
      xmlXPathFreeObject(xpathObj);
      return -1;
    }

    p_end = (char *) p_xmlstr;
    for (i=0; i<p_comp->nVerts; ++i) {
      for (k=0; k<p_comp->p_vertData[j].dim; ++k) {
        indx = k + i*(p_comp->p_vertData[j].dim) +
          p_comp->p_vertData[j].offset*(p_comp->nVerts);
        p_comp->a_scalar0[indx] = strtod(p_end, &p_end);
      }
    }

    xmlFree(p_xmlstr); /* cleanup */
    xmlXPathFreeObject(xpathObj);
    nodes = NULL;
  } /* done for TriData */

  return 0;
}

/**
 * -- trix_fillContent(): driver function to read in vert coordinates,
 *                        connectivity, and tri and vert data
 */
int trix_fillContent(xmlDocPtr doc, const p_tsTriangulation p_config,
                     const int nComps, const int firstComp)
{
  xmlXPathContextPtr xpathCtx;
  xmlXPathObjectPtr  xpathObj, xpathSubset;
  xmlNodeSetPtr      nodes;
  char               xpathExpr[STRING_LEN];
  int                nNodes, icomp;
  p_tsTriangulation  p_comp;

  xpathCtx = xmlXPathNewContext(doc); /* create xpath evaluation context */
  if(xpathCtx == NULL) {
    ERR("trix_fillContent: unable to create new XPath context\n");
    return(-1);
  }

  for (icomp=firstComp; icomp<nComps; ++icomp) {
    p_comp = p_config + icomp;
                                          /* grab xml node of this component */
    snprintf(xpathExpr,sizeof(xpathExpr),
             "//UnstructuredGrid/Piece[@Name='%s']",
             p_comp->geomName);

    xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
    if(xpathObj == NULL) {
      ERR("Unable to evaluate expression \"%s\"\n", xpathExpr);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    nodes  = xpathObj->nodesetval;
    nNodes = (nodes) ? nodes->nodeNr : 0;
    if (1 != nNodes) { /* there can be only one component with this name */
      ERR("Too many nodes %d %s\n", nNodes,xpathExpr);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
                                    /* set xpath root node to this component */
    xpathCtx->node = xpathObj->nodesetval->nodeTab[0];
                                                              /* vert coords */
    snprintf(xpathExpr,sizeof(xpathExpr), ".//Points/DataArray");

    xpathSubset = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
    if(xpathSubset == NULL) {
      ERR("Unable to evaluate expression \"%s\"\n", xpathExpr);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    nodes  = xpathSubset->nodesetval;
    nNodes = (nodes) ? nodes->nodeNr : 0;
    if (1 != nNodes) {
      ERR("Too many nodes %d %s\n", nNodes,xpathExpr);
      xmlXPathFreeObject(xpathSubset);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }

    trix_readVertCoord(p_config, icomp, nodes->nodeTab[0]);

    xmlXPathFreeObject(xpathSubset); /* cleanup */
    nodes = NULL;
                                                        /* read connectivity */
    snprintf(xpathExpr,sizeof(xpathExpr),
             ".//Cells/DataArray[@Name='connectivity']");

    xpathSubset = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
    if(xpathSubset == NULL) {
      ERR("readTriX: Failed to evaluate expression \"%s\"\n", xpathExpr);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }
    nodes  = xpathSubset->nodesetval;
    nNodes = (nodes) ? nodes->nodeNr : 0;
    if (1 != nNodes) {
      ERR("readTriX: Too many nodes %d\n", nNodes);
      xmlXPathFreeObject(xpathSubset);
      xmlXPathFreeObject(xpathObj);
      xmlXPathFreeContext(xpathCtx);
      return(-1);
    }

    trix_readConnectivity(p_config, icomp, nodes->nodeTab[0]);

    xmlXPathFreeObject(xpathSubset); /* cleanup */
    nodes = NULL;
                                                   /* read tri and vert data */
    trix_readTriData(xpathCtx, p_comp);
    trix_readVertData(xpathCtx, p_comp);

    xmlXPathFreeObject(xpathObj); /* cleanup Piece */
  }

  xmlXPathFreeContext(xpathCtx);

  return 0;
}

#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

/**
 * -- trix_writeBasicTri()
 */
int trix_writeBasicTri(xmlTextWriterPtr writer, const p_tsTriangulation p_comp,
                       const int options)
{
  int rc, i;

                     /* ---------- VTK Points (vert coordinates) ----------- */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "Points");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterStartElement(writer, BAD_CAST "DataArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  if (TRIX_DP_VERTS & options && p_comp->a_dpVerts) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                     BAD_CAST "Float64");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }
  else {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                     BAD_CAST "Float32");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "NumberOfComponents",
                                   BAD_CAST "3");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "format",
                                   BAD_CAST "ascii");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  xmlTextWriterWriteFormatString(writer,"\n");

  /* note: 'options' checks seems unnecessary, more straightforward to just
   * write doubles if doubles are present, otherwise write singles
   */
  if (TRIX_DP_VERTS & options && p_comp->a_dpVerts) {
    NOTE("c3dio: double precision write of vert coordinates\n");
    for (i=0; i<p_comp->nVerts; ++i) {
      /* IEEE Standard 754: writing out 17 decimals so when the number is
       * converted back it will match the original. No error checking for
       * speed
       */
      xmlTextWriterWriteFormatString(writer,
                                     " %26.17e %26.17e %26.17e\n",
                                     p_comp->a_dpVerts[i].x[X],
                                     p_comp->a_dpVerts[i].x[Y],
                                     p_comp->a_dpVerts[i].x[Z]);
    }
  }
  else {
    for (i=0; i<p_comp->nVerts; ++i) {
      /* IEEE Standard 754: writing out 9 decimals so when the number is
       * converted back it will match the original. No error checking for
       * speed
       */
      xmlTextWriterWriteFormatString(writer,
                                     " %16.9e %16.9e %16.9e\n",
                                     p_comp->a_Verts[i].x[X],
                                     p_comp->a_Verts[i].x[Y],
                                     p_comp->a_Verts[i].x[Z]);
    }
  }

  rc = xmlTextWriterEndElement(writer); /* end DataArray */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  rc = xmlTextWriterEndElement(writer); /* end Points */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

                      /* ---------- VTK Cells (tri connectivity) ----------- */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "Cells");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterStartElement(writer, BAD_CAST "DataArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                   BAD_CAST "Int32");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name",
                                   BAD_CAST "connectivity");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "format",
                                   BAD_CAST "ascii");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  xmlTextWriterWriteFormatString(writer,"\n");

  for (i=0; i<p_comp->nTris; ++i) { /* no error checking for speed */
    xmlTextWriterWriteFormatString(writer, " %d %d %d\n",
                                   p_comp->a_Tris[i].vtx[0],
                                   p_comp->a_Tris[i].vtx[1],
                                   p_comp->a_Tris[i].vtx[2]);
  }

  rc = xmlTextWriterEndElement(writer); /* end DataArray */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

                  /* --------------- VTK offsets --------------------------- */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "DataArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                   BAD_CAST "Int32");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name",
                                   BAD_CAST "offsets");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "format",
                                   BAD_CAST "ascii");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  xmlTextWriterWriteFormatString(writer,"\n");

  for (i=0; i<p_comp->nTris; ++i) { /* no error checking for speed */
    xmlTextWriterWriteFormatString(writer, " %d", (i+1)*3);
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

                  /* --------------- VTK type ------------------------------ */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "DataArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                   BAD_CAST "UInt32");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name",
                                   BAD_CAST "types");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "format",
                                   BAD_CAST "ascii");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  xmlTextWriterWriteFormatString(writer,"\n");

  for (i=0; i<p_comp->nTris; ++i) { /* no error checking for speed */
    xmlTextWriterWriteFormatString(writer, " %d", 5);
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  rc = xmlTextWriterEndElement(writer); /* end Cells */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  return 0;
}

/**
 * -- trix_writeData()
 */
int trix_writeData(xmlTextWriterPtr writer, const char *const p_name,
                   const int dim, const int offset, const double *const p_data,
                   const int nd, const int type, const int info)
{
  int rc, i, j, indx;
  const double zero = 0.;

  rc = xmlTextWriterStartElement(writer, BAD_CAST "DataArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name", BAD_CAST p_name);
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "NumberOfComponents",
                                         "%d", dim);
  if (rc < 0) {
    ERR("xmlTextWriterWriteFormatAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                   BAD_CAST vtkTypes[type]);
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "format",
                                   BAD_CAST "ascii");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "TRIXtype",
                                   BAD_CAST trixTypes[info]);
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  xmlTextWriterWriteFormatString(writer,"\n");

  if ( 7 > type ) { /* type is integer */
    for (i=0; i<nd; ++i) {
      for (j=0; j<dim; ++j) {
        indx = j+i*dim + nd*offset; /* no error checking for speed */
        xmlTextWriterWriteFormatString(writer," %d",(int) p_data[indx]);
      }
      xmlTextWriterWriteFormatString(writer,"\n");
    }
  }
  else {
    for (i=0; i<nd; ++i) {
      for (j=0; j<dim; ++j) {
        indx = j+i*dim + nd*offset; /* no error checking for speed */
        if (ABS(p_data[indx]) < 1.e-14) { /* MACHINE_ZERO */
          xmlTextWriterWriteFormatString(writer," %22.14e", zero);
        } else {
          xmlTextWriterWriteFormatString(writer," %22.14e", p_data[indx]);
        }
      }
      xmlTextWriterWriteFormatString(writer,"\n");
    }
  }

  rc = xmlTextWriterEndElement(writer); /* end DataArray */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  return 0;
}

/**
 * -- trix_writeComps(): here we write out component tags from
 *                       a_Tris[].Comps into a data-array called
 *                       IntersectComponents.  Historically, this data
 *                       field was only needed by Intersect, but other
 *                       codes, e.g. GMP, use it too
 */
int trix_writeComps(xmlTextWriterPtr writer, const int nTris, const  p_tsTri p_Tris)
{
  int rc, i;

  rc = xmlTextWriterStartElement(writer, BAD_CAST "DataArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name",
                                   BAD_CAST "IntersectComponents");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "NumberOfComponents",
                                         "%d", 1);
  if (rc < 0) {
    ERR("xmlTextWriterWriteFormatAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                   BAD_CAST vtkTypes[3]);
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "format",
                                   BAD_CAST "ascii");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "TRIXtype",
                                   BAD_CAST trixTypes[3]);
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  xmlTextWriterWriteFormatString(writer,"\n");

  for (i=0; i<nTris; ++i) {
    /* no error checking for speed and use one-offset */
    xmlTextWriterWriteFormatString(writer," %d\n", p_Tris[i].Comp+1);
  }

  rc = xmlTextWriterEndElement(writer); /* end DataArray */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  return 0;
}

/**
 * -- io_writeSurfTrix(): Top-Level API routine to write configuration
 *                        in xml format.  Note wr fill an xml buffer
 *                        is used for speed. The XML buffer is dumped
 *                        to disk at the end of this function.
 */
int io_writeSurfTrix(const p_tsTriangulation p_config, const int nComps,
                     const char *const p_fileName, const int options)
{
  int i, j, rc, tritags, cellDataOpen=0;
  xmlTextWriterPtr writer;
  xmlBufferPtr buf;

  p_tsTriangulation p_comp;
  FILE *p_strm;

  LIBXML_TEST_VERSION; /* init libxml */

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
                                                       /* start root element */
  rc = xmlTextWriterStartElement(writer, BAD_CAST "VTKFile");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type",
                                   BAD_CAST "UnstructuredGrid");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "version",
                                   BAD_CAST "0.1");
  if (rc < 0) {
    ERR("xmlTextWriterWriteAttribute failed\n");
    return -1;
  }

  rc = xmlTextWriterStartElement(writer, BAD_CAST "UnstructuredGrid");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return -1;
  }

  for (i=0; i<nComps; ++i) { /* write-out all components */
    p_comp = p_config + i;

    rc = xmlTextWriterStartElement(writer, BAD_CAST "Piece");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return -1;
    }

    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name",
                                     BAD_CAST p_comp->geomName);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return -1;
    }

    rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "NumberOfPoints",
                                           "%d", p_comp->nVerts);
    if (rc < 0) {
      ERR("xmlTextWriterWriteFormatAttribute failed\n");
      return -1;
    }

    rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "NumberOfCells",
                                           "%d", p_comp->nTris);
    if (rc < 0) {
      ERR("xmlTextWriterWriteFormatAttribute failed\n");
      return -1;
    }
                             /* write basic triangulation for this component */
    rc = trix_writeBasicTri(writer, p_comp, options);
    if (rc < 0) {
      ERR("trix_writeBasicTri failed\n");
      return -1;
    }

    /* Write out component tags used by Intersect, taken from a_Tri[].Comp
     * field to a label named IntersectComponents. If Comp contains BAD_INDEX,
     * this is skipped. Also, if triData contains both GMPtags and
     * IntersectComponents, then this is skipped again because they will be
     * written out below.
     */
    tritags = 0;
    for (j=0; j<p_comp->nTriData; ++j) { /* write tri data */
      if ( 0 == strcmp(p_comp->p_triData[j].name,"IntersectComponents") )
        tritags++;
      if ( 0 == strcmp(p_comp->p_triData[j].name,"GMPtags") )
        tritags++;
    }

    if (2 != tritags) {
      if ( BAD_INDEX != p_comp->a_Tris[0].Comp ) {
        /* start CellData */
        rc = xmlTextWriterStartElement(writer, BAD_CAST "CellData");
        if (rc < 0) {
          ERR("xmlTextWriterStartElement failed\n");
          return -1;
        }
        cellDataOpen = 1;

        rc = trix_writeComps(writer, p_comp->nTris, p_comp->a_Tris);
        if (rc < 0) {
          ERR("trix_writeComps failed\n");
          return -1;
        }
      }
    }

    for (j=0; j<p_comp->nTriData; ++j) { /* write tri data */
      /* skip IntersectComponents, these are taken from .Comp field */
      if ((2 != tritags) &&
          (0 == strcmp(p_comp->p_triData[j].name,"IntersectComponents")) )
        continue;

      if (! cellDataOpen ) {                     /* start CellData */
        rc = xmlTextWriterStartElement(writer, BAD_CAST "CellData");
        if (rc < 0) {
          ERR("xmlTextWriterStartElement failed\n");
          return -1;
        }
        cellDataOpen = 1;
      }

      rc = trix_writeData(writer, p_comp->p_triData[j].name,
                          p_comp->p_triData[j].dim,
                          p_comp->p_triData[j].offset,
                          p_comp->a_scalar0_t,
                          p_comp->nTris,
                          (int) p_comp->p_triData[j].type,
                          (int) p_comp->p_triData[j].info);
      if (rc < 0) {
        ERR("trix_writeTriData failed\n");
        return -1;
      }
    }

    if ( cellDataOpen ) {
      rc = xmlTextWriterEndElement(writer);          /* end CellData element */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
      cellDataOpen = 0;
    }

    if (0 != p_comp->nVertData) {                         /* start PointData */
      rc = xmlTextWriterStartElement(writer, BAD_CAST "PointData");
      if (rc < 0) {
        ERR("xmlTextWriterStartElement failed\n");
        return -1;
      }

      for (j=0; j<p_comp->nVertData; ++j) { /* write vert data */
        rc = trix_writeData(writer, p_comp->p_vertData[j].name,
                            p_comp->p_vertData[j].dim,
                            p_comp->p_vertData[j].offset,
                            p_comp->a_scalar0,
                            p_comp->nVerts,
                            (int) p_comp->p_vertData[j].type,
                            (int) p_comp->p_vertData[j].info);
        if (rc < 0) {
          ERR("writeCompXML failed\n");
          return -1;
        }
      }

      rc = xmlTextWriterEndElement(writer);         /* end PointData element */
      if (rc < 0) {
        ERR("xmlTextWriterEndElement failed\n");
        return -1;
      }
    } /* end if nVertData */

    rc = xmlTextWriterEndElement(writer);              /* end Piece element */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return -1;
    }
  } /* end for comps */

  rc = xmlTextWriterEndElement(writer);      /* end UnstructuredGrid element */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  rc = xmlTextWriterEndElement(writer);              /* end VTKFile element */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return -1;
  }

  rc = xmlTextWriterEndDocument(writer);
  if (rc < 0) {
    ERR("trix_writeConfig: Error at xmlTextWriterEndDocument\n");
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

  /* force to disk: experimental to see if it reduces intermittent failures in
   * executables that follow codes that use c3dio
   */
  {
    int fd, rc;
    fd = fileno(p_strm);
    rc = fsync(fd);
    if ( 0 != rc ) {
      WARN("fsync on mesh io failed\n");
    }
  }

  fclose(p_strm);

  xmlBufferFree(buf);
  xmlCleanupParser(); /* shutdown libxml2 */
  xmlMemoryDump(); /* this is to debug memory for regression tests */

  return 0;
}
#else

/**
 * -- io_writeSurfTrix(): Just a stub for io_writeSurfTrix() to enable
 *                        compilation with no functioning writer
 */
int io_writeSurfTrix(const p_tsTriangulation p_config, const int nComps,
                     const char *const p_fileName)
{
  ERR("LIBXML2 writer not compiled\n");
  return -1;
}

#endif

/**
 * -- io_writeTecplot(): for each comp in config, write a tecplot file
 *                       with triData written as cell-centered
 *                       data. Data in a_Tris[].Comp field is written
 *                       in variable Comp
 */
int io_writeTecplot(FILE *p_strm, const p_tsTriangulation p_surf)
{
  const int        nVerts       = p_surf->nVerts;
  const int        nTris        = p_surf->nTris;
  const p_tsVertex p_V          = p_surf->a_Verts;
  const p_tsTri    p_T          = p_surf->a_Tris;
  const int        nVertData    = p_surf->nVertData;
  const int        nTriData     = p_surf->nTriData;
  const double *const p_sc0     = p_surf->a_scalar0;
  const double *const p_sc0_t   = p_surf->a_scalar0_t;
  const int        nVertScalars = p_surf->nVertScalars;
  const int        nTriScalars  = p_surf->nTriScalars;

  int i,j,k,m,indx,lines,remains,oneOffset;
  char cLine[10*STRING_LEN], labels[STRING_LEN];

  const int numbersPerLine = 5;

  if ( 0 == nVerts || 0 == nTris ) return -1;
  /* header info */
  strcpy(cLine,"VARIABLES = X Y Z");
  /* vert scalars are written out node based, tri scalars are written out as
   * cell-centered data
   */
  if (nVertData) {
    for (i = 0; i < nVertData; ++i) {
      /* if data element has dimension greater than one, append subscript to
       * name
       */
      if ( 1 != p_surf->p_vertData[i].dim ) {
        for (j = 0; j < p_surf->p_vertData[i].dim; ++j) {
          snprintf(labels,STRING_LEN*sizeof(char),"  %s_%d",
                   p_surf->p_vertData[i].name,j);
          strcat(cLine, labels);
        }
      }
      else {
        snprintf(labels,STRING_LEN*sizeof(char),
                 "  %s", p_surf->p_vertData[i].name);
        strcat(cLine, labels);
      }
    }
  }
                   /* first write out component tags and then other tri tags */
  strcat(cLine, "  Comp");
  if (nTriData) {
    for (i = 0; i < nTriData; ++i) {
      /* if data element has dimension greater than one, append subscript to
       * name
       */
      if ( 1 != p_surf->p_triData[i].dim ) {
        for (j = 0; j < p_surf->p_triData[i].dim; ++j) {
          snprintf(labels,STRING_LEN*sizeof(char),
                   "  %s_%d", p_surf->p_triData[i].name,j);
          strcat(cLine, labels);
        }
      }
      else {
        snprintf(labels,STRING_LEN*sizeof(char),
                 "  %s", p_surf->p_triData[i].name);
        strcat(cLine, labels);
      }
    }
  }
  fprintf(p_strm, "%s\n", cLine);

  fprintf(p_strm,"ZONE T=\"%s\", N=%d, E=%d,  DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE\n",
          p_surf->geomName, nVerts, nTris);

  if ( 0 == nTriData ) { /* 3 for vert coords and 1 for Comp field */
    fprintf(p_strm,"VARLOCATION=(%d=CELLCENTERED)\n",3+nVertScalars+1);
  }
  else {
    fprintf(p_strm,"VARLOCATION=([%d-%d]=CELLCENTERED)\n",3+nVertScalars+1,
            3+nVertScalars+nTriScalars+1);
  }

  /* tecplot has a limit on line length so we workaround it */
  lines   = (int) floor((double) nVerts/numbersPerLine);
  remains = nVerts - lines*numbersPerLine;

  k = 0;
  for (i=0; i<lines; ++i) { /* x-coord */
    for (j=0; j<numbersPerLine; ++j) {
      fprintf(p_strm," %16.9e ", p_V[k].x[X]);
      k++;
    }
    fprintf(p_strm,"\n");
  }
  for (i=0; i<remains; ++i) {
    fprintf(p_strm," %16.9e ", p_V[k].x[X]);
    k++;
  }
  fprintf(p_strm,"\n");

  k = 0;
  for (i=0; i<lines; ++i) { /* y-coord */
    for (j=0; j<numbersPerLine; ++j) {
      fprintf(p_strm," %16.9e ", p_V[k].x[Y]);
      k++;
    }
    fprintf(p_strm,"\n");
  }
  for (i=0; i<remains; ++i) {
    fprintf(p_strm," %16.9e ", p_V[k].x[Y]);
    k++;
  }
  fprintf(p_strm,"\n");

  k = 0;
  for (i=0; i<lines; ++i) { /* z-coord */
    for (j=0; j<numbersPerLine; ++j) {
      fprintf(p_strm," %16.9e ", p_V[k].x[Z]);
      k++;
    }
    fprintf(p_strm,"\n");
  }
  for (i=0; i<remains; ++i) {
    fprintf(p_strm," %16.9e ", p_V[k].x[Z]);
    k++;
  }
  fprintf(p_strm,"\n");

  if (nVertData) {
    ASSERT(p_sc0);
    for (k = 0; k < nVertData; ++k) {
      for (j = 0; j < p_surf->p_vertData[k].dim; ++j) {

        for (i=0; i<lines; ++i) {
          for (m=0; m<numbersPerLine; ++m) {
            indx = j + (i*numbersPerLine+m)*(p_surf->p_vertData[k].dim) +
              nVerts*(p_surf->p_vertData[k].offset);
            fprintf(p_strm, " %16.9e", p_sc0[indx]);
          }
          fprintf(p_strm,"\n");
        }
        for (i=0; i<remains; ++i) {
          indx = j + (lines*numbersPerLine+i)*(p_surf->p_vertData[k].dim) +
            nVerts*(p_surf->p_vertData[k].offset);
          fprintf(p_strm, " %16.9e", p_sc0[indx]);
        }
        fprintf(p_strm,"\n");
      }
    }
  }

  lines   = (int) floor((double) nTris/numbersPerLine);
  remains = nTris - lines*numbersPerLine;

  k = oneOffset = 0;
  if (BAD_INDEX != p_T[0].Comp) oneOffset = 1;

  for (i = 0; i < lines; ++i) {
    for (j=0; j<numbersPerLine; ++j) {
      fprintf(p_strm, " %d", p_T[k].Comp+oneOffset); /* one-offset */
      k++;
    }
    fprintf(p_strm,"\n");
  }
  for (i=0; i<remains; ++i) {
    fprintf(p_strm, " %d", p_T[k].Comp+oneOffset);
    k++;
  }
  fprintf(p_strm,"\n");

  if (nTriData) {
    ASSERT(p_sc0);
    for (k = 0; k < nTriData; ++k) {
      for (j = 0; j < p_surf->p_triData[k].dim; ++j) {

        for (i = 0; i < lines; ++i) {
          for (m=0; m<numbersPerLine; ++m) {
            indx = j + (i*numbersPerLine+m)*(p_surf->p_triData[k].dim) +
              nTris*(p_surf->p_triData[k].offset);
            fprintf(p_strm, " %16.9e", p_sc0_t[indx]);
          }
          fprintf(p_strm, "\n");
        }
        for (i=0; i<remains; ++i) {
          indx = j + (lines*numbersPerLine+i)*(p_surf->p_triData[k].dim) +
            nTris*(p_surf->p_triData[k].offset);
          fprintf(p_strm, " %16.9e", p_sc0_t[indx]);
        }
        fprintf(p_strm, "\n");
      }
    }
  }
                                             /* connectivity with one offset */
  for (i=0; i<nTris; ++i) {
    fprintf(p_strm,"%d %d %d\n", p_T[i].vtx[X]+1, p_T[i].vtx[Y]+1,
            p_T[i].vtx[Z]+1);
  }

  return 0;
}

/**
 * -- getTRIXtype()
 */
TRIXtype getTRIXtype(const char *type, const int len)
{
  TRIXtype rv;
  if ( 0==strncmp(type, "SHAPE_LINEARIZATION", len*sizeof(int)) ) {
    rv = TRIX_shapeLinearization;
  }
  else if ( 0==strncmp(type, "FLOW_VARIABLE",  len*sizeof(int)) ) {
    rv = TRIX_flowVariable;
  }
  else if ( 0==strncmp(type, "COMPONENT_TAG",  len*sizeof(int)) ) {
    rv = TRIX_componentTag;
  }
  else if ( 0==strncmp(type, "OTHER",          len*sizeof(int)) ) {
    rv = TRIX_other;
  }
  else {
    rv = TRIX_unset;
  }
  return (rv);
}

/**
 * -- getVTKtype()
 */
VTKtype getVTKtype(const char *type, const int len)
{
  VTKtype rv = VTK_unset;
  if ( 0==strncmp(type, "Int8", len*sizeof(int)) ) {
    rv = VTK_Int8;
  }
  else if ( 0==strncmp(type, "UInt8",   len*sizeof(int)) ) {
    rv = VTK_UInt8;
  }
  else if ( 0==strncmp(type, "Int16",   len*sizeof(int)) ) {
    rv = VTK_Int16;
  }
  else if ( 0==strncmp(type, "UInt16",  len*sizeof(int)) ) {
    rv = VTK_UInt16;
  }
  else if ( 0==strncmp(type, "Int32",   len*sizeof(int)) ) {
    rv = VTK_Int32;
  }
  else if ( 0==strncmp(type, "UInt32",  len*sizeof(int)) ) {
    rv = VTK_UInt32;
  }
  else if ( 0==strncmp(type, "Int64",   len*sizeof(int)) ) {
    rv = VTK_Int64;
  }
  else if ( 0==strncmp(type, "UInt64",  len*sizeof(int)) ) {
    rv = VTK_UInt64;
  }
  else if ( 0==strncmp(type, "Float32", len*sizeof(int)) ) {
    rv = VTK_Float32;
  }
  else if ( 0==strncmp(type, "Float64", len*sizeof(int)) ) {
    rv = VTK_Float64;
  }
  return (rv);
}
