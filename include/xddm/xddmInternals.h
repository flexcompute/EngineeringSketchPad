/*
 * $Id: xddmInternals.h,v 1.2 2014/05/07 23:55:32 mnemec Exp $
 */

/**
 * Internal data-structures and functions of XDDM library
 *
 * Marian Nemec, STC, June 2013
 */

#ifndef __XDDM_INTERNALS_H_
#define __XDDM_INTERNALS_H_

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

struct xddmXmlAttributes {
  int n;
  const xmlChar** a_names;
  xmlChar** a_values;
  xmlNodePtr p_node;
};
typedef struct xddmXmlAttributes tsXddmXmlAttr;
typedef tsXddmXmlAttr *p_tsXddmXmlAttr;

int countNodes(const xmlXPathObjectPtr xpathObj);
p_tsXddmXmlAttr * xddmParseXpathObj(const xmlXPathObjectPtr xpathObj);
p_tsXddmXmlAttr xddmParseNode(const xmlNodePtr p_node);
char * fillXmlString(const xmlChar *p_xml);
char * fillString(const char *p_xml);
double fillDouble(xmlChar * p_xml);
void freeAttributes(p_tsXddmXmlAttr p_xmlAttr);
void freeXddmAttr(p_tsXddmXmlAttr *pp_xmlAttr);
void shutdown_libxml(xmlDocPtr doc);

p_tsXddmVar  allocVariable(const int n);
p_tsXddmAPar allocAnalysis(const int n);
p_tsXmTess   allocTessellate(const int n);

int xddm_writeAttr(xmlTextWriterPtr writer, p_tsXmParent p_node);

int xddm_writeVariable(xmlTextWriterPtr writer, const p_tsXddmVar p_v);
int xddm_writeAnalysis(xmlTextWriterPtr writer, const p_tsXddmAPar p_a);

#endif
