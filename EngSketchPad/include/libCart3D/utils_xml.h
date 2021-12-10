/*
 * Utils for XML
 * Marian Nemec, Nov 2005
 */

/*
 * $Id: utils_xml.h,v 1.1 2008/10/14 13:51:00 aftosmis Exp $
 */

#ifndef __UTILS_XML_H_
#define __UTILS_XML_H_

#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include "memory_util.h"

char*    uxml_copyChar(const xmlChar *source);
void     uxml_deleteChar(const char *source);
xmlChar *uxml_convertInput(const char *in, const char *encoding);

#endif
