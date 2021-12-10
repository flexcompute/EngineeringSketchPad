/*
 * $Id: main.c,v 1.2 2014/05/12 22:05:21 mnemec Exp $
 */

/**
 * Example client of libxddm. For XDDM documentation, see
 * $CART3D/doc/xddm/xddm.html. The library uses XML Path Language (XPath) to
 * navigate the elements of XDDM documents. For XPath tutorials, see the web,
 * e.g. www.developer.com/net/net/article.php/3383961/NET-and-XML-XPath-Queries.htm
 *
 * Usage: xddm_test <valid_xddm_filename> <xpath_expression>
 *
 * Dependency: libxml2, www.xmlsoft.org. This library is usually present on
 * most systems, check existence of 'xml2-config' script
 *
 * Marian Nemec, STC, June 2013
 */

#include <stdio.h>
#include <stdlib.h>

#include "xddm.h"

int
main(int argc, char *argv[])
{
  int opts = 0;
  p_tsXddm p_xddm = NULL;

  if( argc != 3 ) {
    ERR("Must have two arguments: xddm_test <xddm_filename> <xpath_expression>\n");
    return 1;
  } 

  opts |= 1; /* be verbose */

  p_xddm = xddm_readFile(argv[1], argv[2], opts);
  if (NULL == p_xddm) {
    ERR("xddm_readFile failed to parse\n");
    return 1;
  }

  xddm_echo(p_xddm);

  xddm_writeFile("tester_out.xml", p_xddm, opts);

  xddm_updateAnalysisParams(argv[1], p_xddm, opts);
  
  xddm_free(p_xddm);

  p_xddm = xddm_readFile("tester_out.xml", argv[2], opts);
  if (NULL == p_xddm) {
    ERR("xddm_readFile failed to parse tester_out.xml\n");
    return 1;
  }
  xddm_free(p_xddm);

  return 0;
}
