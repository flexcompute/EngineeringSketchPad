
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
 * $Id: main.c,v 1.5 2022/11/07 23:01:39 mnemec Exp $
 */

/* open source */

/**
 * Example client of libxddm. For XDDM documentation, see
 * $CART3D/doc/xddm/xddm.html. The library uses XML Path Language
 * (XPath) to search XDDM documents. For XPath
 * tutorials, see the web, for example:
 * https://www.developer.com/microsoft/dotnet/net-and-xml-xpath-queries/
 *
 * Usage: test_xddm <xddm_filename> <xpath_expression>
 *
 * Dependency: libxml2
 * https://gitlab.gnome.org/GNOME/libxml2/-/wikis/homewww.xmlsoft.org
 * This library is usually present on most systems, check via
 * existence of 'xml2-config' script.
 */

#include <stdio.h>
#include <stdlib.h>

#include "xddm.h"

int main(int argc, char *argv[])
{
  unsigned options = 0;
  p_tsXddm p_xddm = NULL;
  unsigned indent = 3;

  if( argc != 3 ) {
    ERR("Needs 2 arguments: xddm_tester <xddm_filename> <xpath_expression>\n");
    return 1;
  }

  options |= XDDM_VERBOSE;

  p_xddm = xddm_readFile(argv[1], argv[2], &options);
  if (NULL == p_xddm) {
    ERR("xddm_readFile failed to parse\n");
    return 1;
  }

  printf(" o Evaluated \'%s\' for file \'%s\'\n", argv[2], argv[1]);
  xddm_echo(p_xddm, indent);

  xddm_writeFile("xtest_echo.xml", p_xddm, options);

  /* use this to update */
  /* xddm_updateAnalysisParams(argv[1], p_xddm, options); */

  xddm_free(p_xddm, options);

  {
    /* Here we define an optimization problem to demo various XDDM functions
     */
    p_tsXddm p_x, p_y, p_z;
    double x = 1.1;
    double min_x = -1;
    double max_x = 1;

    p_x = xddm_new(1);
    p_y = xddm_new(1);
    p_z = xddm_new(1);

    /* set root node */
    xddm_setID(p_x->p_e, "Optimize", NULL, NULL);
    xddm_setID(p_y->p_e, "InletPressRatioBC", NULL, NULL);
    xddm_setID(p_z->p_e, "DesingPoint", "test", "my comment");

    xddm_addAttribute("Geometry", "Airfoil", &p_z->p_e->nAttr, &p_z->p_e->p_attr);
    xddm_addAttribute("Aero", "aero.csh", &p_z->p_e->nAttr, &p_z->p_e->p_attr);

    p_z->a_kids = p_y;
    p_z->nk = 1;

    p_x->a_kids = p_z;
    p_x->nk = 1;

    /* create variables */
    p_x->nv = 1;
    p_x->a_v = xddm_newVariable(p_x->nv);
    xddm_setID(p_x->a_v[0].p_e, "Variable", "AR", "Aspect Ratio");
    xddm_setVariable(&p_x->a_v[0], &x, NULL, &min_x, &max_x, NULL);

    p_y->nv = 1;
    p_y->a_v = xddm_newVariable(p_y->nv);
    xddm_setID(p_y->a_v[0].p_e, "Variable", "BackPressure", NULL);
    xddm_setVariable(&p_y->a_v[0], &x, NULL, &min_x, &max_x, NULL);

    /* create constants */
    p_x->nc = 1;
    p_x->a_c = xddm_newVariable(p_x->nc);
    xddm_setID(p_x->a_c[0].p_e, "Constant", "Span", "Wing span");
    x = 4;
    xddm_setVariable(&p_x->a_c[0], &x, NULL, NULL, NULL, NULL);

    /* objective function */
    p_x->nj = 1;
    p_x->a_j = xddm_newFunction(p_x->nj);
    xddm_setID(p_x->a_j[0].p_e, "Objective", "Drag", "This is drag objective");
    xddm_setElementFlags(p_x->a_j[0].p_e, XDDM_LIN);
    xddm_setFunction(&p_x->a_j[0], "Cd^2", NULL, NULL, NULL);

    xddm_echo(p_x, indent);
    options |= XDDM_LIN; /* set global Configure */
    xddm_writeFile("xtest_opt.xml", p_x, options);

    xddm_free(p_x, options);
  }

  return 0;
}
