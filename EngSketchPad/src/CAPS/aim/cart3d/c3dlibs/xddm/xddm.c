
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
 * $Id: xddm.c,v 1.7 2022/11/07 23:01:39 mnemec Exp $
 */

/* open source */

/**
 * Constructors, destructors and utility functions to manipulate
 * Extensible-Design-Description-Markup (XDDM) elements. This library
 * provides functions to read and write XDDM files, and to handle XDDM
 * data structures.
 *
 * Dependency: libxml2
 * https://gitlab.gnome.org/GNOME/libxml2/-/wikis/homewww.xmlsoft.org
 * This library is usually present on most systems, check via
 * existence of 'xml2-config' script.
 */

#include "xddm.h"
#include "xddmInternals.h"

/**
 * xddm_new(): constructor of XDDM data structure with basic initialization
 */
p_tsXddm
xddm_new(const size_t n)
{
  size_t i;
  p_tsXddm p_xddm = NULL;

  p_xddm = malloc(n * sizeof(*p_xddm));
  if (NULL == p_xddm) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_tsXddm p_x = p_xddm + i;
    p_x->p_e    = xddm_newElement(1);
    p_x->nv     = 0;
    p_x->nc     = 0;
    p_x->na     = 0;
    p_x->naf    = 0;
    p_x->nt     = 0;
    p_x->nf     = 0;
    p_x->ns     = 0;
    p_x->nj     = 0;
    p_x->ncon   = 0;
    p_x->ng     = 0;
    p_x->a_v    = NULL;
    p_x->a_c    = NULL;
    p_x->a_ap   = NULL;
    p_x->a_afun = NULL;
    p_x->a_t    = NULL;
    p_x->a_f    = NULL;
    p_x->a_s    = NULL;
    p_x->a_j    = NULL;
    p_x->a_con  = NULL;
    p_x->a_geo  = NULL;

    p_x->a_kids = NULL;
    p_x->nk     = 0;
  }

  return p_xddm;
}

/**
 * xddm_newAttribute(): constructor for a name-value pairs
 */
p_tsXddmAttr
xddm_newAttribute(const size_t n)
{
  size_t i;
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

/**
 * xddm_newElement(): constructor for XDDM elements
 */
p_tsXddmElem
xddm_newElement(const size_t n)
{
  size_t i;
  p_tsXddmElem p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_nn      = NULL;
    p_v[i].p_id      = NULL;
    p_v[i].p_xn      = NULL;
    p_v[i].p_comment = NULL;
    p_v[i].p_attr    = NULL;
    p_v[i].flags     = 0;
    p_v[i].nAttr     = 0;
  }

  return p_v;
}

/**
 * xddm_setID(): set node name, ID and Comment attributes
 */
void
xddm_setID(const p_tsXddmElem p_e, const char *const p_nn,
           const char *const p_id, const char *const p_comment)
{
  if (! p_e) return;

  if (p_nn) {
    if (!p_e->p_nn) {
      p_e->p_nn = fillString(p_nn);
    }
    else {
      ERR("cannot overwrite name %s with %s\n", p_e->p_nn, p_id);
    }
  }

  if (p_id) {
    if (!p_e->p_id) {
      p_e->p_id = fillString(p_id);
    }
    else {
      ERR("cannot overwrite ID %s with %s\n", p_e->p_id, p_id);
    }
  }

  if (p_comment) {
    if (!p_e->p_comment) {
      p_e->p_comment = fillString(p_comment);
    }
    else {
      ERR("cannot overwrite comment %s with %s\n", p_e->p_comment, p_comment);
    }
  }
}

void
xddm_setElementXN(const p_tsXddmElem p_e, const char *const p_xn)
{
  if (! p_e) return;
  if (p_xn && !p_e->p_xn) p_e->p_xn = fillString(p_xn);
}

void
xddm_setElementFlags(const p_tsXddmElem p_e, const xddmOptions flag)
{
  if (! p_e) return;
  p_e->flags |= flag;
}

/**
 * xddm_newVariable(): constructor for XDDM variables and constants
 */
p_tsXddmVar
xddm_newVariable(const size_t n)
{
  size_t i;
  p_tsXddmVar p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_e         = xddm_newElement(1);
    p_v[i].val         = UNSET;
    p_v[i].minVal      = UNSET;
    p_v[i].maxVal      = UNSET;
    p_v[i].typicalSize = UNSET;
    p_v[i].fdstep      = UNSET;
  }

  return p_v;
}

/**
 * xddm_setVariable(): sets xddm values of design variables and constants
 */
void
xddm_setVariable(const p_tsXddmVar p_v, const double *const p_val,
                 const double *const p_typicalSize,
                 const double *const p_minVal,
                 const double *const p_maxVal,
                 const double *const p_fdstep)
{
  if (! p_v) return;
  if (p_val)         p_v->val         = *p_val;
  if (p_typicalSize) p_v->typicalSize = *p_typicalSize;
  if (p_minVal)      p_v->minVal      = *p_minVal;
  if (p_maxVal)      p_v->maxVal      = *p_maxVal;
  if (p_fdstep)      p_v->fdstep      = *p_fdstep;
}

p_tsXddmAFun
xddm_newAeroFun(const size_t n)
{
  size_t i;
  p_tsXddmAFun p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_e    = xddm_newElement(1);
    p_v[i].p_text = NULL;
  }

  return p_v;
}

p_tsXddmGeom
xddm_newIntersect(const size_t n)
{
  size_t i;
  p_tsXddmGeom p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_e        = xddm_newElement(1);
    p_v[i].p_parts    = NULL;
    p_v[i].p_comp2tri = NULL;
    p_v[i].p_cutout   = NULL;
    p_v[i].p_overlap  = NULL;
    p_v[i].p_ps       = NULL;
  }

  return p_v;
}

p_tsXddmAPar
xddm_newAnalysis(const size_t n)
{
  size_t i;
  p_tsXddmAPar p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_e       = xddm_newElement(1);
    p_v[i].p_afun    = NULL;
    p_v[i].val       = UNSET;
    p_v[i].derr      = UNSET;
    p_v[i].ndvs      = 0;
    p_v[i].a_lin     = NULL;
    p_v[i].pa_dvs    = NULL;
  }

  return p_v;
}

/**
 * xddm_newFunction(): constructor for functions, objectives and constraints
 */
p_tsXddmFun
xddm_newFunction(const size_t n)
{
  size_t i;
  p_tsXddmFun p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_e       = xddm_newElement(1);
    p_v[i].p_expr    = NULL;
    p_v[i].val       = UNSET;
    p_v[i].min       = UNSET;
    p_v[i].max       = UNSET;
    p_v[i].ndvs      = 0;
    p_v[i].a_lin     = NULL;
    p_v[i].pa_dvs    = NULL;
  }

  return p_v;
}

/**
 * xddm_setFunction(): sets expression, values and limits for
 * functions, objectives and constraints
 */
void
xddm_setFunction(p_tsXddmFun p_f, const char *const p_expr,
                 const double *const p_val, const double *const p_min,
                 const double *const p_max)
{
  if (!p_f) return;
  if ( p_expr && !p_f->p_expr) p_f->p_expr = fillString(p_expr);
  if ( p_val ) p_f->val = *p_val;
  if ( p_min ) p_f->min = *p_min;
  if ( p_max ) p_f->max = *p_max;
}

p_tsXddmSum
xddm_newSum(const size_t n)
{
  size_t i;
  p_tsXddmSum p_v=NULL;

  p_v = malloc(n * sizeof(*p_v));
  if (NULL == p_v) {
    ERR("malloc failed\n");
    exit(1);
  }

  for (i=0; i<n; i++) {
    p_v[i].p_e       = xddm_newElement(1);
    p_v[i].p_expr    = NULL;
    p_v[i].val       = UNSET;
    p_v[i].ndvs      = 0;
    p_v[i].a_lin     = NULL;
    p_v[i].pa_dvs    = NULL;
  }

  return p_v;
}

/**
 * Prepends t into s. Assumes s has enough space allocated for the
 * combined string.
 */
void
xddm_prepend(char* s, const char* t)
{
  size_t len = strlen(t);
  size_t i;

  memmove(s + len, s, strlen(s) + 1);

  for (i = 0; i < len; ++i) {
    s[i] = t[i];
  }
}

/**
 * xddm_getName(): Construct full name of xddm node, for example
 * Optimize__Variable__Alpha. This is a concatenation of node names
 * and IDs up to document root.
 */
char *
xddm_getName(xmlNodePtr p_node, const unsigned options)
{
  xmlAttr   *p_attribute;
  xmlChar   *p_xml;
  char      *name=NULL;
  size_t     len = 0;
  xmlNodePtr p_orig;

  if (!p_node) return name;

  p_orig = p_node;
  /* determine name length */
  while (p_node && p_node->name) {
    p_attribute = p_node->properties;
    while (p_attribute && p_attribute->name && p_attribute->children) {
      if ( 0 == xmlStrcasecmp( p_attribute->name, BAD_CAST "id" ) ) {
        p_xml = xmlGetProp(p_node, p_attribute->name);
        if (len) len+=2;
        len += strlen( (char *) p_xml );
        free(p_xml);
      }
      p_attribute = p_attribute->next;
    }
    if (len) len+=2;
    len += strlen( (char *)p_node->name );
    p_node = p_node->parent;
  }

  len += 1; /* add terminating char */

  if ( len > MAX_STR_LEN ) {
    ERR("huge string\n");
    exit(1);
  }

  name = malloc( len * sizeof(*name) );
  if ( NULL == name ) {
    ERR("malloc failed\n");
    exit(1);
  }
  memset((void*) name,'\0',len);

  p_node = p_orig;

  while (p_node && p_node->name) {
    p_attribute = p_node->properties;
    while (p_attribute && p_attribute->name && p_attribute->children) {
      if ( 0 == xmlStrcasecmp( p_attribute->name, BAD_CAST "id" ) ) {
        p_xml = xmlGetProp(p_node, p_attribute->name);
        if (strlen(name)) xddm_prepend(name, "__");
        xddm_prepend(name, (char *) p_xml);
        free(p_xml);
      }
      p_attribute = p_attribute->next;
    }
    if (strlen(name) > 0) xddm_prepend(name, "__");
    xddm_prepend(name, (char *) p_node->name);
    p_node = p_node->parent;
  }

  if (XDDM_DEBUG & options) {
    printf("xddm_getName: name %s len %zu\n", name, len);
  }

  return name;
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

void
xddm_addAttribute(const char * p_name, const char * p_value,
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
  (*pp_attr)[(*p_nAttr)-1].p_name  = fillString(p_name);
  (*pp_attr)[(*p_nAttr)-1].p_value = fillString(p_value);

  // printf("%zu %s %s\n",*p_nAttr,(*pp_attr)[(*p_nAttr)-1].p_name,(*pp_attr)[(*p_nAttr)-1].p_value);
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
    pp_xmlAttr = NULL;
  }

  return;
}

static void
echoAttributes(const size_t n, const p_tsXddmAttr p_at, unsigned indent)
{
  unsigned i;
  indent += 3;
  if ( n < 1 ) return;
  printf("%*cCustom Attributes\n",indent,' ');
  for (i=0; i<n; i++) {
    printf("%*c%2d. '%s' '%s'\n", indent, ' ', i, p_at[i].p_name,
           p_at[i].p_value);
  }
}

static void
echoElement(const p_tsXddmElem p_e, const size_t index, unsigned indent)
{
  if (!p_e) return;

  printf("%*c%zu. Node=%s\n", indent, ' ', index, p_e->p_nn);

  indent += 3;

  if (p_e->p_id) printf("%*cID=%s\n",    indent, ' ', p_e->p_id);
  if (p_e->p_xn) printf("%*cXDDM Name=%s\n", indent, ' ', p_e->p_xn);

  if ( XDDM_LIN   & p_e->flags ) printf("%*cLinearization Required\n", indent, ' ');
  if ( XDDM_NOLIN & p_e->flags ) printf("%*cLinearization Not Required\n", indent, ' ');

  if (p_e->p_comment) printf("%*c%s\n", indent, ' ', p_e->p_comment);

  echoAttributes(p_e->nAttr, p_e->p_attr, indent);
}

static void
echoIntersect(const p_tsXddmGeom p_g, const size_t index, unsigned indent)
{
  if (!p_g) return;

  echoElement(p_g->p_e, index, indent);

  indent += 3;

  if (p_g->p_parts)    printf("%*cParts=\"%s\"\n", indent, ' ', p_g->p_parts);
  if (p_g->p_comp2tri) printf("%*cComp2tri=\"%s\"\n", indent, ' ', p_g->p_comp2tri);
  if (p_g->p_cutout)   printf("%*cCutout=\"%s\"\n", indent, ' ', p_g->p_cutout);
  if (p_g->p_overlap)  printf("%*cOverlap=\"%s\"\n", indent, ' ', p_g->p_overlap);
  if (p_g->p_ps)       printf("%*cPS=\"%s\"\n", indent, ' ', p_g->p_ps);
}

static void
echoVariable(const p_tsXddmVar p_v, const size_t index, unsigned indent)
{
  if (!p_v) return;

  echoElement(p_v->p_e, index, indent);

  indent += 3;

  if (UNSET != p_v->val)    printf("%*cValue=%g\n", indent, ' ', p_v->val);
  if (UNSET != p_v->minVal) printf("%*cMin=%g", indent, ' ', p_v->minVal);
  if (UNSET != p_v->maxVal) printf(" Max=%g", p_v->maxVal);
  if (UNSET != p_v->typicalSize) printf(" TypicalSize=%g",p_v->typicalSize);
  if (UNSET != p_v->fdstep) printf(" FDstep=%g",p_v->fdstep);
  if (UNSET != p_v->minVal || UNSET != p_v->maxVal ||
      UNSET != p_v->typicalSize || UNSET != p_v->fdstep) printf("\n");
}

static void
echoAeroFun(const p_tsXddmAFun p_af, const size_t index, unsigned indent)
{
  if (!p_af) return;

  echoElement(p_af->p_e, index, indent);

  if (p_af->p_text) {
    printf("%*c--- AeroFun TEXT ---\n", indent, ' ');
    printf("%s\n", p_af->p_text);
    printf("%*c--- END AeroFun ---\n", indent, ' ');
  }
}

static void
echoSensitivity(char **const names, const double *const values, const int n,
                unsigned indent)
{
  int i;
  printf("%*cSensitivity array\n", indent,' ');
  indent += 2;
  for (i=0; i<n; i++) {
    printf("%*cDV=%s Value=%g\n", indent, ' ',names[i], values[i]);
  }
}

static void
echoAnalysis(const p_tsXddmAPar p_a, const size_t index, unsigned indent)
{
  if (!p_a) return;

  echoElement(p_a->p_e, index, indent);

  indent += 3;

  if (UNSET != p_a->val)  printf("%*cValue=%g\n", indent,  ' ', p_a->val);
  if (UNSET != p_a->derr) printf("%*cDiscretization Error=%g\n", indent,  ' ',
                                 p_a->derr);

  if ( p_a->p_afun )   echoAeroFun(p_a->p_afun, 0, indent);
  if ( p_a->ndvs > 0 ) echoSensitivity(p_a->pa_dvs, p_a->a_lin, p_a->ndvs,
                                       indent);
}

static void
echoFunction(const p_tsXddmFun p_f, size_t index, unsigned indent)
{
  echoElement(p_f->p_e, index, indent);

  indent += 3;

  if (p_f->p_expr)       printf("%*cExpr=%s\n", indent, ' ', p_f->p_expr);
  if (UNSET != p_f->val) printf("%*cValue=%g\n", indent,  ' ', p_f->val);
  if (UNSET != p_f->min) printf("%*cMin=%g\n", indent, ' ', p_f->min);
  if (UNSET != p_f->max) printf("%*cMax=%g\n", indent, ' ', p_f->max);

  if ( XDDM_BOUND_UPPER & p_f->p_e->flags ) printf("%*cUpper Bound\n", indent, ' ');
  if ( XDDM_BOUND_LOWER & p_f->p_e->flags ) printf("%*cLower Bound\n", indent, ' ');

  if ( p_f->ndvs > 0 ) echoSensitivity(p_f->pa_dvs, p_f->a_lin, p_f->ndvs,
                                       indent);
}

static void
echoSum(const p_tsXddmSum p_s, size_t index, unsigned indent)
{
  echoElement(p_s->p_e, index, indent);

  indent += 3;

  if (p_s->p_expr)       printf("%*cExpr=%s\n", indent, ' ', p_s->p_expr);
  if (UNSET != p_s->val) printf("%*cValue=%g\n", indent,  ' ', p_s->val);

  if ( p_s->ndvs > 0 ) echoSensitivity(p_s->pa_dvs, p_s->a_lin, p_s->ndvs,
                                       indent);
}

static void
echoBasicTypes(const p_tsXddm p_xddm, int indent)
{
  size_t i;

  if (p_xddm->ng > 0) {
    printf("%*co Number of Geometry nodes = %zu\n",indent,' ',p_xddm->ng);
    for (i=0; i<p_xddm->ng; i++) {
      echoIntersect(&p_xddm->a_geo[i], i, indent+3);
    }
  }

  if (p_xddm->nv > 0) {
    printf("%*co Number of variables = %zu\n",indent,' ',p_xddm->nv);
    for (i=0; i<p_xddm->nv; i++) {
      echoVariable(&p_xddm->a_v[i], i, indent+3);
    }
  }

  if (p_xddm->nc > 0) {
    printf("%*co Number of constants = %zu\n",indent,' ',p_xddm->nc);
    for (i=0; i<p_xddm->nc; i++) {
      echoVariable(&p_xddm->a_c[i], i, indent+3);
    }
  }

  if (p_xddm->na > 0) {
    printf("%*co Number of analysis parameters = %zu\n",indent,' ',p_xddm->na);
    for (i=0; i<p_xddm->na; i++) {
      echoAnalysis(&p_xddm->a_ap[i], i, indent+3);
    }
  }

  if (p_xddm->naf > 0) {
    printf("%*co Number of AeroFun elememnts = %zu\n",indent,' ',p_xddm->naf);
    for (i=0; i<p_xddm->naf; i++) {
      echoAeroFun(&p_xddm->a_afun[i], i, indent+3);
    }
  }

  if (p_xddm->nf > 0) {
    printf("%*co Number of functions = %zu\n",indent,' ',p_xddm->nf);
    for (i=0; i<p_xddm->nf; i++) {
      echoFunction(&p_xddm->a_f[i], i, indent+3);
    }
  }

  if (p_xddm->ns > 0) {
    printf("%*co Number of sums = %zu\n",indent,' ',p_xddm->ns);
    for (i=0; i<p_xddm->ns; i++) {
      echoSum(&p_xddm->a_s[i], i, indent+3);
    }
  }

  if (p_xddm->nj > 0) {
    printf("%*co Number of objectives = %zu\n",indent,' ',p_xddm->nj);
    for (i=0; i<p_xddm->nj; i++) {
      echoFunction(&p_xddm->a_j[i], i, indent+3);
    }
  }

  if (p_xddm->ncon > 0) {
    printf("%*co Number of constraints = %zu\n",indent,' ',p_xddm->ncon);
    for (i=0; i<p_xddm->ncon; i++) {
      echoFunction(&p_xddm->a_con[i], i, indent+3);
    }
  }

  if (p_xddm->nt > 0) {
    printf("%*co Number of tessellate elements = %zu\n",indent,' ',p_xddm->nt);
    for (i=0; i<p_xddm->nt; i++) {
      echoElement(&p_xddm->a_t[i], i, indent+2);
    }
  }
}

void
xddm_echo(const p_tsXddm p_xddm, unsigned indent)
{
  size_t i;

  echoElement(p_xddm->p_e, 0, indent);

  echoBasicTypes(p_xddm, indent);

  if (p_xddm->nk > 0) {
    printf("%*co Number of XDDM kids = %zu\n",indent,' ',p_xddm->nk);
    for (i=0; i<p_xddm->nk; i++) {
      xddm_echo(&p_xddm->a_kids[i], indent+3);
    }
  }

  return;
}

/**
 * destructors for all datatypes
 */

void
xddm_freeElement(p_tsXddmElem p_e)
{
  size_t i;
  if (p_e->p_nn)      free(p_e->p_nn);
  if (p_e->p_id)      free(p_e->p_id);
  if (p_e->p_xn)      free(p_e->p_xn);
  if (p_e->p_comment) free(p_e->p_comment);
  for (i=0; i<p_e->nAttr; i++) {
    if (p_e->p_attr[i].p_name ) free(p_e->p_attr[i].p_name);
    if (p_e->p_attr[i].p_value) free(p_e->p_attr[i].p_value);
  }
  if (p_e->p_attr) free(p_e->p_attr);

  p_e->flags = 0;
  p_e->nAttr = 0;

  return;
}

void
xddm_freeAeroFun(p_tsXddmAFun p_f)
{
  if (p_f->p_e) {
    xddm_freeElement(p_f->p_e);
    free(p_f->p_e);
    p_f->p_e = NULL;
  }
  if (p_f->p_text) free(p_f->p_text);
  p_f->p_text = NULL;
  return;
}

void
xddm_freeFunction(p_tsXddmFun p_f)
{
  size_t i;

  if (p_f->p_e) {
    xddm_freeElement(p_f->p_e);
    free(p_f->p_e);
    p_f->p_e = NULL;
  }

  if (p_f->p_expr) free(p_f->p_expr);
  /* derivative arrays */
  if (p_f->a_lin)  free(p_f->a_lin);
  for (i=0; i<p_f->ndvs; i++) {
    if (p_f->pa_dvs[i]) free(p_f->pa_dvs[i]);
  }
  if (p_f->pa_dvs) free(p_f->pa_dvs);
  p_f->ndvs = 0;

  return;
}

void
xddm_freeBasicTypes(p_tsXddm p_x)
{
  size_t i,j;

  if (!p_x) return;

  if (p_x->a_v) {
    for (i=0; i<p_x->nv; i++) {
      if (p_x->a_v[i].p_e) {
        xddm_freeElement(p_x->a_v[i].p_e);
        free(p_x->a_v[i].p_e);
        p_x->a_v[i].p_e = NULL;
      }
    }
    free(p_x->a_v);
    p_x->a_v = NULL;
    p_x->nv = 0;
  }

  if (p_x->a_c) {
    for (i=0; i<p_x->nc; i++) {
      if (p_x->a_c[i].p_e) {
        xddm_freeElement(p_x->a_c[i].p_e);
        free(p_x->a_c[i].p_e);
        p_x->a_c[i].p_e = NULL;
      }
    }
    free(p_x->a_c);
    p_x->a_c = NULL;
    p_x->nc = 0;
  }

  if (p_x->a_ap) {
    for (i=0; i<p_x->na; i++) {
      p_tsXddmAPar p_a = &p_x->a_ap[i];

      if (p_a->p_e) {
        xddm_freeElement(p_a->p_e);
        free(p_a->p_e);
        p_a->p_e = NULL;
      }

      if (p_a->a_lin) free(p_a->a_lin);
      for (j=0; j<p_a->ndvs; j++) {
        if (p_a->pa_dvs[j]) free(p_a->pa_dvs[j]);
      }
      if (p_a->pa_dvs) free(p_a->pa_dvs);
      p_a->ndvs = 0;

      if (p_a->p_afun) {
        xddm_freeAeroFun(p_a->p_afun);
        free(p_a->p_afun);
        p_a->p_afun = NULL;
      }
    }
    free(p_x->a_ap);
    p_x->a_ap = NULL;
    p_x->na = 0;
  }

  if (p_x->a_afun) {
    for (i=0; i<p_x->naf; i++) xddm_freeAeroFun(p_x->a_afun+i);
    free(p_x->a_afun);
    p_x->a_afun = NULL;
    p_x->naf = 0;
  }

  if (p_x->a_t) {
    for (i=0; i<p_x->nt; i++) xddm_freeElement(&p_x->a_t[i]);
    free(p_x->a_t);
    p_x->a_t = NULL;
    p_x->nt = 0;
  }

  if (p_x->a_f) {
    for (i=0; i<p_x->nf; i++) xddm_freeFunction(p_x->a_f+i);
    free(p_x->a_f);
    p_x->a_f = NULL;
    p_x->nf = 0;
  }

  if (p_x->a_j) {
    for (i=0; i<p_x->nj; i++) xddm_freeFunction(p_x->a_j+i);
    free(p_x->a_j);
    p_x->a_j = NULL;
    p_x->nj = 0;
  }

  if (p_x->a_con) {
    for (i=0; i<p_x->nc; i++) xddm_freeFunction(p_x->a_con+i);
    free(p_x->a_con);
    p_x->a_con = NULL;
    p_x->ncon = 0;
  }

  if (p_x->a_s) {
    for (i=0; i<p_x->ns; i++) {
      p_tsXddmSum p_s = &p_x->a_s[i];
      if (p_s->p_e) {
        xddm_freeElement(p_s->p_e);
        free(p_s->p_e);
        p_s->p_e = NULL;
      }
      if (p_s->a_lin) free(p_s->a_lin);
      for (j=0; j<p_s->ndvs; j++) {
        if (p_s->pa_dvs[j]) free(p_s->pa_dvs[j]);
      }
      if (p_s->pa_dvs) free(p_s->pa_dvs);
      p_s->ndvs = 0;
    }
    free(p_x->a_s);
    p_x->a_s = NULL;
    p_x->ns = 0;
  }

  if (p_x->a_geo) {
    for (i=0; i<p_x->ng; i++) {
      p_tsXddmGeom p_g = &p_x->a_geo[i];
      xddm_freeElement(p_g->p_e);
      free(p_g->p_e);
      if (p_g->p_parts)    free(p_g->p_parts);
      if (p_g->p_comp2tri) free(p_g->p_comp2tri);
      if (p_g->p_cutout)   free(p_g->p_cutout);
      if (p_g->p_overlap)  free(p_g->p_overlap);
      if (p_g->p_ps)       free(p_g->p_ps);
    }
    free(p_x->a_geo);
    p_x->a_geo = NULL;
    p_x->ng = 0;
  }

  return;
}

/**
 * xddm_free(): destructor of XDDM datastructure
 */
void
xddm_free(p_tsXddm p_xddm, unsigned flags)
{
  size_t i;

  if (!p_xddm) return;

  if ( XDDM_VERBOSE  & flags )
    printf(" o Freeing XDDM\n");

  if (p_xddm->a_kids) {
    for (i=0; i<p_xddm->nk; i++) {
      xddm_free(p_xddm->a_kids+i, flags);
    }
    p_xddm->nk = 0;
    p_xddm->a_kids = NULL;
  }

  if (p_xddm->p_e) {
    xddm_freeElement(p_xddm->p_e);
    free(p_xddm->p_e);
    p_xddm->p_e = NULL;
  }

  xddm_freeBasicTypes(p_xddm);

  free(p_xddm);
  p_xddm = NULL;

  if ( XDDM_VERBOSE  & flags )
    printf("   done\n");

  return;
}

/**
 * shutdown_libxml(): done with libxml2
 */
void
shutdown_libxml()
{
  xmlCleanupParser(); /* shutdown libxml */
  xmlMemoryDump();    /* this is to debug memory for regression tests */
}
