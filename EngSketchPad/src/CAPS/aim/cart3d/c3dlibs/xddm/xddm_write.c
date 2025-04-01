
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
 * $Id: xddm_write.c,v 1.3 2022/11/07 23:01:39 mnemec Exp $
 */

/* open source */

/**
 * Write XDDM data structure to file. Recursion is used to write
 * nested elements.
 *
 * Dependency: libxml2
 * https://gitlab.gnome.org/GNOME/libxml2/-/wikis/homewww.xmlsoft.org
 * This library is usually present on most systems, check via
 * existence of 'xml2-config' script.
 */

#include "xddm.h"
#include "xddmInternals.h"

#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

#define MY_ENCODING "ISO-8859-1"

/**
 * Local prototypes
 */
static int xddm_write(xmlTextWriterPtr writer, const p_tsXddm p_xddm,
                      unsigned options);
static int xddm_writeElement(xmlTextWriterPtr writer, const p_tsXddmElem p_node,
                             const unsigned options);
static int xddm_writeVariable(xmlTextWriterPtr writer, const p_tsXddmVar p_v,
                              const unsigned options);
static int xddm_writeAnalysis(xmlTextWriterPtr writer, const p_tsXddmAPar p_a,
                              const unsigned options);
static int xddm_writeAeroFun(xmlTextWriterPtr writer, const p_tsXddmAFun p_af,
                             const unsigned options);
static int xddm_writeIntersect(xmlTextWriterPtr writer, const p_tsXddmGeom p_g,
                               const unsigned options);
static int xddm_writeBasicTypes(xmlTextWriterPtr writer, const p_tsXddm p_xddm,
                                const unsigned flags);
static int xddm_writeFunction(xmlTextWriterPtr writer, const p_tsXddmFun p_f,
                              const unsigned options);
static int xddm_writeSensitivity(xmlTextWriterPtr writer, const size_t ndv,
                                 double *a_lin, char **pa_dvs);
static int xddm_writeSum(xmlTextWriterPtr writer, const p_tsXddmSum p_f,
                         const unsigned options);
static int xddm_writeConfigure(xmlTextWriterPtr writer, const unsigned options);

/**
 * xddm_writeFile(): write xddm to file
 */
int
xddm_writeFile(const char *const p_fileName, const p_tsXddm p_xddm,
               const unsigned options)
{
  xmlTextWriterPtr writer;
  xmlBufferPtr buf;
  int fd, rc=0;
  FILE *p_strm;

  if (!p_xddm) return 1;

  { LIBXML_TEST_VERSION };

  if ( XDDM_VERBOSE & options ) printf(" o Writing \'%s\'\n", p_fileName);

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

  rc = xmlTextWriterSetIndent(writer, 1);                   /* pretty indent */
  if (rc < 0) {
    ERR("xmlTextWriterSetIndent failed\n");
    return 1;
  }

  rc = xddm_write(writer, p_xddm, options);
  if (0 != rc) {
    ERR("xddm_write failed, status %d\n", rc);
    return 1;
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
    WARN("fsync failed\n");
  }
  fclose(p_strm);

  if (buf) xmlBufferFree(buf);
  shutdown_libxml();

  if ( XDDM_VERBOSE & options ) {
    printf("   Wrote \'%s\'\n",p_fileName);
  }

  return rc;
}

/**
 * xddm_write(): write all elements and attributes, recursively
 */
static int
xddm_write(xmlTextWriterPtr writer, const p_tsXddm p_xddm, unsigned options)
{
  size_t i;
  int rc;

  rc = xddm_writeElement(writer, p_xddm->p_e, options); /* write root element */
  if (0 != rc) {
    ERR("xddm_writeElement failed, status %d\n", rc);
    return 1;
  }

  if (XDDM_LIN & options || XDDM_NOLIN & options ) {
                                                  /* write Configure element */
    rc = xddm_writeConfigure(writer, options);
    if (0 != rc) {
      ERR("xddm_writeConfigure failed, status %d\n", rc);
      return 1;
    }

    /* Condifure element should be printed only once, clear bits for
     * recursive calls
     */
    options &= ~(1UL << 2); /* clear linearization, 3rd bit */
    options &= ~(1UL << 3); /* clear no linearization, 4th bit */
  }

  rc = xddm_writeBasicTypes(writer, p_xddm, options);
  if (0 != rc) {
    ERR("xddm_writeBasicTypes failed, status %d\n", rc);
    return 1;
  }

  for (i=0; i<p_xddm->nk; i++) {
    rc = xddm_write( writer, p_xddm->a_kids + i, options);
    if (0 != rc) {
      ERR("xddm_write failed, status %d\n", rc);
      return 1;
    }
  }

  rc = xmlTextWriterEndElement(writer);                  /* end root element */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed here\n");
    return rc;
  }

  return 0;
}

/**
 * xddm_writeConfigure(): write Configure node
 */
static int
xddm_writeConfigure(xmlTextWriterPtr writer, const unsigned options)
{
  int rc=0;

  if (!writer) return 1;

  rc = xmlTextWriterStartElement(writer, BAD_CAST "Configure");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return rc;
  }

  if (XDDM_LIN & options) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity",
                                     BAD_CAST "Required");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }
  else if (XDDM_NOLIN & options) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity",
                                     BAD_CAST "None");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return 0;
}

/**
 * xddm_writeBasicTypes(): write XDDM basic types, i.e. Variables,
 * Constants, Functions, ...
 */
static int
xddm_writeBasicTypes(xmlTextWriterPtr writer, const p_tsXddm p_xddm,
                     const unsigned options)
{
  size_t i;
  int rc;

  if (XDDM_VERBOSE & options)
    printf("   Writing BasicTypes\n");

  for (i=0; i<p_xddm->ng; i++) {
    rc = xddm_writeIntersect(writer, p_xddm->a_geo+i, options);
    if (rc != 0) {
      ERR("xddm_writeIntersect failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->naf; i++) {
    rc = xddm_writeAeroFun(writer, p_xddm->a_afun+i, options);
    if (rc) {
      ERR("xddm_writeAeroFun failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->nv; i++) {
    rc = xddm_writeVariable(writer, &p_xddm->a_v[i], options);
    if (rc != 0) {
      ERR("xddm_writeVariable failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->nc; i++) {
    rc = xddm_writeVariable(writer, &p_xddm->a_c[i], options);
    if (rc != 0) {
      ERR("xddm_writeVariable failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->na; i++) {
    rc = xddm_writeAnalysis(writer, &p_xddm->a_ap[i], options);
    if (rc != 0) {
      ERR("xddm_writeAnalysis failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->nf; i++) {
    rc = xddm_writeFunction(writer, &p_xddm->a_f[i], options);
    if (rc != 0) {
      ERR("xddm_writeFunction failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->ns; i++) {
    rc = xddm_writeSum(writer, &p_xddm->a_s[i], options);
    if (rc != 0) {
      ERR("xddm_writeSum failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->nj; i++) {
    rc = xddm_writeFunction(writer, &p_xddm->a_j[i], options);
    if (rc != 0) {
      ERR("xddm_writeFunction failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->ncon; i++) {
    rc = xddm_writeFunction(writer, &p_xddm->a_con[i], options);
    if (rc != 0) {
      ERR("xddm_writeFunction failed, status %d\n", rc);
      return rc;
    }
  }

  for (i=0; i<p_xddm->nt; i++) {
    rc = xddm_writeElement(writer, &p_xddm->a_t[i], options);
    if (rc != 0) {
      ERR("xddm_writeElement failed, status %d\n", rc);
      return rc;
    }
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return rc;
    }
  }

  return 0;
}

/**
 * xddm_writeElement(): write XDDM element
 */
static int
xddm_writeElement(xmlTextWriterPtr writer, const p_tsXddmElem p_node,
                  const unsigned options)
{
  size_t j;
  int rc=0;

  if (!p_node) return 1;

  if ( NULL == p_node->p_nn ) return 1;

  rc = xmlTextWriterStartElement(writer, BAD_CAST p_node->p_nn);
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return rc;
  }

  if (p_node->p_id) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "ID", BAD_CAST p_node->p_id);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  /* only write the local sensitivity attribute if global options are not set */
  if ( ! (XDDM_LIN & options || XDDM_NOLIN & options) ) {
    if (XDDM_LIN & p_node->flags) {
      rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity",
                                       BAD_CAST "Required");
      if (rc < 0) {
        ERR("xmlTextWriterWriteAttribute failed\n");
        return rc;
      }
    }
    else if (XDDM_NOLIN & p_node->flags) {
      rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Sensitivity",
                                       BAD_CAST "None");
      if (rc < 0) {
        ERR("xmlTextWriterWriteAttribute failed\n");
        return rc;
      }
    }
  }

  if (XDDM_BOUND_UPPER & p_node->flags) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Bound",
                                     BAD_CAST "Upper");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  if (XDDM_BOUND_LOWER & p_node->flags) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Bound",
                                     BAD_CAST "Lower");
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  if (p_node->p_comment) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comment",
                                     BAD_CAST p_node->p_comment);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  for (j=0; j<p_node->nAttr; j++) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST p_node->p_attr[j].p_name,
                                     BAD_CAST p_node->p_attr[j].p_value);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  /* do not end element here in case we need to nest another one */

  return 0;
}

/**
 * xddm_writeDouble(): write number with safeguards
 */
static int
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

/**
 * xddm_writeVariable(): write Variable element
 */
static int
xddm_writeVariable(xmlTextWriterPtr writer, const p_tsXddmVar p_v,
                   const unsigned options)
{
  int rc=0;

  if (!writer) return 0;
  if (!p_v)    return 0;

  rc = xddm_writeElement( writer, p_v->p_e, options);
  if (0 != rc) {
    ERR("xddm_writeElement failed, status %d\n", rc);
    return rc;
  }

  rc = xddm_writeDouble(writer, "Value", p_v->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  rc = xddm_writeDouble(writer, "Min", p_v->minVal);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  rc = xddm_writeDouble(writer, "Max", p_v->maxVal);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  rc = xddm_writeDouble(writer, "TypicalSize", p_v->typicalSize);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  rc = xddm_writeDouble(writer, "FDstep", p_v->fdstep);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  rc = xmlTextWriterEndElement(writer);         /* end variable */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return 0;
}

/**
 * xddm_writeAnalysis(): write Analysis parameter
 */
static int
xddm_writeAnalysis(xmlTextWriterPtr writer, const p_tsXddmAPar p_a,
                   const unsigned options)
{
  int rc=0;

  if (!writer) return 0;
  if (!p_a)    return 0;

  rc = xddm_writeElement( writer, p_a->p_e, options);
  if (0 != rc) {
    ERR("xddm_writeElement failed, status %d\n", rc);
    return rc;
  }

  rc = xddm_writeDouble(writer, "Value", p_a->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  rc = xddm_writeDouble(writer, "DiscretizationError", p_a->derr);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return rc;
  }

  if (p_a->p_afun) {
    rc = xddm_writeAeroFun(writer, p_a->p_afun, options);
    if (rc) {
      ERR("xddm_writeAeroFun failed\n");
      return rc;
    }
  }

  if ( p_a->ndvs > 0 ) {
    rc = xddm_writeSensitivity(writer, p_a->ndvs, p_a->a_lin, p_a->pa_dvs);
    if (rc != 0) {
      ERR("xddm_writeSensitivity failed\n");
      return rc;
    }
  }

  rc = xmlTextWriterEndElement(writer);         /* end analysis */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed Analysis\n");
    return rc;
  }

  return 0;
}

/**
 * xddm_writeSensitivity(): write Sensitivity vector
 */
static int
xddm_writeSensitivity(xmlTextWriterPtr writer, const size_t ndv, double *a_lin,
                      char **pa_dvs)
{
  int rc=0;
  size_t i;

  if (!writer) return 0;
  if (!ndv)    return 0;

  rc = xmlTextWriterStartElement(writer, BAD_CAST "SensitivityArray");
  if (rc < 0) {
    ERR("xmlTextWriterStartElement failed\n");
    return rc;
  }

  for (i=0; i<ndv; i++) {
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Sensitivity");
    if (rc < 0) {
      ERR("xmlTextWriterStartElement failed\n");
      return rc;
    }

    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "P", BAD_CAST pa_dvs[i]);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }

    rc = xddm_writeDouble(writer, "Value", a_lin[i]);
    if (rc) {
      ERR("xddm_writeDouble failed\n");
      return rc;
    }

    rc = xmlTextWriterEndElement(writer);         /* end sensitivity */
    if (rc < 0) {
      ERR("xmlTextWriterEndElement failed\n");
      return rc;
    }
  }

  rc = xmlTextWriterEndElement(writer);         /* end sensitivity array */
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return 0;
}

static int
xddm_writeIntersect(xmlTextWriterPtr writer, const p_tsXddmGeom p_g,
                    const unsigned options)
{
  int rc=0;

  if (!writer) return 0;
  if (!p_g)    return 0;

  rc = xddm_writeElement( writer, p_g->p_e, options);
  if (0 != rc) {
    ERR("xddm_writeElement failed, status %d\n", rc);
    return 1;
  }

  if (p_g->p_parts) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Parts",
                                     BAD_CAST p_g->p_parts);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  if (p_g->p_comp2tri) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Comp2tri",
                                     BAD_CAST p_g->p_comp2tri);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  if (p_g->p_cutout) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Cutout",
                                     BAD_CAST p_g->p_cutout);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  if (p_g->p_overlap) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Overlap",
                                     BAD_CAST p_g->p_overlap);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  if (p_g->p_ps) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "PS",
                                     BAD_CAST p_g->p_ps);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return 0;
}

int
xddm_writeFunction(xmlTextWriterPtr writer, const p_tsXddmFun p_f,
                   const unsigned options)
{
  int rc=0;

  if (!writer) return 0;
  if (!p_f)    return 0;

  rc = xddm_writeElement( writer, p_f->p_e, options);
  if (0 != rc) {
    ERR("xddm_writeElement failed, status %d\n", rc);
    return 1;
  }

  if (p_f->p_expr) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Expr",
                                     BAD_CAST p_f->p_expr);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_f->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  rc = xddm_writeDouble(writer, "Min", p_f->min);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  rc = xddm_writeDouble(writer, "Max", p_f->max);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  if (p_f->ndvs > 0) {
    xddm_writeSensitivity(writer, p_f->ndvs, p_f->a_lin, p_f->pa_dvs);
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return 0;
}

static int
xddm_writeSum(xmlTextWriterPtr writer, const p_tsXddmSum p_f,
              const unsigned options)
{
  int rc=0;

  if (!writer) return 0;
  if (!p_f)    return 0;

  rc = xddm_writeElement( writer, p_f->p_e, options);
  if (0 != rc) {
    ERR("xddm_writeElement failed, status %d\n", rc);
    return 1;
  }

  if (p_f->p_expr) {
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Expr",
                                     BAD_CAST p_f->p_expr);
    if (rc < 0) {
      ERR("xmlTextWriterWriteAttribute failed\n");
      return rc;
    }
  }

  rc = xddm_writeDouble(writer, "Value", p_f->val);
  if (rc) {
    ERR("xddm_writeDouble failed\n");
    return 1;
  }

  if (p_f->ndvs > 0) {
    xddm_writeSensitivity(writer, p_f->ndvs, p_f->a_lin, p_f->pa_dvs);
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return 0;
}

/**
 * xddm_writeAeroFun(): write AeroFun functional
 */
static int
xddm_writeAeroFun(xmlTextWriterPtr writer, const p_tsXddmAFun p_af,
                  const unsigned options)
{
  int rc=0;

  if (!p_af) return 1;

  rc = xddm_writeElement( writer, p_af->p_e, options);
  if (rc) {
    ERR("xddm_writeElement failed\n");
    return 1;
  }

  if (p_af->p_text) {
    rc = xmlTextWriterWriteFormatString(writer, "%s", p_af->p_text);
    if (rc < 0) {
      ERR("xmlTextWriterWriteFormatString failed\n");
      return rc;
    }
  }

  rc = xmlTextWriterEndElement(writer);
  if (rc < 0) {
    ERR("xmlTextWriterEndElement failed\n");
    return rc;
  }

  return rc;
}

static p_tsXddmAPar
xddm_getAnalysis(const p_tsXddmAPar a_ap, const size_t n, const char *const p_xname)
{
  p_tsXddmAPar rv = NULL;
  size_t i;

  if (!a_ap)    return rv;
  if (!p_xname) return rv;

  for (i=0; i<n; i++) {
    p_tsXddmAPar p_a = a_ap + i;
    if ( strncmp(p_xname, p_a->p_e->p_xn, strlen(p_xname)) ) {
      rv = p_a;
      break;
    }
  }

  return rv;
}

static p_tsXddmAPar
xddm_findAnalysis(const p_tsXddm p_xddm, const char *const p_xname)
{
  p_tsXddmAPar rv = NULL;
  size_t i;

  if (!p_xddm)  return rv;
  if (!p_xname) return rv;

  rv = xddm_getAnalysis(p_xddm->a_ap, p_xddm->na, p_xname);
  if (rv) {
    return rv;
  }

  if (p_xddm->a_kids) {
    for (i=0; i<p_xddm->nk; i++) {
      rv = xddm_findAnalysis(p_xddm->a_kids+i, p_xname);
      if (rv) {
        return rv;
      }
    }
  }

  return rv;
}

/**
 * xddm_updateAnalysisParams(): update values and sensitivities of analysis parameters
 */
int
xddm_updateAnalysisParams(const char *const p_fileName, const p_tsXddm p_xddm,
                          const unsigned options)
{
  unsigned libxmlOpts=0;

  xmlDocPtr          doc;
  xmlParserCtxtPtr   ctxt;
  xmlXPathContextPtr xpathCtx = NULL;
  xmlXPathObjectPtr  xddmObj  = NULL;
  FILE *p_f;

  /* initialize and check libxml2 version */
  { LIBXML_TEST_VERSION };

  xmlXPathInit();  /* initialize xpath */

  if ( XDDM_VERBOSE & options ) {
    printf("    o  Parsing file \"%s\" with libxml2\n", p_fileName);
  }

  ctxt = xmlNewParserCtxt(); /* create a document parser context */
  if (ctxt == NULL) {
    ERR("failed to allocate parser context\n");
    xmlCleanupParser();
    return 1;
  }

  doc = xmlCtxtReadFile(ctxt, p_fileName, NULL, libxmlOpts);
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
     ERR("xmlXPathNewContext failed to create xpath context\n");
     if (doc) xmlFreeDoc(doc);
     shutdown_libxml();
     return 1;
  }
                                                   /* get all Analysis nodes */
  xddmObj = xmlXPathEvalExpression(BAD_CAST "//Analysis", xpathCtx);
  if (xddmObj) {
    if ( ! xmlXPathNodeSetIsEmpty(xddmObj->nodesetval) ) {
      const xmlNodeSetPtr nodes = xddmObj->nodesetval;
      size_t nNodes = (nodes) ? nodes->nodeNr : 0;
      size_t i,j;
      p_tsXddmAPar p_a = NULL;

      for (i=0; i < nNodes; i++) {
        if ( nodes->nodeTab[i]->type != XML_ELEMENT_NODE ) continue;
        {
          xmlNodePtr node = nodes->nodeTab[i];
          char   *p_xname = xddm_getName(node, options);
          xmlChar *xmlstr;
          xmlNodePtr sarray_node, new_node;
          xmlAttrPtr new_attr;

          p_a = xddm_findAnalysis(p_xddm, p_xname);

          if (p_a && p_a->val != UNSET) {
            xmlstr = xmlXPathCastNumberToString(p_a->val);
            if (xmlHasProp(node, BAD_CAST "Value")) {
              printf("setting to %g\n",p_a->val);
              xmlSetProp(node, BAD_CAST "Value", xmlstr);
            }
            else {
              xmlNewProp(node, BAD_CAST "Value", xmlstr);
            }
            xmlFree(xmlstr);

            if (p_a->ndvs > 0 && p_a->a_lin) {
              /* should be more careful here, check if it exists, etc., next time */
              sarray_node = xmlNewChild(node, NULL, BAD_CAST "SensitivityArray", NULL);
              for (j=0; j<p_a->ndvs; j++) {
                new_node = xmlNewChild(sarray_node, NULL, BAD_CAST "Sensitivity", NULL);
                new_attr = xmlNewProp(new_node, BAD_CAST "P", BAD_CAST p_a->pa_dvs[j]);
                xmlstr   = xmlXPathCastNumberToString(p_a->a_lin[j]);
                new_attr = xmlNewProp(new_node, BAD_CAST "Value", xmlstr);
                if (NULL==new_attr) {
                  ERR("xmlNewProp failed\n");
                  exit(1);
                }
                xmlFree(xmlstr);
              }
            }
          }
          if (p_xname) free(p_xname);
        }
      }
    } /* if not empty */
    xmlXPathFreeObject(xddmObj);
    xddmObj = NULL;
  }

  p_f = fopen("tester_update_1.xml", "w");
  if (p_f == NULL) {
    printf("fopen failed\n");
    return(-1);
  }

  xmlKeepBlanksDefault(0);
  xmlDocFormatDump(p_f, doc, 1);
  xmlSaveFormatFile("tester_update_2.xml", doc, 1);
  fclose(p_f);

  /* all done, clean up and return
   */
  if (doc) xmlFreeDoc(doc);
  shutdown_libxml();

  return 0;
}

#undef MY_ENCODING

#else

int
xddm_writeFile(const char *const p_fileName, p_tsXddm p_xddm, const unsigned options);
{
  ERR("LIBXML2 writer not compiled\n");
  return;
}

#endif
