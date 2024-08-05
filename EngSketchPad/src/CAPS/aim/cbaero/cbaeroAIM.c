/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             CBAero AIM
 *
 *     Written by Dr. Ryan Durscher AFRL/RQVC
 *
 *
 */

/*!\mainpage Introduction
 * \tableofcontents
 * \section overviewCBAERO CBAero AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with NASA Ames's CBAero \cite CBAero. CBAero (Configuration Based Aerodynamics)
 * software package is an engineering level aero-thermodynamics tool for predicting the aerodynamic and
 * aero-thermodynamic environments of general vehicle configurations. Currently only a subset of CBAero's
 * input options have been exposed in the analysis interface module (AIM), but features can easily be included
 * as future needs arise.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsCBAERO and \ref aimOutputsCBAERO, respectively.
 *
 * Geometric attributes recognized by the AIM are provided in \ref attributeCBAERO.
 *
 * The accepted and expected geometric representation are detailed in \ref geomRepIntentCBAERO.
 *
 */

 /*! \page attributeCBAERO Attribution
 *
 * The following list of attributes drives the CBAero geometric definition.
 *
 *  - <b> capsLength</b> This attribute defines the length units that the *.csm file is generated in.  CBAero grids
 *  MUST be in units of meter, as such the geometry is scaled accordingly based on this value.
 *
 *  - <b> capsReferenceArea</b>  [Optional] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the reference area in CBAero's input file with its units assumed to be consistent with
 * the attribute "capsLength". No conversion takes place if "capsLength" isn't set.  This value may be alternatively set
 * through an input value, "ReferenceArea" (see \ref aimInputsCBAERO)
 *
 *  - <b> capsReferenceChord</b> and <b> capsReferenceSpan</b> [Optional] These attribute may exist on any <em> Body</em>. Their
 * value will be used as the reference moment lengths in CBAero's input file with their units assumed to be consistent with
 * the attribute "capsLength". No conversion takes place if "capsLength" isn't set.  These values may be alternatively set
 * through an input value, "Moment_Length" (see \ref aimInputsCBAERO)
 *
 *  - <b> capsReferenceX</b>, <b> capsReferenceY</b>, and <b>capsReferenceZ</b> [Optional]
 * These attribute may exist on any <em> Body</em>. Their
 * value will be used as the center of gravity (CG) location in CBAero's input file with their units assumed to be consistent with
 * the attribute "capsLength". No conversion takes place if "capsLength" isn't set.  These values may be alternatively set
 * through an input value, "Moment_Center" (see \ref aimInputsCBAERO)
 *
 */

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"
#include "miscUtils.h"
#include "jsonUtils.h"

#include "fastWriter.h"

#include "errno.h"


#ifdef WIN32
#define snprintf    _snprintf
#define strcasecmp  stricmp
#define strncasecmp _strnicmp
#define strtok_r    strtok_s
#else
#include <unistd.h>
#include <limits.h>
#endif

enum aimInputs
{
  inProj_Name = 1,                 /* index is 1-based */
  inMach,
  inDynamic_Pressure,
  inAlpha,
  inBeta,
  inReferenceArea,
  inReferenceChord,
  inReferenceSpan,
  inMoment_Center,
  inFlow_Type,
  inCritical_Transition,
  inPlanet,
  inDefault_Body_Method,
  inDefault_Wing_Method,
  inDefault_Low_Speed_Method,
  inLeading_Edge_Suction,
  inAero_Surface,
  inMaterial_Group,
  inNumParallelCase,
  inNumThreadPerCase,
  inMesh_Morph,
  inSurface_Mesh,
  NUMINPUT = inSurface_Mesh       /* Total number of inputs */
};

enum aimOutputs
{
  outBeta = 1,                    /* index is 1-based */
  outAlpha,
  outDynamic_Pressure,
  outMach,
  outPerTrb,
  outCLtot,
  outCDtot,
  outCMYtot,
  outLoDtot,
  outCL_p,
  outCD_p,
  outCL_v,
  outCD_v,
  outStagnation_Temperature,
  outStagnation_Radius,
  outConvective_Flux,
  outRadiative_Flux,
  outCL_Trefftz,
  outCD_Trefftz,
  NUMOUTPUT = outCD_Trefftz       /* Total number of outputs */
};

#define MXCHAR  255

static char cbaeroInput[] = "cbaeroInput.txt";

//#define DEBUG

typedef struct {

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    // Project name
    char *projectName;

    double Sref, Cref, Bref, Xref, Yref, Zref;

    // Mesh reference obtained from meshing AIM
    aimMeshRef *meshRefIn, meshRefObj;

} aimStorage;


static int initialize_aimStorage(aimStorage *cbaeroInstance) {

  int status;

  // Set initial values for cbaeroInstance
  cbaeroInstance->projectName = NULL;

  // Container for attribute to index map
  status = initiate_mapAttrToIndexStruct(&cbaeroInstance->groupMap);
  if (status != CAPS_SUCCESS) return status;

  cbaeroInstance->meshRefIn = NULL;
  aim_initMeshRef(&cbaeroInstance->meshRefObj, aimUnknownMeshType);

  cbaeroInstance->Sref = 0;
  cbaeroInstance->Cref = 0;
  cbaeroInstance->Bref = 0;
  cbaeroInstance->Xref = 0;
  cbaeroInstance->Yref = 0;
  cbaeroInstance->Zref = 0;

  return status;
}


static int destroy_aimStorage(aimStorage *cbaeroInstance) {

  int status;

  // projectName is just a reference
  cbaeroInstance->projectName = NULL;

  // Attribute to index map
  status = destroy_mapAttrToIndexStruct(&cbaeroInstance->groupMap);

  // Surface Mesh
  aim_freeMeshRef(&cbaeroInstance->meshRefObj);
  cbaeroInstance->meshRefIn = NULL;

  cbaeroInstance->Sref = 0;
  cbaeroInstance->Cref = 0;
  cbaeroInstance->Bref = 0;
  cbaeroInstance->Xref = 0;
  cbaeroInstance->Yref = 0;
  cbaeroInstance->Zref = 0;

  return status;
}


static int cbaero_selectFlow(void *aimInfo, char * string) {
  int value;
  enum cbearoFlowTypeEnum {cbFreeTransition=0, cbLaminar=1, cbTurbulent=2, cbInviscid=3};

  if      (strcasecmp(string, "FreeTransition") == 0 || strcasecmp(string, "0") == 0) value = cbFreeTransition;
  else if (strcasecmp(string, "Laminar")        == 0 || strcasecmp(string, "1") == 0) value = cbLaminar;
  else if (strcasecmp(string, "Turbulent")      == 0 || strcasecmp(string, "2") == 0) value = cbTurbulent;
  else if (strcasecmp(string, "Inviscid")       == 0 || strcasecmp(string, "3") == 0) value = cbInviscid;
  else {
    AIM_ERROR(aimInfo, "Invalid flow type, %s, options: FreeTransition, Laminar, Turbulent, Inviscid!", string);
    value = CAPS_NOTFOUND;
  }
  return value;
}

static int cbaero_selectHighSpeed(void *aimInfo, char * string) {
  int value;
  enum cbearoAeroHighMethodEnum {Base=0, ModifiedNewtonian=3, TangentCone=21, TangentConeNormalShock=22, TangentWedge=31, TangentWedgeNormalShock=32, FreeMolecular=99};

  if      (strcasecmp(string, "ModifiedNewtonian")      == 0) value = ModifiedNewtonian;
  else if (strcasecmp(string, "TangentCone")            == 0) value = TangentCone;
  else if (strcasecmp(string, "TangentConeNormalShock") == 0) value = TangentConeNormalShock;
  else if (strcasecmp(string, "TangentWedge")           == 0) value = TangentWedge;
  else if (strcasecmp(string, "TangentWedgeNormalShock")== 0) value = TangentWedgeNormalShock;
  else if (strcasecmp(string, "FreeMolecular")          == 0) value = FreeMolecular;
  else if (strcasecmp(string, "Base")                   == 0) value = Base;
  else {
    AIM_ERROR(aimInfo, "Invalid high speed method, %s, options: ModifiedNewtonian, TangentCone, TangentConeNormalShock, TangentWedge, TangentWedgeNormalShock, FreeMolecular!\n", string);
    value = CAPS_NOTFOUND;
  }
  return value;
}

static int cbaero_selectLowSpeed(void *aimInfo, char * string) {
  int value;
  enum cbearoAeroLowMethodEnum {FastPanel=1, LowAR=2};

  if      (strcasecmp(string, "FastPanel") == 0) value = FastPanel;
  else if (strcasecmp(string, "LowAR")     == 0) value = LowAR;
  else  {
    AIM_ERROR(aimInfo, "Invalid low speed method, %s, options: FastPanel, TangLowARentCone!", string);
    value = CAPS_NOTFOUND;
  }
  return value;
}

static int cbearo_writeTag(void *aimInfo, const char *projectName,
                           const mapAttrToIndexStruct *groupMap, aimMeshRef *meshRef) {
  int status; // Function status return

  int i, j; // Indexing

  char filename[PATH_MAX];

  int state, nglobal;
  ego body;

  int itri;
  int iface, nFace;
  ego *efaces=NULL;

  int plen, tlen;
  const double *points, *uv;
  const int *ptypes, *pindexs, *tris, *tric;
  const char *groupName = NULL;
  int groupIndex, elementID = 1;

  char tagFolder[] = "TaggedRegions";
  char tagExt[] = ".tag";
  char tagAll[] = ".ALL.taglist";

  FILE *fp = NULL;

  int numTagElement = 0;
  int *tagElement=NULL;
  FILE **tagFiles=NULL;

  printf("Writing CBAero tagged regions - %s\n", tagFolder);

  // Determine number of elements that have a given index
  numTagElement = groupMap->numAttribute;
  AIM_ALLOC(tagElement, numTagElement, int, aimInfo, status);
  for (i = 0; i < numTagElement; i++) tagElement[i] = 0;

  // Figure out how many elements for each tag
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &efaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(efaces, aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      status = retrieve_CAPSGroupAttr(efaces[iface], &groupName);
      AIM_STATUS (aimInfo, status );
      AIM_NOTNULL(groupName, aimInfo, status);
      status = get_mapAttrToIndexIndex(groupMap, groupName, &groupIndex);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
        AIM_ADDLINE(aimInfo, "Available capsGroup attributes:");
        for (j = 0; j < groupMap->numAttribute; j++) {
            AIM_ADDLINE(aimInfo, "\t%s", groupMap->attributeName[j]);
        }
        goto cleanup;
      }

      tagElement[groupIndex-1] += tlen;
    }
    AIM_FREE(efaces);
  }

  // Write tagged files
  // Folder - TaggedRegions
  // File - "ProjectName".ALL.taglist
  //  -- contains a list of all files in folder
  //  --- #numElement projectname.tagName
  /// --- ......
  //  File - "projectname.tagName.tag"
  // -- #numElement
  // -- elementID
  // -- ....

  // TaggedRegions folder
  status = aim_mkDir(aimInfo, tagFolder);
  AIM_STATUS(aimInfo, status);

  // ALL.taglist
#ifdef WIN32
  snprintf(filename, PATH_MAX, "%s\\%s%s", tagFolder, projectName, tagAll);
#else
  snprintf(filename, PATH_MAX, "%s/%s%s", tagFolder, projectName, tagAll);
#endif

  fp = aim_fopen(aimInfo, filename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file: %s", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  for (i = 0; i < groupMap->numAttribute; i++) {
    fprintf(fp, "%d %s.%s", tagElement[i], projectName, groupMap->attributeName[i]);
    if (i-1 != groupMap->numAttribute) fprintf(fp, "\n");
  }

  if (fp != NULL) fclose(fp);
  fp = NULL;

  // Write tag files
  /*@+voidabstract@*/
  AIM_ALLOC(tagFiles, numTagElement, FILE*, aimInfo, status);
  /*@-voidabstract@*/
  for (i = 0; i < numTagElement; i++) tagFiles[i] = NULL;

  for (i = 0; i < numTagElement; i++ ) {

#ifdef WIN32
    snprintf(filename, PATH_MAX, "%s\\%s.%s%s", tagFolder, projectName, groupMap->attributeName[i], tagExt);
#else
    snprintf(filename, PATH_MAX, "%s/%s.%s%s", tagFolder, projectName, groupMap->attributeName[i], tagExt);
#endif

    tagFiles[i] = aim_fopen(aimInfo, filename, "w");
    if (tagFiles[i] == NULL) {
      AIM_ERROR(aimInfo, "Unable to open file: %s", filename);
      status = CAPS_IOERR;
      goto cleanup;
    }
  }

  // Write total number of elements in tag files
  for (i = 0; i < numTagElement; i++ ) {
    fprintf(tagFiles[i], "%d\n", tagElement[i]);
  }

  // Write element indexes to tag files
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &efaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(efaces, aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      status = retrieve_CAPSGroupAttr(efaces[iface], &groupName);
      AIM_STATUS (aimInfo, status );
      AIM_NOTNULL(groupName, aimInfo, status);
      status = get_mapAttrToIndexIndex(groupMap, groupName, &groupIndex);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
        goto cleanup;
      }

      for (itri = 0; itri < tlen; itri++, elementID++)
        fprintf(tagFiles[groupIndex-1], "%d\n", elementID);
    }
    AIM_FREE(efaces);
  }

  status  = CAPS_SUCCESS;

cleanup:

  AIM_FREE(tagElement);
  AIM_FREE(efaces);

  if (tagFiles != NULL) {
    for (i = 0; i < numTagElement; i++) if (tagFiles[i] != NULL) fclose(tagFiles[i]);
    AIM_FREE(tagFiles);
  }

  return status;
}

static void write_InputArray(enum aimInputs varName, FILE *fp, capsValue *aimInputs) {

  int i; // Indexing

  double value;

  fprintf(fp,"%d\n", aimInputs[varName-1].length);
  for (i = 0; i < aimInputs[varName-1].length; i++) {

    if (aimInputs[varName-1].length == 1) {

      value = aimInputs[varName-1].vals.real;

    } else {
      value = aimInputs[varName-1].vals.reals[i];
    }

    if (i+1 == aimInputs[varName-1].length){
      fprintf(fp, "%.18e", value); // No space afterwards

    } else {
      fprintf(fp, "%.18e ", value);

    }
  }

  fprintf(fp, "\n");
}

static int cbaero_writeInput(void *aimInfo, const aimStorage *cbaeroInstance, capsValue *aimInputs)
{

  int status; // Function return status
  char filename[PATH_MAX];

  FILE *fp = NULL;
  char fileExt[] = ".cbaero";

  printf("Writing CBAero input file - %s%s\n", cbaeroInstance->projectName, fileExt);

  snprintf(filename, PATH_MAX, "%s%s", cbaeroInstance->projectName, fileExt);

  fp = aim_fopen(aimInfo, filename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file: %s", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  // Mesh format
  fprintf(fp, "FileType =  fast\n");

  // Reference values
  fprintf(fp, "Sref = %.18e\n", cbaeroInstance->Sref);
  fprintf(fp, "Cref = %.18e\n", cbaeroInstance->Cref);
  fprintf(fp, "Bref = %.18e\n", cbaeroInstance->Bref);

  fprintf(fp, "X_cg = %.18e\n", cbaeroInstance->Xref);
  fprintf(fp, "Y_cg = %.18e\n", cbaeroInstance->Yref);
  fprintf(fp, "Z_cg = %.18e\n", cbaeroInstance->Zref);

  fprintf(fp, "scale = %.18e\n", 1.0);

  status = cbaero_selectFlow(aimInfo, aimInputs[inFlow_Type-1].vals.string);
  if (status < CAPS_SUCCESS) AIM_STATUS(aimInfo, status);
  fprintf(fp, "flotyp = %d\n", status);

  fprintf(fp, "retm_c = %.18e\n", aimInputs[inCritical_Transition-1].vals.real);
  fprintf(fp, "retm_t = %.18e\n", 0.0);
  fprintf(fp, "strm_line_dt = %.18e\n", 0.25);

  fprintf(fp, "Planet = %s\n", aimInputs[inPlanet-1].vals.string);

  fprintf(fp, "Mach Number\n");
  (void) write_InputArray(inMach, fp, aimInputs);

  fprintf(fp, "Dynamic Pressure (Bars)\n");
  (void) write_InputArray(inDynamic_Pressure, fp, aimInputs);

  fprintf(fp, "Angle of Attack\n");
  (void) write_InputArray(inAlpha, fp, aimInputs);

  fprintf(fp, "Angle of SideSlip\n");
  (void) write_InputArray(inBeta, fp, aimInputs);

  //fprintf(fp, "Control Surfaces\n");
  //fprintf(fp, "%d\n", 0);

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);

  return status;
}

static int cbaero_appendSetupControl(FILE *fp, /*@unused@*/ aimMeshRef *meshRef, /*@unused@*/ const mapAttrToIndexStruct *groupMap)
{//, int numTuple, capsTuple controlTuple[]) {

  int status;

  fprintf(fp,"Control Surfaces:\n");
  fprintf(fp,"0\n");

  status = CAPS_SUCCESS;

//cleanup:
  return status;
}

static int cbaero_appendSetupAero(void *aimInfo, FILE *fp, aimMeshRef *meshRef,
                                  const mapAttrToIndexStruct *groupMap,
                                  int numTuple, capsTuple *panelTuple,
                                  char *bodyMethod, char *wingMethod)
{
  int status; // Status return

  int i, j; // Indexing

  int index;
  int value, value2;

  int state, nglobal, nTri = 0;
  ego body;

  int itri;
  int iface, nFace;
  ego *efaces=NULL;

  int plen, tlen;
  const double *points, *uv;
  const int *ptypes, *pindexs, *tris, *tric;
  const char *groupName = NULL;
  char *tupleValue = NULL;
  int groupIndex, elementID = 1;

  enum cbearoAeroSurfaceIndexEnum {Body=1010, Base=1030, Wing=1020, Inlet=9000, Cowl=11000, Nozzle=10000};

  int numBaseTri = 0;
  int *baseTri = NULL;

  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      nTri += tlen;
    }
  }


  AIM_ALLOC(baseTri, nTri, int, aimInfo, status);
  for (i = 0; i < nTri; i++) baseTri[i] = 0;

  fprintf(fp, "Independent Panel Method flags:\n");

  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &efaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(efaces, aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      status = retrieve_CAPSGroupAttr(efaces[iface], &groupName);
      AIM_STATUS (aimInfo, status );
      AIM_NOTNULL(groupName, aimInfo, status);
      status = get_mapAttrToIndexIndex(groupMap, groupName, &groupIndex);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
        goto cleanup;
      }

      index = -1;
      for (j = 0; j < numTuple; j++) {
        status = get_mapAttrToIndexIndex(groupMap, panelTuple[j].name, &index);
        if (status != CAPS_SUCCESS) {
          AIM_ERROR(aimInfo, "Aero_Surface name '%s' not found in capsGroup map!", panelTuple[j].name);
          AIM_ADDLINE(aimInfo, "Available capsGroup attributes:");
          for (j = 0; j < groupMap->numAttribute; j++) {
            AIM_ADDLINE(aimInfo, "\t%s", groupMap->attributeName[j]);
          }
          goto cleanup;
        }
        if (groupIndex == index) {
          tupleValue = string_removeQuotation(panelTuple[j].value);
          break;
        }
      }

      if  (tupleValue == NULL)
        value = cbaero_selectHighSpeed(aimInfo, bodyMethod); // Default to Body
      else {
        if (strcasecmp(tupleValue, "Body")  == 0) value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
        else if (strcasecmp(tupleValue, "Base")  == 0) {
          for (itri = 0; itri < tlen; itri++) {
            baseTri[numBaseTri+itri] = elementID + itri;
          }
          numBaseTri += tlen;
          value = cbaero_selectHighSpeed(aimInfo, "Base");
        }
        else if (strcasecmp(tupleValue, "Wing")  == 0) value = cbaero_selectHighSpeed(aimInfo, wingMethod);
        else if (strcasecmp(tupleValue, "Inlet") == 0) value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
        else if (strcasecmp(tupleValue, "Cowl")  == 0) value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
        else if (strcasecmp(tupleValue, "Nozzle")== 0) value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
        else {
          AIM_ERROR(aimInfo, "Invalid Aero_Surface '%s' value '%s'. Options: Body, Base, Wing, Inlet, Cowl, Nozzle!", panelTuple[j].name, panelTuple[j].value);
          status = CAPS_NOTFOUND;
          goto cleanup;
        }
        AIM_FREE(tupleValue);
      }

      if (value < 0) {
        status = value;
        AIM_STATUS(aimInfo, status);
      }

      for (itri = 0; itri < tlen; itri++) {
        fprintf(fp,"%d\n", value);
      }
      elementID += tlen;
    }
    AIM_FREE(efaces);
  }


  // Groups
  fprintf(fp,"%d\n", groupMap->numAttribute); // Number of groups
  for (i = 0; i < groupMap->numAttribute; i++) {

    index = -1;
    for (j = 0; j < numTuple; j++) {

      status = get_mapAttrToIndexIndex(groupMap, panelTuple[j].name, &index);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "Attribute name '%s' not found in capsGroup map!\n", panelTuple[j].name);
        AIM_ADDLINE(aimInfo, "Available capsGroup attributes:");
        for (j = 0; j < groupMap->numAttribute; j++) {
          AIM_ADDLINE(aimInfo, "\t%s", groupMap->attributeName[j]);
        }
        goto cleanup;
      }
      if (groupMap->attributeIndex[i] == index) {
        tupleValue = string_removeQuotation(panelTuple[j].value);
        break;
      }
    }

    if  (tupleValue == NULL)  {
      value2 = Body;
      value = cbaero_selectHighSpeed(aimInfo, bodyMethod); // Default to Body
    } else {
      AIM_NOTNULL(tupleValue, aimInfo, status);
      if (strcasecmp(tupleValue, "Body")  == 0) {
        value2 = Body;
        value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
      } else if (strcasecmp(tupleValue, "Base")  == 0) {
        value2 = Base;
        value = cbaero_selectHighSpeed(aimInfo, "Base");
      } else if (strcasecmp(tupleValue, "Wing")  == 0) {
        value2 = Wing;
        value = cbaero_selectHighSpeed(aimInfo, wingMethod);
      } else if (strcasecmp(tupleValue, "Inlet") == 0) {
        value2 = Inlet;
        value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
      } else if (strcasecmp(tupleValue, "Cowl")  == 0) {
        value2 = Cowl;
        value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
      } else if (strcasecmp(tupleValue, "Nozzle")== 0) {
        value2 = Nozzle;
        value = cbaero_selectHighSpeed(aimInfo, bodyMethod);
      } else {
        AIM_ERROR(aimInfo, "Invalid Aero_Surface '%s' value '%s'. Options: Body, Base, Wing, Inlet, Cowl, Nozzle!", panelTuple[j].name, panelTuple[j].value);
        status = CAPS_NOTFOUND;
        goto cleanup;
      }
      AIM_FREE(tupleValue);
    }

    if (value < 0) {
      status = value;
      AIM_STATUS(aimInfo, status);
    }

    fprintf(fp, "%d\n", value2); // cbearoAeroSurfaceIndexEnum
    fprintf(fp, "%s\n", groupMap->attributeName[i]);
    fprintf(fp, "%d\n", value);
  }

  // Surface IDs - Which group
  fprintf(fp,"Surface IDs\n");
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &efaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(efaces, aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      status = retrieve_CAPSGroupAttr(efaces[iface], &groupName);
      AIM_STATUS (aimInfo, status );
      AIM_NOTNULL(groupName, aimInfo, status);
      status = get_mapAttrToIndexIndex(groupMap, groupName, &groupIndex);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
        goto cleanup;
      }

      for (itri = 0; itri < tlen; itri++)
        fprintf(fp, "%d\n", groupIndex);
    }
    AIM_FREE(efaces);
  }

  // Base Triangles
  fprintf(fp,"Base Triangles:\n");
  fprintf(fp, "%d\n", numBaseTri);
  for (i = 0; i < numBaseTri; i++) {
    fprintf(fp, "%d\n", baseTri[i]);
  }

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(efaces);
  AIM_FREE(baseTri);
  AIM_FREE(tupleValue);
  return status;
}

static int cbaero_optimization(void *aimInfo, FILE *fp, aimMeshRef *meshRef)
{
  int status; // Status return

  int i; // Indexing

  int state, nglobal, nTri = 0;
  ego body;

  int itri;
  int iface, nFace;

  int plen, tlen;
  const double *points, *uv;
  const int *ptypes, *pindexs, *tris, *tric;

  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      nTri += tlen;
    }
  }

  fprintf(fp,"Optimization Data:\n");
  for (itri = 0; itri < nTri; itri++) {
    fprintf(fp,"0\n");
  }

  fprintf(fp,"250\n");
  fprintf(fp,"1.000000\n");
  for (i = 0; i < 58; i++) {
    fprintf(fp,"0\n");
  }

  status = CAPS_SUCCESS;

cleanup:
  return status;
}

static int cbaero_material(void *aimInfo, FILE *fp, aimMeshRef *meshRef,
                           const mapAttrToIndexStruct *groupMap,
                           int numTuple, capsTuple *materialTuple)
{
  int status; // Status return

  int i, j, k; // Indexing

  char *materialName=NULL;
  int surfaceType;
  double emissivity;

  int state, nglobal;
  ego body;

  int itri;
  int iface, nFace;
  ego *efaces=NULL;

  int plen, tlen, groupIndex, index;
  const double *points, *uv;
  const int *ptypes, *pindexs, *tris, *tric;
  const char *groupName;

  char *keyValue = NULL;
  int *numGroupName=NULL;
  char ***groupNames=NULL;

  const char *materials[17] =
  {"Non Catalytic",
   "Fully Catalytic",
   "RCG",
   "TUFI",
   "ORCC Coated ACC",
   "PCC Coated NEXTEL 440",
   "RLV Design Goal",
   "SIRCA",
   "TABI",
   "Grey C-9 Coated NEXTEL 440",
   "SiC – SiC",
   "C-CAT",
   "LVP Coated ACC",
   "SiC Coated Carbon – Russian",
   "INCONEL 617 PreOxidized",
   "Mars Fully Catalytic, No O2 Recombination",
   "Venus Fully Catalytic, No CO Oxidation"
  };


  if (numTuple == 0) {
    // Default material if non specified
    fprintf(fp,"Material Group Data:\n");
    fprintf(fp,"1\n");
    fprintf(fp,"0 1 1.0\n");
    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(body, NULL, FACE, &nFace,  NULL);
      AIM_STATUS(aimInfo, status);

      for (iface = 0; iface < nFace; iface++) {
        status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        for (itri = 0; itri < tlen; itri++)
          fprintf(fp, "1\n");
      }
    }
    status = CAPS_SUCCESS;
    goto cleanup;
  }

  AIM_ALLOC(numGroupName, numTuple, int, aimInfo, status);
  for (i = 0; i < numTuple; i++) numGroupName[i] = 0;
  AIM_ALLOC(groupNames, numTuple, char**, aimInfo, status);
  for (i = 0; i < numTuple; i++) groupNames[i] = NULL;

  /*!\page cbaeroMaterialGroup CBAero Material Group
   * Structure for the Material_Group tuple  = ("Material Name", "Value").
   * The "Value" must be a JSON String dictionary.
   *
   * If no Material_Group is specified, a default material of "Non Catalytic" with emissity of "1.0" is applied to all capsGroups.
   *
   * \section jsonMaterialGroup Material Group JSON String Dictionary
   *
   * For the JSON string "Value" dictionary
   *  (e.g. "Value" = {"surfaceType": 1, "emissivity": 0.8})
   *  the following keywords ( = default values) may be used:
   */

  fprintf(fp,"Material Group Data:\n");
  fprintf(fp,"%d\n", numTuple);

  for (i = 0; i < numTuple; i++) {

    /*!\page cbaeroMaterialGroup
     * <ul>
     * <li> <B>surfaceType</B> </li> <br>
     *    The surface type must be a string integer, e.g. "0" or "10, <br>
     *    or a case insensitive partial match to one of:
     * - "Non Catalytic" <br>
     * - "Fully Catalytic" <br>
     * - "RCG" <br>
     * - "TUFI" <br>
     * - "ORCC Coated ACC" <br>
     * - "PCC Coated NEXTEL 440" <br>
     * - "RLV Design Goal" <br>
     * - "SIRCA" <br>
     * - "TABI" <br>
     * - "Grey C-9 Coated NEXTEL 440" <br>
     * - "SiC – SiC" <br>
     * - "C-CAT" <br>
     * - "LVP Coated ACC" <br>
     * - "SiC Coated Carbon – Russian" <br>
     * - "INCONEL 617 PreOxidized" <br>
     * - "Mars Fully Catalytic, No O2 Recombination" <br>
     * - "Venus Fully Catalytic, No CO Oxidation" <br>
     * </ul>
     */

    status = json_getString(materialTuple[i].value, "surfaceType", &materialName);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"surfaceType\" in 'Material_Group' input");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    AIM_NOTNULL(materialName, aimInfo, status);

    surfaceType = -1;
    for (j = 0; j < sizeof(materials)/sizeof(char*); j++) {
      if (strcasecmp(materialName, materials[j]) == 0) {
        surfaceType = j;
        break;
      }
    }
    if (surfaceType == -1) {
      if(strspn(materialName, "0123456789") == strlen(materialName))
        surfaceType = atoi(materialTuple[i].name);
    }
    if (surfaceType == -1) {
      AIM_ERROR(aimInfo, "Material_Group '%s' surfaceType '%s' not found!", materialTuple[i].name, materialName);
      AIM_ADDLINE(aimInfo, "Available material names:");
      for (j = 0; j < sizeof(materials)/sizeof(char*); j++) {
        AIM_ADDLINE(aimInfo, "\t%s", materials[j]);
      }
      status = CAPS_NOTFOUND;
      goto cleanup;
    }
    AIM_FREE(materialName);


    /*!\page cbaeroMaterialGroup
     * <ul>
     * <li> <B>emissivity</B> </li> <br>
     *    Emissivity of the material [0 - 1]
     * </ul>
     */
    status = json_getDouble(materialTuple[i].value, "emissivity", &emissivity);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "missing required entry \"emissivity\" in 'Material_Group' input");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    /*! \page cbaeroMaterialGroup
     *
     * <ul>
     *  <li> <B>groupName = "(no default)"</B> </li> <br>
     *  Single or list of <c>capsGroup</c> names on which to apply the material
     *  (e.g. "Name1" or ["Name1","Name2",...]).
     * </ul>
     */
    status = search_jsonDictionary( materialTuple[i].value, "groupName", &keyValue);
    if (status == CAPS_SUCCESS) {
      status = string_toStringDynamicArray(keyValue, &numGroupName[i], &groupNames[i]);
      AIM_FREE(keyValue);
      AIM_STATUS(aimInfo, status);
    } else {
      AIM_ERROR(aimInfo, "missing required entry \"groupName\" in 'Material_Group' input");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    fprintf(fp,"%d 1 %.8f\n", surfaceType, emissivity);
  }

  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &efaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(efaces, aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      status = retrieve_CAPSGroupAttr(efaces[iface], &groupName);
      AIM_STATUS (aimInfo, status );
      AIM_NOTNULL(groupName, aimInfo, status);
      status = get_mapAttrToIndexIndex(groupMap, groupName, &groupIndex);
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
        goto cleanup;
      }

      index = -1;
      for (j = 0; j < numTuple; j++) {
        for (k = 0; k < numGroupName[j]; k++) {
          get_mapAttrToIndexIndex(groupMap, groupNames[j][k], &index);
          if (groupIndex == index) {
            for (itri = 0; itri < tlen; itri++) {
              fprintf(fp,"%d\n", j+1);
            }
            break;
          }
        }
        if (groupIndex == index) break;
      }

      if (j == numTuple) {
        AIM_ERROR(aimInfo, "capsGroup '%s' not found in Material_Group groupName!", groupName);
        AIM_ADDLINE(aimInfo, "Available Material_Group groupName:");
        for (j = 0; j < numTuple; j++) {
          for (k = 0; k < numGroupName[j]; k++) {
            AIM_ADDLINE(aimInfo, "\t%s",  groupNames[j][k]);
          }
        }
        status = CAPS_NOTFOUND;
        goto cleanup;
      }

    }
    AIM_FREE(efaces);

  }

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(materialName);
  AIM_FREE(efaces);
  if (groupNames != NULL)
    for (i = 0; i < numTuple; i++)
      AIM_FREE(groupNames[i]);

  AIM_FREE(numGroupName);
  AIM_FREE(groupNames);

  return status;
}

static int cbaero_structure(void *aimInfo, FILE *fp, aimMeshRef *meshRef)
{
  int status; // Status return

  int i; // Indexing

  int state, nglobal, nTri = 0;
  ego body;

  int itri;
  int iface, nFace;

  int plen, tlen;
  const double *points, *uv;
  const int *ptypes, *pindexs, *tris, *tric;

  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      nTri += tlen;
    }
  }

  fprintf(fp,"Structure Zones Data (v.10):\n");
  fprintf(fp,"0\n");
  for (itri = 0; itri < nTri; itri++) {
    fprintf(fp,"0\n");
  }

  status = CAPS_SUCCESS;

cleanup:
  return status;
}

static int cbaero_writeSetup(void *aimInfo, const aimStorage *cbaeroInstance, capsValue *aimInputs)
{
  int status; // Function return status

  int i; // Indexing

  char filename[PATH_MAX];

  int state, nglobal;
  ego body;

  int itri;
  int iface, nFace;

  int plen, tlen;
  const double *points, *uv;
  const int *ptypes, *pindexs, *tris, *tric;

  FILE *fp = NULL;
  char fileExt[] = ".stp";

  aimMeshRef *meshRef;
  const mapAttrToIndexStruct *groupMap;

  printf("Writing CBAero setup file - %s%s\n", cbaeroInstance->projectName, fileExt);

  meshRef = cbaeroInstance->meshRefIn; // Get pointer to simplify writing
  groupMap = &cbaeroInstance->groupMap;

  snprintf(filename, PATH_MAX, "%s%s", cbaeroInstance->projectName, fileExt);

  fp = aim_fopen(aimInfo, filename, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file: %s", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  status = cbaero_selectHighSpeed(aimInfo, aimInputs[inDefault_Body_Method-1].vals.string);
  if (status < CAPS_SUCCESS) AIM_STATUS(aimInfo, status);
  fprintf(fp, "HypersonicBodyMethod: %d\n", status);

  status = cbaero_selectHighSpeed(aimInfo, aimInputs[inDefault_Wing_Method-1].vals.string);
  if (status < CAPS_SUCCESS) AIM_STATUS(aimInfo, status);
  fprintf(fp, "HypersonicWingMethod: %d\n", status);

  status = cbaero_selectLowSpeed(aimInfo, aimInputs[inDefault_Low_Speed_Method-1].vals.string);
  if (status < CAPS_SUCCESS) AIM_STATUS(aimInfo, status);
  fprintf(fp, "LowSpeedMethod: %d\n", status);

  fprintf(fp, "Leading Edge Suction:  %.6f\n", aimInputs[inLeading_Edge_Suction-1].vals.real);
  fprintf(fp, "Mangler Setting: 2\n");

  status = cbaero_appendSetupAero(aimInfo, fp, meshRef,
                                  groupMap,
                                  aimInputs[inAero_Surface-1].length,
                                  aimInputs[inAero_Surface-1].vals.tuple,
                                  aimInputs[inDefault_Body_Method-1].vals.string,
                                  aimInputs[inDefault_Wing_Method-1].vals.string);
  AIM_STATUS(aimInfo, status);

  // TPS - Needs to be enhanced
  fprintf(fp,"TPS Zones.v.3.4:\n");
  fprintf(fp,"1\n");
  fprintf(fp,"0\n");
  fprintf(fp,"0\n");
  fprintf(fp,"5\n");
  fprintf(fp,"DefaultZone\n");
  fprintf(fp,"StackUpFile\n");
  fprintf(fp,"StructuresStackUpFile\n");
  fprintf(fp,"MarginsFile\n");
  fprintf(fp,"MarginsFile\n");
  fprintf(fp,"None\n");

  // Write element indexes to tag files
  for (i = 0; i < meshRef->nmap; i++) {
    status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nglobal);
    AIM_STATUS(aimInfo, status);

    status = EG_getBodyTopos(body, NULL, FACE, &nFace, NULL);
    AIM_STATUS(aimInfo, status);

    for (iface = 0; iface < nFace; iface++) {
      status = EG_getTessFace(meshRef->maps[i].tess, iface + 1, &plen, &points, &uv, &ptypes, &pindexs,
                              &tlen, &tris, &tric);
      AIM_STATUS(aimInfo, status);

      for (itri = 0; itri < tlen; itri++)
        fprintf(fp,"1\n");
    }
  }

  fprintf(fp,"Trajectories:\n");
  fprintf(fp,"0\n");
  fprintf(fp,"Wake Edges:\n");
  fprintf(fp,"0\n");
  fprintf(fp,"0\n");

  //if (aimInputs[aim_getIndex(aimInfo, "Control_Surface", ANALYSISIN)-1].nullVal ==  NotNull) {
  status = cbaero_appendSetupControl(fp, meshRef, groupMap);//,
  //aimInputs[aim_getIndex(aimInfo, "Control_Surface", ANALYSISIN)-1].length,
  //aimInputs[aim_getIndex(aimInfo, "Control_Surface", ANALYSISIN)-1].vals.tuple);
  AIM_STATUS(aimInfo, status);


  fprintf(fp,"Trajectory Constraint Triangles:\n");
  fprintf(fp,"0\n");

  status = cbaero_optimization(aimInfo, fp, meshRef);
  AIM_STATUS(aimInfo, status);

  fprintf(fp,"Wake Carry Thrus:\n");
  fprintf(fp,"0\n");
  fprintf(fp,"Wake Nodes:\n");
  fprintf(fp,"0\n");
  fprintf(fp,"Wake Tris:\n");
  fprintf(fp,"0\n");

  status = cbaero_material(aimInfo, fp, meshRef,
                           groupMap,
                           aimInputs[inMaterial_Group-1].length,
                           aimInputs[inMaterial_Group-1].vals.tuple);
  AIM_STATUS(aimInfo, status);

  status = cbaero_structure(aimInfo, fp, meshRef);
  AIM_STATUS(aimInfo, status);

  fprintf(fp,"0\n"); // not sure what these are for...
  fprintf(fp,"0\n"); // not sure what these are for...
  fprintf(fp,"Not_Specified\n"); // not sure what these are for...
  fprintf(fp,"0\n"); // not sure what these are for...
  fprintf(fp,"0\n"); // not sure what these are for...
  fprintf(fp,"0.0\n"); // not sure what these are for...

  fprintf(fp,"Propulsion Data\n"); // colon missing on purpose (why cbaero?)
  fprintf(fp,"0\n");

  status = CAPS_SUCCESS;

cleanup:
  if (fp != NULL) fclose(fp);

  return status;
}

/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{

  int status = CAPS_SUCCESS; // Function status return

  aimStorage *cbaeroInstance = NULL;

#ifdef DEBUG
  printf("\n cbaeroAIM/aimInitialize   inst = %d!\n", inst);
#endif

  /* specify the number of analysis input and out "parameters" */
  *nIn     = NUMINPUT;
  *nOut    = NUMOUTPUT;
  if (inst == -1) return CAPS_SUCCESS;

  /*! \page geomRepIntentCBAERO Geometry Representation
   * The geometric representation for the CBAero AIM requires that the body be either a solid
   * body (SOLIDBODY) or a manifold sheet body (SHEETBODY).
   */

  /* specify the field variables this analysis can generate */
  *nFields = 0;
  *franks  = NULL;
  *fnames  = NULL;
  *fInOut  = NULL;

  // Allocate cbaeroInstance
  AIM_ALLOC(cbaeroInstance, 1, aimStorage, aimInfo, status);

  // Set initial values for cbaeroInstance
  status = initialize_aimStorage(cbaeroInstance);
  AIM_STATUS(aimInfo, status);

  *instStore = cbaeroInstance;

  status = CAPS_SUCCESS;

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
  /*! \page aimInputsCBAERO AIM Inputs
   * The following list outlines the CBAero inputs along with their default values available
   * through the AIM interface.
   */
  int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" cbaeroAIM/aimInputs instance = %d  index = %d!\n", iIndex, index);
#endif

  *ainame = NULL;

  // CBAero Inputs
  if (index == inProj_Name) {
    *ainame              = EG_strdup("Proj_Name");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("cbaero_CAPS");
    defval->lfixed       = Change;

    /*! \page aimInputsCBAERO
     * - <B> Proj_Name = "cbaero_CAPS"</B> <br>
     * This corresponds to the project "root" name.
     */
  } else if (index == inMach) {
    *ainame              = EG_strdup("Mach"); // Mach number
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = NULL;
    defval->lfixed       = Change;
    defval->sfixed       = Fixed;
    defval->dim          = 1;
    defval->vals.real    = 0.0;

    /*! \page aimInputsCBAERO
     * - <B> Mach = 0.0 (default) or [0.0, ... , 0.0] </B> <br>
     *  Mach number (can be a single or array of values).
     *
     */
  } else if (index == inDynamic_Pressure) {
    *ainame              = EG_strdup("Dynamic_Pressure"); // Dynamic pressure
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = NULL;
    defval->lfixed       = Change;
    defval->sfixed       = Fixed;
    defval->dim          = 1;
    defval->units        = EG_strdup("bar");
    defval->vals.real    = 0.0;

    /*! \page aimInputsCBAERO
     * - <B> Dynamic_Pressure = 0.0 (default) or [0.0, ... , 0.0] </B> <br>
     *  Dynamic pressure [bar] value (can be a single or array of values).
     *
     */
  } else if (index == inAlpha) {
    *ainame              = EG_strdup("Alpha");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = EG_strdup("degree");
    defval->lfixed       = Change;
    defval->sfixed       = Fixed;
    defval->dim          = 1;
    defval->vals.real    = 0.0;

    /*! \page aimInputsCBAERO
     * - <B> Alpha = 0.0 (default) or [0.0, ... , 0.0] </B> <br>
     *  Angle of attack [degree] (can be a single or array of values).
     */
  } else if (index == inBeta) {
    *ainame              = EG_strdup("Beta");
    defval->type         = Double;
    defval->nullVal      = NotNull;
    defval->units        = EG_strdup("degree");
    defval->lfixed       = Change;
    defval->sfixed       = Fixed;
    defval->dim          = 1;
    defval->vals.real    = 0.0;
    /*! \page aimInputsCBAERO
     * - <B> Beta = 0.0 (default) or [0.0, ... , 0.0] </B> <br>
     *  Sideslip angle  (can be a single or array of values).
     */
  } else if (index == inReferenceArea) {
    *ainame              = EG_strdup("ReferenceArea");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = EG_strdup("meter^2");
    defval->lfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 0.0;

    /*! \page aimInputsCBAERO
     * - <B>ReferenceArea = NULL </B> <br>
     * This sets the reference area for used in force and moment calculations.
     *  Alternatively, the geometry (body) attribute (see \ref attributeCBAERO) "capsReferenceArea" maybe used to specify this variable
     * (note: values set through the AIM input will supersede the attribution value).
     */
  } else if (index == inReferenceChord) {
    *ainame              = EG_strdup("ReferenceChord");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = EG_strdup("meter");
    defval->lfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 0.0;

    /*! \page aimInputsCBAERO
     * - <B>ReferenceChord = NULL </B> <br>
     * This sets the reference chord for used in force and moment calculations.
     *  Alternatively, the geometry (body) attribute (see \ref attributeCBAERO) "capsReferenceChord" maybe used to specify this variable
     * (note: values set through the AIM input will supersede the attribution value).
     */

  } else if (index == inReferenceSpan) {
    *ainame              = EG_strdup("ReferenceSpan");
    defval->type         = Double;
    defval->nullVal      = IsNull;
    defval->units        = EG_strdup("meter");
    defval->lfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 0.0;

    /*! \page aimInputsCBAERO
     * - <B>ReferenceSpan = NULL </B> <br>
     * This sets the reference span for used in force and moment calculations.
     *  Alternatively, the geometry (body) attribute (see \ref attributeCBAERO) "capsReferenceSpan" maybe used to specify this variable
     * (note: values set through the AIM input will supersede the attribution value).
     */

  } else if (index == inMoment_Center) {
    *ainame              = EG_strdup("Moment_Center");
    defval->type          = Double;
    defval->dim           = 1;
    defval->length        = 3;
    defval->nrow          = 3;
    defval->ncol          = 1;
    defval->units         =  EG_strdup("meter");
    AIM_ALLOC(defval->vals.reals, defval->length, double, aimInfo, status);
    defval->vals.reals[0] = 0.0;
    defval->vals.reals[1] = 0.0;
    defval->vals.reals[2] = 0.0;
    defval->nullVal       = IsNull;
    defval->lfixed        = Fixed;

    /*! \page aimInputsCBAERO
     * - <B>Moment_Center = [0.0, 0.0, 0.0] (NULL)</B> <br>
     * Array values correspond to the x, y, and z center of gravity (CG) locations [meter].
     * Alternatively, the geometry (body) attributes (see \ref attributeCBAERO) "capsReferenceX", "capsReferenceY",
     * and "capsReferenceZ" may be used to specify the center of gravity, respectively
     * (note: values set through the AIM input will supersede the attribution values).
     */
  }  else if (index == inFlow_Type) {
    *ainame              = EG_strdup("Flow_Type");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("Inviscid");
    defval->lfixed       = Change;
    /*! \page aimInputsCBAERO
     * - <B>Flow_Type = "Inviscid" </B> <br>
     * Type of flow to consider. Options (=corresponding integer code): FreeTransition(=0), Laminar(=1), Turbulent(=2), Inviscid(=3).
     */
  } else if (index == inCritical_Transition) {
    *ainame              = EG_strdup("Critical_Transition");
    defval->type         = Double;
    defval->nullVal      = NotNull;
    defval->units        = NULL;
    defval->lfixed       = Fixed;
    defval->sfixed       = Fixed;
    defval->dim          = 0;
    defval->vals.real    = 220.0;
    /*! \page aimInputsCBAERO
     * - <B>Critical_Transition = 220.0 </B> <br>
     * Critical ratio of Re-theta (Reynolds based on momentum thickness) and Ma (Mach number) for transition.
     */
  } else if (index == inPlanet) {
    *ainame              = EG_strdup("Planet");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("EARTH");
    defval->lfixed       = Change;

    /*! \page aimInputsCBAERO
     * - <B> Planet = "EARTH"</B> <br>
     * Planet type. Options include “MERCURY”, “VENUS”, “EARTH”, “MARS”, “JUPITER”, “SATURN”,
     * “URANUS”, “NEPTUNE”, and “PLUTO”.
     *
     */
  } else if (index == inDefault_Body_Method) {
    *ainame              = EG_strdup("Default_Body_Method");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("ModifiedNewtonian");
    defval->lfixed       = Change;

    /*! \page aimInputsCBAERO
     * - <B> Default_Body_Method = "ModifiedNewtonian"</B> <br>
     * Default hypersonic base method. Options (=corresponding integer code): ModifiedNewtonian(=3),
     * TangentCone(=21), TangentConeNormalShock(=22), TangentWedge(=31), TangentWedgeNormalShock(=32), FreeMolecular(=99).
     */
  } else if (index == inDefault_Wing_Method) {
    *ainame              = EG_strdup("Default_Wing_Method");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("ModifiedNewtonian");
    defval->lfixed       = Change;

    /*! \page aimInputsCBAERO
     * - <B> Default_Wing_Method = "ModifiedNewtonian"</B> <br>
     * Default hypersonic aerodynamic wing method. Options (=corresponding integer code): ModifiedNewtonian(=3),
     * TangentCone(=21), TangentConeNormalShock(=22), TangentWedge(=31), TangentWedgeNormalShock(=32), FreeMolecular(=99).
     */
  } else if (index == inDefault_Low_Speed_Method) {
    *ainame              = EG_strdup("Default_Low_Speed_Method");
    defval->type         = String;
    defval->nullVal      = NotNull;
    defval->vals.string  = EG_strdup("FastPanel");
    defval->lfixed       = Change;

    /*! \page aimInputsCBAERO
     * - <B> Default_Low_Speed_Method = "FastPanel"</B> <br>
     * Default low speed method. Options (=corresponding integer code): FastPanel(=1), LowAR(=2).
     *
     */
  } else if (index == inLeading_Edge_Suction) {
    *ainame              = EG_strdup("Leading_Edge_Suction");
    defval->type         = Double;
    defval->nullVal      = NotNull;
    defval->units        = NULL;
    defval->lfixed       = Fixed;
    defval->sfixed       = Fixed;
    defval->dim          = Scalar;
    defval->vals.real    = 1.0;
    defval->limits.dlims[0] = -1.0; // Limit of accepted values
    defval->limits.dlims[1] =  1.0;

    /*! \page aimInputsCBAERO
     * - <B> Leading_Edge_Suction = 1.00</B> <br>
     * Default low speed method integer tag. Range [-1.0, 1.0]
     */

  } else if (index == inAero_Surface) {
    *ainame              = EG_strdup("Aero_Surface");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->dim          = Vector;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;

    /*! \page aimInputsCBAERO
     * - <B> Aero_Surface = NULL</B> <br>
     * Defines the type of aero. surface by associating a "capsGroups" attribute name with a particular panel method - ("capsGroup Name", "Value"),
     * where "Value" can either be "Body", "Base", "Wing", "Inlet", "Cowl", or "Nozzle". If a capsGroup panel method is not defined it will be assumed
     * to be a "Body".
     *
     */

  } else if (index == inMaterial_Group) {
    *ainame              = EG_strdup("Material_Group");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->dim          = Vector;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;

    /*! \page aimInputsCBAERO
     * - <B> Material_Group = NULL</B> <br>
     * Defines the type of aero. surface by associating a "capsGroups" attributes with a particular material group - ("Material Name", "Value"),
     * where "Value" must be a JSON String dictionary, see \ref cbaeroMaterialGroup for additional details.
    */

  } else if (index == inNumParallelCase) {
    *ainame              = EG_strdup("NumParallelCase");
    defval->type         = Integer;
    defval->vals.integer = 1;

    /*! \page aimInputsCBAERO
     * - <B> NumParallelCase = 1</B> <br>
     * Set CBAero -mp to define number of Mach, dynamic pressure, and angle of attack cases to solve simultaneously.<br>
     * May be used in conjunction with NumThreadPerCase.
     */

  } else if (index == inNumThreadPerCase) {
    *ainame              = EG_strdup("NumThreadPerCase");
    defval->type         = Integer;
    defval->vals.integer = 1;

    /*! \page aimInputsCBAERO
     * - <B> NumThreadPerCase = 1</B> <br>
     * Set CBAero -omp to define number of threads to solve a each Mach, dynamic pressure, and angle of attack cases.<br>
     * May be used in conjunction with NumParallelCase.
     */

  } else if (index == inMesh_Morph) {
    *ainame              = EG_strdup("Mesh_Morph");
    defval->type         = Boolean;
    defval->lfixed       = Fixed;
    defval->vals.integer = (int) false;
    defval->dim          = Scalar;
    defval->nullVal      = NotNull;

    /*! \page aimInputsCBAERO
     * - <B> Mesh_Morph = False</B> <br>
     * Project previous surface mesh onto new geometry.
     */

  } else if (index == inSurface_Mesh) {
    *ainame             = EG_strdup("Surface_Mesh");
    defval->type        = PointerMesh;
    defval->dim         = Vector;
    defval->lfixed      = Change;
    defval->sfixed      = Change;
    defval->vals.AIMptr = NULL;
    defval->nullVal     = IsNull;
    AIM_STRDUP(defval->units     , "meter"   , aimInfo, status);
    AIM_STRDUP(defval->meshWriter, MESHWRITER, aimInfo, status);

    /*! \page aimInputsCBAERO
     * - <B>Surface_Mesh = NULL</B> <br>
     * A Surface_Mesh link.
     */
  } else {
    AIM_ERROR(aimInfo, "Unknown input index $%d", index);
    status = CAPS_RANGEERR;
    goto cleanup;
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
  // Function return flag
  int status = CAPS_SUCCESS;

  int i;
  int foundSref=(int)false, foundCref=(int)false, foundBref=(int)false;
  int foundXref=(int)false, foundYref=(int)false, foundZref=(int)false;

  // AIM input bodies
  const char *intent;
  int  numBody;
  ego *bodies = NULL;

  // EGADS return values
  int          atype, alen;
  const int    *ints;
  const char   *string;
  const double *reals;

  const char *lengthUnits=NULL;
  double scaleFactor = 1.0;

  aimStorage *cbaeroInstance = (aimStorage *)instStore;

  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Free our meshRef
  (void) aim_freeMeshRef(&cbaeroInstance->meshRefObj);

  if (aimInputs[inSurface_Mesh-1].nullVal == IsNull &&
      aimInputs[inMesh_Morph-1].vals.integer == (int) false) {
    AIM_ANALYSISIN_ERROR(aimInfo, inSurface_Mesh, "'Surface_Mesh' input must be linked to an output 'Surface_Mesh'");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intent, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

#ifdef DEBUG
  printf(" cbaeroAIM/aimUpdateState instance = %d  numBody = %d!\n", iIndex, numBody);
#endif

  if ((numBody <= 0) || (bodies == NULL)) {
    AIM_ERROR(aimInfo, "No body!\n");
    return CAPS_SOURCEERR;
  }

  if (numBody > 1) {
    AIM_ERROR(aimInfo, "CBAero can only accept a single body! numBody = %d\n", numBody);
    return CAPS_SOURCEERR;
  }

  status = aim_capsLength(aimInfo, &lengthUnits);
  AIM_NOTFOUND(aimInfo, status);
  if (status == CAPS_NOTFOUND) {
    AIM_ERROR(aimInfo, "capsLength attribute must be specified for CBAero");
    goto cleanup;
  }
  AIM_NOTNULL(lengthUnits, aimInfo, status);

  status = aim_convert(aimInfo, 1, lengthUnits, &scaleFactor, "meter", &scaleFactor);
  AIM_STATUS(aimInfo, status);

  // Loop over bodies and look for reference quantity attributes
  for (i=0; i < numBody; i++) {
    status = EG_attributeRet(bodies[i], "capsReferenceArea",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      if (atype == ATTRREAL && alen == 1) {
        cbaeroInstance->Sref = reals[0] * scaleFactor * scaleFactor;
        foundSref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceArea should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceChord",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      if (atype == ATTRREAL && alen == 1) {
        cbaeroInstance->Cref = reals[0] * scaleFactor;
        foundCref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceChord should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceSpan",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      if (atype == ATTRREAL && alen == 1) {
        cbaeroInstance->Bref = reals[0] * scaleFactor;
        foundBref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceSpan should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceX",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {

      if (atype == ATTRREAL && alen == 1) {
        cbaeroInstance->Xref = reals[0] * scaleFactor;
        foundXref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceX should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceY",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {

      if (atype == ATTRREAL && alen == 1) {
        cbaeroInstance->Yref = reals[0] * scaleFactor;
        foundYref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceY should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceZ",
                             &atype, &alen, &ints, &reals, &string);
    if (status == EGADS_SUCCESS){

      if (atype == ATTRREAL && alen == 1) {
        cbaeroInstance->Zref = reals[0] * scaleFactor;
        foundZref = (int)true;
      } else {
        AIM_ERROR(aimInfo, "capsReferenceZ should be followed by a single real value!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
      }
    }
  }

  if (aimInputs[inReferenceArea-1].nullVal == NotNull) {
    cbaeroInstance->Sref = aimInputs[inReferenceArea-1].vals.real;
    foundSref = (int)true;
  }
  if (aimInputs[inReferenceChord-1].nullVal == NotNull) {
    cbaeroInstance->Cref = aimInputs[inReferenceChord-1].vals.real;
    foundCref = (int)true;
  }
  if (aimInputs[inReferenceSpan-1].nullVal == NotNull) {
    cbaeroInstance->Bref = aimInputs[inReferenceSpan-1].vals.real;
    foundBref = (int)true;
  }

  if (foundSref == (int)false) {
    AIM_ERROR(aimInfo, "capsReferenceArea is not set on any body and 'ReferenceArea' input not set!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (foundCref == (int)false) {
    AIM_ERROR(aimInfo, "capsReferenceChord is not set on any body and 'ReferenceChord' input not set!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }
  if (foundBref == (int)false) {
    AIM_ERROR(aimInfo, "capsReferenceSpan is not set on any body and 'ReferenceSpan' input not set!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Check for moment reference overwrites
  if (aimInputs[inMoment_Center-1].nullVal == NotNull) {

    cbaeroInstance->Xref = aimInputs[inMoment_Center-1].vals.reals[0];
    cbaeroInstance->Yref = aimInputs[inMoment_Center-1].vals.reals[1];
    cbaeroInstance->Zref = aimInputs[inMoment_Center-1].vals.reals[2];
  } else {
    if (foundXref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceX is not set on any body and 'Moment_Center' input not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    if (foundYref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceY is not set on any body and 'Moment_Center' input not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    if (foundZref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceZ is not set on any body and 'Moment_Center' input not set!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  }

  // Get project name
  cbaeroInstance->projectName = aimInputs[inProj_Name-1].vals.string;

  // Get mesh
  cbaeroInstance->meshRefIn = (aimMeshRef *) aimInputs[inSurface_Mesh-1].vals.AIMptr;

  if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true &&
      cbaeroInstance->meshRefIn == NULL) { // If we are mighty morphing

    // Lets "load" the meshRef now since it's not linked
    status = aim_loadMeshRef(aimInfo, &cbaeroInstance->meshRefObj);
    AIM_STATUS(aimInfo, status);

    // Mighty Morph the mesh
    status = aim_morphMeshUpdate(aimInfo, &cbaeroInstance->meshRefObj, numBody, bodies);
    AIM_STATUS(aimInfo, status);
    /*@-immediatetrans@*/
    cbaeroInstance->meshRefIn = &cbaeroInstance->meshRefObj;
    /*@+immediatetrans@*/
  }
  AIM_NOTNULL(cbaeroInstance->meshRefIn, aimInfo, status);

  // Get attribute to index mapping
  status = create_MeshRefToIndexMap(aimInfo, cbaeroInstance->meshRefIn, &cbaeroInstance->groupMap);
  AIM_STATUS(aimInfo, status);

  status = CAPS_SUCCESS;

cleanup:

  return status;
}

// ********************** AIM Function Break *****************************
int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{
  // Function return flag
  int status = CAPS_SUCCESS;

  int i;

  char filename[PATH_MAX];
  char gridfile[PATH_MAX];
  FILE *fp = NULL;

  const aimStorage *cbaeroInstance = (const aimStorage *)instStore;

  AIM_NOTNULL(aimInputs, aimInfo, status);

  /* symlink the mesh file */
  snprintf(gridfile, PATH_MAX, "%s%s", cbaeroInstance->meshRefIn->fileName, MESHEXTENSION);
  snprintf(filename, PATH_MAX, "%s%s", cbaeroInstance->projectName, ".msh");

  status = aim_symLink(aimInfo, gridfile, filename);
  AIM_STATUS(aimInfo, status);

  for (i = 1; i <= NUMINPUT; i++) {
    if (aim_newAnalysisIn(aimInfo, i) == CAPS_SUCCESS) {
      status = cbaero_writeInput(aimInfo, cbaeroInstance, aimInputs);
      AIM_STATUS(aimInfo, status);
      break;
    }
  }

  // Write Tag files
  status = cbearo_writeTag(aimInfo,
                           cbaeroInstance->projectName,
                           &cbaeroInstance->groupMap,
                           cbaeroInstance->meshRefIn);
  AIM_STATUS(aimInfo, status);

  status = cbaero_writeSetup(aimInfo, cbaeroInstance, aimInputs);
  AIM_STATUS(aimInfo, status);


  // write input file
  fp = aim_fopen(aimInfo, cbaeroInput, "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open %s", cbaeroInput);
    status = CAPS_IOERR;
    goto cleanup;
  }

  // execute cbaero command
  fprintf(fp,
           " -mp %d -omp %d %s",
           aimInputs[inNumParallelCase-1].vals.integer,
           aimInputs[inNumThreadPerCase-1].vals.integer,
           aimInputs[inProj_Name-1].vals.string);

  fclose(fp); fp = NULL;

  status = CAPS_SUCCESS;

cleanup:
  if (fp != NULL) fclose(fp);
  return status;
}

// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ const void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteRefine AIM Execution
   *
   * If auto execution is enabled when creating an CBAero AIM,
   * the AIM will execute refine just-in-time with the command line:
   *
   * \code{.sh}
   * cbaero $(cat cbaeroInput.txt) > cbaeroOutput.txt
   * \endcode
   *
   * where preAnalysis generated the file "cbaeroInput.txt" which contains commandline arguments for cbaero.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an refine AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or refine can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * cbaero.preAnalysis()
   *
   * print ("\n\nRunning......")
   * cbaero.system("cbaero $(cat cbaeroInput.txt) > cbaeroOutput.txt"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * cbaero.postAnalysis()
   * \endcode
   */


  int status = CAPS_SUCCESS;
  char command[PATH_MAX];

  *state = 0;

  // execute cbaero command
  snprintf(command, PATH_MAX,
           "cbaero $(cat %s) > cbaeroOutput.txt",
           cbaeroInput);

  status = aim_system(aimInfo, NULL, command);
  AIM_STATUS(aimInfo, status, "Failed to execute: %s", command);

  status = CAPS_SUCCESS;
cleanup:
  return status;
}

// ********************** AIM Function Break *****************************
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  return CAPS_SUCCESS;
}

// ********************** AIM Function Break *****************************
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{
  /*! \page aimOutputsCBAERO AIM Outputs
   * The following list outlines the CBAero outputs available through the AIM interface. All variables currently
   * correspond to values found in the *.plt file
   */
 int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" cbaeroAIM/aimOutputs instance = %d  index = %d!\n", iIndex, index);
#endif

  form->type   = Double;
  form->lfixed = Change;
  form->sfixed = Fixed;
  form->dim    = 1;
  form->length = 1;
  form->nrow   = 1;
  form->ncol   = 1;
  form->vals.real = 0.0;
  form->vals.reals = NULL;

  if (index == outBeta) {
    *aoname = EG_strdup("Beta");
    form->units = EG_strdup("degree");
  }
  else if (index == outAlpha) {
    *aoname = EG_strdup("Alpha");
    form->units = EG_strdup("degree");
  }
  else if (index == outDynamic_Pressure) {
    *aoname = EG_strdup("Dynamic_Pressure");
    form->units = EG_strdup("bar");
  }
  else if (index == outMach) {
    *aoname = EG_strdup("Mach");
  }
  else if (index == outPerTrb) {
    *aoname = EG_strdup("PerTrb");
  }

  /*! \page aimOutputsCBAERO
   * Reiterate inputs (based on cases):
   * - <B>Beta</B> = Sideslip [degree].
   * - <B>Alpha</B> = Angle of attack [degree].
   * - <B>Dynamic_Pressure</B> = Dynamic pressure [bar].
   * - <B>Mach</B> = Mach number.
   */

  else if (index == outPerTrb) *aoname = EG_strdup("PerTrb");
  /*! \page aimOutputsCBAERO
   * Per-Trb:
   * - <B>PerTrb</B> = PerTrb.
   */

  // Total Forces - Pressure + Viscous
  else if (index == outCLtot) *aoname = EG_strdup("CLtot");
  else if (index == outCDtot) *aoname = EG_strdup("CDtot");
  else if (index == outCMYtot) *aoname = EG_strdup("CMYtot");
  else if (index == outLoDtot) *aoname = EG_strdup("LoDtot");

  /*! \page aimOutputsCBAERO
   * Net Forces - Pressure + Viscous:
   * - <B>CLtot</B> = The lift coefficient.
   * - <B>CDtot</B> = The drag coefficient.
   * - <B>CMYtot</B> = The moment coefficient about the y-axis.
   * - <B>LoDtot</B> = Lift to drag ratio.
   */

  // Pressure Forces
  else if (index == outCL_p) *aoname = EG_strdup("CL_p");
  else if (index == outCD_p) *aoname = EG_strdup("CD_p");

  /*! \page aimOutputsCBAERO
   * Pressure Forces:
   * - <B>CL_p</B> = The lift coefficient - pressure contribution only.
   * - <B>CD_p</B> = The drag coefficient - pressure contribution only.
   */

  // Viscous Forces
  else if (index == outCL_v) *aoname = EG_strdup("CL_v");
  else if (index == outCD_v) *aoname = EG_strdup("CD_v");

  /*! \page aimOutputsCBAERO
   * Viscous Forces:
   * - <B>CL_v</B> = The lift coefficient - viscous contribution only.
   * - <B>CD_v</B> = The drag coefficient - viscous contribution only.
   */

  // Heating
  else if (index == outStagnation_Temperature) {
    *aoname = EG_strdup("Stagnation_Temperature");
    form->units = EG_strdup("kelvin");
  }
  else if (index == outStagnation_Radius) {
    *aoname = EG_strdup("Stagnation_Radius");
    form->units = EG_strdup("meter");
  }
  else if (index == outConvective_Flux) {
    *aoname = EG_strdup("Convective_Flux");
    form->units = EG_strdup("watt per centimeter^2");
  }
  else if (index == outRadiative_Flux) {
    *aoname = EG_strdup("Radiative_Flux");
    form->units = EG_strdup("watt per centimeter^2");
  }

  /*! \page aimOutputsCBAERO
   * Aero-thermal:
   * - <B>Stagnation_Temperature</B> = Stagnation temperature [K].
   * - <B>Stagnation_Radius</B> = = Stagnation radius [m].
   * - <B>Convective_Flux</B> = Convective heat flux [W/cm<sup>2</sup>].
   * - <B>Radiative_Flux</B> = Radiation heat flux [W/cm<sup>2</sup>].
   */

  // Trefftz
  else if (index == outCL_Trefftz) *aoname = EG_strdup("CL_Trefftz");
  else if (index == outCD_Trefftz) *aoname = EG_strdup("CD_Trefftz");

  /*! \page aimOutputsCBAERO
   * Trefftz:
   * - <B>CL_Trefftz</B> = Trefftz lift coefficient.
   * - <B>CD_Trefftz</B> = Trefftz drag coefficient.
   */

  else {

    printf(" cbaeroAIM/aimOutputs index = %d NOT Found!\n", index);
    return CAPS_NOTFOUND;
  }

  AIM_NOTNULL(*aoname, aimInfo, status);

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
// Calculate CBAero output
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, /*@unused@*/ int index,
                  /*@unused@*/ capsValue *val)
{
  int status; // Function status;

  // excess data entries in case future versions of xfoil change
#define MAX_DATA_ENTRY 25
  int numDataEntry = 0;
  double dataLine[MAX_DATA_ENTRY];
  char  headers[MAX_DATA_ENTRY][40];
  const char *valHeader;
  int valIndex;

  int i; // Indexing

  size_t linecap = 0;
  char *line = NULL, *rest = NULL, *token = NULL; // Temporary line holder

  char filename[PATH_MAX]; // File to open
  char fileExt[] = ".plt";

  FILE *fp = NULL; // File pointer

  int numCase = 0;
  double *tempVal = NULL;

  int cbVersion = 6;

  aimStorage *cbaeroInstance = (aimStorage *)instStore;

  #ifdef DEBUG
    printf(" cbaeroAIM/aimCalcOutput instance = %d  index = %d!\n", iIndex, index);
  #endif

  if (val->length > 1) {
    AIM_FREE(val->vals.reals);
  } else {
    val->vals.real = 0.0;
  }

  val->nrow = 1;
  val->ncol = 1;
  val->length = val->nrow*val->ncol;

  // Open cbaero plt file
  snprintf(filename, PATH_MAX, "%s%s", cbaeroInstance->projectName, fileExt);

  fp = aim_fopen(aimInfo, filename, "r");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file: %s\n", filename);
    status = CAPS_IOERR;
    goto cleanup;
  }

  status = getline(&line, &linecap, fp);
  if (status <= 0) AIM_STATUS(aimInfo, (status = CAPS_IOERR));
  AIM_NOTNULL(line, aimInfo, status);

  // check if the BLOCK is first, or the header line is first
  if (strncmp(line, "BLOCK", 5) == 0) {
    // block is first (CBAero version 5), the next line is the header
    cbVersion = 5;

    status = getline(&line, &linecap, fp);
    if (status <= 0) AIM_STATUS(aimInfo, (status = CAPS_IOERR));
  }

  // Parse the header information,
  rest = line;

  numDataEntry = 0;
  while ((token = strtok_r(rest, " ", &rest))) {
    if (token[0] == '\n') continue;
    strcpy(headers[numDataEntry], token);
    //printf("'%s'\n", headers[numDataEntry]);
    numDataEntry++;
    if (numDataEntry > MAX_DATA_ENTRY) {
      AIM_ERROR(aimInfo, "More than %d columns in %s is not expected!",
                MAX_DATA_ENTRY, filename);
      status = CAPS_IOERR;
      goto cleanup;
    }
  }

  // Headers expected in output file that correspond to the AIM output names

  if (index == outBeta)
    valHeader = "Beta";
  else if (index == outAlpha)
    valHeader = "Alpha";
  else if (index == outDynamic_Pressure)
    valHeader = "Bars";
  else if (index == outMach)
    valHeader = "Mach";
  else if (index == outPerTrb)
    valHeader = "PerTrb";
  else if (index == outCLtot)
    valHeader = "CL";
  else if (index == outCDtot)
    valHeader = "CD";
  else if (index == outCMYtot)
    valHeader = "Cm";
  else if (index == outLoDtot)
    valHeader = "LoD";
  else if (index == outCL_p)
    valHeader = "CLp";
  else if (index == outCD_p)
    valHeader = "CDp";
  else if (index == outCL_v)
    valHeader = "CLf";
  else if (index == outCD_v)
    valHeader = "CDf";
  else if (index == outStagnation_Temperature)
    valHeader = "StagTemp";
  else if (index == outConvective_Flux)
    valHeader = "QdotConv";
  else if (index == outRadiative_Flux)
    valHeader = "QdotRad";
  else if (index == outStagnation_Radius)
    valHeader = "StagRadius";
  else if (index == outCL_Trefftz)
    valHeader = "CL_Trefftz";
  else if (index == outCD_Trefftz)
    valHeader = "CD_Trefftz";
  else {
    AIM_ERROR(aimInfo, "Developer error: Unknown variable index %d", index);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  // Find which column contains the requested data
  valIndex = 0;
  while( valIndex < numDataEntry &&
      strncasecmp(headers[valIndex], valHeader, strlen(valHeader)) != 0 ) valIndex++;

  if (valIndex == numDataEntry) {
    //AIM_WARNING(aimInfo, "Could not find '%s' header in %s", valHeader, filename);
    val->nullVal = IsNull;
    status = CAPS_SUCCESS;
    goto cleanup;
  }

  rewind(fp);

  // Scan the file for the BLOCK string
  while (getline(&line, &linecap, fp) >= 0) {

    if (line == NULL) continue;

    if (strncmp(line, "BLOCK", 5) == 0) {

      if (cbVersion == 5) {
        // skip the header
        status = getline(&line, &linecap, fp);
        if (status <= 0) AIM_STATUS(aimInfo, (status = CAPS_IOERR));
      }

      while (getline(&line, &linecap, fp) >= 0) {
        if (strlen(line) == 1) break;
        numCase += 1;
      }
      continue;
    }
  }

  // Did we find any cases?
  if (numCase <= 0) {
    AIM_ERROR(aimInfo, "No BLOCK Case found in CBAero outout!");
    status = CAPS_NOTFOUND;
    goto cleanup;
  }

  AIM_ALLOC(tempVal, numCase, double, aimInfo, status);

  rewind(fp);

  numCase = 0;
  // Scan the file for the data
  while (getline(&line, &linecap, fp) >= 0) {
    AIM_NOTNULL(line, aimInfo, status);

    if (strncmp(line, "BLOCK", 5) == 0) {

      if (cbVersion == 5) {
        // skip the header
        status = getline(&line, &linecap, fp);
        if (status <= 0) AIM_STATUS(aimInfo, (status = CAPS_IOERR));
      }

      while (getline(&line, &linecap, fp) >= 0) {

        if (strlen(line) == 1) break;

        rest = line;
        for (i = 0; i < numDataEntry; i++) {
          token = strtok_r(rest, " ", &rest);
          status = sscanf(token, "%lf", &dataLine[i]);
          if (status <= 0) AIM_STATUS(aimInfo, (status = CAPS_IOERR));
        }

        tempVal[numCase] = dataLine[valIndex];
        numCase += 1;
      }
    }
  }

  // Transfer value(s)
  if (numCase > 1) {

    AIM_ALLOC(val->vals.reals, numCase, double, aimInfo, status);

    for (i = 0; i < numCase; i++) {
      val->vals.reals[i] = tempVal[i];
    }

  } else {

    val->vals.real = tempVal[0];
  }

  val->nrow = numCase;
  val->ncol = 1;
  val->length = val->nrow*val->ncol;

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  if (line != NULL) free(line);

  AIM_FREE(tempVal);

  return status;
}

// ********************** AIM Function Break *****************************
void aimCleanup(void *instStore)
{
  int status; // Returning status

#ifdef DEBUG
  printf(" egadsTessAIM/aimClenup!\n");
#endif

  aimStorage *cbaeroInstance = (aimStorage *)instStore;

  // Clean up cbaeroInstance data
  status  = destroy_aimStorage(cbaeroInstance);
  if (status != CAPS_SUCCESS) printf("Status %d during destroy_aimStorage", status);

  AIM_FREE(cbaeroInstance);
}

