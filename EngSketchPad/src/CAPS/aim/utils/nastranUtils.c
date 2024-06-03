// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.


#ifdef HAVE_PYTHON
#include "Python.h" // Bring in Python API
#endif


#include <string.h>
#include <math.h>

#include "aimUtil.h"

#include "feaTypes.h"  // Bring in FEA structures
#include "capsTypes.h" // Bring in CAPS types

#include "feaUtils.h" // Bring in FEA utility functions
#include "vlmUtils.h" // Bring in VLM utility functions
#include "miscUtils.h" //Bring in misc. utility functions
#include "nastranUtils.h" // Bring in nastran utility header
#include "nastranCards.h" // Bring in nastran cards

#include "cardUtils.h"

#ifdef HAVE_PYTHON
#include "nastranOP2Reader.h" // Bring in Cython generated header file

#ifndef CYTHON_PEP489_MULTI_PHASE_INIT
#define CYTHON_PEP489_MULTI_PHASE_INIT (PY_VERSION_HEX >= 0x03050000)
#endif

#if CYTHON_PEP489_MULTI_PHASE_INIT
static int nastranOP2Reader_Initialized = (int)false;
#endif

#endif

#ifdef WIN32
#define strcasecmp  stricmp
#else
#include <unistd.h>
#endif



#define PI        3.1415926535897931159979635



static int _getDesignVariableIDSet(const feaProblemStruct *feaProblem,
                                   int numDesignVariableNames,
                                   char **designVariableNames,
                                   int *numDesignVariableID,
                                   int **designVariableIDSet) {
    int i, status;

    int numDesignVariables;
    int *designVariableIDs;
    feaDesignVariableStruct **designVariables = NULL;

    if (numDesignVariableNames == 0) {
        *numDesignVariableID = 0;
        *designVariableIDSet = NULL;
        return CAPS_SUCCESS;
    }

    status = fea_findDesignVariablesByNames(
        feaProblem,
        numDesignVariableNames, designVariableNames,
        &numDesignVariables, &designVariables
    );

    if (status == CAPS_NOTFOUND) {
        PRINT_WARNING("Only %d of %d design variables found",
                      numDesignVariables, numDesignVariableNames);
    }
    else if (status != CAPS_SUCCESS) goto cleanup;

    designVariableIDs = EG_alloc(sizeof(int) * numDesignVariables);

    for (i = 0; i < numDesignVariables; i++) {

        designVariableIDs[i] = designVariables[i]->designVariableID;
    }

    *numDesignVariableID = numDesignVariables;
    *designVariableIDSet = designVariableIDs;

    cleanup:

        if (designVariables != NULL) {
            EG_free(designVariables);
        }
        return status;
}

static int _getDesignResponseIDSet(const feaProblemStruct *feaProblem,
                                   int numDesignResponseNames,
                                   char **designResponseNames,
                                   int *numDesignResponseID,
                                   int **designResponseIDSet) {
    int i, status;

    int numDesignResponses;
    int *designResponseIDs;
    feaDesignResponseStruct **designResponses = NULL;

    if (numDesignResponseNames == 0) {
        *numDesignResponseID = 0;
        *designResponseIDSet = NULL;
        return CAPS_SUCCESS;
    }

    status = fea_findDesignResponsesByNames(
        feaProblem,
        numDesignResponseNames, designResponseNames,
        &numDesignResponses, &designResponses
    );

    if (status == CAPS_NOTFOUND) {
        PRINT_WARNING("Only %d of %d design responses found",
                      numDesignResponses, numDesignResponseNames);
    }
    else if (status != CAPS_SUCCESS) goto cleanup;

    designResponseIDs = EG_alloc(sizeof(int) * numDesignResponses);

    for (i = 0; i < numDesignResponses; i++) {

        designResponseIDs[i] = 100000 + designResponses[i]->responseID;
    }

    *numDesignResponseID = numDesignResponses;
    *designResponseIDSet = designResponseIDs;

    cleanup:

        if (designResponses != NULL) {
            EG_free(designResponses);
        }
        return status;
}


static int _getEquationResponseIDSet(const feaProblemStruct *feaProblem,
                                     int numEquationResponseNames,
                                     char *const*equationResponseNames,
                                     int *numEquationResponseID,
                                     int **equationResponseIDSet) {
    int i, status;

    int numEquationResponses;
    int *equationResponseIDs;
    feaDesignEquationResponseStruct **equationResponses = NULL;

    if (numEquationResponseNames == 0) {
        *numEquationResponseID = 0;
        *equationResponseIDSet = NULL;
        return CAPS_SUCCESS;
    }
    status = fea_findEquationResponsesByNames(
        feaProblem,
        numEquationResponseNames, equationResponseNames,
        &numEquationResponses, &equationResponses
    );
    if (status == CAPS_NOTFOUND) {
        PRINT_WARNING("Only %d of %d design equation responses found",
                      numEquationResponses, numEquationResponseNames);
    }
    else if (status != CAPS_SUCCESS) goto cleanup;

    equationResponseIDs = EG_alloc(sizeof(int) * numEquationResponses);

    for (i = 0; i < numEquationResponses; i++) {

        equationResponseIDs[i] = 200000 + equationResponses[i]->equationResponseID;
    }

    *numEquationResponseID = numEquationResponses;
    *equationResponseIDSet = equationResponseIDs;

    cleanup:

        if (equationResponses != NULL) {
            EG_free(equationResponses);
        }
        return status;
}

static int _getEquationID(const feaProblemStruct *feaProblem, char *equationName, int *equationID) {

    int status;

    feaDesignEquationStruct *equation;

    status = fea_findEquationByName(feaProblem, equationName, &equation);
    if (status != CAPS_SUCCESS) return status;

    *equationID = equation->equationID;

    return status;
}


// Write SET case control card
int nastran_writeSetCard(FILE *fp, int n, int numSetID, int *setID) {

    int i, maxCharPerID = 10; // 8 per field, 2 for command and space

    int bufferLength = 0, addLength, lineLength = 0;
    char *buffer = NULL, *continuation = "\n\t        ";

    if (numSetID == 0) {
        PRINT_ERROR("Empty case control set, n = %d", n);
    }
    else if (numSetID == 1) {
        fprintf(fp, "\tSET %d = %d\n", n, setID[0]);
    }
    else {

        buffer = EG_alloc(sizeof(char) * (maxCharPerID * numSetID
                                          + 100 * strlen(continuation) // enough for 100 continuations
                                          + 1));
        if (buffer == NULL) {
            return EGADS_MALLOC;
        }

        for (i = 0; i < numSetID-1; i++) {

            // dry run, do we pass the 72-char limit ?
            addLength = snprintf(NULL, 0, "%d, ", setID[i]);

            if (lineLength + addLength >= 72 ) {
                addLength = snprintf(buffer + bufferLength, strlen(continuation) + 1, "%s", continuation);
                lineLength = addLength - 1; // -1 because dont count newline
                bufferLength += addLength;
            }

            addLength = snprintf(buffer + bufferLength, maxCharPerID, "%d, ", setID[i]);
            lineLength += addLength;
            bufferLength += addLength;
        }

        snprintf(buffer + bufferLength, maxCharPerID, "%d", setID[numSetID-1]);

        fprintf(fp, "\tSET %d = %s\n", n, buffer);

        if (buffer != NULL) EG_free(buffer);
    }

    return CAPS_SUCCESS;
}


// Write a Nastran element cards not supported by mesh_writeNastran in meshUtils.c
// Elements are ordered by property ID with the property name label above
int nastran_writeSubElementCard(void *aimInfo,
                                FILE *fp,
                                const meshStruct *feaMesh,
                                int numProperty,
                                const feaPropertyStruct *feaProperty,
                                const feaFileFormatStruct *feaFileFormat) {

    int status;

    int i, j; // Indexing
    int iprop;

    int *mcid;
    double *theta, zoff;

    int gID;
    meshNodeStruct *nodeA, *nodeB;
    double vectorX[3];

    feaMeshDataStruct *feaData;

    if (numProperty > 0 && feaProperty == NULL) return CAPS_NULLVALUE;
    if (feaMesh == NULL) return CAPS_NULLVALUE;

    if (feaMesh->meshType == VolumeMesh) return CAPS_SUCCESS;

    // loop over elements in order of property ID
    for (iprop = 0; iprop < numProperty; iprop++) {

        fprintf(fp, "$ Femap Property : %s\n", feaProperty[iprop].name);

        // Write all elements for the current property ID
        for (i = 0; i < feaMesh->numElement; i++) {

            if (feaMesh->element[i].analysisType != MeshStructure) continue;

            feaData = (feaMeshDataStruct *) feaMesh->element[i].analysisData;

            // only consider the current property
            if (feaData->propertyID != feaProperty[iprop].propertyID) continue;

            if (feaData->coordID != 0){
                mcid = &feaData->coordID;
                theta = NULL;
            } else {
                mcid = NULL;
                theta = NULL;
            }

            zoff = feaProperty[iprop].membraneThickness * feaProperty[iprop].zOffsetRel / 100.0;

            if (feaMesh->element[i].elementType == Node &&
                feaData->elementSubType == ConcentratedMassElement) {

//                if (found == (int) false) {
//                    printf("No property information found for element %d of type \"ConcentratedMass\"!", feaMesh->element[i].elementID);
//                    continue;
//                }

                status = nastranCard_conm2(
                    fp,
                    &feaMesh->element[i].elementID, // eid
                    feaMesh->element[i].connectivity, // g
                    &feaData->coordID, // cid
                    &feaProperty[iprop].mass, // m
                    feaProperty[iprop].massOffset, // x
                    feaProperty[iprop].massInertia, // i
                    feaFileFormat->gridFileType);
                AIM_STATUS(aimInfo, status);
            }

            if (feaMesh->element[i].elementType == Line) {

                if (feaData->elementSubType == BarElement) {

//                    if (found == (int) false) {
//                        printf("No property information found for element %d of type \"Bar\"!", feaMesh->element[i].elementID);
//                        continue;
//                    }

                    gID = feaMesh->element[i].connectivity[0];
                    nodeA = &feaMesh->node[gID];

                    gID = feaMesh->element[i].connectivity[1];
                    nodeB = &feaMesh->node[gID];

                    double vec[3];

                    vec[0] = nodeB->xyz[0] - nodeA->xyz[0];
                    vec[1] = nodeB->xyz[1] - nodeA->xyz[1];
                    vec[2] = nodeB->xyz[2] - nodeA->xyz[2];

                    double flat_vec[3] = {vec[0], vec[1], 0.0};
                    double vec_dot;
                    double flat_norm;
                    double vec_proj[3];

                    const double pi = 4.0 * atan(1.0);
                    double ang = -pi/2;
                    double rot_z[9] = {cos(ang), -sin(ang), 0.,
                                       sin(ang), cos(ang), 0.,
                                       1., 0., 0.};

                    vec_dot = dot_DoubleVal(vec, flat_vec);
                    flat_norm = sqrt(pow(flat_vec[0], 2) + pow(flat_vec[1], 2) + pow(flat_vec[2], 2));

                    for (j = 0; j < 3; j++) {
                        vec_proj[j] = (vec_dot/flat_norm) * (flat_vec[j]/flat_norm);
                        vectorX[j] = rot_z[3*j]*vec_proj[j] + rot_z[(3*j)+1]*vec_proj[j] + rot_z[(3*j)+2]*vec_proj[j];
                    }

                    status = nastranCard_cbar(
                        fp,
                        &feaMesh->element[i].elementID, // eid
                        &feaData->propertyID, // pid
                        feaMesh->element[i].connectivity, // g
                        vectorX,
    //                    feaProperty[j].orientationVec, // x
                        NULL, // g0
                        NULL, // pa
                        NULL, // pb
                        NULL, // wa
                        NULL, // wb
                        feaFileFormat->gridFileType);
                    AIM_STATUS(aimInfo, status);
                }

                if (feaData->elementSubType == BeamElement) {
                    printf("Beam elements not supported yet - Sorry !\n");
                    return CAPS_NOTIMPLEMENT;
                }
            }

            if ( feaMesh->element[i].elementType == Triangle) {

                if (feaData->elementSubType == ShellElement) {

//                    if (found == (int) false) {
//                        printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
//                        continue;
//                    }

                    status = nastranCard_ctria3(
                        fp,
                        &feaMesh->element[i].elementID, // eid
                        &feaData->propertyID, //pid
                        feaMesh->element[i].connectivity, // g
                        theta, // theta
                        mcid, // mcid
                        &zoff, // zoffs
                        NULL, // t
                        feaFileFormat->gridFileType);
                    AIM_STATUS(aimInfo, status);
                }
            }


            if ( feaMesh->element[i].elementType == Triangle_6) {

                if (feaData->elementSubType == ShellElement) {

//                    if (found == (int) false) {
//                        printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
//                        continue;
//                    }

                    status = nastranCard_ctria6(
                        fp,
                        &feaMesh->element[i].elementID, // eid
                        &feaData->propertyID, //pid
                        feaMesh->element[i].connectivity, // g
                        theta, // theta
                        mcid, // mcid
                        &zoff, // zoffs
                        NULL, // t
                        feaFileFormat->gridFileType);
                    AIM_STATUS(aimInfo, status);
                }
            }

            if ( feaMesh->element[i].elementType == Quadrilateral) {

                if (feaData->elementSubType == ShearElement) {

//                    if (found == (int) false) {
//                        printf("No property information found for element %d of type \"ShearElement\"!", feaMesh->element[i].elementID);
//                        continue;
//                    }

                    status = nastranCard_cshear(
                        fp,
                        &feaMesh->element[i].elementID, // eid
                        &feaData->propertyID, //pid
                        feaMesh->element[i].connectivity, // g
                        feaFileFormat->gridFileType);
                    AIM_STATUS(aimInfo, status);
                }


                if (feaData->elementSubType == ShellElement) {

//                    if (found == (int) false) {
//                        printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
//                        continue;
//                    }

                    status = nastranCard_cquad4(
                        fp,
                        &feaMesh->element[i].elementID, // eid
                        &feaData->propertyID, //pid
                        feaMesh->element[i].connectivity, // g
                        theta, // theta
                        mcid, // mcid
                        &zoff, // zoffs
                        NULL, // t
                        feaFileFormat->gridFileType);
                    AIM_STATUS(aimInfo, status);
                }
            }

            if ( feaMesh->element[i].elementType == Quadrilateral_8) {

                if (feaData->elementSubType == ShellElement) {

//                    if (found == (int) false) {
//                        printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
//                        continue;
//                    }

                    status = nastranCard_cquad8(
                        fp,
                        &feaMesh->element[i].elementID, // eid
                        &feaData->propertyID, //pid
                        feaMesh->element[i].connectivity, // g
                        theta, // theta
                        mcid, // mcid
                        &zoff, // zoffs
                        NULL, // t
                        feaFileFormat->gridFileType);
                    AIM_STATUS(aimInfo, status);
                }
            }
        }
    }

    status = CAPS_SUCCESS;
cleanup:
    return status;
}

// Write a Nastran connections card from a feaConnection structure
int nastran_writeConnectionCard(FILE *fp, const feaConnectionStruct *feaConnect, const feaFileFormatStruct *feaFileFormat) {

    int status;

    if (fp == NULL) return CAPS_IOERR;
    if (feaConnect == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // Mass
    if (feaConnect->connectionType == Mass) {

        status = nastranCard_cmass2(
            fp,
            &feaConnect->elementID, // eid
            &feaConnect->mass, // m
            &feaConnect->connectivity[0], // g1
            &feaConnect->connectivity[1], // g2
            &feaConnect->componentNumberStart, // c1
            &feaConnect->componentNumberEnd, // c2
            feaFileFormat->gridFileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    // Spring
    if (feaConnect->connectionType == Spring) {

        status = nastranCard_celas2(
            fp,
            &feaConnect->elementID, // eid
            &feaConnect->stiffnessConst, // k
            &feaConnect->connectivity[0], // g1
            &feaConnect->connectivity[1], // g2
            &feaConnect->componentNumberStart, // c1
            &feaConnect->componentNumberEnd, // c2
            &feaConnect->dampingConst, // ge
            &feaConnect->stressCoeff, // s
            feaFileFormat->gridFileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    // Damper
    if (feaConnect->connectionType == Damper) {

        status = nastranCard_cdamp2(
            fp,
            &feaConnect->elementID, // eid
            &feaConnect->dampingConst, // b
            &feaConnect->connectivity[0], // g1
            &feaConnect->connectivity[1], // g2
            &feaConnect->componentNumberStart, // c1
            &feaConnect->componentNumberEnd, // c2
            feaFileFormat->gridFileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    // Rigid Body
    if (feaConnect->connectionType == RigidBody) {

        status = nastranCard_rbe2(
            fp,
            &feaConnect->elementID, // eid
            &feaConnect->connectivity[0], // gn
            &feaConnect->dofDependent, // cm
            1, &feaConnect->connectivity[1], // gm
            feaFileFormat->gridFileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    // Rigid Body Interpolate
    if (feaConnect->connectionType == RigidBodyInterpolate) {

        status = nastranCard_rbe3(
            fp,
            &feaConnect->elementID, // eid
            &feaConnect->connectivity[1], // refgrid
            &feaConnect->dofDependent, // refc
            feaConnect->numMaster,
            feaConnect->masterWeighting, // wt
            feaConnect->masterComponent, // c
            feaConnect->masterIDSet, // g
            0,
            NULL, // gm
            NULL, // cm
            feaFileFormat->gridFileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write a Nastran AERO card from a feaAeroRef structure
int nastran_writeAEROCard(FILE *fp, const feaAeroRefStruct *feaAeroRef, const feaFileFormatStruct *feaFileFormat) {

    if (fp == NULL) return CAPS_IOERR;
    if (feaAeroRef == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    return nastranCard_aero(
        fp,
        &feaAeroRef->coordSystemID, // acsid
        &feaAeroRef->refVelocity, // velocity
        &feaAeroRef->refChord, // refc
        &feaAeroRef->refDensity, // rhoref
        NULL, // symxz
        NULL, // symxy
        feaFileFormat->fileType
    );
}

// Write a Nastran AEROS card from a feaAeroRef structure
int nastran_writeAEROSCard(FILE *fp, const feaAeroRefStruct *feaAeroRef, const feaFileFormatStruct *feaFileFormat) {

    if (fp == NULL) return CAPS_IOERR;
    if (feaAeroRef == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;


    return nastranCard_aeros(
        fp,
        &feaAeroRef->coordSystemID, // acsid
        &feaAeroRef->rigidMotionCoordSystemID, // rcsid
        &feaAeroRef->refChord, // refc
        &feaAeroRef->refSpan, // refb
        &feaAeroRef->refArea, // refs
        &feaAeroRef->symmetryXZ, // symxz
        &feaAeroRef->symmetryXY, // symxy
        feaFileFormat->fileType
    );
}

// Write Nastran SET1 card from a feaAeroStruct
int nastran_writeSet1Card(FILE *fp, const feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat) {

    if (fp == NULL) return CAPS_IOERR;
    if (feaAero == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    return nastranCard_set1(
        fp,
        &feaAero->surfaceID, // sid
        feaAero->numGridID, feaAero->gridIDSet, // g
        feaFileFormat->fileType
    );
}

// Write Nastran Spline1 cards from a feaAeroStruct
int nastran_writeAeroSplineCard(FILE *fp, const feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat) {

    int numSpanWise;
    int boxBegin, boxEnd;

    if (fp == NULL) return CAPS_IOERR;
    if (feaAero == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (feaAero->vlmSurface.NspanTotal > 0)
        numSpanWise = feaAero->vlmSurface.NspanTotal;
    else if (feaAero->vlmSurface.NspanSection > 0)
        numSpanWise = (feaAero->vlmSurface.numSection-1)*feaAero->vlmSurface.NspanSection;
    else {
        printf("Error: Only one of numSpanTotal and numSpanPerSection can be non-zero!\n");
        printf("       numSpanTotal      = %d\n", feaAero->vlmSurface.NspanTotal);
        printf("       numSpanPerSection = %d\n", feaAero->vlmSurface.NspanSection);
        return CAPS_BADVALUE;
    }

    boxBegin = feaAero->surfaceID;
    boxEnd = feaAero->surfaceID + numSpanWise*feaAero->vlmSurface.Nchord - 1;

    return nastranCard_spline1(
        fp,
        &feaAero->surfaceID, // eid
        &feaAero->surfaceID, // caero
        &boxBegin, // box1
        &boxEnd, // box2
        &feaAero->surfaceID, // setg
        NULL, // dz
        feaFileFormat->fileType
    );
}

static inline double _getSectionChordLength(vlmSectionStruct *section) {

    return sqrt(pow(section->xyzTE[0] - section->xyzLE[0], 2) +
                pow(section->xyzTE[1] - section->xyzLE[1], 2) +
                pow(section->xyzTE[2] - section->xyzLE[2], 2));
}

// Get divisions as equal fractions from 0.0 to 1.0
static inline int _getDivisions(int numDivs, double **divisionsOut) {

    int i;

    double *divisions = NULL;

    divisions = EG_alloc(numDivs * sizeof(double));
    if (divisions == NULL) return EGADS_MALLOC;

    divisions[0] = 0.0;
    for (i = 1; i < numDivs-1; i++) {
        divisions[i] = divisions[i-1] + 1. / numDivs;
    }
    divisions[numDivs-1] = 1.0;

    *divisionsOut = divisions;

    return CAPS_SUCCESS;
}

// Determine index of closest division percent to control surface percent chord
static inline int _getClosestDivisionIndex(int numDivs, double *divs, double percentChord, int *closestDivIndexOut) {

    int i, closestDivIndex = 0;
    double closestDivDist = 1.0, divDist;

    for (i = 0; i < numDivs; i++) {
        divDist = fabs(percentChord - divs[i]);
        if (divDist < closestDivDist) {
            closestDivDist = divDist;
            closestDivIndex = i;
        }
    }

    if ((closestDivIndex == 0) || (closestDivIndex == numDivs-1) || (closestDivDist == 1.0)) {
        return CAPS_BADVALUE;
    }

    *closestDivIndexOut = closestDivIndex;

    return CAPS_SUCCESS;
}

// Get set of box IDs corresponding to control surface
static int _getControlSurfaceBoxIDs(int boxBeginID, int numChordDivs,
                       /*@unused@*/ double *chordDivs,
                                    int numSpanDivs,
                       /*@unused@*/ double *spanDivs,
                                    int hingelineIndex,
                                    int isTrailing, int *numBoxIDsOut, int **boxIDsOut) {

    int ichord, ispan, csBoxIndex, boxCount, chordDivIndex;

    int boxID, numBoxIDs, *boxIDs = NULL;

    numBoxIDs = (numChordDivs-1) * (numSpanDivs-1);

    // conservative allocate
    boxIDs = EG_alloc(numBoxIDs * sizeof(int));
    if (boxIDs == NULL) return EGADS_MALLOC;

    boxCount = 0;
    csBoxIndex = 0;

    for (ispan = 0; ispan < numSpanDivs-1; ispan++) {

        for (ichord = 0; ichord < numChordDivs-1; ichord++) {

            boxID = boxBeginID + boxCount++;

            chordDivIndex = ichord + 1;

            if (!isTrailing && chordDivIndex <= hingelineIndex) {
                boxIDs[csBoxIndex++] = boxID;
            }
            else if (isTrailing && chordDivIndex > hingelineIndex) {
                boxIDs[csBoxIndex++] = boxID;
            }
        }
    }

    numBoxIDs = csBoxIndex;

    *numBoxIDsOut = numBoxIDs;
    *boxIDsOut = boxIDs;

    return CAPS_SUCCESS;
}


static int _writeAesurfCard(FILE *fp,
                            const feaAeroStruct *feaAero,
                            vlmSectionStruct *rootSection,
                            vlmSectionStruct *tipSection,
                            const feaFileFormatStruct *feaFileFormat) {

    int i, j, status;

    int controlID, coordSystemID, aelistID;
    int numChordDivs, numSpanDivs, numBoxes;
    double *chordDivs = NULL, *spanDivs = NULL;
    int *boxIDs = NULL, hingelineDivIndex;
    double xyzHingeVec[3] = {0.0, 0.0, 0.0};
    double xyzHingeVecNorm[3] = {0.0, 0.0, 0.0};
    double pointA[3] = {0.0, 0.0, 0.0};
    double pointB[3] = {0.0, 0.0, 0.0};
    double pointC[3] = {0.0, 0.0, 0.0};

    int found;

    vlmControlStruct *controlSurface, *controlSurface2 = NULL;

    for (i = 0; i < rootSection->numControl; i++) {

        controlSurface = &rootSection->vlmControl[i];

        // find matching control surface info on tip section
        found = (int) false;
        controlSurface2 = NULL;
        for (j = 0; j < tipSection->numControl; j++) {
            if (strcmp(rootSection->vlmControl[i].name, tipSection->vlmControl[j].name) == 0) {
                controlSurface2 = &tipSection->vlmControl[j];
                found = (int) true;
                break;
            }
        }
        if (!found) continue;

//        // get hingeline vector
//        for (j = 0; j < 3; j++) {
//            xyzHingeVec[j] = controlSurface2->xyzHinge[j] - controlSurface->xyzHinge[j];
//        }
//
//        controlID = feaAero->surfaceID + i;
//
//        // get control surface coordinate system points
//        pointA[0] = controlSurface->xyzHinge[0];
//        pointA[1] = controlSurface->xyzHinge[1];
//        pointA[2] = controlSurface->xyzHinge[2];
//
//        pointB[0] = pointA[0];
//        pointB[1] = pointA[1];
//        pointB[2] = pointA[2] + 1;
//
//        pointC[0] = pointA[0] + 1;
//        pointC[1] = xyzHingeVec[0]/xyzHingeVec[1] * pointC[0];
//        pointC[2] = pointA[2] + 0.5;


//        printf("X: %f: %f\n", controlSurface2->xyzHinge[0], controlSurface->xyzHinge[0]);
//        printf("Y: %f: %f\n", controlSurface2->xyzHinge[1], controlSurface->xyzHinge[1]);
//        printf("Z: %f: %f\n", controlSurface2->xyzHinge[2], controlSurface->xyzHinge[2]);

        // get hingeline vector : only in x-y, assum Z is zero
        for (j = 0; j < 2; j++) {
            xyzHingeVec[j] = controlSurface2->xyzHinge[j] - controlSurface->xyzHinge[j];
        }

        for (j = 0; j < 2; j++) {
            xyzHingeVecNorm[j] = xyzHingeVec[j] / sqrt(pow(xyzHingeVec[0],2) + pow(xyzHingeVec[1],2) + pow(xyzHingeVec[2],2));
        }

        controlID = feaAero->surfaceID + i;
        // get control surface coordinate system points, average the z location
        pointA[0] = controlSurface->xyzHinge[0];
        pointA[1] = controlSurface->xyzHinge[1];
        pointA[2] = (controlSurface->xyzHinge[2]+controlSurface2->xyzHinge[2])/2;

        pointB[0] = pointA[0];
        pointB[1] = pointA[1];
        pointB[2] = pointA[2] + 10.0;

        pointC[0] = pointA[0] + xyzHingeVecNorm[1]*10;
        pointC[1] = pointA[1] - xyzHingeVecNorm[0]*10;
        pointC[2] = pointA[2] + 1.0;

        // get chordwise division fractions
        numChordDivs = feaAero->vlmSurface.Nchord + 1;
        status = _getDivisions(numChordDivs, &chordDivs);
        if (status != CAPS_SUCCESS) goto cleanup;

        // get spanwise division fractions
        numSpanDivs = feaAero->vlmSurface.NspanTotal + 1;
        status = _getDivisions(numSpanDivs, &spanDivs);
        if (status != CAPS_SUCCESS) goto cleanup;

        // hingeline is the closest chord div to percent chord
        status = _getClosestDivisionIndex(
            numChordDivs, chordDivs, controlSurface->percentChord, &hingelineDivIndex);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = _getControlSurfaceBoxIDs(feaAero->surfaceID,
                                            numChordDivs, chordDivs,
                                            numSpanDivs, spanDivs,
                                            hingelineDivIndex, controlSurface->leOrTe,
                                            &numBoxes, &boxIDs);
        if (status != CAPS_SUCCESS) goto cleanup;

        coordSystemID = controlID;
        status = nastranCard_cord2r(
            fp,
            &coordSystemID, // cid
            NULL, // rid
            pointA, // a
            pointB, // b
            pointC, // c
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) goto cleanup;

        aelistID = controlID;
        status = nastranCard_aelist(
            fp,
            &aelistID,
            numBoxes,
            boxIDs,
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) goto cleanup;

        status = nastranCard_aesurf(
            fp,
            &controlID, // id, ignored
            controlSurface->name, // label
            &coordSystemID, // cid
            &aelistID, // alid
            NULL, // eff
            "LDW", // ldw
            NULL, // crefc
            NULL, // crefs
            NULL, // pllim
            NULL, // pulim
            NULL, // hmllim
            NULL, // hmulim
            NULL, // tqllim
            NULL, // tqulim
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) goto cleanup;

        if (chordDivs != NULL) {
            EG_free(chordDivs);
            chordDivs = NULL;
        }
        if (spanDivs != NULL) {
            EG_free(spanDivs);
            spanDivs = NULL;
        }
        if (boxIDs != NULL) {
            EG_free(boxIDs);
            boxIDs = NULL;
        }
    }

    status = CAPS_SUCCESS;

    cleanup:

        if (chordDivs != NULL) EG_free(chordDivs);
        if (spanDivs != NULL) EG_free(spanDivs);
        if (boxIDs != NULL) EG_free(boxIDs);

        return status;
}

// Write Nastran CAERO1 cards from a feaAeroStruct
int nastran_writeCAeroCard(FILE *fp, const feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat) {

    int status;

    int i, sectionIndex; // Indexing

    const int *nspan, *nchord, *lspan, *lchord;
    int defaultIGroupID = 1;

    double chordLength12, chordLength43;
    double *xyz1, *xyz4;

    vlmSectionStruct *rootSection, *tipSection;

    if (fp == NULL) return CAPS_IOERR;
    if (feaAero == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    //printf("NUMBER OF SECTIONS = %d\n", feaAero->vlmSurface.numSection);
    for (i = 0; i < feaAero->vlmSurface.numSection-1; i++) {

        // If Cspace and/or Sspace is something (to be defined) lets write a AEFact card instead with our own distributions
        if (feaAero->vlmSurface.Sspace == 0.0) {
            nspan = &feaAero->vlmSurface.NspanTotal;
            lspan = NULL;
        } else {
            nspan = NULL;
            lspan = 0; // TODO: aefact card id
            PRINT_WARNING("Definition of spanwise boxes via LSPAN not implemented yet!\n");
        }

        if (feaAero->vlmSurface.Cspace == 0.0) {
            nchord = &feaAero->vlmSurface.Nchord;
            lchord = NULL;
        } else {
            nchord = NULL;
            lchord = 0; // TODO: aefact card id
            PRINT_WARNING("Definition of chordwise boxes via LCHORD not implemented yet!\n");
        }

        sectionIndex = feaAero->vlmSurface.vlmSection[i].sectionIndex;
        rootSection = &feaAero->vlmSurface.vlmSection[sectionIndex];
        xyz1 = rootSection->xyzLE;
        chordLength12 = _getSectionChordLength(rootSection);

        sectionIndex = feaAero->vlmSurface.vlmSection[i+1].sectionIndex;
        tipSection = &feaAero->vlmSurface.vlmSection[sectionIndex];
        xyz4 = tipSection->xyzLE;
        chordLength43 = _getSectionChordLength(tipSection);

        // Write necessary PAER0 card
        status = nastranCard_paero1(
            fp,
            &feaAero->surfaceID, // pid
            0, NULL, // b
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        status = nastranCard_caero1(
            fp,
            &feaAero->surfaceID, // eid
            &feaAero->surfaceID, // pid
            &feaAero->coordSystemID, // cp
            nspan, // nspan
            nchord, // nchord
            lspan, // lspan
            lchord, // lchord
            &defaultIGroupID, // igid
            xyz1, // xyz1
            xyz4, // xyz4
            &chordLength12, // x12
            &chordLength43, // x43
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        if (rootSection->numControl > 0) {

            status = _writeAesurfCard(fp, feaAero, rootSection, tipSection, feaFileFormat);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    return CAPS_SUCCESS;
}

// Write Nastran coordinate system card from a feaCoordSystemStruct structure
int nastran_writeCoordinateSystemCard(FILE *fp, const feaCoordSystemStruct *feaCoordSystem, const feaFileFormatStruct *feaFileFormat) {

    double pointA[3] = {0, 0, 0};
    double pointB[3] = {0, 0, 0};
    double pointC[3] = {0, 0, 0};

    if (fp == NULL) return CAPS_IOERR;
    if (feaCoordSystem == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    pointA[0] = feaCoordSystem->origin[0];
    pointA[1] = feaCoordSystem->origin[1];
    pointA[2] = feaCoordSystem->origin[2];

    pointB[0] = feaCoordSystem->normal3[0] + pointA[0];
    pointB[1] = feaCoordSystem->normal3[1] + pointA[1];
    pointB[2] = feaCoordSystem->normal3[2] + pointA[2];

    pointC[0] = feaCoordSystem->normal1[0] + pointB[0];
    pointC[1] = feaCoordSystem->normal1[1] + pointB[1];
    pointC[2] = feaCoordSystem->normal1[2] + pointB[2];

    // Rectangular
    if (feaCoordSystem->coordSystemType == RectangularCoordSystem) {

        return nastranCard_cord2r(
            fp,
            &feaCoordSystem->coordSystemID, // cid
            &feaCoordSystem->refCoordSystemID, //rid
            pointA, // a
            pointB, // b
            pointC, // c
            feaFileFormat->fileType
        );
    }
    // Spherical
    else if (feaCoordSystem->coordSystemType == SphericalCoordSystem) {

        return nastranCard_cord2s(
            fp,
            &feaCoordSystem->coordSystemID, // cid
            &feaCoordSystem->refCoordSystemID, //rid
            pointA, // a
            pointB, // b
            pointC, // c
            feaFileFormat->fileType
        );
    }
    // Cylindrical
    else if (feaCoordSystem->coordSystemType == CylindricalCoordSystem) {

        return nastranCard_cord2c(
            fp,
            &feaCoordSystem->coordSystemID, // cid
            &feaCoordSystem->refCoordSystemID, //rid
            pointA, // a
            pointB, // b
            pointC, // c
            feaFileFormat->fileType
        );
    }
    else {

        PRINT_ERROR("Unrecognized coordinate system type !!\n");
        return CAPS_BADVALUE;
    }
}

// Write combined Nastran constraint card from a set of constraint IDs.
// 	The combined constraint ID is set through the constraintID variable.
int nastran_writeConstraintADDCard(FILE *fp, int constraintID, int numSetID, int constraintSetID[], const feaFileFormatStruct *feaFileFormat) {

    if (fp == NULL) return CAPS_IOERR;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    return nastranCard_spcadd(
        fp, &constraintID, numSetID, constraintSetID, feaFileFormat->fileType);
}

// Write Nastran constraint card from a feaConstraint structure
int nastran_writeConstraintCard(FILE *fp, const feaConstraintStruct *feaConstraint, const feaFileFormatStruct *feaFileFormat) {

    int status;
    int i; // Indexing;

    if (fp == NULL) return CAPS_IOERR;
    if (feaConstraint == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (feaConstraint->constraintType == Displacement) {

        for (i = 0; i < feaConstraint->numGridID; i++) {

            status = nastranCard_spc(
                fp,
                &feaConstraint->constraintID, // sid
                1, // currently always 1 single tuple of Gi Ci Di
                &feaConstraint->gridIDSet[i], // g
                &feaConstraint->dofConstraint, // c
                &feaConstraint->gridDisplacement, // d
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    if (feaConstraint->constraintType == ZeroDisplacement) {

        for (i = 0; i < feaConstraint->numGridID; i++) {

            status = nastranCard_spc1(
                fp,
                &feaConstraint->constraintID, // sid
                &feaConstraint->dofConstraint, // c
                1, // currently always 1 single Gi value
                &feaConstraint->gridIDSet[i], // g
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    return CAPS_SUCCESS;
}

// Write Nastran support card from a feaSupport structure - withID = NULL or false SUPORT, if withID = true SUPORT1
int nastran_writeSupportCard(FILE *fp, const feaSupportStruct *feaSupport, const feaFileFormatStruct *feaFileFormat, const int *withID) {

    int status;

    int i; // Indexing;

    if (fp == NULL) return CAPS_IOERR;
    if (feaSupport == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    for (i = 0; i < feaSupport->numGridID; i++) {

        if (withID != NULL && *withID == (int) true) {

            status = nastranCard_suport1(
                fp,
                &feaSupport->supportID, // sid
                1,
                &feaSupport->gridIDSet[i], // id
                &feaSupport->dofSupport, // c
                feaFileFormat->fileType
            );
        }
        else {

            status = nastranCard_suport(
                fp,
                1,
                &feaSupport->gridIDSet[i], // id
                &feaSupport->dofSupport, // c
                feaFileFormat->fileType
            );
        }

        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write a Nastran Material card from a feaMaterial structure
int nastran_writeMaterialCard(FILE *fp, const feaMaterialStruct *feaMaterial, const feaFileFormatStruct *feaFileFormat) {

    double strainAllowable = 1.0;
    int status = 0;

    const double *g1z = NULL, *g2z = NULL, *xt = NULL, *xc = NULL;
    const double *yt = NULL, *yc = NULL, *s = NULL, *strn = NULL;
    const double *g = NULL, *a = NULL, *tref = NULL, *ge = NULL;

    const double *cs=NULL, *ec=NULL, *gc=NULL;
    const double *alpha0=NULL, *sb=NULL, *ef1=NULL;
    const double *nuf12=NULL, *msmf=NULL;
    const double *pnpt=NULL, *pnpc=NULL;
    const char   *ft=NULL;
    const double *nb=NULL, *e3=NULL;
    const double *nu23=NULL, *nu31=NULL;
    const double *e1rsf=NULL, *e2rsf=NULL;
    const double *g12rsf=NULL, *g1zrsf=NULL, *g2zrsf=NULL;
    const int    *te1rsf=NULL, *te2rsf=NULL;
    const int    *tg12rsf=NULL, *tg1zrsf=NULL, *tg2zrsf=NULL;

    if (fp == NULL) return CAPS_IOERR;
    if (feaMaterial == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (feaMaterial->thermalExpCoeff != 0.0) {
        a = &feaMaterial->thermalExpCoeff;
    }
    if (feaMaterial->temperatureRef != 0.0) {
        tref = &feaMaterial->temperatureRef;
    }
    if (feaMaterial->dampingCoeff != 0.0) {
        ge = &feaMaterial->dampingCoeff;
    }

    // Isotropic
    if (feaMaterial->materialType == Isotropic) {

        status = nastranCard_mat1(
            fp,
            &feaMaterial->materialID, // mid
            &feaMaterial->youngModulus, // e
            g, // g
            &feaMaterial->poissonRatio, // nu
            &feaMaterial->density, // rho
            a, // a
            tref, // tref
            ge, // ge
            &feaMaterial->tensionAllow, // st
            &feaMaterial->compressAllow, // sc
            &feaMaterial->shearAllow, // ss
            NULL, // mcsid
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        // thermal properties
        if (feaMaterial->kappa != 0.0) {
            status = nastranCard_mat4(fp,
                &feaMaterial->materialID, // mid, shifted for pyTACS to read
                &feaMaterial->kappa, // k
                &feaMaterial->specificHeat, // cp
                &feaMaterial->density, // rho
                NULL, // H
                NULL, // HGEN
                feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }

    }

    // Orthotropic
    if (feaMaterial->materialType == Orthotropic) {

        // thermal properties
        if (feaMaterial->K[0] != 0.0 ||
            feaMaterial->K[1] != 0.0 ||
            feaMaterial->K[2] != 0.0 ||
            feaMaterial->K[3] != 0.0 ||
            feaMaterial->K[4] != 0.0 ||
            feaMaterial->K[5] != 0.0) {
            status = nastranCard_mat5(fp,
                &feaMaterial->materialID,
                feaMaterial->K, // K
                &feaMaterial->specificHeat, // cp
                &feaMaterial->density, // rho
                NULL, // HGEN
                feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }

        if (feaMaterial->shearModulusTrans1Z != 0)
            g1z = &feaMaterial->shearModulusTrans1Z;

        if (feaMaterial->shearModulusTrans2Z != 0)
            g2z = &feaMaterial->shearModulusTrans2Z;

        if (feaMaterial->tensionAllow != 0)
            xt = &feaMaterial->tensionAllow;

        if (feaMaterial->compressAllow != 0)
            xc = &feaMaterial->compressAllow;

        if (feaMaterial->tensionAllowLateral != 0)
            yt = &feaMaterial->tensionAllowLateral;

        if (feaMaterial->compressAllowLateral != 0)
            yc = &feaMaterial->compressAllowLateral;

        if (feaMaterial->shearAllow != 0)
            s = &feaMaterial->shearAllow;

        if (feaMaterial->allowType != 0)
            strn = &strainAllowable;

        if (feaMaterial->honeycombCellSize >= 0)
            cs = &feaMaterial->honeycombCellSize;

        if (feaMaterial->honeycombYoungModulus >= 0)
            ec = &feaMaterial->honeycombYoungModulus;

        if (feaMaterial->honeycombShearModulus >= 0)
            gc = &feaMaterial->honeycombShearModulus;

        if (feaMaterial->fractureAngle > 0)
            alpha0 = &feaMaterial->fractureAngle;

        if (feaMaterial->interlaminarShearAllow >= 0)
            sb = &feaMaterial->interlaminarShearAllow;

        if (feaMaterial->fiberYoungModulus > 0)
            ef1 = &feaMaterial->fiberYoungModulus;

        if (feaMaterial->fiberPoissonRatio >= 0)
            nuf12 = &feaMaterial->fiberPoissonRatio;

        if (feaMaterial->meanStressFactor >= 0)
            msmf = &feaMaterial->meanStressFactor;

        if (feaMaterial->transTensionSlope >= 0)
            pnpt = &feaMaterial->transTensionSlope;

        if (feaMaterial->transCompressionSlope >= 0)
            pnpc = &feaMaterial->transCompressionSlope;

        if (feaMaterial->compositeFailureTheory != NULL)
            ft = feaMaterial->compositeFailureTheory;

        if (feaMaterial->interlaminarNormalStressAllow >= 0)
            nb = &feaMaterial->interlaminarNormalStressAllow;

        if (feaMaterial->youngModulusThick >= 0)
            e3 = &feaMaterial->youngModulusThick;

        if (feaMaterial->poissonRatio23 != 0)
            nu23 = &feaMaterial->poissonRatio23;

        if (feaMaterial->poissonRatio31 != 0)
            nu31 = &feaMaterial->poissonRatio31;

        if (feaMaterial->youngModulusFactor >= 0)
            e1rsf = &feaMaterial->youngModulusFactor;

        if (feaMaterial->youngModulusLateralFactor >= 0)
            e2rsf = &feaMaterial->youngModulusLateralFactor;

        if (feaMaterial->shearModulusFactor >= 0)
            g12rsf = &feaMaterial->shearModulusFactor;

        if (feaMaterial->shearModulusTrans1ZFactor >= 0)
            g1zrsf = &feaMaterial->shearModulusTrans1ZFactor;

        if (feaMaterial->shearModulusTrans2ZFactor >= 0)
            g2zrsf = &feaMaterial->shearModulusTrans2ZFactor;

        return nastranCard_mat8(
            fp,
            &feaMaterial->materialID, // mid
            &feaMaterial->youngModulus, // e1
            &feaMaterial->youngModulusLateral, // e2
            &feaMaterial->poissonRatio, // nu12
            &feaMaterial->shearModulus, // g12
            g1z, // g1z
            g2z, // g2z
            &feaMaterial->density, // rho
            &feaMaterial->thermalExpCoeff, // a1
            &feaMaterial->thermalExpCoeffLateral, // a2
            &feaMaterial->temperatureRef, // tref
            xt, // xt
            xc, // xc
            yt, // yt
            yc, // yc
            s, // s
            &feaMaterial->dampingCoeff, // ge
            NULL, // f12
            strn, // strn
            cs, ec, gc,
            alpha0, sb, ef1,
            nuf12, msmf,
            pnpt, pnpc,
            ft,
            nb, e3,
            nu23, nu31,
            e1rsf, e2rsf,
            g12rsf, g1zrsf, g2zrsf,
            te1rsf, te2rsf,
            tg12rsf, tg1zrsf, tg2zrsf,
            feaFileFormat->fileType
        );
    }

    // Anisothotropic
    if (feaMaterial->materialType == Anisothotropic) {

        return nastranCard_mat2(
            fp,
            &feaMaterial->materialID, // mid
            &feaMaterial->gmat[0], // g11
            &feaMaterial->gmat[1], // g12
            &feaMaterial->gmat[2], // g13
            &feaMaterial->gmat[3], // g22
            &feaMaterial->gmat[4], // g23
            &feaMaterial->gmat[5], // g33
            &feaMaterial->density, // rho
            NULL, // a1
            NULL, // a2
            NULL, // a3
            NULL, // &feaMaterial->temperatureRef, // tref
            NULL, // &feaMaterial->dampingCoeff, // ge
            NULL, // st
            NULL, // sc
            NULL, // ss
            NULL, // mcsid
            NULL, // ge11
            NULL, // ge12
            NULL, // ge13
            NULL, // ge22
            NULL, // ge23
            NULL, // ge33
            feaFileFormat->fileType
        );
    }

    return CAPS_SUCCESS;
}

// Write a Nastran Property card from a feaProperty structure
int nastran_writePropertyCard(FILE *fp, const feaPropertyStruct *feaProperty, const feaFileFormatStruct *feaFileFormat) {

    const double *nsm = NULL; // massPerLength, massPerArea field

    // pshell
    const int *mid2 = NULL, *mid3 = NULL;
    const double *i12t3 = NULL, *tst = NULL;

    // pcomp
    char *lam = NULL;

    if (fp == NULL) return CAPS_IOERR;
    if (feaProperty == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    fprintf(fp, "$ Femap Property  : %s\n", feaProperty->name);

    //          1D Elements        //

    // Rod
    if (feaProperty->propertyType ==  Rod) {

        return nastranCard_prod(
            fp,
            &feaProperty->propertyID, // pid
            &feaProperty->materialID, // mid
            &feaProperty->crossSecArea, // a
            &feaProperty->torsionalConst, // j
            &feaProperty->torsionalStressReCoeff, // c
            &feaProperty->massPerLength, // nsm
            feaFileFormat->fileType
        );
    }

    // Bar
    if (feaProperty->propertyType ==  Bar) {

        if (feaProperty->crossSecType != NULL) {

            return nastranCard_pbarl(
                fp,
                &feaProperty->propertyID, // pid
                &feaProperty->materialID, // mid
                feaProperty->crossSecType, // type
                NULL, // f0
                10, feaProperty->crossSecDimension,
                &feaProperty->massPerLength, // nsm
                feaFileFormat->fileType
            );

        } else {

            return nastranCard_pbar(
                fp,
                &feaProperty->propertyID, // pid
                &feaProperty->materialID, // mid
                &feaProperty->crossSecArea, // a
                &feaProperty->zAxisInertia, // i1
                &feaProperty->yAxisInertia, // i2
                NULL, // i12
                &feaProperty->torsionalConst, // j
                &feaProperty->massPerLength, // nsm
                NULL, // c
                NULL, // d
                NULL, // e
                NULL, // f
                NULL, // k1
                NULL, // k2
                feaFileFormat->fileType
            );
        }
    }

    //          2D Elements        //

    // Shell
    if (feaProperty->propertyType == Shell) {

        if (feaProperty->materialBendingID != 0) {
            mid2 = &feaProperty->materialBendingID;
            i12t3 = &feaProperty->bendingInertiaRatio;
        }

        if (feaProperty->materialShearID != 0) {
            mid3 = &feaProperty->materialShearID;
            tst = &feaProperty->shearMembraneRatio;
        }

        if (feaProperty->massPerArea != 0) {
            nsm = &feaProperty->massPerArea;
        }

        return nastranCard_pshell(
            fp,
            &feaProperty->propertyID, // pid
            &feaProperty->materialID, // mid
            &feaProperty->membraneThickness, // t
            mid2, // mid2
            i12t3, // i12t3
            mid3, // mid3
            tst, // tst
            nsm, // nsm
            NULL, // z1
            NULL, // z2
            NULL, // mid4
            feaFileFormat->fileType
        );
    }

    // Shear
    if (feaProperty->propertyType == Shear) {

        if (feaProperty->massPerArea != 0) {
            nsm = &feaProperty->massPerArea;
        }

        return nastranCard_pshear(
            fp,
            &feaProperty->propertyID, // pid
            &feaProperty->materialID, // mid
            &feaProperty->membraneThickness, // t
            nsm, // nsm
            NULL, // f1
            NULL, // f2
            feaFileFormat->fileType
        );
    }

    // Composite
    if (feaProperty->propertyType == Composite) {

        if (feaProperty->massPerArea != 0) {
            nsm = &feaProperty->massPerArea;
        }

        if (feaProperty->compositeSymmetricLaminate == (int) true) {
            lam = "SYM";
        }

        return nastranCard_pcomp(
            fp,
            &feaProperty->propertyID, // pid
            NULL, // z0
            nsm, // nsm
            &feaProperty->compositeShearBondAllowable, // sb
            feaProperty->compositeFailureTheory,
            NULL, // tref
            NULL, // ge
            lam, // lam
            feaProperty->numPly, // number of layers
            feaProperty->compositeMaterialID, // mid
            feaProperty->compositeThickness, // t
            feaProperty->compositeOrientation, // theta
            NULL, // sout
            feaFileFormat->fileType
        );
    }

    //          3D Elements        //

    // Solid
    if (feaProperty->propertyType == Solid) {

        return nastranCard_psolid(
            fp,
            &feaProperty->propertyID, // pid
            &feaProperty->materialID, // mid
            NULL, // cordm
            NULL, // in
            NULL, // stress
            NULL, // isop
            NULL, // fctn
            feaFileFormat->fileType
        );
    }

    return CAPS_SUCCESS;
}

// Write a combined Nastran load card from a set of load IDs. Uses the feaLoad structure to get the local scale factor
// 	for the load. The overall load scale factor is set to 1. The combined load ID is set through the loadID variable.
int nastran_writeLoadADDCard(FILE *fp, int loadID, int numSetID, int loadSetID[], feaLoadStruct feaLoad[], const feaFileFormatStruct *feaFileFormat) {

    int status;

    int i; // Indexing

    double overallScale = 1.0, *loadScaleFactors = NULL;

    if (fp == NULL) return CAPS_IOERR;
    if (numSetID > 0 && feaLoad == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    loadScaleFactors = (double *) EG_alloc(sizeof(double) * numSetID);
    if (loadScaleFactors == NULL) {
        return EGADS_MALLOC;
    }

    for (i = 0; i < numSetID; i++) {
        // ID's are 1 bias
        loadScaleFactors[i] = feaLoad[loadSetID[i]-1].loadScaleFactor;
    }

    status = nastranCard_load(
        fp,
        &loadID,
        &overallScale,
        numSetID,
        loadScaleFactors,
        loadSetID,
        feaFileFormat->fileType
    );

    EG_free(loadScaleFactors);

    return status;
}

// Write a Nastran Load card from a feaLoad structure
int nastran_writeLoadCard(FILE *fp, const feaLoadStruct *feaLoad, const feaFileFormatStruct *feaFileFormat)
{
    int status;

    int i; // Indexing

    if (fp == NULL) return CAPS_IOERR;
    if (feaLoad == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // Concentrated force at a grid point
    if (feaLoad->loadType ==  GridForce) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = nastranCard_force(
                fp,
                &feaLoad->loadID, // sid
                &feaLoad->gridIDSet[i], // g
                &feaLoad->coordSystemID, // cid
                &feaLoad->forceScaleFactor, // f
                feaLoad->directionVector, // n
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Concentrated moment at a grid point
    if (feaLoad->loadType ==  GridMoment) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = nastranCard_moment(
                fp,
                &feaLoad->loadID, // sid
                &feaLoad->gridIDSet[i], // g
                &feaLoad->coordSystemID, // cid
                &feaLoad->momentScaleFactor, // m
                feaLoad->directionVector, // n
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Gravitational load
    if (feaLoad->loadType ==  Gravity) {

        status = nastranCard_grav(
            fp,
            &feaLoad->loadID, // sid
            &feaLoad->coordSystemID, // cid
            &feaLoad->gravityAcceleration, // g
            feaLoad->directionVector, // n
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    // Pressure load
    if (feaLoad->loadType ==  Pressure) {

        for (i = 0; i < feaLoad->numElementID; i++) {

            status = nastranCard_pload2(
                fp,
                &feaLoad->loadID, // sid
                &feaLoad->pressureForce, // p
                1, &feaLoad->elementIDSet[i], // eid
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Pressure load at element Nodes
    if (feaLoad->loadType ==  PressureDistribute) {

        for (i = 0; i < feaLoad->numElementID; i++) {

            status = nastranCard_pload4(
                fp,
                &feaLoad->loadID, // sid
                &feaLoad->elementIDSet[i], // eid
                feaLoad->pressureDistributeForce, // j
                NULL, // g1
                NULL, // g3
                NULL, // cid
                NULL, // n
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Pressure load at element Nodes - different distribution per element
    if (feaLoad->loadType ==  PressureExternal) {

        for (i = 0; i < feaLoad->numElementID; i++) {

            status = nastranCard_pload4(
                fp,
                &feaLoad->loadID, // sid
                &feaLoad->elementIDSet[i], // eid
                &feaLoad->pressureMultiDistributeForce[4*i], // j
                NULL, // g1
                NULL, // g3
                NULL, // cid
                NULL, // n
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Rotational velocity
    if (feaLoad->loadType ==  Rotational) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = nastranCard_rforce(
                fp,
                &feaLoad->loadID, // sid
                &feaLoad->gridIDSet[i], // g
                &feaLoad->coordSystemID, // cid
                &feaLoad->angularVelScaleFactor, // a
                feaLoad->directionVector, // r
                NULL, // method
                &feaLoad->angularAccScaleFactor, // racc
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Default Thermal load at all unspecified grid point
    if (feaLoad->loadType == Thermal ||
        feaLoad->loadType == ThermalExternal) {

        status = nastranCard_tempd(
            fp,
            1, &feaLoad->loadID, // sid
            &feaLoad->temperatureDefault, // t
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    // Thermal load at a grid point
    if (feaLoad->loadType == Thermal) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = nastranCard_temp(
                fp,
                &feaLoad->loadID, // sid
                1, &feaLoad->gridIDSet[i], // g
                &feaLoad->temperature, // t
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Thermal load distributed at grid coordinates
    if (feaLoad->loadType == ThermalExternal) {

        status = nastranCard_temp(
            fp,
            &feaLoad->loadID, // sid
            feaLoad->numGridID, feaLoad->gridIDSet, // g
            feaLoad->temperatureMultiDistribute, // t
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write a Nastran Analysis card from a feaAnalysis structure
int nastran_writeAnalysisCard(FILE *fp,
                              feaAnalysisStruct *feaAnalysis,
                              const feaFileFormatStruct *feaFileFormat) {

    int status;

    int i; // Indexing
    int analysisID, densityID, machID, velocityID;

    double *mach = NULL;

    // trim
    int numParams, paramIndex;
    double *paramMags;
    char **paramLabels;

    if (fp == NULL) return CAPS_IOERR;
    if (feaAnalysis == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // Eigenvalue
    if (feaAnalysis->analysisType == Modal || feaAnalysis->analysisType == AeroelasticFlutter || feaAnalysis->analysisType == Gust ) {

        if (feaAnalysis->extractionMethod != NULL &&
            strcasecmp(feaAnalysis->extractionMethod, "Lanczos") == 0) {

            status = nastranCard_eigrl(
                fp,
                &feaAnalysis->analysisID, // sid
                &feaAnalysis->frequencyRange[0], // v1
                &feaAnalysis->frequencyRange[1], // v2
                &feaAnalysis->numDesiredEigenvalue, // nd
                NULL, // msglvl
                NULL, // maxset
                NULL, // shfscl
                feaAnalysis->eigenNormalization,
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;

        } else {

            status = nastranCard_eigr(
                fp,
                &feaAnalysis->analysisID, // sid
                feaAnalysis->extractionMethod, // method
                &feaAnalysis->frequencyRange[0], // v1
                &feaAnalysis->frequencyRange[1], // v2
                &feaAnalysis->numEstEigenvalue, // ne
                &feaAnalysis->numDesiredEigenvalue, // nd
                feaAnalysis->eigenNormalization, // norm
                &feaAnalysis->gridNormalization, // g
                &feaAnalysis->componentNormalization, //c
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    if (feaAnalysis->analysisType == AeroelasticTrim) {

        numParams = (feaAnalysis->numRigidConstraint
                     + feaAnalysis->numControlConstraint);

        paramLabels = EG_alloc(sizeof(char *) * numParams);
        if (paramLabels == NULL) {
            return EGADS_MALLOC;
        }

        paramMags = EG_alloc(sizeof(double) * numParams);
        if (paramMags == NULL) {
            AIM_FREE(paramLabels);
            return EGADS_MALLOC;
        }

        paramIndex = 0;

        // Rigid Constraints
        for (i = 0; i < feaAnalysis->numRigidConstraint; i++) {
            paramLabels[paramIndex] = feaAnalysis->rigidConstraint[i];
            paramMags[paramIndex++] = feaAnalysis->magRigidConstraint[i];
        }

        // Control Constraints
        for (i = 0; i < feaAnalysis->numControlConstraint; i++) {
            paramLabels[paramIndex] = feaAnalysis->controlConstraint[i];
            paramMags[paramIndex++] = feaAnalysis->magControlConstraint[i];
        }

        // Mach number if set
        if (feaAnalysis->machNumber != NULL && feaAnalysis->machNumber[0] > 0.0) {
            mach = feaAnalysis->machNumber;
        }

        status = nastranCard_trim(
            fp,
            &feaAnalysis->analysisID, // id
            mach, // mach
            &feaAnalysis->dynamicPressure, // q
            numParams,
            paramLabels, // label
            paramMags, // ux
            feaFileFormat->fileType
        );

        AIM_FREE(paramLabels);
        AIM_FREE(paramMags);

        if (status != CAPS_SUCCESS) return status;
    }

    if (feaAnalysis->analysisType == AeroelasticFlutter) {

        fprintf(fp,"%s\n","$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|");
        // MKAERO1, FLUTTER, FLFACT for density, mach and velocity

        // Write MKAERO1 INPUT
        status = nastranCard_mkaero1(
            fp,
            feaAnalysis->numMachNumber, feaAnalysis->machNumber, // m
            feaAnalysis->numReducedFreq, feaAnalysis->reducedFreq, // k
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        //fprintf(fp,"%s","$LUTTER SID     METHOD  DENS    MACH    RFREQ     IMETH   NVAL/OMAX   EPS   CONT\n");

        analysisID = 100 + feaAnalysis->analysisID;
        densityID = 10 * feaAnalysis->analysisID + 1;
        machID = 10 * feaAnalysis->analysisID + 2;
        velocityID = 10 * feaAnalysis->analysisID + 3;

        // Write FLUTTER INPUT
        status = nastranCard_flutter(
            fp,
            &analysisID, // sid
            "PK", // method
            &densityID, // dens
            &machID, // mach
            &velocityID, // rfreq
            "L", // imeth
            &feaAnalysis->numDesiredEigenvalue, // nvalue
            NULL, // eps
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        fprintf(fp,"$ DENSITY\n");
        status = nastranCard_flfact(
            fp,
            &densityID,
            1, &feaAnalysis->density,
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        fprintf(fp,"$ MACH\n");
        status = nastranCard_flfact(
            fp,
            &machID,
            feaAnalysis->numMachNumber, feaAnalysis->machNumber,
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        fprintf(fp,"$ DYANMIC PRESSURE=%f\n", feaAnalysis->dynamicPressure);

        // Setup the default flutter velocities  if not specified
        status = fea_defaultFlutterVelocity(feaAnalysis);
        if (status != CAPS_SUCCESS) return status;

        if (feaAnalysis->visualFlutter == (int) true) {
            for ( i = 0; i < feaAnalysis->numFlutterVel; i++) {
                feaAnalysis->flutterVel[i] = feaAnalysis->flutterVel[i] * -1;
            }
        }

        fprintf(fp,"$ VELOCITY\n");
        status = nastranCard_flfact(fp,
                                    &velocityID,
                                    feaAnalysis->numFlutterVel,
                                    feaAnalysis->flutterVel,
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    if (feaAnalysis->analysisType == Gust) {

        fprintf(fp,"%s\n","$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|");
        // MKAERO1, FLUTTER, FLFACT for density, mach and velocity

        // Write MKAERO1 INPUT
        status = nastranCard_mkaero1(
            fp,
            feaAnalysis->numMachNumber, feaAnalysis->machNumber, // m
            feaAnalysis->numReducedFreq, feaAnalysis->reducedFreq, // k
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        //fprintf(fp,"%s","$LUTTER SID     METHOD  DENS    MACH    RFREQ     IMETH   NVAL/OMAX   EPS   CONT\n");

        analysisID = 100 + feaAnalysis->analysisID;
        densityID = 10 * feaAnalysis->analysisID + 1;
        machID = 10 * feaAnalysis->analysisID + 2;
        velocityID = 10 * feaAnalysis->analysisID + 3;

//        // Write FLUTTER INPUT
//        status = nastranCard_flutter(
//            fp,
//            &analysisID, // sid
//            "PK", // method
//            &densityID, // dens
//            &machID, // mach
//            &velocityID, // rfreq
//            "L", // imeth
//            &feaAnalysis->numDesiredEigenvalue, // nvalue
//            NULL, // eps
//            feaFileFormat->fileType
//        );
//        if (status != CAPS_SUCCESS) return status;
//
//        fprintf(fp,"$ DENSITY\n");
//        status = nastranCard_flfact(
//            fp,
//            &densityID,
//            1, &feaAnalysis->density,
//            feaFileFormat->fileType
//        );
//        if (status != CAPS_SUCCESS) return status;

        fprintf(fp,"$ MACH\n");
        status = nastranCard_flfact(
            fp,
            &machID,
            feaAnalysis->numMachNumber, feaAnalysis->machNumber,
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;
//
//        fprintf(fp,"$ DYANMIC PRESSURE=%f\n", feaAnalysis->dynamicPressure);
//
//        // Setup the default flutter velocities  if not specified
//        status = fea_defaultFlutterVelocity(feaAnalysis);
//        if (status != CAPS_SUCCESS) return status;
//
//        if (feaAnalysis->visualFlutter == (int) true) {
//            for ( i = 0; i < feaAnalysis->numFlutterVel; i++) {
//                feaAnalysis->flutterVel[i] = feaAnalysis->flutterVel[i] * -1;
//            }
//        }
//
//        fprintf(fp,"$ VELOCITY\n");
//        status = nastranCard_flfact(fp,
//                                    &velocityID,
//                                    feaAnalysis->numFlutterVel,
//                                    feaAnalysis->flutterVel,
//                                    feaFileFormat->fileType);
//        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write a combined Nastran design constraint card from a set of constraint IDs.
//  The combined design constraint ID is set through the constraint ID variable.
int nastran_writeDesignConstraintADDCard(FILE *fp, int constraintID, int numSetID, const int designConstraintSetID[], const feaFileFormatStruct *feaFileFormat) {

    if (fp == NULL) return CAPS_IOERR;
    if (numSetID > 0 && designConstraintSetID == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (numSetID > 0) {
        return nastranCard_dconadd(
            fp,
            &constraintID, // dcid
            numSetID, designConstraintSetID, // dc
            feaFileFormat->fileType
        );
    }
    else {
        return CAPS_SUCCESS;
    }
}

// Write design constraint/optimization information from a feaDesignConstraint structure
int nastran_writeDesignConstraintCard(FILE *fp,
                                      const feaDesignConstraintStruct *feaDesignConstraint,
                                      const feaProblemStruct *feaProblem,
                                      const feaFileFormatStruct *feaFileFormat) {

    int status;

    int  i, j, found; // Index

    int elementStressLocation[4], drespID, dresp2ID, responseAttrB;
    int modeID, densityID, machID, velocityID, *atti_array[4];
    char label[9];

    int axialStressCode = 2; //torsionalStressCode = 4;
    int failureCriterionCode = 5;
    int numEquationResponseID = 0, *equationResponseIDSet = NULL;

    char tempString[16];
    char *valstr;

    if (fp == NULL) return CAPS_IOERR;
    if (feaDesignConstraint == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (feaDesignConstraint->designConstraintType == PropertyDesignCon) {

        for (i = 0; i < feaDesignConstraint->numPropertyID; i++) {

            if (feaDesignConstraint->propertySetType[i] == Rod) {

                drespID = feaDesignConstraint->designConstraintID + 10000;

                status = nastranCard_dconstr(
                    fp,
                    &feaDesignConstraint->designConstraintID,
                    &drespID,
                    &feaDesignConstraint->lowerBound,
                    &feaDesignConstraint->upperBound,
                    feaFileFormat->fileType
                );
                if (status != CAPS_SUCCESS) return status;

                //$...STRUCTURAL RESPONSE IDENTIFICATION
                //$DRESP1 ID      LABEL   RTYPE   PTYPE   REGION  ATTA    ATTB    ATT1    +
                //$+      ATT2    ...

                // ------------- STRESS RESPONSE ---------------------------------------------------------------------
                if (strcmp(feaDesignConstraint->responseType,"STRESS") == 0) {

                    // make label field value
                    status = convert_integerToString(drespID, 6, 0, tempString);
                    if (status != CAPS_SUCCESS) return status;
                    snprintf(label, 9, "R%s", tempString);

                    status = nastranCard_dresp1(
                        fp,
                        &drespID, // id
                        label, // label
                        feaDesignConstraint->responseType, // rtype
                        "PROD", // ptype
                        NULL, // region
                        Integer, &axialStressCode, // atta
                        Integer, NULL, // attb
                        Integer, 1, &feaDesignConstraint->propertySetID[i], // att
                        (int) false, // blankIntegerMode
                        feaFileFormat->fileType
                    );
                    if (status != CAPS_SUCCESS) return status;
                }
            } else if (feaDesignConstraint->propertySetType[i] == Bar) {
                // Nothing set yet
            } else if (feaDesignConstraint->propertySetType[i] == Shell) {

                // ------------- STRESS RESPONSE ---------------------------------------------------------------------
                if (strcmp(feaDesignConstraint->responseType,"STRESS") == 0) {
                    // elementStressLocation[0] = 7; // Major principal at Z1
                    elementStressLocation[0] = 9; // Von Mises (or max shear) at Z1
                    // elementStressLocation[2] = 16; // Minor principal at Z2
                    // elementStressLocation[3] = 17; // von Mises (or max shear) at Z2

                    for (j = 0; j < 1; j++) {

                        // DCONSTR
                        drespID = feaDesignConstraint->designConstraintID + 10000 + j*1000;

                        status = nastranCard_dconstr(
                            fp,
                            &feaDesignConstraint->designConstraintID,
                            &drespID,
                            &feaDesignConstraint->lowerBound,
                            &feaDesignConstraint->upperBound,
                            feaFileFormat->fileType
                        );
                        if (status != CAPS_SUCCESS) return status;

                        // DRESP1

                        // make label field value
                        status = convert_integerToString(drespID, 6, 0, tempString);
                        if (status != CAPS_SUCCESS) return status;
                        snprintf(label, 9, "R%s", tempString);

                        status = nastranCard_dresp1(
                            fp,
                            &drespID, // id
                            label, // label
                            feaDesignConstraint->responseType, // rtype
                            "PSHELL", // ptype
                            NULL, // region
                            Integer, &elementStressLocation[j], // atta
                            Integer, NULL, // attb
                            Integer, 1, &feaDesignConstraint->propertySetID[i], // att
                            (int) false, // blankIntegerMode
                            feaFileFormat->fileType
                        );
                        if (status != CAPS_SUCCESS) return status;
                    }
                }

            } else if (feaDesignConstraint->propertySetType[i] == Composite) {

                // ------------- CFAILURE RESPONSE ---------------------------------------------------------------------
                if (strcmp(feaDesignConstraint->responseType,"CFAILURE") == 0) {

                    // DCONSTR

                    drespID = feaDesignConstraint->designConstraintID + 10000;

                    status = nastranCard_dconstr(
                        fp,
                        &feaDesignConstraint->designConstraintID,
                        &drespID,
                        &feaDesignConstraint->lowerBound,
                        &feaDesignConstraint->upperBound,
                        feaFileFormat->fileType
                    );
                    if (status != CAPS_SUCCESS) return status;

                    // DRESP1

                    // make label field value
                    status = convert_integerToString(drespID, 6, 0, tempString);
                    if (status != CAPS_SUCCESS) return status;
                    snprintf(label, 9, "L%s", tempString);

                     if (feaDesignConstraint->fieldPosition == 0) {
                        found = 0;
                        // OPTIONS ARE "Ti", "THETAi", "LAMINAi" all = Integer i
                        valstr = NULL;
                        // skan for T, THETA or LAMINA

                        valstr = strstr(feaDesignConstraint->fieldName, "THETA");
                        if (valstr != NULL) {
                            // little trick to get integer value of character digit, i.e. '1' - '0' = 1
                            responseAttrB = feaDesignConstraint->fieldName[5] - '0';
                            found = 1;
                        }

                        if (found == 0) {
                            valstr = strstr(feaDesignConstraint->fieldName, "LAMINA");
                            if (valstr != NULL) {
                                responseAttrB = feaDesignConstraint->fieldName[6] - '0';
                                found = 1;
                            }
                        }

                        if (found == 0) {
                            valstr = strstr(feaDesignConstraint->fieldName, "T");
                            if (valstr != NULL) {
                                responseAttrB = feaDesignConstraint->fieldName[1] - '0';
                                found = 1;
                            }
                        }

                        if (found == 0) {
                            printf("  WARNING: Could not determine what Lamina to apply constraint too, using default = 1\n");
                            printf("  String Entered: %s\n", feaDesignConstraint->fieldName);
                            responseAttrB = 1;
                        }
                    } else {
                        responseAttrB = feaDesignConstraint->fieldPosition;
                    }

                    status = nastranCard_dresp1(
                        fp,
                        &drespID, // id
                        label, // label
                        feaDesignConstraint->responseType, // rtype
                        "PCOMP", // ptype
                        NULL, // region
                        Integer, &failureCriterionCode, // atta
                        Integer, &responseAttrB, // attb
                        Integer, 1, &feaDesignConstraint->propertySetID[i], // att
                        (int) false, // blankIntegerMode
                        feaFileFormat->fileType
                    );
                    if (status != CAPS_SUCCESS) return status;
                }

            } else if (feaDesignConstraint->propertySetType[i] == Solid) {
                // Nothing set yet
            }
        }

    } else if (feaDesignConstraint->designConstraintType == FlutterDesignCon) {

        drespID = feaDesignConstraint->designConstraintID + 10000;
        dresp2ID = feaDesignConstraint->designConstraintID + 20000;

        if (feaDesignConstraint->numModes > 0) {
            modeID = 1000 + (10 * feaDesignConstraint->designConstraintID + 1);
        } else {
            modeID = 0;
        }

        densityID = 1000 + (10 * feaDesignConstraint->designConstraintID + 2);
        machID = 1000 + (10 * feaDesignConstraint->designConstraintID + 3);
        velocityID = 1000 + (10 * feaDesignConstraint->designConstraintID + 4);

        if (modeID == 0) {
            atti_array[0] = NULL;

        } else {
            atti_array[0] = &modeID;
        }

        atti_array[1] = &densityID;
        atti_array[2] = &machID;
        atti_array[3] = &velocityID;

        status = nastranCard_dconstr(fp,
                                     &feaDesignConstraint->designConstraintID,
                                     &dresp2ID,
                                     &feaDesignConstraint->lowerBound,
                                     &feaDesignConstraint->upperBound,
                                     feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        // Flutter design response
        status = nastranCard_dresp1(fp,
                                    &drespID, // id
                                    "FLUTTER", // label
                                    "FLUTTER", // rtype
                                    NULL, // ptype
                                    NULL, // region
                                    Integer, NULL, // atta
                                    Integer, NULL, // attb
                                    Integer, 4, atti_array, // atti
                                    (int) true, // blankIntegerMode
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        char *equationLines = "F(A) = (A-0.03)/0.1";
        int numEquationLines = 1;
        int equationID = 1000 + (10 * feaDesignConstraint->designConstraintID + 5);

        status = nastranCard_deqatn(fp,
                                    &equationID,
                                    numEquationLines,
                                    &equationLines);
        if (status != CAPS_SUCCESS) return status;

        status = nastranCard_dresp2(fp,
                                    &dresp2ID, // id
                                    "GDAMP", // label
                                    &equationID, // eqid
                                    NULL, //region
                                    0, NULL, // dvid
                                    0, NULL, // labl
                                    1, &drespID, // nr
                                    0, NULL, // g
                                    NULL, // c
                                    0, NULL, // nrr
                                    feaFileFormat->fileType);

        if (status != CAPS_SUCCESS) return status;

        // Modes for flutter constraint
        if (feaDesignConstraint->numModes > 0) {
            status = nastranCard_set1(fp,
                                      &modeID,
                                      feaDesignConstraint->numModes,
                                      feaDesignConstraint->modes,
                                      feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }

        // Density for flutter constraint
        status = nastranCard_flfact(fp,
                                    &densityID,
                                    feaDesignConstraint->numDensity,
                                    feaDesignConstraint->density,
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        // Mach number for flutter constraint
        status = nastranCard_flfact(fp,
                                    &machID,
                                    feaDesignConstraint->numMach,
                                    feaDesignConstraint->Mach,
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        // Velocities for flutter constraint
        status = nastranCard_flfact(fp,
                                    &velocityID,
                                    feaDesignConstraint->numVelocity,
                                    feaDesignConstraint->velocity,
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;


    } else if (feaDesignConstraint->designConstraintType == EquationResponseDesignCon) {
        // DRESP2
        status = _getEquationResponseIDSet(feaProblem,
                                           1,
                                           &feaDesignConstraint->responseType,
                                           &numEquationResponseID,
                                           &equationResponseIDSet);
        if (status != CAPS_SUCCESS) return status;

        status = nastranCard_dconstr(fp,
                                     &feaDesignConstraint->designConstraintID,
                                     equationResponseIDSet,
                                     &feaDesignConstraint->lowerBound,
                                     &feaDesignConstraint->upperBound,
                                     feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    if (equationResponseIDSet != NULL) {
        EG_free(equationResponseIDSet);
    }

    return CAPS_SUCCESS;
}


// get element type value for TYPE field in DVCREL* card
/*@null@*/
static char * _getElementTypeIdentifier(int elementType, int elementSubType) {

    char *identifier = NULL;

    if (elementType == Node) {

        if (elementSubType == ConcentratedMassElement) {
            identifier = "CONM2";
        }
        // else {
        //     identifier = "CONM2"; // default
        // }
    }

    if (elementType == Line) {

        if (elementSubType == BarElement) {
            identifier = "CBAR";
        }
        else if (elementSubType == BeamElement) {
            // beam elements not supported yet
            identifier = NULL;
        }
        else {
            identifier = "CROD"; // default
        }
    }

    if (elementType == Triangle) {

        if (elementSubType == ShellElement) {
            identifier = "CTRIA3";
        }
        else {
            identifier = "CTRIA3"; // default
        }
    }

    if (elementType == Triangle_6) {

        if (elementSubType == ShellElement) {
            identifier = "CTRIA6";
        }
        else {
            identifier = "CTRIA6"; // default
        }
    }

    if (elementType == Quadrilateral) {

        if (elementSubType == ShearElement) {
            identifier = "CSHEAR";
        }
        else if (elementSubType == ShellElement) {
            identifier = "CQUAD4";
        }
        else {
            identifier = "CQUAD4"; // default
        }
    }

    if (elementType == Quadrilateral_8) {

        if (elementSubType == ShellElement) {
            identifier = "CQUAD8";
        }
        else {
            identifier = "CQUAD8"; // default
        }
    }

    if (identifier != NULL) {
        return EG_strdup(identifier);
    }
    else {
        return NULL;
    }
}

// Write design variable/optimization information from a feaDesignVariable structure
int nastran_writeDesignVariableCard(void *aimInfo, FILE *fp, const feaDesignVariableStruct *feaDesignVariable, const feaFileFormatStruct *feaFileFormat) {

    // int  i;

    int status; // Function return status

    const int *ddval = NULL;
    const double *xlb = NULL, *xub = NULL, *delxv = NULL;
    // char *type = NULL;

    int dlinkID;//, uniqueID;

    // char *pname = NULL, *cname = NULL;

    if (fp == NULL) return CAPS_IOERR;
    if (feaDesignVariable == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (feaDesignVariable->numDiscreteValue == 0) {
        xlb = &feaDesignVariable->lowerBound;
        xub = &feaDesignVariable->upperBound;
        delxv = &feaDesignVariable->maxDelta;
    }
    else {
        ddval = &feaDesignVariable->designVariableID;
    }

    status = nastranCard_desvar(
        fp,
        &feaDesignVariable->designVariableID, // id
        feaDesignVariable->name, // label
        &feaDesignVariable->initialValue, // xinit
        xlb, // xlb
        xub, // xub
        delxv, // delxv
        ddval, // ddval
        feaFileFormat->fileType
    );
    AIM_STATUS(aimInfo, status);

    if (ddval != NULL) {

        // Write DDVAL card
        status = nastranCard_ddval(
            fp,
            ddval, // id
            feaDesignVariable->numDiscreteValue,
            feaDesignVariable->discreteValue, // dval
            feaFileFormat->fileType
        );
        AIM_STATUS(aimInfo, status);
    }

    if (feaDesignVariable->numIndependVariable > 0) {

        dlinkID = feaDesignVariable->designVariableID + 10000;

        status = nastranCard_dlink(
            fp,
            &dlinkID, // id
            &feaDesignVariable->designVariableID, // ddvid
            &feaDesignVariable->variableWeight[0], // c0
            &feaDesignVariable->variableWeight[1], // cmult
            feaDesignVariable->numIndependVariable,
            feaDesignVariable->independVariableID, // idv
            feaDesignVariable->independVariableWeight, // c
            feaFileFormat->fileType
        );
        AIM_STATUS(aimInfo, status);
    }

cleanup:
    return status;
}

// Write design variable relation information from a feaDesignVariableRelation structure
int nastran_writeDesignVariableRelationCard(void *aimInfo,
                                            FILE *fp,
                                            const feaDesignVariableRelationStruct *feaDesignVariableRelation,
                               /*@unused@*/ const feaProblemStruct *feaProblem,
                                            const feaFileFormatStruct *feaFileFormat)
{
    int status = CAPS_SUCCESS;
    int i, j, uniqueID;

    int numDesignVariable, *designVariableSetID = NULL;
    feaDesignVariableStruct **designVariableSet = NULL, *designVariable=NULL;

    int numRelation, relationIndex;
    int *relationSetID = NULL, *relationSetType = NULL, *relationSetSubType = NULL;

    char *type = NULL, *fieldName = NULL;

    if (fp == NULL) return CAPS_IOERR;
    if (feaDesignVariableRelation == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    numDesignVariable = feaDesignVariableRelation->numDesignVariable;
    designVariableSet = feaDesignVariableRelation->designVariableSet;

    designVariableSetID = EG_alloc(sizeof(int) * numDesignVariable);
    if (designVariableSetID == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    numRelation = 0;
    relationIndex = 0;
    for (i = 0; i < numDesignVariable; i++) {

        designVariable = designVariableSet[i];

        designVariableSetID[i] = designVariable->designVariableID;

        if (feaDesignVariableRelation->componentType == MaterialDesignVar) {

            if (numRelation == 0) {
                numRelation = feaDesignVariableRelation->numMaterialID;
                relationSetID = EG_alloc(numRelation * sizeof(int));
                relationSetType = EG_alloc(numRelation * sizeof(int));
            }
            else {
                numRelation += feaDesignVariableRelation->numMaterialID;
                relationSetID = EG_reall(relationSetID, numRelation * sizeof(int));
                relationSetType = EG_reall(relationSetType, numRelation * sizeof(int));
            }

            if (relationSetID == NULL || relationSetType == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (j = 0; j < feaDesignVariableRelation->numMaterialID; j++) {

                relationSetID[relationIndex] = feaDesignVariableRelation->materialSetID[j];
                relationSetType[relationIndex] = feaDesignVariableRelation->materialSetType[j];
                relationIndex++;
            }
        }

        else if (feaDesignVariableRelation->componentType == PropertyDesignVar) {

            if (numRelation == 0) {
                numRelation = feaDesignVariableRelation->numPropertyID;
                relationSetID = EG_alloc(numRelation * sizeof(int));
                relationSetType = EG_alloc(numRelation * sizeof(int));
            }
            else {
                numRelation += feaDesignVariableRelation->numPropertyID;
                relationSetID = EG_reall(relationSetID, numRelation * sizeof(int));
                relationSetType = EG_reall(relationSetType, numRelation * sizeof(int));
            }

            if (relationSetID == NULL || relationSetType == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (j = 0; j < feaDesignVariableRelation->numPropertyID; j++) {
                relationSetID[relationIndex] = feaDesignVariableRelation->propertySetID[j];
                relationSetType[relationIndex] = feaDesignVariableRelation->propertySetType[j];
                relationIndex++;
            }
        }

        else if (feaDesignVariableRelation->componentType == ElementDesignVar) {

            if (numRelation == 0) {
                numRelation = feaDesignVariableRelation->numElementID;
                relationSetID = EG_alloc(numRelation * sizeof(int));
                relationSetType = EG_alloc(numRelation * sizeof(int));
                relationSetSubType = EG_alloc(numRelation * sizeof(int));
            }
            else {
                numRelation += feaDesignVariableRelation->numElementID;
                relationSetID = EG_reall(relationSetID, numRelation * sizeof(int));
                relationSetType = EG_reall(relationSetType, numRelation * sizeof(int));
                relationSetSubType = EG_reall(relationSetSubType, numRelation * sizeof(int));
            }

            if (relationSetID == NULL || relationSetType == NULL || relationSetSubType == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (j = 0; j < feaDesignVariableRelation->numElementID; j++) {

                relationSetID[relationIndex] = feaDesignVariableRelation->elementSetID[j];
                relationSetType[relationIndex] = feaDesignVariableRelation->elementSetType[j];
                relationSetSubType[relationIndex] = feaDesignVariableRelation->elementSetSubType[j];
                relationIndex++;
            }
        }
    }

    // get *PNAME field value
    if (feaDesignVariableRelation->fieldPosition == 0) {
        AIM_STRDUP(fieldName, feaDesignVariableRelation->fieldName, aimInfo, status);
    } else {
        AIM_ALLOC(fieldName, 8, char, aimInfo, status);
        status = convert_integerToString(feaDesignVariableRelation->fieldPosition, 7, 1, fieldName);
        AIM_STATUS(aimInfo, status);
    }

    if (feaDesignVariableRelation->componentType == MaterialDesignVar) {

        for (i = 0; i < numRelation; i++) {

            uniqueID = feaDesignVariableRelation->relationID * 100 + i;

            // UnknownMaterial, Isotropic, Anisothotropic, Orthotropic, Anisotropic
            if (relationSetType[i] == Isotropic) {
                type = "MAT1";
            } else if (relationSetType[i] == Anisothotropic) {
                type = "MAT2";
            } else if (relationSetType[i] == Orthotropic) {
                type = "MAT8";
            } else if (relationSetType[i] == Anisotropic) {
                type = "MAT9";
            } else {
                AIM_ERROR(aimInfo, "Unknown material type: %d", relationSetType[i]);
                status = CAPS_NOTIMPLEMENT;
                goto cleanup;
            }

            if (feaDesignVariableRelation->isEquation == (int)true)
            {
              status = nastranCard_dvmrel1(
                  fp,
                  &uniqueID, // id
                  type, // type
                  &relationSetID[i], // mid
                  fieldName, // mpname
                  NULL, // mpmin
                  NULL, // mpmax
                  &feaDesignVariableRelation->constantRelationCoeff, // c0
                  numDesignVariable,
                  designVariableSetID, // dvid
                  feaDesignVariableRelation->linearRelationCoeff, // coeff
                  feaFileFormat->fileType
              );
              AIM_STATUS(aimInfo, status);
            }
        }
    }

    else if (feaDesignVariableRelation->componentType == PropertyDesignVar) {

        for (i = 0; i < numRelation; i++) {

            uniqueID = feaDesignVariableRelation->relationID * 100 + i;

            // UnknownProperty, Rod, Bar, Shear, Shell, Composite, Solid
            if (relationSetType[i] == Rod) {
                type = "PROD";
            } else if (relationSetType[i] == Bar) {
                type = "PBAR";
            } else if (relationSetType[i] == Shell) {
                type = "PSHELL";
            } else if (relationSetType[i] == Shear) {
                type = "PSHEAR";
            } else if (relationSetType[i] == Composite) {
                type = "PCOMP";
            } else if (relationSetType[i] == Solid) {
                type = "PSOLID";
            } else {
                AIM_ERROR(aimInfo, "Unknown property type: %d", relationSetType[i]);
                status = CAPS_NOTIMPLEMENT;
                goto cleanup;
            }

            status = nastranCard_dvprel1(
                fp,
                &uniqueID, // id
                type, // type
                &relationSetID[i], // pid
                NULL, // fid
                fieldName, // pname
                NULL, // pmin
                NULL, // pmax
                &feaDesignVariableRelation->constantRelationCoeff, // c0
                numDesignVariable,
                designVariableSetID, // dvid
                feaDesignVariableRelation->linearRelationCoeff, // coeff
                feaFileFormat->fileType
            );

            AIM_STATUS(aimInfo, status);
        }
    }

    else if (feaDesignVariableRelation->componentType == ElementDesignVar) {

        for (i = 0; i < numRelation; i++) {

            uniqueID = feaDesignVariableRelation->relationID * 10000 + i;

            // Element types:    UnknownMeshElement, Node, Line, Triangle, Triangle_6, Quadrilateral, Quadrilateral_8, Tetrahedral, Tetrahedral_10, Pyramid, Prism, Hexahedral
            // Element subtypes: UnknownMeshSubElement, ConcentratedMassElement, BarElement, BeamElement, ShellElement, ShearElement
            type = _getElementTypeIdentifier(relationSetType[i], relationSetSubType[i]);
            if (type == NULL) {

                AIM_ERROR(aimInfo, "Unknown element type and/or subtype: %d %d",
                                   relationSetType[i],
                                   relationSetSubType[i]);
                status = CAPS_BADVALUE;
                goto cleanup;

            }

            status = nastranCard_dvcrel1(
                fp,
                &uniqueID, // id
                type, // type
                &relationSetID[i], // eid
                fieldName, // cpname
                NULL, // cpmin
                NULL, // cpmax
                &feaDesignVariableRelation->constantRelationCoeff, // c0
                numDesignVariable,
                designVariableSetID, // dvid
                feaDesignVariableRelation->linearRelationCoeff, // coeff
                feaFileFormat->fileType
            );

            AIM_FREE(type);

            AIM_STATUS(aimInfo, status);
        }
    }

    else {

        AIM_ERROR(aimInfo, "Unknown design variable relation type: %d",
                    feaDesignVariableRelation->componentType);
        status = CAPS_BADVALUE;
        goto cleanup;
    }


cleanup:

    AIM_FREE(fieldName);
    // AIM_FREE(designVariableSet);
    AIM_FREE(designVariableSetID);
    AIM_FREE(relationSetID);
    AIM_FREE(relationSetType);
    AIM_FREE(relationSetSubType);

    return status;
}

static int _getNextEquationLine(char **equationLines, const int lineIndex,
                                const char *equationString, const int lineMaxChar) {

    int numPrint, equationLength;

    char lineBuffer[80];

    equationLength = strlen(equationString);

    if (equationLength < lineMaxChar) {
        // full equation fits in first line
        numPrint = snprintf(lineBuffer, lineMaxChar, "%s;", equationString);
    }
    else {
        numPrint = snprintf(lineBuffer, lineMaxChar, "%s", equationString);
    }
    if (numPrint >= lineMaxChar) {
        numPrint = lineMaxChar - 1; // -1 due to NULL terminator
    }
    equationLines[lineIndex] = EG_strdup(lineBuffer);

    return numPrint;
}


static int _getEquationLines(const feaDesignEquationStruct *feaEquation,
                             int *numEquationLines,
                             char ***equationLines) {

    int i;

    int numLines, lineIndex, numPrint, equationLength;
    char *equationString, **lines;

    // conservative estimate number of lines required
    numLines = 0;
    for (i = 0; i < feaEquation->equationArraySize; i++) {
        equationString = feaEquation->equationArray[i];
        if (i == 0)
            numLines += 1 + ((strlen(equationString)+1) / 56);
        else
            numLines += 1 + ((strlen(equationString)+1) / 64);
    }

    if (numLines == 0) {
        PRINT_WARNING("Empty equation: %s", feaEquation->name);
        return CAPS_SUCCESS;
    }

    // assume equationLines is uninitialized
    // alloc enough lines for double number of chars to be conservative
    // 64 char per line
    lines = EG_alloc(sizeof(char *) * numLines);
    if (lines == NULL) {
        return EGADS_MALLOC;
    }
    lineIndex = 0;

    // first equation string

    // first line is 56 char
    equationString = feaEquation->equationArray[0];
    equationLength = strlen(equationString);

    numPrint = _getNextEquationLine(
        lines, lineIndex, equationString, 56);
    lineIndex++;

    // each cont line is 64 char
    // TODO: for readability, should each cont line start with 8 blank spaces?
    while (numPrint < equationLength) {

        numPrint += _getNextEquationLine(
            lines, lineIndex, equationString + numPrint, 64);
        lineIndex++;
    }

    // for each remaining equation string
    for (i = 1; i < feaEquation->equationArraySize; i++) {

        equationString = feaEquation->equationArray[i];
        equationLength = strlen(equationString);
        numPrint = 0;

        while (numPrint < equationLength) {

            numPrint += _getNextEquationLine(
                lines, lineIndex, equationString + numPrint, 64);
            lineIndex++;
        }
    }

    *numEquationLines = lineIndex;
    *equationLines = lines;

    return CAPS_SUCCESS;
}

// Write equation information from a feaDesignEquation structure
int nastran_writeDesignEquationCard(FILE *fp,
                                    const feaDesignEquationStruct *feaEquation,
                       /*@unused@*/ const feaFileFormatStruct *fileFormat) {

    int status;

    int numEquationLines = 0;
    char **equationLines = NULL;

    if (fp == NULL) return CAPS_IOERR;

    status = _getEquationLines(feaEquation, &numEquationLines, &equationLines);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = nastranCard_deqatn(fp, &feaEquation->equationID, numEquationLines, equationLines);

    cleanup:

        if (equationLines != NULL) {
            string_freeArray(numEquationLines, &equationLines);
        }

        return status;
}

// Write design table constants information from a feaDesignTable structure
int nastran_writeDesignTableCard(FILE *fp, const feaDesignTableStruct *feaDesignTable, const feaFileFormatStruct *fileFormat) {

    if (feaDesignTable->numConstant > 0) {

        return nastranCard_dtable(
            fp,
            feaDesignTable->numConstant,
            feaDesignTable->constantLabel, // labl
            feaDesignTable->constantValue, // valu
            fileFormat->fileType
        );
    }
    else {
        return CAPS_SUCCESS;
    }

}

int nastran_writeMasssetCard(FILE *fp,
                             int incrementID,
                             const feaFileFormatStruct *fileFormat) {

    int status;
    int masssetID;

    masssetID = 100 + incrementID;
    status = nastranCard_massset(fp,
                                 &masssetID, // eid
                                 &incrementID,        // id2
                                 fileFormat->gridFileType);
    return status;
}

int nastran_writeMassIncrementSet(FILE *fp,
                                  int incrementID,
                                  int numIncrementID,
                                  const int massIncrementID[],
                                  const feaMassIncrementStruct *feaMassIncrement,
                                  const feaFileFormatStruct *fileFormat) {

    int eid, i, status = CAPS_SUCCESS;

    fprintf(fp,"%s\n","$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|");
    fprintf(fp,"begin massid=%d label='subcase %d'\n", incrementID, incrementID);

    for (i = 0; i < numIncrementID; i++) {
//        printf("ID: %d %d\n", massIncrementID[i], feaMassIncrement[massIncrementID[i]-1].propertySetID[0]);
        eid = incrementID*10000 + massIncrementID[i];
        status = nastranCard_conm2(fp,
                                   &eid, // eid
                                   &feaMassIncrement[massIncrementID[i]-1].propertySetID[0], // g
                                   NULL, // cid
                                   &feaMassIncrement[massIncrementID[i]-1].increment, // m
                                   NULL, // x
                                   NULL, // i
                                   fileFormat->gridFileType);
    }

    return status;
}

// Write design response type "DISP"
static int _writeDesignResponseDISP(FILE *fp, const feaDesignResponseStruct *feaDesignResponse, const feaFileFormatStruct *fileFormat) {

    int status;

    int drespID = 100000 + feaDesignResponse->responseID;

    status = nastranCard_dresp1(
        fp,
        &drespID,
        feaDesignResponse->name, // label
        feaDesignResponse->responseType, // rtype
        NULL, // ptype
        NULL, // region
        Integer, &feaDesignResponse->component, // atta
        Integer, &feaDesignResponse->attb, // attb
        Integer, 1, &feaDesignResponse->gridID, // atti
        (int) false, // blankIntegerMode
        fileFormat->fileType
    );

    return status;
}

//// Write design response type "FREQ"
//static int _writeDesignResponseFREQ(FILE *fp, const feaDesignResponseStruct *feaDesignResponse, const feaFileFormatStruct *fileFormat) {
//
//    int status;
//
//    int drespID = 100000 + feaDesignResponse->responseID;
//
//    status = nastranCard_dresp1(
//        fp,
//        &drespID,
//        feaDesignResponse->name, // label
//        feaDesignResponse->responseType, // rtype
//        NULL, // ptype
//        NULL, // region
//        Integer, &feaDesignResponse->component, // atta
//        Integer, NULL, // attb
//        Integer, 1, NULL, // atti
//        fileFormat->fileType
//    );
//
//    return status;
//}

// Write design response information from a feaDesignResponse structure
int nastran_writeDesignResponseCard(FILE *fp, const feaDesignResponseStruct *feaDesignResponse, const feaFileFormatStruct *fileFormat) {

    const char *responseType = feaDesignResponse->responseType;

    if (strcmp(responseType, "DISP") == 0) {
        return _writeDesignResponseDISP(fp, feaDesignResponse, fileFormat);
    }
//    else if (strcmp(responseType, "FREQ") == 0) {
//        return _writeDesignResponseFREQ(fp, feaDesignResponse, fileFormat);
//    }
    else {
        PRINT_ERROR("Unknown responseType: %s", responseType);
        return CAPS_BADVALUE;
    }
}

// Write design equation response information from a feaDesignEquationResponse structure
int nastran_writeDesignEquationResponseCard(FILE *fp,
                                            const feaDesignEquationResponseStruct *feaEquationResponse,
                                            const feaProblemStruct *feaProblem,
                                            const feaFileFormatStruct *fileFormat) {

    int status;

    int drespID, equationID;
    int numDesignVariableID=0, *designVariableIDSet = NULL;
    int numConstant = 0;
    char **constantLabelSet = NULL;
    int numResponseID = 0, *responseIDSet = NULL;
    int numGrid = 0, *gridIDSet = NULL, *dofNumberSet = NULL;
    int numEquationResponseID = 0, *equationResponseIDSet = NULL;

    status = _getEquationID(feaProblem, feaEquationResponse->equationName, &equationID);
    if (status != CAPS_SUCCESS) {
        PRINT_ERROR("Unable to get equation ID for name: %s - status: %d",
                    feaEquationResponse->equationName,
                    status);
        goto cleanup;
    }

    // DESVAR
    status = _getDesignVariableIDSet(
        feaProblem,
        feaEquationResponse->numDesignVariable, feaEquationResponse->designVariableNameSet,
        &numDesignVariableID, &designVariableIDSet
    );
    if (status != CAPS_SUCCESS) goto cleanup;

    // DTABLE
    numConstant = feaEquationResponse->numConstant;
    constantLabelSet = feaEquationResponse->constantLabelSet;

    // DRESP1
    status = _getDesignResponseIDSet(
        feaProblem,
        feaEquationResponse->numResponse, feaEquationResponse->responseNameSet,
        &numResponseID, &responseIDSet
    );
    if (status != CAPS_SUCCESS) goto cleanup;

    // DNODE
    numGrid = 0;
    gridIDSet = NULL;
    dofNumberSet = NULL;

    // DRESP2
    status = _getEquationResponseIDSet(
        feaProblem,
        feaEquationResponse->numEquationResponse, feaEquationResponse->equationResponseNameSet,
        &numEquationResponseID, &equationResponseIDSet
    );
    if (status != CAPS_SUCCESS) goto cleanup;

    drespID = 200000 + feaEquationResponse->equationResponseID;

    status = nastranCard_dresp2(
        fp,
        &drespID, // id
        feaEquationResponse->name, // label
        &equationID, // eqid
        NULL, //region
        numDesignVariableID, designVariableIDSet, // dvid
        numConstant, constantLabelSet, // labl
        numResponseID, responseIDSet, // nr
        numGrid, gridIDSet, dofNumberSet, // g, c
        numEquationResponseID, equationResponseIDSet, // nrr
        fileFormat->fileType
    );

    cleanup:

        if (designVariableIDSet != NULL) {
            EG_free(designVariableIDSet);
        }
        if (responseIDSet != NULL) {
            EG_free(responseIDSet);
        }
        if (equationResponseIDSet != NULL) {
            EG_free(equationResponseIDSet);
        }

        return status;
}

// Write design optimization parameter information from a feaDesignOptParam structure
int nastran_writeDesignOptParamCard(FILE *fp, const feaDesignOptParamStruct *feaDesignOptParam, const feaFileFormatStruct *fileFormat) {

    if (feaDesignOptParam->numParam > 0) {

        return nastranCard_doptprm(
            fp,
            feaDesignOptParam->numParam,
            feaDesignOptParam->paramLabel, // param
            feaDesignOptParam->paramType,
            feaDesignOptParam->paramValue, // val
            fileFormat->fileType
        );
    }
    else {
        return CAPS_SUCCESS;
    }

}

// Read data from a Nastran F06 file to determine the number of eignevalues
int nastran_readF06NumEigenValue(FILE *fp, int *numEigenVector) {

    int status; // Function return

    char *beginEigenLine = "                                              R E A L   E I G E N V A L U E S";
    char *endEigenLine = "1";
    int keepCollecting = (int) true;

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    int tempInt[2];
    double tempDouble[5];

    *numEigenVector = 0;

    if (fp == NULL) return CAPS_IOERR;

    while (*numEigenVector == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if (status < 0) break;

        // See how many Eigen-Values we have
        if (strncmp(beginEigenLine, line,strlen(beginEigenLine)) == 0) {

            // Skip ahead 2 lines
            status = getline(&line, &linecap, fp);
            if (status < 0) break;
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            while (keepCollecting == (int) true) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;

                if (strncmp(endEigenLine, line, strlen(endEigenLine)) == 0) {
                    keepCollecting = (int) false;
                    break;
                }

                status = sscanf(line, "%d%d%lf%lf%lf%lf%lf", &tempInt[0],
                                                             &tempInt[1],
                                                             &tempDouble[0],
                                                             &tempDouble[1],
                                                             &tempDouble[2],
                                                             &tempDouble[3],
                                                             &tempDouble[4]);

                if (status != 7) {
                    keepCollecting = (int) false;
                    break;
                }

                *numEigenVector += 1;
            }
        }
    }

    if (line != NULL) free(line);

    // Rewind the file
    rewind(fp);

    if (*numEigenVector == 0) return CAPS_NOTFOUND;
    else return CAPS_SUCCESS;
}

// Read data from a Nastran F06 file and load it into a dataMatrix[numEigenVector][numGridPoint*8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int nastran_readF06EigenVector(FILE *fp, int *numEigenVector, int *numGridPoint, double ***dataMatrix) {

    int status = 0; // Function return

    int i, j, eigenValue = 0; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *beginEigenLine = "      EIGENVALUE =";
    char *endEigenLine = "1";
    char tempString[30];

    int numVariable = 8; // Grid Id, Coord Id, T1, T2, T3, R1, R2, R3

    printf("Reading Nastran FO6 file - extracting Eigen-Vectors!\n");

    *numEigenVector = 0;
    *numGridPoint = 0;

    if (fp == NULL) return CAPS_IOERR;

    // See how many Eigen-Values we have
    status = nastran_readF06NumEigenValue(fp, numEigenVector);
    // *numEigenVector = 10;
    printf("\tNumber of Eigen-Vectors = %d\n", *numEigenVector);
    if (status != CAPS_SUCCESS) return status;

    // Loop through file line by line until we have determined how many Eigen-Values and grid points we have
    while (*numGridPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if (status < 0) break;

        // Look for start of Eigen-Vector 1
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            // Fast forward 3 lines
            for (i = 0; i < 3; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through lines counting the number of grid points
            while (getline(&line, &linecap, fp) >= 0) {

                if (strncmp(endEigenLine, line, strlen(endEigenLine)) == 0) break;

                *numGridPoint +=1;
            }
        }
    }

    if (line != NULL) EG_free(line);
    line = NULL;

    printf("\tNumber of Grid Points = %d for each Eigen-Vector\n", *numGridPoint);
    if (*numGridPoint == 0) return CAPS_NOTFOUND;

    // Rewind the file
    rewind(fp);

    // Allocate dataMatrix array
    if (*dataMatrix != NULL) EG_free(*dataMatrix);

    *dataMatrix = (double **) EG_alloc(*numEigenVector *sizeof(double *));
    if (*dataMatrix == NULL) return EGADS_MALLOC; // If allocation failed ....

    for (i = 0; i < *numEigenVector; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc((*numGridPoint)*numVariable*sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

                if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
            }

            if ((*dataMatrix) != NULL) EG_free((*dataMatrix));

            return EGADS_MALLOC;
        }
    }

    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {

        // Look for start of Eigen-Vector
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            printf("\tLoading Eigen-Vector = %d\n", eigenValue+1);

            // Fast forward 3 lines
            for (i = 0; i < 3; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            if (line != NULL) EG_free(line);
            line = NULL;

            // Loop through the file and fill up the data matrix
            for (i = 0; i < *numGridPoint; i++) {
                for (j = 0; j < numVariable; j++) {

                    if (j == 0) {
                        fscanf(fp, "%lf", &(*dataMatrix)[eigenValue][j+numVariable*i]);
                        //printf("Data (i,j) = (%d, %d) = %f\n", i,j, (*dataMatrix)[eigenValue][j+numVariable*i]);
                        fscanf(fp, "%s", tempString);
                        //printf("String = %s\n", tempString);
                        j = j + 1;
                        (*dataMatrix)[eigenValue][j+numVariable*i] = 0.0;
                    } else fscanf(fp, "%lf", &(*dataMatrix)[eigenValue][j+numVariable*i]);

                    //printf("Data (i,j) = (%d, %d) = %f\n", i,j, (*dataMatrix)[eigenValue][j+numVariable*i]);
                }
            }

            eigenValue += 1;
        }

        if (eigenValue == *numEigenVector) break;
    }

    if (line != NULL) EG_free(line);

    return CAPS_SUCCESS;
}

// Read data from a Nastran F06 file and load it into a dataMatrix[numEigenVector][5]
// where variables are eigenValue, eigenValue(radians), eigenValue(cycles), generalized mass, and generalized stiffness.                                                                   MASS              STIFFNESS
int nastran_readF06EigenValue(FILE *fp, int *numEigenVector, double ***dataMatrix) {

    int status; // Function return

    int i, j;// Indexing

    int tempInt, eigenValue =0;

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *beginEigenLine = "                                              R E A L   E I G E N V A L U E S";
    //char *endEigenLine = "1";

    int numVariable = 5; // EigenValue, eigenValue(radians), eigenValue(cycles), generalized mass, and generalized stiffness.

    printf("Reading Nastran FO6 file - extracting Eigen-Values!\n");

    *numEigenVector = 0;

    if (fp == NULL) return CAPS_IOERR;

    // See how many Eigen-Values we have
    status = nastran_readF06NumEigenValue(fp, numEigenVector);
    printf("\tNumber of Eigen-Values = %d\n", *numEigenVector);
    if (status != CAPS_SUCCESS) return status;

    // Allocate dataMatrix array
    AIM_FREE(*dataMatrix);

    *dataMatrix = (double **) EG_alloc(*numEigenVector *sizeof(double *));
    if (*dataMatrix == NULL) return EGADS_MALLOC; // If allocation failed ....

    for (i = 0; i < *numEigenVector; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc(numVariable*sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

              AIM_FREE((*dataMatrix)[j]);
            }

            AIM_FREE(*dataMatrix);

            return EGADS_MALLOC;
        }
    }

    // Loop through the file again and pull out data
    while (eigenValue != *numEigenVector) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if (status < 0) break;

        // Look for start of Eigen-Vector
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            // Fast forward 2 lines
            for (i = 0; i < 2; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            for (i = 0; i < *numEigenVector; i++) {
                fscanf(fp, "%d", &eigenValue);
                printf("\tLoading Eigen-Value = %d\n", eigenValue);

                fscanf(fp, "%d", &tempInt);

                // Loop through the file and fill up the data matrix
                for (j = 0; j < numVariable; j++) {

                    fscanf(fp, "%lf", &(*dataMatrix)[i][j]);

                }
            }
        }
    }

    if (line != NULL) free(line);

    return CAPS_SUCCESS;
}

// Read data from a Nastran F06 file and load it into a dataMatrix[numGridPoint][8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int nastran_readF06Displacement(FILE *fp, int subcaseId, int *numGridPoint, double ***dataMatrix) {

    int status; // Function return

    int i, j; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *outputSubcaseLine = "0                                                                                                            SUBCASE ";
    char *displacementLine ="                                             D I S P L A C E M E N T   V E C T O R";
    char *beginSubcaseLine=NULL;
    char *endSubcaseLine = "1";
    char tempString[30];

    int numVariable = 8; // Grid Id, coord Id T1, T2, T3, R1, R2, R3
    int intLength;
    int lineFastForward = 0;

    printf("Reading Nastran FO6 file - extracting Displacements!\n");

    *numGridPoint = 0;

    if (fp == NULL) return CAPS_IOERR;

    // Rewind the file
    rewind(fp);

    if      (subcaseId >= 1000) intLength = 4;
    else if (subcaseId >= 100) intLength = 3;
    else if (subcaseId >= 10) intLength = 2;
    else intLength = 1;

    if (subcaseId > 0) {
        beginSubcaseLine = (char *) EG_alloc((strlen(outputSubcaseLine)+intLength+1)*sizeof(char));
        if (beginSubcaseLine == NULL) return EGADS_MALLOC;

        snprintf(beginSubcaseLine,strlen(outputSubcaseLine)+intLength+1,"%s%d",outputSubcaseLine, subcaseId);

        beginSubcaseLine[strlen(outputSubcaseLine)+intLength] = '\0';

        lineFastForward = 4;

    } else {

        beginSubcaseLine = (char *) EG_alloc((strlen(displacementLine)+1)*sizeof(char));
        if (beginSubcaseLine == NULL) return EGADS_MALLOC;
        snprintf(beginSubcaseLine,strlen(displacementLine)+1,"%s",displacementLine);

        //beginSubcaseLine[strlen(displacementLine)] = '\0';

        lineFastForward = 2;
    }

    // Loop through file line by line until we have determined how many grid points we have
    while (*numGridPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if (status < 0) break;

        // Look for start of subcaseId
        if (strncmp(beginSubcaseLine, line, strlen(beginSubcaseLine)) == 0) {

            // Fast forward lines
            for (i = 0; i < lineFastForward; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through lines counting the number of grid points
            while (getline(&line, &linecap, fp) >= 0) {

                if (strncmp(endSubcaseLine, line,strlen(endSubcaseLine)) == 0)break;
                *numGridPoint +=1;
            }
        }
    }
    if (line != NULL) EG_free(line);
    line = NULL;
    linecap = 0;

    printf("Number of Grid Points = %d\n", *numGridPoint);

    if (*numGridPoint == 0) {
        printf("Either data points  = 0 and/or subcase wasn't found\n");

        if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);
        return CAPS_NOTFOUND;
    }

    // Rewind the file
    rewind(fp);

    // Allocate dataMatrix array
    //if (*dataMatrix != NULL) EG_free(*dataMatrix);

    *dataMatrix = (double **) EG_alloc(*numGridPoint *sizeof(double *));
    if (*dataMatrix == NULL) {
        if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);
        return EGADS_MALLOC; // If allocation failed ....
    }

    for (i = 0; i < *numGridPoint; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc(numVariable*sizeof(double));

        if ( (*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

                if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
            }

            if ((*dataMatrix) != NULL) EG_free((*dataMatrix));

            if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);

            return EGADS_MALLOC;
        }
    }

    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {

        // Look for start of Displacement
        if (strncmp(beginSubcaseLine, line, strlen(beginSubcaseLine)) == 0) {

            printf("Loading displacements for Subcase = %d\n", subcaseId);

            // Fast forward lines
            for (i = 0; i < lineFastForward; i++) {
                status = getline(&line, &linecap, fp);
                //printf("Line - %s\n", line);
                if (status < 0) {
                    printf("Unable to fast forward through file- status %d\n", status);
                    break;
                }
            }


            // Loop through the file and fill up the data matrix
            for (i = 0; i < (*numGridPoint); i++) {
                for (j = 0; j < numVariable; j++) {

                    if (j == 0){// || j % numVariable+1 == 0) {
                        fscanf(fp, "%lf", &(*dataMatrix)[i][j]);
                        fscanf(fp, "%s", tempString);
                        j = j + 1;
                        (*dataMatrix)[i][j] = 0.0;
                    } else fscanf(fp, "%lf", &(*dataMatrix)[i][j]);
                }
            }

            break;
        }
    }

    if (beginSubcaseLine != NULL) EG_free(beginSubcaseLine);
    if (line != NULL) EG_free(line);

    printf("Done reading displacements for Subcase = %d\n", subcaseId);
    return CAPS_SUCCESS;
}

// Read weight value from a Nastran F06 file
int nastran_readF06GridPointWeightGeneratorOutput(FILE *fp, double *mass, double cg[3],
                                                  double is[3], double iq[3], double q[9])
{
    int i, j, status;

    char labelArray[3];
    double massArray[3];
    double cgMatrix[3][3], isMatrix[3][3], iqArray[3], qMatrix[3][3];

    char *beginGRDPNTLine = "                           O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R";

    size_t linecap = 0;
    char *line = NULL; // Temporary line holder

    int foundGRDPNT = (int) false;

    if (fp == NULL) return CAPS_IOERR;

    while (!foundGRDPNT) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if (status < 0) break;

        // If we enter the Grid Point Generator section, get weight
        if (strncmp(beginGRDPNTLine, line, strlen(beginGRDPNTLine)) == 0) {

            // Skip to line containing MASS (15 lines down)
            for (i = 0; i < 15; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Get XYZ label, MASS, and C.G. values in next 3 lines
            for (i = 0; i < 3; i++) {

                sscanf(line, " %c%lf%lf%lf%lf", &labelArray[i],
                                                &massArray[i],
                                                &cgMatrix[i][0],
                                                &cgMatrix[i][1],
                                                &cgMatrix[i][2]);

                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // All mass values should be equal
            if (massArray[0] != massArray[1] ||
                massArray[1] != massArray[2] ||
                massArray[2] != massArray[0]) {
                printf("Not all extracted MASS values are equal\n");
                if (line != NULL) free(line);
                return CAPS_BADVALUE;
            }

            // Set mass
            if (mass != NULL) {
                *mass = massArray[0];
                // printf("\tMASS: %lf\n", *mass);
            }

            // Set C.G.
            if (cg != NULL) {
                cg[0] = cgMatrix[1][0];
                cg[1] = cgMatrix[0][1];
                cg[2] = cgMatrix[0][2];
                // printf("\tC.G.: [%lf, %lf, %lf]\n", cg[0], cg[1], cg[2]);
            }

            // Skip to I(S) matrix, (1 line down)
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            // Get I(S) values in next 3 lines
            for (i = 0; i < 3; i++) {

                sscanf(line, " * %lf%lf%lf *", &isMatrix[i][0],
                                               &isMatrix[i][1],
                                               &isMatrix[i][2]);

                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // set I(S)
            if (is != NULL) {
                is[I11] = isMatrix[0][0]; // ixx
                is[I21] = isMatrix[0][1]; // ixy
                is[I22] = isMatrix[1][1]; // iyy
                is[I31] = isMatrix[0][2]; // ixz
                is[I32] = isMatrix[1][2]; // iyz
                is[I33] = isMatrix[2][2]; // izz
                // printf("\tI(S): [%lf, %lf, %lf, %lf, %lf, %lf]\n",
                //         is[0], is[1], is[2], is[3], is[4], is[5]);
            }

            // Skip to I(Q) matrix, (1 line down)
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            // Get I(Q) values in next 3 lines
            for (i = 0; i < 3; i++) {
                // Get I(Q) value
                sscanf(line, " * %lf *", &iqArray[i]);

                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // set I(Q)
            if (iq != NULL) {
                iq[0] = iqArray[0];
                iq[1] = iqArray[1];
                iq[2] = iqArray[2];
                // printf("\tI(Q): [%lf, %lf, %lf]\n", iq[0], iq[1], iq[2]);
            }

            // Skip to Q matrix, (1 line down)
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            // Get Q values in next 3 lines
            for (i = 0; i < 3; i++) {

                sscanf(line, " * %lf%lf%lf *", &qMatrix[i][0],
                                               &qMatrix[i][1],
                                               &qMatrix[i][2]);

                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // set Q
            if (q != NULL) {
                for (i = 0; i < 3; i++) {
                    for (j = 0; j < 3; j++) {
                        q[i + 3 * j] = qMatrix[i][j];
                    }
                }
                // printf("\tQ   : [[%lf, %lf, %lf],\n"
                //        "\t       [%lf, %lf, %lf],\n"
                //        "\t       [%lf, %lf, %lf]]\n", q[0], q[1], q[2],
                //                                       q[3], q[4], q[5],
                //                                       q[6], q[7], q[9]);
            }

            foundGRDPNT = (int) true;
        }

    }

    if (line != NULL) free(line);

    // Rewind the file
    rewind(fp);

    if (!foundGRDPNT) return CAPS_NOTFOUND;

    return CAPS_SUCCESS;
}

// Read objective values for a Nastran OP2 file  and liad it into a dataMatrix[numPoint]
int nastran_readOP2Objective(/*@unused@*/char *filename, int *numData,  double **dataMatrix) {

    int status;

    *numData = 0;
    *dataMatrix = NULL;

#ifdef HAVE_PYTHON
    PyObject* mobj = NULL;
#if CYTHON_PEP489_MULTI_PHASE_INIT
    PyModuleDef *mdef = NULL;
    PyObject *modname = NULL;
#endif
#endif

#ifdef HAVE_PYTHON
        {
            // Flag to see if Python was already initialized - i.e. the AIM was called from within Python
            int initPy = (int) false;

            printf("\nUsing Python to read OP2 file\n");

            // Initialize python
            if (Py_IsInitialized() == 0) {
                printf("\tInitializing Python within AIM\n\n");
                Py_InitializeEx(0);
                initPy = (int) true;
            } else {
                initPy = (int) false;
            }

            PyGILState_STATE gstate;
            gstate = PyGILState_Ensure();

            // Taken from "main" by running cython with --embed
            #if PY_MAJOR_VERSION < 3
                initnastranOP2Reader();
            #elif CYTHON_PEP489_MULTI_PHASE_INIT
                if (nastranOP2Reader_Initialized == (int)false || initPy == (int)true) {
                    nastranOP2Reader_Initialized = (int)true;
                    mobj = PyInit_nastranOP2Reader();
                    if (!PyModule_Check(mobj)) {
                        mdef = (PyModuleDef *) mobj;
                        modname = PyUnicode_FromString("nastranOP2reader");
                        mobj = NULL;
                        if (modname) {
                            mobj = PyModule_NewObject(modname);
                            Py_DECREF(modname);
                            if (mobj) PyModule_ExecDef(mobj, mdef);
                        }
                    }
                }
            #else
                mobj = PyInit_nastranOP2Reader();
            #endif

            if (PyErr_Occurred()) {
                PyErr_Print();
                #if PY_MAJOR_VERSION < 3
                    if (Py_FlushLine()) PyErr_Clear();
                #endif
                /* Release the thread. No Python API allowed beyond this point. */
                PyGILState_Release(gstate);
                status = CAPS_BADOBJECT;
                goto cleanup;
            }

            status = nastran_getObjective((const char *) filename, numData, dataMatrix);
            if (status == -1) {
                printf("\tError: Python error occurred while reading OP2 file\n");
            } else {
                printf("\tDone reading OP2 file with Python\n");
            }

            if (PyErr_Occurred()) {
                PyErr_Print();
                #if PY_MAJOR_VERSION < 3
                    if (Py_FlushLine()) PyErr_Clear();
                #endif
            }

            /* Release the thread. No Python API allowed beyond this point. */
            PyGILState_Release(gstate);

            // Close down python
            if (initPy == (int) false) {
                printf("\n");
            } else {
                printf("\tClosing Python\n");
                Py_Finalize(); // Do not finalize if the AIM was called from Python
            }

            if (status == -1) status = CAPS_NOTFOUND;
            else {
                if (numData == 0) status = CAPS_BADVALUE;
                else              status = CAPS_SUCCESS;
            }
        }
#else

    status = CAPS_NOTIMPLEMENT;
    goto cleanup;
#endif

cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Status %d during nastran_readOP2Objective\n", status);

    return status;
}

// lagrange interpolation derivative
static double _dL(double x, double x0, double x1, double x2) {

    return ((x-x2) + (x-x1)) / ((x0-x1) * (x0-x2));
}

// Get interpolated z coordinate, using 3 bracketing points
// from xi, zi to define interpolating function
static double _dzdx(double x, int n, double *xi, double *zi) {

    int i, j, firstBracketIndex = 0;
    double dz, xbracket[3]={0,0,0}, zbracket[3]={0,0,0};

    // get 3 bracketing points
    for (i = 0; i < n; i++) {

        // if xi coord is greater than target x, found bracketing point
        if (xi[i] > x) {

            // first bracketing point is previous point
            if (i != n-1) {
                firstBracketIndex = (i-1);
            }
            // unless this is last point, then first is two points before
            else {
                firstBracketIndex = (i-2);
            }

            for (j = 0; j < 3; j++) {
                xbracket[j] = xi[firstBracketIndex + j];
                zbracket[j] = zi[firstBracketIndex + j];
            }

            break;
        }

        // if last iteration and not found, error
        if (i == n-1) {
            PRINT_ERROR("Could not find bracketing point in dzdx: %f!", x);
            return 0.0;
        }
    }

    dz = ( zbracket[0] * _dL(x, xbracket[0], xbracket[1], xbracket[2])
         + zbracket[1] * _dL(x, xbracket[1], xbracket[0], xbracket[2])
         + zbracket[2] * _dL(x, xbracket[2], xbracket[0], xbracket[1]));

    return dz;
}

static double _getEndDownwash(double x, int n, double *xi, double *zi) {

    return atan(_dzdx(x, n, xi, zi));
}

static double _getEndGlobalDownwash(vlmSectionStruct *section) {

    double xVec, zVec;

    xVec = section->xyzTE[0] - section->xyzLE[0];
    zVec = section->xyzTE[2] - section->xyzLE[2];

    return atan(zVec/xVec);
}

static double _getPanelDownwash(double wroot, double wtip, double yroot, double ytip, double yj) {

    return wroot + (wtip - wroot) * ((yj - yroot) / (ytip - yroot));
}

static int _getSectionCamberTwist(void *aimInfo,
                                  vlmSectionStruct *sectionRoot, vlmSectionStruct *sectionTip,
                                  int numChord, int numSpan, int *numPanelOut, double **downwashOut) {

    int status;

    int i, ichord, ispan, imid; // Indexing

    int numPanel, numChordDiv, numSpanDiv;
    double *xCoordRoot = NULL, *xCoordTip = NULL, *yCoord = NULL;
    double *zCamberRoot = NULL, *zCamberTip=0;
    double xmid, ymid, wroot, wtip, wij, yroot, ytip;
    double *downwash = NULL;

    if (sectionRoot == NULL) return CAPS_NULLVALUE;
    if (sectionTip == NULL) return CAPS_NULLVALUE;

    numChordDiv = numChord + 1;
    numSpanDiv = numSpan + 1;

    // get normalized chordwise coordinates and camber line

    status = vlm_getSectionCamberLine(aimInfo, sectionRoot,
                                      1.0, // Cosine distribution
                                      (int) true, numChordDiv,
                                      &xCoordRoot, &zCamberRoot);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = vlm_getSectionCamberLine(aimInfo, sectionTip,
                                      1.0, // Cosine distribution
                                      (int) true, numChordDiv,
                                      &xCoordTip, &zCamberTip);
    if (status != CAPS_SUCCESS) goto cleanup;

    // get normalized spanwise coordinates

    yCoord = EG_alloc(numSpanDiv * sizeof(double));

    for (i = 0; i < numSpanDiv; i++) {
        yCoord[i] = i / (numSpanDiv-1.0);
    }

    // get panel downwashes

    numPanel = numSpan * numChord;

    downwash = EG_alloc(numPanel * sizeof(double));
    if (downwash == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    imid = 0;
    for (ispan = 0; ispan < numSpanDiv-1; ++ispan) {

        for (ichord = 0; ichord < numChordDiv-1; ichord++) {

            // mid panel coordinates
            xmid = (xCoordRoot[ichord] + xCoordRoot[ichord+1]) / 2; // xCoordRoot and xCoordTip should be the same since normalized
            ymid  = (yCoord[ispan] + yCoord[ispan+1]) / 2;

            // wroot
            wroot = _getEndDownwash(xmid, numChordDiv, xCoordRoot, zCamberRoot);
            wroot += _getEndGlobalDownwash(sectionRoot);

            // wtip
            wtip = _getEndDownwash(xmid, numChordDiv, xCoordTip, zCamberTip);
            wtip += _getEndGlobalDownwash(sectionTip);

            // yroot, ytip
            yroot = yCoord[0];
            ytip = yCoord[numSpanDiv-1];

            // wij
            wij = _getPanelDownwash(wroot, wtip, yroot, ytip, ymid);

            downwash[imid] = wij;

            imid++;
        }
    }

    // printf("downwash: ");
    // for (i = 0; i < numPanel; i++) {
    //     printf("%f: %d, ", downwash[i], ((downwash[i] > -0.23) && (downwash[i] < 0.05)));
    // }
    // printf("\n");

    status = CAPS_SUCCESS;

    cleanup:

        if (status == CAPS_SUCCESS) {
            *numPanelOut = numPanel;
            *downwashOut = downwash;
        }
        else {
            if (downwash != NULL) EG_free(downwash);
        }
        if (xCoordRoot != NULL) EG_free(xCoordRoot);
        if (xCoordTip != NULL) EG_free(xCoordTip);
        if (yCoord != NULL) EG_free(yCoord);
        if (zCamberRoot != NULL) EG_free(zCamberRoot);
        if (zCamberTip != NULL) EG_free(zCamberTip);

        return status;
}

// Write Nastran DMI cards for downwash matrix from collection of feaAeroStructs
int nastran_writeAeroCamberTwist(void *aimInfo, FILE *fp, int numAero, feaAeroStruct *feaAero, const feaFileFormatStruct *feaFileFormat) {

    int i, j, iAero, status;

    int numPanel, numSectionPanel;
    int form, tin, tout;
    double *downwash = NULL, *sectionDownwash = NULL;

    feaAeroStruct *aero;

    if (fp == NULL) return CAPS_IOERR;
    if (feaAero == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    numPanel = 0;
    numSectionPanel = 0;

    for (iAero = 0; iAero < numAero; iAero++) {

        aero = &feaAero[iAero];

        for (i = 0; i < aero->vlmSurface.numSection-1; i++) {

            status = _getSectionCamberTwist(aimInfo,
                                            &aero->vlmSurface.vlmSection[i],
                                            &aero->vlmSurface.vlmSection[i+1],
                                            aero->vlmSurface.Nchord,
                                            aero->vlmSurface.NspanTotal,
                                            &numSectionPanel, &sectionDownwash);
            if (status != CAPS_SUCCESS) goto cleanup;

            if (downwash == NULL) {
                downwash = EG_alloc(numSectionPanel * sizeof(double));
                if (downwash == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }
            }
            else {
                downwash = EG_reall(downwash, (numPanel + numSectionPanel) * sizeof(double));
                if (downwash == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }
            }

            for (j = 0; j < numSectionPanel; j++) {
                downwash[numPanel++] = sectionDownwash[j];
            }

            EG_free(sectionDownwash);
            sectionDownwash = NULL;
        }
    }

    form = 2;
    tin = 1;
    tout = 0;

    // Write DIM card
    status = nastranCard_dmi(
        fp,
        "W2GJ",
        &form, // form
        &tin, // tin
        &tout, // tout
        numPanel, // m
        1, // n
        downwash, // a
        NULL, // b
        feaFileFormat->fileType
    );
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;

    cleanup:

        if (sectionDownwash != NULL) EG_free(sectionDownwash);
        if (downwash != NULL) EG_free(downwash);

        return status;
}
