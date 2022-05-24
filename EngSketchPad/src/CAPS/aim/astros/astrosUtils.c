// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>
#include <math.h>
#include <ctype.h>

#include "aimUtil.h"

#include "feaTypes.h"    // Bring in FEA structures
#include "capsTypes.h"   // Bring in CAPS types

#include "miscUtils.h"   // Bring in misc. utility functions
#include "vlmUtils.h"    // Bring in vlm. utility functions
#include "meshUtils.h"   // Bring in mesh utility functions
#include "cardUtils.h"   // Bring in card utility functions
#include "astrosUtils.h" // Bring in astros utility header
#include "astrosCards.h"

#ifdef WIN32
#define strcasecmp  stricmp
#endif


// Write a Astros connections card from a feaConnection structure
int astros_writeConnectionCard(FILE *fp, feaConnectionStruct *feaConnect,
                               feaFileFormatStruct *feaFileFormat)
{
    int status = CAPS_SUCCESS;
    int compNumbers[2];

    if (fp == NULL) return CAPS_IOERR;
    if (feaConnect == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    compNumbers[0] = feaConnect->componentNumberStart;
    compNumbers[1] = feaConnect->componentNumberEnd;

    // Mass
    if (feaConnect->connectionType == Mass) {

        status = astrosCard_cmass2(
            fp,
            &feaConnect->elementID, // eid
            &feaConnect->mass, // m
            feaConnect->connectivity, // Gi
            compNumbers, // Ci
            NULL, // tmin
            NULL, // tmax
            feaFileFormat->gridFileType);
        if (status != CAPS_SUCCESS) return status;

    }

    // Spring
    if (feaConnect->connectionType == Spring) {

        status = astrosCard_celas2(fp,
                                   &feaConnect->elementID, // eid
                                   &feaConnect->stiffnessConst, // k
                                   feaConnect->connectivity, // Gi
                                   compNumbers, // Ci
                                   &feaConnect->dampingConst, // ge
                                   &feaConnect->stressCoeff, // s
                                   NULL, // tmin
                                   NULL, // tmax
                                   feaFileFormat->gridFileType);
        if (status != CAPS_SUCCESS) return status;

    }

    // Damper
    if (feaConnect->connectionType == Damper) {

        printf("Astros doesn't appear to support CDAMP2 cards!\n");
        return CAPS_NOTFOUND;
    }

    // Rigid Body
    if (feaConnect->connectionType == RigidBody) {

        status = astrosCard_rbe2(fp,
                                 &feaConnect->connectionID, // setid
                                 &feaConnect->elementID, // eid
                                 &feaConnect->connectivity[0], // gn
                                 &feaConnect->dofDependent, // cm
                                 1, &feaConnect->connectivity[1], // GMi
                                 feaFileFormat->gridFileType);
        if (status != CAPS_SUCCESS) return status;

    }

    // Rigid Body Interpolate
    if (feaConnect->connectionType == RigidBodyInterpolate) {

        status = astrosCard_rbe3(fp,
                                 &feaConnect->connectionID, // setid
                                 &feaConnect->elementID, // eid
                                 &feaConnect->connectivity[1], // refg
                                 &feaConnect->dofDependent, // refc
                                 feaConnect->numMaster,
                                 feaConnect->masterWeighting, // wt
                                 feaConnect->masterComponent, // c
                                 feaConnect->masterIDSet, // g
                                 0,
                                 NULL, // gm
                                 NULL, // cm
                                 feaFileFormat->gridFileType);
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write a Nastran Load card from a feaLoad structure
int astros_writeLoadCard(FILE *fp, meshStruct *mesh, feaLoadStruct *feaLoad,
                         feaFileFormatStruct *feaFileFormat)
{

    int status = CAPS_SUCCESS;

    int i, j; // Indexing
    int elemIDindex; // used in PLOAD element ID check logic
    int numGridPoints;
    double avgPressure;

    cardStruct card;

#ifdef ASTROS_11_DOES_NOT_HAVE_PLOAD4
    double *pressureDistributeForce = NULL;
#endif

    // Concentrated force at a grid point
    if (feaLoad->loadType ==  GridForce) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = astrosCard_force(fp,
                                      &feaLoad->loadID, // sid
                                      &feaLoad->gridIDSet[i], // g
                                      &feaLoad->coordSystemID, // cid
                                      &feaLoad->forceScaleFactor, // f
                                      feaLoad->directionVector, // Ni
                                      feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Concentrated moment at a grid point
    if (feaLoad->loadType ==  GridMoment) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = astrosCard_moment(fp,
                                       &feaLoad->loadID, // sid
                                       &feaLoad->gridIDSet[i], // g
                                       &feaLoad->coordSystemID, // cid
                                       &feaLoad->momentScaleFactor, // m
                                       feaLoad->directionVector,
                                       feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Gravitational load
    if (feaLoad->loadType ==  Gravity) {

        status = astrosCard_grav(fp,
                                 &feaLoad->loadID, // sid
                                 &feaLoad->coordSystemID, // cid
                                 &feaLoad->gravityAcceleration, // g
                                 feaLoad->directionVector, // Ni
                                 feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    // Pressure load
    if (feaLoad->loadType ==  Pressure) {

        for (i = 0; i < feaLoad->numElementID; i++) {

            elemIDindex = -1;
            // Check to make sure element ID desired matches the id in the mesh input
            if (mesh->element[i].elementID == feaLoad->elementIDSet[i]) {
                elemIDindex = i;
            } else {
                for (j = 0; j < mesh->numElement; j++) {

                    if (mesh->element[j].elementID != feaLoad->elementIDSet[i]) continue;

                    elemIDindex = j;
                    break;

                }
            }

            if (elemIDindex < 0 ) {
                printf("Error: Element index wasn't found!\n");
                return CAPS_BADVALUE;
            }

            if (mesh->element[elemIDindex].elementType == Quadrilateral ||
                mesh->element[elemIDindex].elementType == Triangle) {

                numGridPoints = mesh_numMeshConnectivity(mesh->element[elemIDindex].elementType);

                status = astrosCard_pload(
                    fp,
                    &feaLoad->loadID, // sid
                    &feaLoad->pressureForce, // P
                    numGridPoints,
                    mesh->element[elemIDindex].connectivity, // G
                    feaFileFormat->fileType
                );
                if (status != CAPS_SUCCESS) return status;

            } else {
                printf("Unsupported element type for PLOAD in astrosAIM!\n");
                return CAPS_BADVALUE;
            }

            /*
                fprintf(fp,"%-8s", "PLOAD2");

            tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
            fprintf(fp, "%s%s",delimiter, tempString);
            EG_free(tempString);

            tempString = convert_doubleToString(feaLoad->pressureForce, fieldWidth, 1);
            fprintf(fp, "%s%s",delimiter, tempString);
            EG_free(tempString);

            tempString = convert_integerToString(feaLoad->elementIDSet[i], fieldWidth, 1);
            fprintf(fp, "%s%s\n",delimiter, tempString);
            EG_free(tempString);
             */
        }
    }

    // Pressure load at element Nodes
    if (feaLoad->loadType ==  PressureDistribute) {

        for (i = 0; i < feaLoad->numElementID; i++) {

            status = astrosCard_pload4(fp,
                                       &feaLoad->loadID, // lid
                                       &feaLoad->elementIDSet[i], // eid
                                       4, feaLoad->pressureDistributeForce, // P
                                       NULL, // cid
                                       NULL, // V
                                       feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    // Pressure load at element Nodes - different distribution per element
    if (feaLoad->loadType ==  PressureExternal) {

#ifdef ASTROS_11_DOES_NOT_HAVE_PLOAD4
        for (i = 0; i < feaLoad->numElementID; i++) {
            // fprintf(fp,"%-8s", "PLOAD4");

            // tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
            // fprintf(fp, "%s%s",delimiter, tempString);
            // EG_free(tempString);

            // tempString = convert_integerToString(feaLoad->elementIDSet[i], fieldWidth, 1);
            // fprintf(fp, "%s%s",delimiter, tempString);
            // EG_free(tempString);

            // for (j = 0; j < 4; j++) {
            //     tempString = convert_doubleToString(feaLoad->pressureMultiDistributeForce[4*i + j], fieldWidth, 1);
            //     fprintf(fp, "%s%s",delimiter, tempString);
            //     EG_free(tempString);
            // }

            // fprintf(fp,"\n");
            pressureDistributeForce = EG_alloc(sizeof(double) * 4);
            for (j = 0; j < 4; j++) {
                pressureDistributeForce[j] = feaLoad->pressureMultiDistributeForce[4 * i + j];
            }

            status = astrosCard_pload4(fp,
                                       &feaLoad->loadID, // lid
                                       &feaLoad->elementIDSet[i], // eid
                                       4, pressureDistributeForce, // P
                                       NULL, // cid
                                       NULL, // V
                                       feaFileFormat->fileType);

            EG_free(pressureDistributeForce);

            if (status != CAPS_SUCCESS) return status;
        }
#endif


        for (i = 0; i < feaLoad->numElementID; i++) {

            elemIDindex = -1;
            // Check to make sure element ID desired matches the id in the mesh input
            if (mesh->element[i].elementID == feaLoad->elementIDSet[i]) {
                elemIDindex = i;
            } else {
                for (j = 0; j < mesh->numElement; j++) {

                    if (mesh->element[j].elementID != feaLoad->elementIDSet[i]) continue;

                    elemIDindex = j;
                    break;

                }
            }

            if (elemIDindex < 0 ) {
                printf("Error: Element index wasn't found!\n");
                return CAPS_BADVALUE;
            }

            numGridPoints = mesh_numMeshConnectivity(mesh->element[elemIDindex].elementType);
            avgPressure = 0;
            for (j = 0; j < numGridPoints; j++) {
                avgPressure += feaLoad->pressureMultiDistributeForce[4*i + j];
            }
            avgPressure /= numGridPoints;

            if (mesh->element[elemIDindex].elementType == Quadrilateral ||
                mesh->element[elemIDindex].elementType == Triangle) {

                for (j = 0; j < numGridPoints; j++) {

                    status = astrosCard_pload(
                        fp,
                        &feaLoad->loadID, // sid
                        &avgPressure, // P
                        numGridPoints,
                        mesh->element[elemIDindex].connectivity, // G
                        feaFileFormat->fileType
                    );
                    if (status != CAPS_SUCCESS) return status;
                }

            } else {
                printf("Unsupported element type for PLOAD in astrosAIM!\n");
                return CAPS_BADVALUE;
            }


            // fprintf(fp,"\n");

            /*
                fprintf(fp,"%-8s", "PLOAD2");

            tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
            fprintf(fp, "%s%s",delimiter, tempString);
            EG_free(tempString);

            tempString = convert_doubleToString(feaLoad->pressureForce, fieldWidth, 1);
            fprintf(fp, "%s%s",delimiter, tempString);
            EG_free(tempString);

            tempString = convert_integerToString(feaLoad->elementIDSet[i], fieldWidth, 1);
            fprintf(fp, "%s%s\n",delimiter, tempString);
            EG_free(tempString);
             */
        }

    }

    // Rotational velocity
    if (feaLoad->loadType ==  Rotational) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = card_initiate(&card, "RFORCE", feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;

            card_addInteger(&card, feaLoad->loadID);

            card_addInteger(&card, feaLoad->gridIDSet[i]);

            card_addInteger(&card, feaLoad->coordSystemID);

            card_addDouble(&card, feaLoad->angularVelScaleFactor);

            card_addDoubleArray(&card, 3, feaLoad->directionVector);

            card_addBlank(&card);

            card_addDouble(&card, feaLoad->angularAccScaleFactor);

            card_write(&card, fp);

            card_destroy(&card);
        }
    }
    // Thermal load at a grid point
    if (feaLoad->loadType ==  Thermal) {

        status = astrosCard_tempd(
            fp,
            1,
            &feaLoad->loadID, // SID
            &feaLoad->temperatureDefault, // T
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        for (i = 0; i < feaLoad->numGridID; i++) {

            status = astrosCard_temp(fp,
                                     &feaLoad->loadID, // sid
                                     1,
                                     &feaLoad->gridIDSet[i], // G
                                     &feaLoad->temperature, // T
                                     feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    return CAPS_SUCCESS;
}

/*
 * Flow of ASTROS aero cards:
 *
 * for each vlm panel { (one per WING, CANARD, FIN) {
 *
 *     write vlm chord cuts AEFACT 0.0-100.0 pct
 *     write vlm span cuts AEFACT in PHYSICAL Y-COORDS, MUST *EXACTLY* ALIGN WITH AIRFOIL CARDS
 *     write CAERO6 card (one per WING, CANARD, FIN)
 *
 *     for each airfoil section in panel {
 *         write airfoil chord station AEFACT 0.0-100.0 pct (may share or be separate per airfoil)
 *         write airfoil upper surf AEFACT in pct chord (0.01 = 1%)
 *         write airfoil lower surf AEFACT in pct chord
 *         write AIRFOIL card (ref chord/upper/lower AEFACTS and vlm CAERO6)
 *     }
 * }
 */

// Write a Astros AEROS card from a feaAeroRef structure
int astros_writeAEROSCard(FILE *fp, feaAeroRefStruct *feaAeroRef,
                          feaFileFormatStruct *feaFileFormat)
{
    int *rcsid;

    if (feaAeroRef == NULL) return CAPS_NULLVALUE;

    rcsid = &feaAeroRef->rigidMotionCoordSystemID;
    // if (*rcsid == 0) rcsid = NULL; //TODO: temporarily commenting this out to exactly match previous logic

    return astrosCard_aeros(fp, &feaAeroRef->coordSystemID, rcsid,
                            &feaAeroRef->refChord, &feaAeroRef->refSpan,
                            &feaAeroRef->refArea, NULL, NULL, NULL,
                            feaFileFormat->fileType);
}

// Write a Astros AERO card from a feaAeroRef structure
int astros_writeAEROCard(FILE *fp, feaAeroRefStruct *feaAeroRef,
                         feaFileFormatStruct *feaFileFormat)
{
    double defaultRhoRef = 1.0;

    if (feaAeroRef == NULL) return CAPS_NULLVALUE;

    return astrosCard_aero(fp, &feaAeroRef->coordSystemID, &feaAeroRef->refChord,
                            &defaultRhoRef, feaFileFormat->fileType);
}

// Check to make for the bodies' topology are acceptable for airfoil shapes
// Return: CAPS_SUCCESS if everything is ok
//                 CAPS_SOURCEERR if not acceptable
//                 CAPS_* if something else went wrong
int astros_checkAirfoil(void *aimInfo, feaAeroStruct *feaAero)
{
    int status; // Function return status

    int i;

    if (aimInfo == NULL) return CAPS_NULLVALUE;
    if (feaAero == NULL) return CAPS_NULLVALUE;

    // Loop through sections in surface
    for (i = 0; i < feaAero->vlmSurface.numSection; i++) {
        if ( feaAero->vlmSurface.vlmSection[i].teClass != NODE ) {
            //printf("\tError in astros_checkAirfoil, body has %d nodes and %d edges!\n", numNode, numEdge);
            status = CAPS_SOURCEERR;
            goto cleanup;
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS && status != CAPS_SOURCEERR)
        printf("\tError: Premature exit in astros_checkAirfoil, status = %d\n",
               status);

    return status;
}

// Write out all the Aero cards necessary to define the VLM surface
int astros_writeAeroData(void *aimInfo,
                         FILE *fp,
                         int useAirfoilShape,
                         feaAeroStruct *feaAero,
                         feaFileFormatStruct *feaFileFormat)
{
    int status; // Function return status

    int i, j; // Indexing

    int numPoint = 30;
    double *xCoord=NULL, *yUpper=NULL, *yLower=NULL;
    int sid;

    if (aimInfo == NULL) return CAPS_NULLVALUE;
    if (fp == NULL) return CAPS_IOERR;
    if (feaAero == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // Loop through sections in surface
    for (i = 0; i < feaAero->vlmSurface.numSection; i++) {

        if (useAirfoilShape == (int) true) { // Using the airfoil upper and lower surfaces or panels?

          status = vlm_getSectionCoordX(&feaAero->vlmSurface.vlmSection[i],
                                        1.0, // Cosine distribution
                                        (int) true, (int) true,
                                        numPoint,
                                        &xCoord, &yUpper, &yLower);
          if (status != CAPS_SUCCESS) return status;
          if ((xCoord == NULL) || (yUpper == NULL) || (yLower == NULL)) {
              EG_free(xCoord); xCoord = NULL;
              EG_free(yUpper); yUpper = NULL;
              EG_free(yLower); yLower = NULL;
              return EGADS_MALLOC;
          }

          for (j = 0; j < numPoint; j++) {
              xCoord[j] *= 100.0;
              yUpper[j] *= 1.0;
              yLower[j] *= 1.0;
          }

           fprintf(fp, "$ Upper surface - Airfoil %d (of %d) \n",i+1,
                   feaAero->vlmSurface.numSection);
           sid = feaAero->surfaceID +
                 100*(feaAero->vlmSurface.vlmSection[i].sectionIndex+1) + 1;
           status = astrosCard_aefact(fp, &sid, numPoint, yUpper,
                                      feaFileFormat->fileType);
           if (status != CAPS_SUCCESS) goto cleanup;

           fprintf(fp, "$ Lower surface - Airfoil %d (of %d) \n",
                   i+1, feaAero->vlmSurface.numSection);
           sid++;
           status = astrosCard_aefact(fp, &sid, numPoint, yLower,
                                      feaFileFormat->fileType);
           if (status != CAPS_SUCCESS) goto cleanup;

        } else {

            xCoord = (double*)EG_alloc(numPoint*sizeof(double));
            if (xCoord == NULL) { status = EGADS_MALLOC; goto cleanup; }

            for (j = 0; j < numPoint; j++) {
                xCoord[j] = ( (double) j / ( (double) numPoint - 1.0 ) ) * 100.0;
            }
        }

        // Write chord range
        xCoord[0] = 0.0;
        xCoord[numPoint-1] = 100.0;

        fprintf(fp, "$ Chord - Airfoil %d (of %d) \n",
                i+1, feaAero->vlmSurface.numSection);
        sid = feaAero->surfaceID +
              100*(feaAero->vlmSurface.vlmSection[i].sectionIndex+1);
        status = astrosCard_aefact(fp, &sid, numPoint, xCoord,
                                   feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) goto cleanup;

        EG_free(xCoord); xCoord = NULL;
        EG_free(yUpper); yUpper = NULL;
        EG_free(yLower); yLower = NULL;
    }

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in astros_writeAero, status = %d\n",
               status);

    EG_free(xCoord);
    EG_free(yUpper);
    EG_free(yLower);

    return status;
}

// Write Astros CAERO6 cards from a feaAeroStruct
int astros_writeCAeroCard(FILE *fp, feaAeroStruct *feaAero,
                          feaFileFormatStruct *feaFileFormat)
{
    int status; // Function return status

    int i, j; // Indexing

    double xmin, xmax, result[6];

    int lengthTemp, offSet = 0;
    double *temp = NULL, pi, x;
    int chordID, spanID, defaultIGRP = 1;

    if (feaAero == NULL) return CAPS_NULLVALUE;

    //printf("WARNING CAER06 card may not be correct - need to confirm\n");
    chordID = feaAero->surfaceID + 10*feaAero->surfaceID + 1; // Chord AEFact ID - Coordinate with AEFACT card
    spanID  = feaAero->surfaceID + 10*feaAero->surfaceID + 2; //  Span AEFact ID - Coordinate with AEFACT card

    status = astrosCard_caero6(fp, &feaAero->surfaceID,
                               feaAero->vlmSurface.surfaceType, NULL,
                               &defaultIGRP, &chordID, &spanID,
                               feaFileFormat->fileType);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Write Chord AEFact
    lengthTemp = feaAero->vlmSurface.Nchord+1; // One more point than boxes for spline

    temp = (double *) EG_alloc(lengthTemp*sizeof(double));
    if (temp == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    pi = 3.0*atan(sqrt(3));

    // Set bounds
    temp[0] = 0.0;
    temp[lengthTemp-1] = 100.0;

    for (j = 1; j < lengthTemp-1; j++) {

        x = j*100.0/(lengthTemp-1);

        // Cosine
        if (fabs(feaAero->vlmSurface.Cspace) == 1.0) {

            x = (x - 50.0)/50.0;

            temp[j] = 0.5 * (1.0 +  x + (1.0/pi)*sin(x*pi));

            // Equal spacing
        } else {
            temp[j] = (double) 0.0 + x;
        }
    }

    fprintf(fp, "$ Chord\n");
    status = astrosCard_aefact(fp, &chordID, lengthTemp, temp, feaFileFormat->fileType);
    if (status != CAPS_SUCCESS) goto cleanup;

    if (temp != NULL) EG_free(temp);

    lengthTemp = feaAero->vlmSurface.NspanTotal + 1; //*(feaAero->vlmSurface.numSection-1) + 1;

    temp = (double *) EG_alloc(lengthTemp*sizeof(double));
    if (temp == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // Write Span AEFact
    for (i = 0; i < feaAero->vlmSurface.numSection-1; i++) {

        result[0] = result[3] = feaAero->vlmSurface.vlmSection[i].xyzLE[0];
        result[1] = result[4] = feaAero->vlmSurface.vlmSection[i].xyzLE[1];
        result[2] = result[5] = feaAero->vlmSurface.vlmSection[i].xyzLE[2];

        if (feaAero->vlmSurface.vlmSection[i+1].xyzLE[0] < result[0]) result[0] = feaAero->vlmSurface.vlmSection[i+1].xyzLE[0];
        if (feaAero->vlmSurface.vlmSection[i+1].xyzLE[0] > result[3]) result[3] = feaAero->vlmSurface.vlmSection[i+1].xyzLE[0];
        if (feaAero->vlmSurface.vlmSection[i+1].xyzLE[1] < result[1]) result[1] = feaAero->vlmSurface.vlmSection[i+1].xyzLE[1];
        if (feaAero->vlmSurface.vlmSection[i+1].xyzLE[1] > result[4]) result[4] = feaAero->vlmSurface.vlmSection[i+1].xyzLE[1];
        if (feaAero->vlmSurface.vlmSection[i+1].xyzLE[2] < result[2]) result[2] = feaAero->vlmSurface.vlmSection[i+1].xyzLE[2];
        if (feaAero->vlmSurface.vlmSection[i+1].xyzLE[2] > result[5]) result[5] = feaAero->vlmSurface.vlmSection[i+1].xyzLE[2];

        xmin = result[3] - result[0];
        if (result[4] - result[1] > xmin) xmin = result[4] - result[1];
        if (result[5] - result[2] > xmin) xmin = result[5] - result[2];

        if ((result[4]-result[1])/xmin > 1.e-5) { // Y-ordering

            xmin = feaAero->vlmSurface.vlmSection[i].xyzLE[1];
            xmax = feaAero->vlmSurface.vlmSection[i+1].xyzLE[1];

        } else { // Z-ordering

            xmin = feaAero->vlmSurface.vlmSection[i].xyzLE[2];
            xmax = feaAero->vlmSurface.vlmSection[i+1].xyzLE[2];
        }

        for (j = 0; j < feaAero->vlmSurface.vlmSection[i].Nspan+1; j++) { // One more point than boxes for spline

            temp[j +offSet] = xmin + j*(xmax -xmin)/(feaAero->vlmSurface.vlmSection[i].Nspan);

        }

        // offset so the first point of the section overwrites the last point of the previous section
        offSet += feaAero->vlmSurface.vlmSection[i].Nspan;
    }

    fprintf(fp, "$ Span\n");
    status = astrosCard_aefact(fp, &spanID, lengthTemp, temp, feaFileFormat->fileType);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;

cleanup:

    EG_free(temp);

    return status;

}

// Write out all the Airfoil cards for each each of a surface
int astros_writeAirfoilCard(FILE *fp,
                            /*@unused@*/ int useAirfoilShape, // = true use the airfoils shape, = false panel
                            feaAeroStruct *feaAero,
                            feaFileFormatStruct *feaFileFormat)
{
    int i, status = CAPS_SUCCESS; // Indexing
    int compID, chordID, lowerID, upperID;

    // This assumes the sections are in order
    for (i = 0; i < feaAero->vlmSurface.numSection; i++) {

        compID = feaAero->surfaceID;
        chordID = feaAero->surfaceID + 100*(feaAero->vlmSurface.vlmSection[i].sectionIndex+1);
        upperID = chordID + 1;
        lowerID = chordID + 2;

        status = astrosCard_airfoil(fp, &compID, feaAero->vlmSurface.surfaceType,
                                    NULL, &chordID, &upperID, &lowerID,
                                    NULL, NULL, feaAero->vlmSurface.vlmSection[i].xyzLE,
                                    &feaAero->vlmSurface.vlmSection[i].chord, NULL,
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write Astros Spline1 cards from a feaAeroStruct
int astros_writeAeroSplineCard(FILE *fp, feaAeroStruct *feaAero,
                               feaFileFormatStruct *feaFileFormat)
{
    int firstBoxID, lastBoxID, numSpanWise;

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

    firstBoxID = feaAero->surfaceID;
    lastBoxID = firstBoxID + numSpanWise * feaAero->vlmSurface.Nchord - 1;

    return astrosCard_spline1(fp,
                              &feaAero->surfaceID, // eid
                              NULL, // cp = NULL, CAERO card defines spline plane
                              &feaAero->surfaceID, // macroid
                              &firstBoxID, // box1
                              &lastBoxID, // box2
                              &feaAero->surfaceID, // setg
                              NULL, // dz
                              feaFileFormat->fileType);
}

// Write Astros constraint card from a feaConstraint structure
int astros_writeConstraintCard(FILE *fp, int feaConstraintSetID,
                               feaConstraintStruct *feaConstraint,
                               feaFileFormatStruct *feaFileFormat)
{
    int status;
    int i; // Indexing;

    if (feaConstraint->constraintType == Displacement) {

        for (i = 0; i < feaConstraint->numGridID; i++) {

            status = astrosCard_spc(fp,
                                    &feaConstraintSetID, // sid
                                    1,
                                    &feaConstraint->gridIDSet[i], // G
                                    &feaConstraint->dofConstraint, // C
                                    &feaConstraint->gridDisplacement, // D
                                    feaFileFormat->fileType);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    if (feaConstraint->constraintType == ZeroDisplacement) {

        for (i = 0; i < feaConstraint->numGridID; i++) {

            status = astrosCard_spc1(
                fp,
                &feaConstraintSetID, // sid
                &feaConstraint->dofConstraint, // c
                1, &feaConstraint->gridIDSet[i], // G
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) return status;
        }
    }

    return CAPS_SUCCESS;
}


// Write Astros support card from a feaSupport structure
int astros_writeSupportCard(FILE *fp, feaSupportStruct *feaSupport,
                            feaFileFormatStruct *feaFileFormat)
{
    int status;
    int i; // Indexing;

    if (fp == NULL) return CAPS_IOERR;
    if (feaSupport == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    for (i = 0; i < feaSupport->numGridID; i++) {

        status = astrosCard_suport(fp,
                                   &feaSupport->supportID, // setid
                                   1,
                                   &feaSupport->gridIDSet[i], // ID
                                   &feaSupport->dofSupport, // C
                                   feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write a Astros Property card from a feaProperty structure w/ design parameters
int astros_writePropertyCard(FILE *fp, feaPropertyStruct *feaProperty,
                             feaFileFormatStruct *feaFileFormat,
                             int numDesignVariable,
                             feaDesignVariableStruct feaDesignVariable[])
{
    int status;
    int j, designIndex; // Indexing

    // composite
    char *failureTheory = NULL;
    // int numCompositeLayers;
    // int *compositeMaterialID = NULL;
    // double *compositeThickness = NULL, *compositeOrientation = NULL, *compositeStressOut = NULL;

    // shell
    int *materialBendingID, *materialShearID;
    double *bendingInertiaRatio, *shearMembraneRatio, *massPerArea;

    double *designMin = NULL;

    int found = (int) false;

    if (fp == NULL) return CAPS_IOERR;
    if (feaProperty == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // Check for design minimum
    designMin = NULL;
    found = (int) false;
    for (designIndex = 0; designIndex < numDesignVariable; designIndex++) {

        for (j = 0; j < feaDesignVariable[designIndex].numPropertyID; j++) {

            if (feaDesignVariable[designIndex].propertySetID[j] ==
                feaProperty->propertyID) {
                found = (int) true;
                designMin = &feaDesignVariable[designIndex].lowerBound;
                break;
            }
        }

        if (found == (int) true) break;
    }

    if (feaProperty->massPerArea != 0.0) {
        massPerArea = &feaProperty->massPerArea;

    } else {
        massPerArea = NULL;
    }

    //          1D Elements        //

    // Rod
    if (feaProperty->propertyType ==  Rod) {

        status = astrosCard_prod(fp,
                                 &feaProperty->propertyID, // pid
                                 &feaProperty->materialID, // mid
                                 &feaProperty->crossSecArea, // a
                                 &feaProperty->torsionalConst, // j
                                 &feaProperty->torsionalStressReCoeff, // c
                                 &feaProperty->massPerLength, // nsm
                                 designMin, // tmin
                                 feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    // Bar
    if (feaProperty->propertyType ==  Bar) {

        status = astrosCard_pbar(fp,
                                 &feaProperty->propertyID, // pid
                                 &feaProperty->materialID, // mid
                                 &feaProperty->crossSecArea, // a
                                 &feaProperty->zAxisInertia, // i1
                                 &feaProperty->yAxisInertia, // i2
                                 &feaProperty->torsionalConst, // j
                                 &feaProperty->massPerLength, // nsm
                                 designMin, // tmin
                                 NULL, // k1
                                 NULL, // k2
                                 NULL, // C
                                 NULL, // D
                                 NULL, // E
                                 NULL, // F
                                 NULL, // r12
                                 NULL, // r22
                                 NULL, // alpha
                                 feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    //          2D Elements        //

    if (feaProperty->propertyType == Shear) {

        status = astrosCard_pshear(fp,
                                   &feaProperty->propertyID, // pid
                                   &feaProperty->materialID, // mid
                                   &feaProperty->membraneThickness, // t
                                   massPerArea, // nsm
                                   designMin, // tmin
                                   feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

    }

    if (feaProperty->propertyType == Membrane) {

        status = astrosCard_pqdmem1(fp,
                                    &feaProperty->propertyID, // pid
                                    &feaProperty->materialID, // mid
                                    &feaProperty->membraneThickness, // t
                                    massPerArea, // nsm
                                    designMin, // tmin
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

    }

    // Shell
    if (feaProperty->propertyType == Shell) {

        // Check for design minimum thickness
        found = (int) false;
        for (designIndex = 0; designIndex < numDesignVariable; designIndex++) {
            for (j = 0; j < feaDesignVariable[designIndex].numPropertyID; j++) {

                if (feaDesignVariable[designIndex].propertySetID[j] == feaProperty->propertyID) {
                    //found = (int) true; eja - let design variable upper/lower bounds handle thisMA
                    break;
                }
            }

            if (found == (int) true) break;
        }

        if (feaProperty->materialBendingID != 0) {

            materialBendingID = &feaProperty->materialBendingID;
            bendingInertiaRatio = &feaProperty->bendingInertiaRatio;

        } else { // Print a blank

            if (found == (int) true ||
                    feaProperty->materialShearID != 0 ||
                    feaProperty->massPerArea != 0) {
            }

            materialBendingID = NULL;
            bendingInertiaRatio = NULL;
        }

        if (feaProperty->materialShearID != 0) {

            materialShearID = &feaProperty->materialShearID;
            shearMembraneRatio = &feaProperty->shearMembraneRatio;

        } else { // Print a blank

            // if (found == (int) true || feaProperty->massPerArea != 0) {
            // }

            materialShearID = NULL;
            shearMembraneRatio = NULL;
        }

        status = astrosCard_pshell(fp,
                                   &feaProperty->propertyID, // pid
                                   &feaProperty->materialID, // mid
                                   &feaProperty->membraneThickness, // t
                                   materialBendingID, // mid2
                                   bendingInertiaRatio, // i12t3
                                   materialShearID, // mid3
                                   shearMembraneRatio, // tst
                                   massPerArea, // nsm
                                   NULL, // z1
                                   NULL, // z2
                                   NULL, // mid4
                                   Integer, NULL, // mscid
                                   Integer, NULL, // scsid
                                   NULL, // zoff
                                   designMin, // tmin
                                   feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    // Composite
    if (feaProperty->propertyType == Composite) {

        // FAILURE THEORY
        // HILL, HOFF, TSAI, STRESS, or STRAIN. ASTROS
        // HILL, HOFF, TSAI, STRN. NASTRAN
        if (feaProperty->compositeFailureTheory != NULL) {

            if (strcmp("STRN", feaProperty->compositeFailureTheory) == 0) {
                failureTheory = "STRAIN";
            } else {
                failureTheory = feaProperty->compositeFailureTheory;
            }
        }

        status = astrosCard_pcomp(fp,
                                  &feaProperty->propertyID, // pid
                                  NULL, // z0
                                  massPerArea, // nsm
                                  &feaProperty->compositeShearBondAllowable, // sbond
                                  failureTheory, // ft
                                  designMin, // tmin
                                  NULL, // lopt
                                  feaProperty->numPly,
                                  feaProperty->compositeMaterialID, // MID
                                  feaProperty->compositeThickness, // T
                                  feaProperty->compositeOrientation, // TH
                                  NULL, // SOUT
                                  feaProperty->compositeSymmetricLaminate,
                                  feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    //          3D Elements        //

    // Solid
    if (feaProperty->propertyType == Solid) {

        status = astrosCard_pihex(fp,
                                  &feaProperty->propertyID, // pid
                                  &feaProperty->materialID, // mid
                                  NULL, // cid
                                  NULL, // nip
                                  NULL, // ar
                                  NULL, // alpha
                                  NULL, // beta
                                  feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Write ASTROS element cards not supported by astros_writeMesh
int astros_writeSubElementCard(FILE *fp, meshStruct *feaMesh, int numProperty,
                               feaPropertyStruct *feaProperty,
                               feaFileFormatStruct *feaFileFormat)
{
    int status;
    int i, j; // Indexing
    int  fieldWidth;
    char *tempString = NULL;
    char *delimiter;

    double zoff;
    int *tm;

    int found;
    feaMeshDataStruct *feaData;

    if (numProperty > 0 && feaProperty == NULL) return CAPS_NULLVALUE;
    if (feaMesh == NULL) return CAPS_NULLVALUE;

    if (feaFileFormat->gridFileType == FreeField) {
        delimiter = ",";
        fieldWidth = 7;
    } else {
        delimiter = " ";
        fieldWidth = 7;
    }

    if (feaMesh->meshType == VolumeMesh) return CAPS_SUCCESS;

    for (i = 0; i < feaMesh->numElement; i++) {

        if (feaMesh->element[i].analysisType != MeshStructure) continue;

        feaData = (feaMeshDataStruct *) feaMesh->element[i].analysisData;

        found = (int) false;
        for (j = 0; j < numProperty; j++) {
            if (feaData->propertyID == feaProperty[j].propertyID) {
                found = (int) true;
                break;
            }
        }

        if (feaMesh->element[i].elementType == Node &&
            feaData->elementSubType == ConcentratedMassElement) {

            if (found == (int) false) {
                printf("No property information found for element %d of type \"ConcentratedMass\"!", feaMesh->element[i].elementID);
                continue;
            }

            status = astrosCard_conm2(fp,
                                      &feaMesh->element[i].elementID, // eid
                                      feaMesh->element[i].connectivity, // g
                                      &feaData->coordID, // cid
                                      &feaProperty[j].mass, // m
                                      feaProperty[j].massOffset, // X
                                      feaProperty[j].massInertia, // I
                                      NULL, // tmin
                                      NULL, // tmax
                                      feaFileFormat->gridFileType);
            if (status != CAPS_SUCCESS) return status;
        }

        if (feaMesh->element[i].elementType == Line)  {

            if (feaData->elementSubType == BarElement) {
                printf("Bar elements not supported yet - Sorry !\n");
                return CAPS_NOTIMPLEMENT;
            }

            if (feaData->elementSubType == BeamElement) {
                printf("Beam elements not supported yet - Sorry !\n");
                return CAPS_NOTIMPLEMENT;
            }

        }

        if ( feaMesh->element[i].elementType == Triangle) {

            if (feaData->elementSubType == ShellElement) {

                if (found == (int) false) {
                    printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
                    continue;
                }

                if (feaData->coordID != 0){
                    tm = &feaData->coordID;
                } else {
                    tm = NULL;
                }

                zoff = feaProperty[j].membraneThickness * feaProperty[j].zOffsetRel / 100.0;

                status = astrosCard_ctria3(
                    fp,
                    &feaMesh->element[i].elementID, // eid
                    &feaData->propertyID, // pid
                    feaMesh->element[i].connectivity, // Gi
                    Integer, tm, // tm,
                    &zoff, // zoff
                    NULL, // tmax
                    NULL, // Ti
                    feaFileFormat->gridFileType
                );
                if (status != CAPS_SUCCESS) return status;
            }
        }

        if ( feaMesh->element[i].elementType == Triangle_6) {

            if (feaData->elementSubType == ShellElement) {

                if (found == (int) false) {
                    printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
                    continue;
                }

                fprintf(fp,"%-8s%s%7d%s%7d%s%7d%s%7d%s%7d%s%7d%s%7d%s%7d%s%-8s\n", "CTRIA6",
                        delimiter, feaMesh->element[i].elementID,
                        delimiter, feaData->propertyID,
                        delimiter, feaMesh->element[i].connectivity[0],
                        delimiter, feaMesh->element[i].connectivity[1],
                        delimiter, feaMesh->element[i].connectivity[2],
                        delimiter, feaMesh->element[i].connectivity[3],
                        delimiter, feaMesh->element[i].connectivity[4],
                        delimiter, feaMesh->element[i].connectivity[5],
                        delimiter, "+CT");

                if (feaData->coordID != 0){
                    fprintf(fp, "%-8s%s%7d\n", "+CT", delimiter, feaData->coordID);
                } else {
                    fprintf(fp, "%-8s%s%7s\n", "+CT", delimiter, "");
                }

                // ZoffSet
                tempString = convert_doubleToString( feaProperty[j].membraneThickness*feaProperty[j].zOffsetRel/100.0, fieldWidth, 1);
                if (tempString == NULL) { status = EGADS_MALLOC; return status; }
                fprintf(fp, "%s%s",delimiter, tempString);
                AIM_FREE(tempString);

                fprintf(fp,"\n");
            }
        }

        if ( feaMesh->element[i].elementType == Quadrilateral) {

            if (feaData->elementSubType == ShearElement) {

                if (found == (int) false) {
                    printf("No property information found for element %d of type \"ShearElement\"!",
                           feaMesh->element[i].elementID);
                    continue;
                }

                status = astrosCard_cshear(fp,
                                           &feaMesh->element[i].elementID, // eid
                                           &feaData->propertyID, // pid
                                           feaMesh->element[i].connectivity, // Gi
                                           NULL, // tmax
                                           feaFileFormat->gridFileType);
                if (status != CAPS_SUCCESS) return status;

            }

            if (feaData->elementSubType == ShellElement) {

                if (found == (int) false) {
                    printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
                    continue;
                }

                if (feaData->coordID != 0) {
                    tm = &feaData->coordID;
                } else {
                    tm = NULL;
                }

                zoff = feaProperty[j].membraneThickness * feaProperty[j].zOffsetRel / 100.0;

                status = astrosCard_cquad4(fp,
                                           &feaMesh->element[i].elementID, // eid
                                           &feaData->propertyID, // pid
                                           feaMesh->element[i].connectivity, // Gi
                                           Integer, tm, // tm,
                                           &zoff, // zoff
                                           NULL, // tmax
                                           NULL, // Ti
                                           feaFileFormat->gridFileType);
                if (status != CAPS_SUCCESS) return status;
            }

            if (feaData->elementSubType == MembraneElement) {

                if (found == (int) false) {
                    printf("No property information found for element %d of type \"MembraneElement\"!", feaMesh->element[i].elementID);
                    continue;
                }

                if (feaData->coordID != 0) {
                    tm = &feaData->coordID;
                } else {
                    tm = NULL;
                }

                status = astrosCard_cqdmem1(fp,
                                            &feaMesh->element[i].elementID, // eid
                                            &feaData->propertyID, // pid
                                            feaMesh->element[i].connectivity, // Gi
                                            Integer, tm, // tm,
                                            NULL, // tmax
                                            feaFileFormat->gridFileType);
                if (status != CAPS_SUCCESS) return status;
            }
        }

        if ( feaMesh->element[i].elementType == Quadrilateral_8) {

            if (feaData->elementSubType == ShellElement) {

                if (found == (int) false) {
                    printf("No property information found for element %d of type \"ShellElement\"!", feaMesh->element[i].elementID);
                    continue;
                }

                fprintf(fp,"%-8s%s%7d%s%7d%s%7d%s%7d%s%7d%s%7d%s%7d%s%7d%s%-8s\n", "CQUAD8",
                        delimiter, feaMesh->element[i].elementID,
                        delimiter, feaData->propertyID,
                        delimiter, feaMesh->element[i].connectivity[0],
                        delimiter, feaMesh->element[i].connectivity[1],
                        delimiter, feaMesh->element[i].connectivity[2],
                        delimiter, feaMesh->element[i].connectivity[3],
                        delimiter, feaMesh->element[i].connectivity[4],
                        delimiter, feaMesh->element[i].connectivity[5],
                        delimiter, "+CQ");

                fprintf(fp,"%-8s%s%7d%s%7d", "+CQ",
                        delimiter, feaMesh->element[i].connectivity[6],
                        delimiter, feaMesh->element[i].connectivity[7]);

                fprintf(fp, "%s%7s%s%7s%s%7s%s%7s", delimiter, "",
                                                    delimiter, "",
                                                    delimiter, "",
                                                    delimiter, "");

                // Write coordinate id
                if (feaData->coordID != 0){
                    fprintf(fp, "%s%7d", delimiter, feaData->coordID);
                } else {
                    fprintf(fp, "%s%7s", delimiter, "");
                }

                // ZoffSet
                tempString = convert_doubleToString( feaProperty[j].membraneThickness*feaProperty[j].zOffsetRel/100.0,
                                                    fieldWidth, 1);
                if (tempString == NULL) { status = EGADS_MALLOC; return status; }
                fprintf(fp, "%s%s",delimiter, tempString);
                AIM_FREE(tempString);

                fprintf(fp,"\n");
            }
        }
    }

    return CAPS_SUCCESS;
}

static int _getTrimType(feaAnalysisStruct *feaAnalysis, int *trimType)
{
    // SYMMETRY
    // SYM (0)
    // ANTISYM (-1)
    // ASYM (1)
    if (feaAnalysis->aeroSymmetryXY == NULL ) {
        printf("\t*** Warning *** aeroSymmetryXY Input to AeroelasticTrim Analysis in astrosAIM not defined! Using ASYMMETRIC\n");
        *trimType = 1;
    } else {
        if(strcmp("SYM",feaAnalysis->aeroSymmetryXY) == 0) {
            *trimType = 0;
        } else if(strcmp("ANTISYM", feaAnalysis->aeroSymmetryXY) == 0) {
            *trimType = -1;
        } else if(strcmp("ASYM", feaAnalysis->aeroSymmetryXY) == 0) {
            *trimType = 1;
        } else if(strcmp("SYMMETRIC", feaAnalysis->aeroSymmetryXY) == 0) {
            *trimType = 0;
        } else if(strcmp("ANTISYMMETRIC", feaAnalysis->aeroSymmetryXY) == 0) {
            *trimType = -1;
        } else if(strcmp("ASYMMETRIC", feaAnalysis->aeroSymmetryXY) == 0) {
            *trimType = 1;
        } else {
            printf("\t*** Warning *** aeroSymmetryXY Input %s to astrosAIM not understood! Using ASYMMETRIC\n", feaAnalysis->aeroSymmetryXY );
            *trimType = 1;
        }
    }

    return CAPS_SUCCESS;
}

// get the ASTROS label equivalent of trim parameter `label` and corresponding `symmetryFlag`
static int _getAstrosTrimParamLabelAndSymmetry(const char *label,
                                               char **astrosLabel,
                                               int *symmetryFlag)
{
    // ASTROS VARIABLES
    //SYM: 'NX','NZ','QACCEL','ALPHA','QRATE','THKCAM'
    //ANTISYM: 'NY','PACCEL','RACCEL','BETA','PRATE','RRATE'

    // NASTRAN VARIABLES
    //SYM: URDD1, URDD3, URDD5, ANGLEA, PITCH
    //ANTISYM: URDD2, URDD4, URDD6, SIDES, ROLL, YAW

    // SYMMETRY
    // SYM (0)
    // ANTISYM (-1)

    if(strcmp("URDD1", label) == 0) {
         *astrosLabel = EG_strdup("NX");
         *symmetryFlag = 0;
    } else if(strcmp("URDD2", label) == 0) {
         *astrosLabel = EG_strdup("NY");
         *symmetryFlag = -1;
    } else if(strcmp("URDD3", label) == 0) {
         *astrosLabel = EG_strdup("NZ");
         *symmetryFlag = 0;
    } else if(strcmp("URDD4", label) == 0) {
         *astrosLabel = EG_strdup("PACCEL");
         *symmetryFlag = -1;
    } else if(strcmp("URDD5", label) == 0) {
         *astrosLabel = EG_strdup("QACCEL");
         *symmetryFlag = 0;
    } else if(strcmp("URDD6", label) == 0) {
         *astrosLabel = EG_strdup("RACCEL");
         *symmetryFlag = -1;
    } else if(strcmp("ANGLEA", label) == 0) {
         *astrosLabel = EG_strdup("ALPHA");
         *symmetryFlag = 0;
    } else if(strcmp("PITCH", label) == 0) {
         *astrosLabel = EG_strdup("QRATE");
         *symmetryFlag = 0;
    } else if(strcmp("SIDES", label) == 0) {
         *astrosLabel = EG_strdup("BETA");
         *symmetryFlag = -1;
    } else if(strcmp("ROLL", label) == 0) {
         *astrosLabel = EG_strdup("PRATE");
         *symmetryFlag = -1;
    } else if(strcmp("YAW", label) == 0) {
         *astrosLabel = EG_strdup("RRATE");
         *symmetryFlag = -1;
    } else {
         *astrosLabel = EG_strdup(label);
        // printf("\t*** Warning *** rigidConstraint Input %s to astrosAIM not understood!\n", label);
         *symmetryFlag = 1;
         return CAPS_NOTFOUND;
    }

    return CAPS_SUCCESS;
}

// Write a Astros Analysis card from a feaAnalysis structure
int astros_writeAnalysisCard(FILE *fp, feaAnalysisStruct *feaAnalysis,
                             feaFileFormatStruct *feaFileFormat)
{
    int i; // Indexing
    double velocity, dv, vmin, vmax, velocityArray[21], *vo;

    int symmetryFlag, trimType, symxy, symxz, status, velID, denID;

    char *tempString;

    int numParams, paramIndex;
    astrosCardTrimParamStruct *trimParams = NULL;

    if (fp == NULL) return CAPS_IOERR;
    if (feaAnalysis == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // Eigenvalue
    if (feaAnalysis->analysisType == Modal ||
        feaAnalysis->analysisType == AeroelasticFlutter) {

        if (feaAnalysis->eigenNormaliztion == NULL) {
            PRINT_ERROR("Undefined eigen normalization method");
            return CAPS_NULLVALUE;
        }

        status = astrosCard_eigr(fp,
                                 &feaAnalysis->analysisID, // sid
                                 feaAnalysis->extractionMethod, // method
                                 &feaAnalysis->frequencyRange[0], // fl
                                 &feaAnalysis->frequencyRange[1], // fu
                                 &feaAnalysis->numEstEigenvalue, // ne
                                 &feaAnalysis->numDesiredEigenvalue, // nd
                                 NULL, // e
                                 feaAnalysis->eigenNormaliztion, // norm
                                 &feaAnalysis->gridNormaliztion, // gid
                                 &feaAnalysis->componentNormaliztion, // dof
                                 feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;
    }

    if (feaAnalysis->analysisType == AeroelasticTrim) {

        _getTrimType(feaAnalysis, &trimType);

        // ASTROS VARIABLES
        //SYM: 'NX','NZ','QACCEL','ALPHA','QRATE','THKCAM'
        //ANTISYM: 'NY','PACCEL','RACCEL','BETA','PRATE','RRATE'

        // NASTRAN VARIABLES
        //SYM: URDD1, URDD3, URDD5, ANGLEA, PITCH
        //ANTISYM: URDD2, URDD4, URDD6, SIDES, ROLL, YAW
        fprintf(fp,"%s\n","$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|");
        fprintf(fp,"%s\n","$TRIM   TRIMID  MACH    QDP     TRMTYP  EFFID   VO                      CONT");

        if(feaAnalysis->density > 0.0) {
            velocity = sqrt(2*feaAnalysis->dynamicPressure/feaAnalysis->density);
            vo = &velocity;
        }
        else {
            vo = NULL;
        }

        numParams = (feaAnalysis->numRigidConstraint + 1
                     + feaAnalysis->numRigidVariable
                     + feaAnalysis->numControlConstraint);

        trimParams = EG_alloc(sizeof(astrosCardTrimParamStruct) * numParams);
        if (trimParams == NULL) {
            return EGADS_MALLOC;
        }

        paramIndex = 0;

        // Rigid Constraints
        for (i = 0; i < feaAnalysis->numRigidConstraint; i++) {

            status = _getAstrosTrimParamLabelAndSymmetry(
                feaAnalysis->rigidConstraint[i], &tempString, &symmetryFlag);
            if (status == CAPS_NOTFOUND) {
                printf("\t*** Warning *** rigidConstraint Input %s to astrosAIM not understood!\n", tempString);
            }

            if (trimType == 1 || trimType == symmetryFlag) {

                astrosCard_initiateTrimParam(&trimParams[paramIndex]);
                trimParams[paramIndex].label = tempString;
                trimParams[paramIndex++].value = feaAnalysis->magRigidConstraint[i];
            }
            tempString = NULL;
        }

        // Add THKCAM for ASYM and SYM Cases
        if (trimType != -1) {

            astrosCard_initiateTrimParam(&trimParams[paramIndex]);
            trimParams[paramIndex].label = EG_strdup("THKCAM");
            trimParams[paramIndex++].value = 1.0;
        }
        else {
            numParams--;
        }

        // Rigid Variables
        for (i = 0; i < feaAnalysis->numRigidVariable; i++) {

            status = _getAstrosTrimParamLabelAndSymmetry(
                     feaAnalysis->rigidVariable[i], &tempString, &symmetryFlag);
            if (status == CAPS_NOTFOUND) {
/*@-nullpass@*/
                printf("\t*** Warning *** rigidVariable Input %s to astrosAIM not understood!\n",
                       tempString);
/*@+nullpass@*/
            }

            if (trimType == 1 || trimType == symmetryFlag) {

                astrosCard_initiateTrimParam(&trimParams[paramIndex]);
                trimParams[paramIndex].label = tempString;
                trimParams[paramIndex++].isFree = (int) true;
            }
            tempString = NULL;
        }

        // Control Constraints
        for (i = 0; i < feaAnalysis->numControlConstraint; i++) {

            astrosCard_initiateTrimParam(&trimParams[paramIndex]);
            trimParams[paramIndex].label = EG_strdup(feaAnalysis->controlConstraint[i]);
            trimParams[paramIndex++].value = feaAnalysis->magControlConstraint[i];
        }

        status = astrosCard_trim(fp,
                                 &feaAnalysis->analysisID, // trimid
                                 feaAnalysis->machNumber, // mach
                                 &feaAnalysis->dynamicPressure, // qdp
                                 NULL, // trmtyp
                                 NULL, // effid
                                 vo, // vo
                                 numParams, trimParams, // trim parameters
                                 feaFileFormat->fileType);

        for (i = 0; i < numParams; i++) {
            astrosCard_destroyTrimParam(&trimParams[i]);
        }
        EG_free(trimParams);

        if (status != CAPS_SUCCESS) return status;
    }

    if (feaAnalysis->analysisType == AeroelasticFlutter) {
        if (feaAnalysis->aeroSymmetryXY == NULL ) {
            printf("\t*** Warning *** aeroSymmetryXY Input to AeroelasticFlutter Analysis in astrosAIM not defined! Using ASYMMETRIC\n");
            symxy = 0;
        } else {
            if(strcmp("SYM",feaAnalysis->aeroSymmetryXY) == 0) {
                symxy = 0;
            } else if(strcmp("ANTISYM",feaAnalysis->aeroSymmetryXY) == 0) {
                symxy = -1;
            } else if(strcmp("ASYM",feaAnalysis->aeroSymmetryXY) == 0) {
                symxy = 1;
            } else if(strcmp("SYMMETRIC",feaAnalysis->aeroSymmetryXY) == 0) {
                symxy = 0;
            } else if(strcmp("ANTISYMMETRIC",feaAnalysis->aeroSymmetryXY) == 0) {
                symxy = -1;
            } else if(strcmp("ASYMMETRIC",feaAnalysis->aeroSymmetryXY) == 0) {
                symxy = 1;
            } else {
                printf("\t*** Warning *** aeroSymmetryXY Input %s to astrosAIM not understood! Using ASYMMETRIC\n",feaAnalysis->aeroSymmetryXY );
                symxy = 0;
            }
        }

        if (feaAnalysis->aeroSymmetryXZ == NULL ) {
            printf("\t*** Warning *** aeroSymmetryXY Input to AeroelasticFlutter Analysis in astrosAIM not defined! Using ASYMMETRIC\n");
            symxz = 0;
        } else {
            if(strcmp("SYM",feaAnalysis->aeroSymmetryXZ) == 0) {
                symxz = 0;
            } else if(strcmp("ANTISYM",feaAnalysis->aeroSymmetryXZ) == 0) {
                symxz = -1;
            } else if(strcmp("ASYM",feaAnalysis->aeroSymmetryXZ) == 0) {
                symxz = 1;
            } else if(strcmp("SYMMETRIC",feaAnalysis->aeroSymmetryXZ) == 0) {
                symxz = 0;
            } else if(strcmp("ANTISYMMETRIC",feaAnalysis->aeroSymmetryXZ) == 0) {
                symxz = -1;
            } else if(strcmp("ASYMMETRIC",feaAnalysis->aeroSymmetryXZ) == 0) {
                symxz = 1;
            } else {
                printf("\t*** Warning *** aeroSymmetryXZ Input %s to astrosAIM not understood! Using ASYMMETRIC\n",feaAnalysis->aeroSymmetryXZ );
                symxz = 0;
            }
        }

        if (feaAnalysis->numMachNumber != 0) {
            if (feaAnalysis->numMachNumber > 6) {
                printf("\t*** Warning *** Mach number input for AeroelasticFlutter in astrosAIM must be less than six\n");
            }
        }

        if (feaAnalysis->numReducedFreq != 0) {
            if (feaAnalysis->numReducedFreq > 8) {
                printf("\t*** Warning *** Reduced freq. input for AeroelasticFlutter in astrosAIM must be less than eight\n");
            }
        }

        // Write MKAERO1 INPUT
        status = astrosCard_mkaero1(
            fp,
            &symxz, // symxz
            &symxy, // symxy
            feaAnalysis->numMachNumber, feaAnalysis->machNumber, // M
            feaAnalysis->numReducedFreq, feaAnalysis->reducedFreq, // K
            feaFileFormat->fileType
        );
        if (status != CAPS_SUCCESS) return status;

        fprintf(fp,"%s\n","$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|");
        fprintf(fp,"%s",  "$LUTTER SID     METHOD  DENS    MACH    VEL     MLIST   KLIST   EFFID   CONT\n");
        fprintf(fp,"%s\n","$CONT   SYMXZ   SYMXY   EPS     CURFIT");

        denID = 10 * feaAnalysis->analysisID + 1;
        velID = 10 * feaAnalysis->analysisID + 2;

        status = astrosCard_flutter(fp,
                                    &feaAnalysis->analysisID, // sid
                                    "PK", // method
                                    &denID, // dens
                                    &feaAnalysis->machNumber[0], // mach
                                    &velID, // vel
                                    NULL, // mlist
                                    NULL, // klist
                                    NULL, // effid
                                    &symxz, // symxz
                                    &symxy, // symxy
                                    NULL, // eps
                                    NULL, // curfit
                                    NULL, // nroot
                                    NULL, // vtype
                                    NULL, // gflut
                                    NULL, // gfilter
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        fprintf(fp, "$\n");

        // DENS FLFACT

        status = astrosCard_flfact(fp, &denID, 1, &feaAnalysis->density, feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        // VEL FLFACT
        velocity = sqrt(2*feaAnalysis->dynamicPressure/feaAnalysis->density);
        vmin = velocity / 2.0;
        vmax = 2.0 * velocity;
        dv = (vmax - vmin) / 20.0;

        for (i = 0; i < 21; i++) {
            velocityArray[i] = vmin + (double) i * dv;
        }

        status = astrosCard_flfact(fp, &velID, 21, velocityArray, feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) return status;

        fprintf(fp, "$\n");

    }

    return CAPS_SUCCESS;
}

// Write design variable/optimization information from a feaDesignVariable structure
int astros_writeDesignVariableCard(FILE *fp,
                                   feaDesignVariableStruct *feaDesignVariable,
                                   int numProperty,
                                   feaPropertyStruct feaProperty[],
                                   feaFileFormatStruct *feaFileFormat)
{
    int  i; //Indexing
    int len, composite, numPly = 0;
    long intVal = 0;
    int *layers = NULL, *layrlst = NULL;
    double vmin = 0.0, vmax = 1.0, vinit;

    int status; // Function return status

    char *label = NULL, *propertyType = NULL;
    char *copy;

    composite = (int) false;

    // if (feaDesignVariable->designVariableType != PropertyDesignVar) {
    //     PRINT_ERROR("For ASTROS Optimization designVariableType must be a property not a material");
    //     return CAPS_BADVALUE;
    // }

    if (feaDesignVariable->numDiscreteValue == 0) {
        // DESVARP BID, LINKID, VMIN, VMAX, VINIT, LAYERNUM, LAYRLST, LABEL

        if ( feaDesignVariable->initialValue != 0.0) {
            vmin = feaDesignVariable->lowerBound / feaDesignVariable->initialValue;
            vmax = feaDesignVariable->upperBound / feaDesignVariable->initialValue;
            vinit = feaDesignVariable->initialValue / feaDesignVariable->initialValue;
        } else {
            vmin = 0.0;
            vmax = 1.0;
            vinit = feaDesignVariable->initialValue;
        }

        if (feaDesignVariable->propertySetType != NULL) {
            if (feaDesignVariable->propertySetType[0] == Composite) {
                layrlst = &feaDesignVariable->designVariableID;
            }
        }

        printf("*** WARNING *** For ASTROS Optimization design variable linking has not been implemented yet\n");

        len = strlen(feaDesignVariable->name);
        if (len > 7) {
            label = EG_alloc(sizeof(char) * 9);
            if (label != NULL)
              snprintf(label, 8, "VAR%d", feaDesignVariable->designVariableID);
        } else {
            label = EG_strdup(feaDesignVariable->name);
        }

        status = astrosCard_desvarp(fp,
                                    &feaDesignVariable->designVariableID, // dvid
                                    &feaDesignVariable->designVariableID, // linkid
                                    &vmin, &vmax, &vinit, // vmin, vmax, vinit
                                    NULL, // layrnum
                                    layrlst, // layrlst
                                    label, // label
                                    feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) goto cleanup;

    } else {
        PRINT_ERROR("For ASTROS Optimization designVariables can not be Discrete Values");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    for (i = 0; i < feaDesignVariable->numPropertyID; i++) {
        // PLIST, LINKID, PTYPE, PID1, ...
        // PTYPE =  PROD, PSHEAR, PCOMP, PCOMP1, PCOMP2, PELAS, PSHELL, PMASS, PTRMEM, PQDMEM1, PBAR

        if (feaDesignVariable->propertySetType == NULL) {
            printf("*** WARNING *** For ASTROS Optimization designVariable name \"%s\", propertySetType not set. PLIST entries not written\n",
                   feaDesignVariable->name);
            continue;
        }

        // UnknownProperty, Rod, Bar, Shear, Shell, Composite, Solid
        if (feaDesignVariable->propertySetType[i] == Rod) {
            propertyType = "PROD";
        }
        else if (feaDesignVariable->propertySetType[i] == Bar) {
            propertyType = "PBAR";
        }
        else if (feaDesignVariable->propertySetType[i] == Shell) {
            propertyType = "PSHELL";
        }
        else if (feaDesignVariable->propertySetType[i] == Membrane) {
            propertyType = "PQDMEM1";
        }
        else if (feaDesignVariable->propertySetType[i] == Shear) {
            propertyType = "PSHEAR";
        }
        else if (feaDesignVariable->propertySetType[i] == Composite) {
            propertyType = "PCOMP";
            composite = (int) true;
        }
        else if (feaDesignVariable->propertySetType[i] == Solid) {
            PRINT_ERROR("For ASTROS Optimization designVariables can not relate to PSOLID property types");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        status = astrosCard_plist(fp,
                                  &feaDesignVariable->designVariableID,
                                  propertyType,
                                  1, &feaDesignVariable->propertySetID[i],
                                  feaFileFormat->fileType);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    if (composite == (int) true) {
        // Check the field input
        if (feaDesignVariable->fieldName == NULL) {
            PRINT_ERROR("For ASTROS Optimization designVariables must have fieldName defined");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        // Check if angle is input (i.e. not lamina thickness)
        if (strncmp(feaDesignVariable->fieldName, "THETA", 5) == 0) {
            PRINT_ERROR("For ASTROS Optimization designVariables, fieldName can not be an angle (i.e. THETAi)");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        // Search all properties to determine the number of layers in the composite
        for (i = 0; i < numProperty; i++) {
            if (feaDesignVariable->propertySetID == NULL) {
                printf("*** WARNING *** For ASTROS Optimization designVariable name \"%s\", propertySetID not set.\n",
                       feaDesignVariable->name);
                continue;
            }

            if (feaDesignVariable->propertySetID[0] == feaProperty[i].propertyID) {
                numPly = feaProperty[i].numPly;
                if (feaProperty->compositeSymmetricLaminate == (int) true) {
                    numPly = 2 * numPly;
                }
                break;
            }
        }

        if (strcmp(feaDesignVariable->fieldName,"TALL") == 0) {

            layers = (int *) EG_alloc((numPly)*sizeof(int));
            if (layers == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }
            for (i = 0; i < numPly; i++) {
                layers[i] = i+1;
            }

            status = astrosCard_plylist(
                fp,
                &feaDesignVariable->designVariableID, // sid
                numPly, layers, // P
                feaFileFormat->fileType
            );
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        if (strncmp(feaDesignVariable->fieldName, "T", 1) == 0) {

            if (strncmp(feaDesignVariable->fieldName, "THETA", 5) == 0 ||
                strncmp(feaDesignVariable->fieldName, "TALL",  4) == 0) {
                // do nothing
            } else {
                // Input is T1, T2, etc or T11 ... Need to print the integer part

                copy = feaDesignVariable->fieldName;

                while (*copy) { // While there are more characters to process...
                    if (isdigit(*copy)) { // Upon finding a digit, ...
                        intVal = strtol(copy, &copy, 10); // Read a number, ...
                        //fprintf(fp,"***** OUTPUT %d\n", (int) intVal); // and print it.
                    } else { // Otherwise, move on to the next character.
                        copy++;
                    }
                }

                //len = strlen(feaDesignVariable->fieldName);

                if (layers != NULL) EG_free(layers);
                if (feaProperty->compositeSymmetricLaminate == (int) true) {
                    // need to add sym laminate layer to the PLYLIST
                    // numPly - total including sym muliplyer
                    // intVal - selected play for 1/2 the stack
                    // otherside = numPly - intVal + 1

                    numPly = 2;
                    layers = (int *) EG_alloc(numPly * sizeof(int));
                    if (layers == NULL) {
                        status = EGADS_MALLOC;
                        goto cleanup;
                    }
                    layers[0] = (int) intVal;
                    layers[1] = numPly + 1 - (int) intVal;
                } else {
                    numPly = 1;
                    layers = (int *) EG_alloc(numPly * sizeof(int));
                    if (layers == NULL) {
                        status = EGADS_MALLOC;
                        goto cleanup;
                    }
                    layers[0] = (int) intVal;
                }

                status = astrosCard_plylist(fp,
                                            &feaDesignVariable->designVariableID, // sid
                                            numPly, layers, // P
                                            feaFileFormat->fileType);
                if (status != CAPS_SUCCESS) goto cleanup;

                // if(copy != NULL) EG_free(copy);
            }
        }

    }

    if (feaDesignVariable->numIndependVariable > 0) {
        printf("*** WARNING *** For ASTROS Optimization designVariable name \"%s\", propertySetID not set.\n",
               feaDesignVariable->name);
    }

    /*
    if (feaDesignVariable->numIndependVariable > 0) {
        c = 1;
                // DLINK, ID,  DDVID,  CO,  CMULT,  IDV1,  C1,  IDV2,  C2
                fprintf(fp,"%-8s", "DLINK");

                tempString = convert_integerToString(feaDesignVariable->designVariableID+10000, fieldWidth, 1);
                fprintf(fp, "%s%s",delimiter, tempString);
                EG_free(tempString);

                tempString = convert_integerToString(feaDesignVariable->designVariableID, fieldWidth, 1);
                fprintf(fp, "%s%s",delimiter, tempString);
                EG_free(tempString);

                tempString = convert_doubleToString(feaDesignVariable->variableWeight[0], fieldWidth, 1);
                fprintf(fp, "%s%s",delimiter, tempString);
                EG_free(tempString);

                tempString = convert_doubleToString(feaDesignVariable->variableWeight[1], fieldWidth, 1);
                fprintf(fp, "%s%s",delimiter, tempString);
                EG_free(tempString);

                for (i = 0; i < feaDesignVariable->numIndependVariable; i++) {

                        tempString = convert_integerToString(feaDesignVariable->independVariableID[i], fieldWidth, 1);
                        fprintf(fp, "%s%s",delimiter, tempString);
                        EG_free(tempString);

                        tempString = convert_doubleToString(feaDesignVariable->independVariableWeight[i], fieldWidth, 1);
                        fprintf(fp, "%s%s",delimiter, tempString);
                        EG_free(tempString);

                        if (i > feaDesignVariable->numIndependVariable-1 && i == c) {
                                c = c + 4;
                                if (feaFileFormat->fileType == FreeField) {
                                        fprintf(fp, ",+C\n+, ");
                                } else {
                                        fprintf(fp, "+C%6s\n%8s", "", "");
                                }
                        } else if (i == feaDesignVariable->numIndependVariable-1) {
                                fprintf(fp,"\n");
                        }

                }
    }
     */

    status = CAPS_SUCCESS;

cleanup:

    if (label  != NULL) EG_free(label);
    if (layers != NULL) EG_free(layers);

    return status;
}

// Write design constraint/optimization information from a feaDesignConstraint structure
int astros_writeDesignConstraintCard(FILE *fp, int feaDesignConstraintSetID,
                                     feaDesignConstraintStruct *feaDesignConstraint,
                                     int numMaterial, feaMaterialStruct feaMaterial[],
                                     int numProperty, feaPropertyStruct feaProperty[],
                                     feaFileFormatStruct *feaFileFormat)
{
    int  i, j; // Index
    int  iPID = 0, iMID = 0;
    long intVal = 0;
    double defaultF12 = 0.0, shearStress;

    char *copy;

    for (i = 0; i < feaDesignConstraint->numPropertyID; i++) {

        if (feaDesignConstraint->propertySetType[i] == Rod) {
            // DCONVMP, SID, ST, SC, SS, PTYPE, LAYRNUM, PID1, PID2, CONT
            // CONT, PID2, PID4, ETC ...
            /* SID Stress constraint set identification (Integer > 0).
             * ST Tensile stress limit (Real > 0.0 or blank)
             * Sc Compressive stress limit (Real, Default = ST)
             * ss Shear stress limit (Real > 0.0 or blank)
             * PTYPE Property type (Text)
             * LAYPNUM The layer number of a composite element (Integer > 0 or blank)
             * PIDi Property identification numbers (Integer > 0)
             */

            return astrosCard_dconvmp(fp,
                                      &feaDesignConstraintSetID, // sid
                                      &feaDesignConstraint->upperBound, // st
                                      NULL, // sc
                                      NULL, // ss
                                      "PROD", // ptype
                                      NULL, // layrnum
                                      1,
                                      &feaDesignConstraint->propertySetID[i], // PIDi (only one right now ?)
                                      feaFileFormat->fileType);

        } else if (feaDesignConstraint->propertySetType[i] == Bar) {
            // Nothing set yet

        } else if (feaDesignConstraint->propertySetType[i] == Shell) {
            // DCONVMP, SID, ST, SC, SS, PTYPE, LAYRNUM, PID1, PID2, CONT
            // CONT, PID2, PID4, ETC ...

            shearStress = feaDesignConstraint->upperBound / 2.0;

            return astrosCard_dconvmp(fp,
                                      &feaDesignConstraintSetID, // sid
                                      &feaDesignConstraint->upperBound, // st
                                      NULL, // sc
                                      &shearStress, // ss
                                      "PSHELL", // ptype
                                      NULL, // layrnum
                                      1,
                                      &feaDesignConstraint->propertySetID[i], // PIDi (only one right now ?)
                                      feaFileFormat->fileType);

        } else if (feaDesignConstraint->propertySetType[i] == Composite) {
            /*
             * DCONTWP SID XT XC YT YC SS F12 PTYPE ICONT
             * CONT LAYRNUM PIDI ID2 PID3 -etc-
             * SID Stress constraint set identification (Integer > 0)
             * XT Tensile stress limit in the longitudinal direction (Real > 0.0)
             * XC Compressive stress limit in the longitudinal direction (Real, Default = XT)
             * YT Tensile stress limit in the transverse direction (Real > 0.0)
             * YC Compressive stress limit in the transverse direction (Real, Default = YT)
             * ss Shear stress limit for in-plane stress (Real > 0.0)
             * F12 Tsai-Wu interaction term (Real)
             * PTYPE Property type (Text)
             * LAYRNUM The layer number of a composite element (Integer > 0 or blank)
             * PIDi Property identification numbers (Integer > 0)
             */

            // Property ID for the constraint: feaDesignConstraint->propertySetID[i]
            // LIST of properties and their materials
            // numProperty, feaProperty[].propertyID
            //                          feaProperty[].materialID
            // LIST of materials
            // numMaterial, feaMaterial[].materialID
            //                          feaMaterial[].tensionAllow
            //                          feaMaterial[].compressAllow
            //                          feaMaterial[].tensionAllowLateral
            //                          feaMaterial[].compressAllowLateral
            //                          feaMaterial[].shearAllow

            for (j = 0; j < numProperty; j++) {
                if (feaDesignConstraint->propertySetID[i] == feaProperty[j].propertyID) {
                    iPID = j;
                    break;
                }
            }
            for (j = 0; j < numMaterial; j++) {
                if (feaProperty[iPID].materialID == feaMaterial[j].materialID) {
                    iMID = j;
                    break;
                }
            }

            // get layrnum
            copy = feaDesignConstraint->fieldName;

            if (copy == NULL) {
                return CAPS_NULLVALUE;
            }

            while (*copy) { // While there are more characters to process...
                if (isdigit(*copy)) { // Upon finding a digit, ...
                    intVal = strtol(copy, &copy, 10); // Read a number, ...
                    //fprintf(fp,"***** OUTPUT %d\n", (int) intVal); // and print it.
                } else { // Otherwise, move on to the next character.
                    copy++;
                }
            }

            return astrosCard_dcontwp(fp,
                                      &feaDesignConstraintSetID, // sid
                                      &feaMaterial[iMID].tensionAllow, // xt
                                      &feaMaterial[iMID].compressAllow, // xc
                                      &feaMaterial[iMID].tensionAllowLateral, // yt
                                      &feaMaterial[iMID].compressAllowLateral, // yc
                                      &feaMaterial[iMID].shearAllow, // ss
                                      &defaultF12, // f12
                                      "PCOMP", // ptype
                                      (int *) &intVal, // layrnum
                                      1, &feaDesignConstraint->propertySetID[i], // PIDi (only one right now ?)
                                      feaFileFormat->fileType);

        } else if (feaDesignConstraint->propertySetType[i] == Solid) {
            // Nothing set yet
        }
    }

    return CAPS_SUCCESS;
}

// Read data from a Astros OUT file to determine the number of eignevalues
int astros_readOUTNumEigenValue(FILE *fp, int *numEigenVector)
{
    int status; // Function return

    char *beginEigenLine = "                                   S U M M A R Y   O F   R E A L   E I G E N   A N A L Y S I S";

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    int tempInt[2];

    while (*numEigenVector == 0) {
        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // See how many Eigen-Values we have
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            // Skip ahead 2 lines
            status = getline(&line, &linecap, fp);
            if (status < 0) break;
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            // Grab summary line
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            sscanf(line, "%d EIGENVALUES AND %d EIGENVECTORS",
                   &tempInt[0], &tempInt[1]);
            *numEigenVector = tempInt[1];
        }
    }

    if (line != NULL) EG_free(line);

    // Rewind the file
    rewind(fp);

    if (*numEigenVector == 0) return CAPS_NOTFOUND;
    else return CAPS_SUCCESS;
}

// Read data from a Astros OUT file to determine the number of grid points
int astros_readOUTNumGridPoint(FILE *fp, int *numGridPoint)
{
    int status; // Function return status

    int i, j;

    size_t linecap;

    char *line = NULL; // Temporary line holder

    char *beginEigenLine = "            EIGENVALUE       =";
    char *endEigenLine   = "            EIGENVALUE       =";

    int stop = (int) false;

    *numGridPoint = 0;

    // Loop through file line by line until we have determined how many grid points we have
    while (*numGridPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // Look for start of Eigen-Vector 1
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            // Fast forward 3 lines
            for (i = 0; i < 3; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through lines counting the number of grid points
            while (stop == (int) false) {

                status = getline(&line, &linecap, fp);
                if (status < 0)  break;

                // If we have a new page - skip ahead 8 lines and continue
                if (strncmp("1", line, 1) == 0) {

                    for (j = 0; j < 7; j++) {
                        status = getline(&line, &linecap, fp);
                        if (status < 0) break;

                        if (strncmp(endEigenLine, line, strlen(endEigenLine)) == 0)
                            stop = (int) true;

                    }
                    continue;
                }

                if (strncmp(endEigenLine, line, strlen(endEigenLine)) == 0 ||
                    strlen(line) == 1) break;

                *numGridPoint +=1;
            }
        }
    }

    if (line != NULL) EG_free(line);

    // Rewind the file
    rewind(fp);

    if (*numGridPoint == 0) return CAPS_NOTFOUND;
    else return CAPS_SUCCESS;
}


// Read data from a Astros OUT file and load it into a dataMatrix[numEigenVector][numGridPoint*8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int astros_readOUTEigenVector(FILE *fp, int *numEigenVector, int *numGridPoint,
                              double ***dataMatrix)
{
    int status = CAPS_SUCCESS; // Function return

    int i, j, eigenValue = 0; // Indexing

    size_t linecap;

    char *line = NULL; // Temporary line holder

    char *beginEigenLine = "            EIGENVALUE       =";

    int tempInt;
    char tempString[2];

    int numVariable = 8; // Grid Id, Coord Id, T1, T2, T3, R1, R2, R3

    printf("Reading Astros OUT file - extracting Eigen-Vectors!\n");

    *numEigenVector = 0;
    *numGridPoint = 0;

    // See how many Eigen-Values we have
    status = astros_readOUTNumEigenValue(fp, numEigenVector);
    printf("\tNumber of Eigen-Vectors = %d\n", *numEigenVector);
    if (status != CAPS_SUCCESS) return status;

    status = astros_readOUTNumGridPoint(fp, numGridPoint);
    printf("\tNumber of Grid Points = %d for each Eigen-Vector\n",
           *numGridPoint);
    if (status != CAPS_SUCCESS) return status;

    // Allocate dataMatrix array
    if (*dataMatrix != NULL) EG_free(*dataMatrix);

    *dataMatrix = (double **) EG_alloc((*numEigenVector) *sizeof(double *));
    if (*dataMatrix == NULL) {
        status = EGADS_MALLOC; // If allocation failed ....
        goto cleanup;
    }

    for (i = 0; i < *numEigenVector; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc((*numGridPoint)*
                                               numVariable*sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

                if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
            }

            if ((*dataMatrix) != NULL) EG_free((*dataMatrix));

            status = EGADS_MALLOC;
            goto cleanup;
        }
    }

    eigenValue = 0;
    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {
        if (line == NULL) break;

        // Look for start of Eigen-Vector
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            printf("\tLoading Eigen-Vector = %d\n", eigenValue+1);

            // Fast forward 3 lines
            for (i = 0; i < 3; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            i = 0;
            while (i != *numGridPoint) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;

                // If we have a new page - skip ahead 7 lines and continue
                if (strncmp("1", line, 1) == 0) {

                    for (j = 0; j < 7; j++) {
                        status = getline(&line, &linecap, fp);
                        if (status < 0) break;
                    }
                    continue;
                }

                sscanf(line, "%d%s%lf%lf%lf%lf%lf%lf", &tempInt,
                       tempString,
                       &(*dataMatrix)[eigenValue][2+numVariable*i],
                       &(*dataMatrix)[eigenValue][3+numVariable*i],
                       &(*dataMatrix)[eigenValue][4+numVariable*i],
                       &(*dataMatrix)[eigenValue][5+numVariable*i],
                       &(*dataMatrix)[eigenValue][6+numVariable*i],
                       &(*dataMatrix)[eigenValue][7+numVariable*i]);

                (*dataMatrix)[eigenValue][0+numVariable*i] = (double) i+1;
                (*dataMatrix)[eigenValue][1+numVariable*i] = 0.0;

                /*
                printf("Data = %f %f %f %f %f %f %f %f\n", (*dataMatrix)[eigenValue][0+numVariable*i],
                                                           (*dataMatrix)[eigenValue][1+numVariable*i],
                                                           (*dataMatrix)[eigenValue][2+numVariable*i],
                                                           (*dataMatrix)[eigenValue][3+numVariable*i],
                                                           (*dataMatrix)[eigenValue][4+numVariable*i],
                                                           (*dataMatrix)[eigenValue][5+numVariable*i],
                                                           (*dataMatrix)[eigenValue][6+numVariable*i],
                                                           (*dataMatrix)[eigenValue][7+numVariable*i]);
                 */
                i = i+1;
            }

            eigenValue += 1;

            // Skip ahead 6 lines after reading an eigenvector
            for (j = 0; j < 6; j++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }
        }

        if (eigenValue == *numEigenVector) break;
    }

    if (eigenValue != *numEigenVector) {
        printf("\tOnly %d of %d Eigen-Vectors read!", eigenValue, *numEigenVector);
        status = CAPS_NOTFOUND;
    } else {
        status = CAPS_SUCCESS;
    }

cleanup:
    if (line != NULL) EG_free(line);

    return status;
}

// Read data from a Astros OUT file and load it into a dataMatrix[numEigenVector][5]
// where variables are eigenValue, eigenValue(radians), eigenValue(cycles),
// generalized mass, and generalized stiffness.      MASS              STIFFNESS
int astros_readOUTEigenValue(FILE *fp, int *numEigenVector, double ***dataMatrix)
{
    int status; // Function return

    int i, j;// Indexing

    int tempInt, eigenValue =0;

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *beginEigenLine = "                              ORDER         (RAD/S)**2         (RAD/S)           (HZ)            MASS           STIFFNESS";
    //char *endEigenLine = "1";

    int numVariable = 5; // EigenValue, eigenValue(radians), eigenValue(cycles), generalized mass, and generalized stiffness.

    printf("Reading Astros OUT file - extracting Eigen-Values!\n");

    *numEigenVector = 0;

    // See how many Eigen-Values we have
    status = astros_readOUTNumEigenValue(fp, numEigenVector);
    printf("\tNumber of Eigen-Values = %d\n", *numEigenVector);
    if (status != CAPS_SUCCESS) return status;

    // Allocate dataMatrix array
    if (*dataMatrix != NULL) EG_free(*dataMatrix);

    *dataMatrix = (double **) EG_alloc(*numEigenVector *sizeof(double *));
    if (*dataMatrix == NULL) return EGADS_MALLOC; // If allocation failed ....

    for (i = 0; i < *numEigenVector; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc(numVariable*sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < i; j++) {

                if ((*dataMatrix)[j] != NULL ) EG_free((*dataMatrix)[j]);
            }

            if ((*dataMatrix) != NULL) EG_free((*dataMatrix));

            return EGADS_MALLOC;
        }
    }

    // Loop through the file again and pull out data
    while (eigenValue != *numEigenVector) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // Look for start of Eigen-Vector
        if (strncmp(beginEigenLine, line, strlen(beginEigenLine)) == 0) {

            // Fast forward 1 lines
            for (i = 0; i < 1; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            i = 0;
            while (eigenValue != *numEigenVector) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;

                // If we have a new page - skip ahead 8 lines and continue
                if (strncmp("1", line, 1) == 0) {

                    for (j = 0; j < 8; j++) {
                        status = getline(&line, &linecap, fp);
                        if (status < 0) break;
                    }
                    continue;
                }

                // Loop through the file and fill up the data matrix
                sscanf(line, "%d%d%lf%lf%lf%lf%lf", &eigenValue,
                        &tempInt,
                        &(*dataMatrix)[i][0],
                        &(*dataMatrix)[i][1],
                        &(*dataMatrix)[i][2],
                        &(*dataMatrix)[i][3],
                        &(*dataMatrix)[i][4]);
                printf("\tLoading Eigen-Value = %d\n", eigenValue);
                i += 1;
            }
        }
    }

    if (line != NULL) EG_free(line);

    return CAPS_SUCCESS;
}

// Read data from a Astros OUT file and load it into a dataMatrix[numGridPoint][8]
// where variables are Grid Id, Coord Id, T1, T2, T3, R1, R2, R3
int astros_readOUTDisplacement(FILE *fp, int subcaseId, int *numGridPoint,
                               double ***dataMatrix)
{
    int status; // Function return

    int i, j; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder

    char *outputSubcaseLine = "0                                                                                                            SUBCASE ";
    char *displacementLine ="                                             D I S P L A C E M E N T   V E C T O R";
    char *beginSubcaseLine=NULL;
    char *endSubcaseLine = "1";
    char tempString[2];

    int numVariable = 8; // Grid Id, coord Id T1, T2, T3, R1, R2, R3
    int intLength;
    int lineFastForward = 0;

    if (subcaseId != -1) {
        printf("Reading Astros FO6 file - extracting Displacements!\n");
    }

    *numGridPoint = 0;

    if (*dataMatrix != NULL) {
      printf("Developer error: dataMatrix should be NULL!\n");
      return CAPS_NULLVALUE;
    }


    if      (subcaseId >= 1000) intLength = 4;
    else if (subcaseId >= 100) intLength = 3;
    else if (subcaseId >= 10) intLength = 2;
    else intLength = 1;

    if (subcaseId == -1) {

        /* count number of grid points */
        rewind(fp);

/*@-nullpass@*/
        while (1) {
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            if (strstr(line, "D I S P L A C E M E N T   V E C T O R") != NULL) break;
        }

        (void) getline(&line, &linecap, fp);
        (void) getline(&line, &linecap, fp);

        while (1) {
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            if (strstr(line, "S T R E S S E S") != NULL) break;

            if (strstr(line, "   G   ") != NULL) {
                (*numGridPoint)++;
            }
        }
/*@+nullpass@*/

        /* allocate space for the grid points */
        *dataMatrix = (double **) EG_alloc((*numGridPoint)*sizeof(double *));
        if (*dataMatrix == NULL) { status = EGADS_MALLOC; goto cleanup; }

        for (i = 0; i < (*numGridPoint); i++) {
            (*dataMatrix)[i] = (double *) EG_alloc(8*sizeof(double));
            if ((*dataMatrix)[i] == NULL) { status = EGADS_MALLOC; goto cleanup; }
        }

        /* read the grid points */
        rewind(fp);

/*@-nullpass@*/
        while (1) {
            status = getline(&line, &linecap, fp);
            if (status < 0) break;

            if (strstr(line, "D I S P L A C E M E N T   V E C T O R") != NULL) break;
        }

        (void) getline(&line, &linecap, fp);
        (void) getline(&line, &linecap, fp);

        i = 0;
        while (i < *numGridPoint) {
            int igid;

            status = getline(&line, &linecap, fp);
            if (status == -1) {
              status = CAPS_IOERR;
              goto cleanup;
            }
            if (strstr(line, "   G   ") != NULL) {
                sscanf(line, "%d G %lf %lf %lf %lf %lf %lf",
                       &igid, &(*dataMatrix)[i][2], &(*dataMatrix)[i][3],
                              &(*dataMatrix)[i][4], &(*dataMatrix)[i][5],
                              &(*dataMatrix)[i][6], &(*dataMatrix)[i][7]);
                (*dataMatrix)[i][0] = igid;
                (*dataMatrix)[i][1] = 0;
                i++;
            }
        }
/*@+nullpass@*/

        status = CAPS_SUCCESS;
        goto cleanup;
    }

    if (subcaseId > 0) {
        beginSubcaseLine = (char *) EG_alloc((strlen(outputSubcaseLine)+
                                              intLength+1)*sizeof(char));
        if (beginSubcaseLine == NULL) { status = EGADS_MALLOC; goto cleanup; }

        sprintf(beginSubcaseLine,"%s%d",outputSubcaseLine, subcaseId);

        lineFastForward = 4;

    } else {

        intLength = 0;
        beginSubcaseLine = (char *) EG_alloc((strlen(displacementLine)+1)*
                                             sizeof(char));
        if (beginSubcaseLine == NULL) { status = EGADS_MALLOC; goto cleanup; }
        sprintf(beginSubcaseLine,"%s",displacementLine);

        lineFastForward = 2;
    }
    printf("beginSubcaseLine=%s, lineFastForward=%d\n",
           beginSubcaseLine, lineFastForward);

    beginSubcaseLine[strlen(outputSubcaseLine)+intLength] = '\0';

    // Loop through file line by line until we have determined how many grid points we have
    while (*numGridPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;
        printf("line=%s\n", line);

        // Look for start of subcaseId
        if (strncmp(beginSubcaseLine, line, strlen(beginSubcaseLine)) == 0) {

            // Fast forward lines
            for (i = 0; i < lineFastForward; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through lines counting the number of grid points
            while (getline(&line, &linecap, fp) >= 0) {
                if (strncmp(endSubcaseLine, line, strlen(endSubcaseLine)) == 0) break;
                *numGridPoint +=1;
            }
        }
    }

    printf("Number of Grid Points = %d\n", *numGridPoint);

    if (*numGridPoint == 0) {
        printf("Either data points  = 0 and/or subcase wasn't found\n");

        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Rewind the file
    rewind(fp);

    // Allocate dataMatrix array
    *dataMatrix = (double **) EG_alloc(*numGridPoint *sizeof(double *));
    if (*dataMatrix == NULL) { status = EGADS_MALLOC; goto cleanup; }

    for (i = 0; i < *numGridPoint; i++) {

        (*dataMatrix)[i] = (double *) EG_alloc(numVariable*sizeof(double));

        if ((*dataMatrix)[i] == NULL) { // If allocation failed ....
            for (j = 0; j < *numGridPoint; j++) {
                EG_free((*dataMatrix)[j]);
            }
            EG_free((*dataMatrix));

            status = EGADS_MALLOC;
            goto cleanup;
        }
    }

    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {
        if (line == NULL) break;

        // Look for start of Displacement
        if (strncmp(beginSubcaseLine, line, strlen(beginSubcaseLine)) == 0) {

            printf("Loading displacements for Subcase = %d\n", subcaseId);

            // Fast forward lines
            for (i = 0; i < lineFastForward; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through the file and fill up the data matrix
            for (i = 0; i < (*numGridPoint); i++) {
                for (j = 0; j < numVariable; j++) {

                    if (j == 0 || j % numVariable+1 == 0) {
                        fscanf(fp, "%lf", &(*dataMatrix)[i][j]);
                        fscanf(fp, "%s", tempString);
                        j = j + 1;
                        (*dataMatrix)[i][j] = 0.0;
                    } else {
                        fscanf(fp, "%lf", &(*dataMatrix)[i][j]);
                    }
                }
            }

            break;
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    EG_free(beginSubcaseLine);
    EG_free(line);

    return status;
}

// Write out a DVGRID entry
static int astros_writeDVGRIDCard(FILE *fp, int dvID, meshNodeStruct node,
                                  double scaleCoeff, double designVec[3],
                                  feaFileFormatStruct *feaFileFormat)
{
    int status = CAPS_SUCCESS;
    cardStruct card;

    if (fp == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    status = card_initiate(&card, "DVGRID", feaFileFormat->fileType);
    if (status != CAPS_SUCCESS) goto cleanup;

    card_addInteger(&card, dvID);

    card_addInteger(&card, node.nodeID);

    card_addBlank(&card); // CID blank field

    card_addDouble(&card, scaleCoeff);

    card_addDoubleArray(&card, 3, designVec);

    card_write(&card, fp);

    card_destroy(&card);

cleanup:

    card_destroy(&card);

    return status;
}

static int astros_writeSNORMCard(FILE *fp, meshNodeStruct node, double snorm[3],
                                 int patchID, /*@unused@*/ int cAxis,
                                 feaFileFormatStruct *feaFileFormat)
{
    int status = CAPS_SUCCESS;
    int coordID = 0;
    feaMeshDataStruct *feaData;

    cardStruct card;

    if (fp == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    if (node.analysisType == MeshStructure) {
        feaData = (feaMeshDataStruct *) node.analysisData;
        coordID = feaData->coordID;
    } else {
        printf("Incorrect analysis type for node %d", node.nodeID);
        return CAPS_BADVALUE;
    }

    status = card_initiate(&card, "SNORM", feaFileFormat->fileType);
    if (status != CAPS_SUCCESS) goto cleanup;

    card_addInteger(&card, node.nodeID);

    card_addInteger(&card, coordID);

    card_addDoubleArray(&card, 3, snorm);

    card_addInteger(&card, patchID);

    // need control from pyCAPS script over cAxis value---blank (omit last field) for default
    // card_addInteger(&card, cAxis);

    card_write(&card, fp);

    cleanup:

        card_destroy(&card);

        return status;
}

static int astros_writeSNORMDTCard(FILE *fp, int dvID, meshNodeStruct node,
                                   double snormdt[3], int patchID,
                                   feaFileFormatStruct *feaFileFormat)
{
    // int fieldWidth;
    // char *tempString=NULL;
    // char *delimiter;
    int status = CAPS_SUCCESS;
    int coordID = 0;
    feaMeshDataStruct *feaData;
    double mag, vec[3];

    cardStruct card;

    if (fp == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat == NULL) return CAPS_NULLVALUE;

    // if (feaFileFormat->fileType == FreeField) {
    //     delimiter = ",";
    //     fieldWidth = 8;
    // } else {
    //     delimiter = " ";
    //     fieldWidth = 7;

    // }

    if (node.analysisType == MeshStructure) {
        feaData = (feaMeshDataStruct *) node.analysisData;
        coordID = feaData->coordID;
    } else {
        printf("Incorrect analysis type for node %d", node.nodeID);
        return CAPS_BADVALUE;
    }

    // Get vector length
    mag = dot_DoubleVal(snormdt, snormdt);
    mag = sqrt(mag);

    vec[0] = snormdt[0];
    vec[1] = snormdt[1];
    vec[2] = snormdt[2];

    if (mag != 0.0) {
        vec[0] /= mag;
        vec[1] /= mag;
        vec[2] /= mag;
    }

    status = card_initiate(&card, "SNORMDT", feaFileFormat->fileType);
    if (status != CAPS_SUCCESS) goto cleanup;

    card_addInteger(&card, dvID);

    card_addInteger(&card, node.nodeID);

    card_addInteger(&card, coordID);

    card_addDouble(&card, mag);

    card_addDoubleArray(&card, 3, vec);

    card_addInteger(&card, patchID);

    card_write(&card, fp);

cleanup:

    card_destroy(&card);

    return status;
}

static int check_edgeInCoplanarFace(ego edge, double t, int currentFaceIndex,
                                    ego body, int *coplanarFlag)
{
    int status; // Function return status

    int faceIndex, loopIndex, edgeIndex; // Indexing

    int numBodyFace;
    ego *bodyFace = NULL;

    int numChildrenLoop, numChildrenEdge;
    ego *childrenLoop, *childrenEdge;

    ego geomRef;
    int oclass,  mtype, *senses;
    double data[18], params[2];

    double mag;
    double residual = 1E-6; // Tolerence residual for normal comparison
    double normal[3], normal2[3];

    *coplanarFlag = (int) false;

    status = EG_getBodyTopos(body, NULL, FACE, &numBodyFace, &bodyFace);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (bodyFace == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    // Get uv on face at edge - need to check sense setting currently 0;
    status = EG_getEdgeUV(bodyFace[currentFaceIndex], edge, 0, t, data);
    if (status != EGADS_SUCCESS) goto cleanup;

    params[0] = data[0];
    params[1] = data[2];

    // Get derivative along face
    status = EG_evaluate(bodyFace[currentFaceIndex], params, data);
    if (status != EGADS_SUCCESS) goto cleanup;

    // Get face normal
    (void) cross_DoubleVal(&data[3], &data[6], normal);

    mag = dot_DoubleVal(normal, normal);
    mag = sqrt(mag);

    normal[0] = fabs(normal[0])/mag;
    normal[1] = fabs(normal[1])/mag;
    normal[2] = fabs(normal[2])/mag;

    for (faceIndex = 0; faceIndex < numBodyFace; faceIndex++) {

        if (faceIndex == currentFaceIndex) continue;

        status = EG_getTopology(bodyFace[faceIndex],
                                &geomRef, &oclass, &mtype, data,
                                &numChildrenLoop, &childrenLoop, &senses);
        if (status != EGADS_SUCCESS) goto cleanup;

        for (loopIndex = 0; loopIndex < numChildrenLoop; loopIndex++) {

            status = EG_getTopology(childrenLoop[loopIndex],
                                    &geomRef, &oclass, &mtype, data,
                                    &numChildrenEdge, &childrenEdge, &senses);
            if (status != EGADS_SUCCESS) goto cleanup;

            for (edgeIndex = 0; edgeIndex < numChildrenEdge; edgeIndex++) {

                if (edge != childrenEdge[edgeIndex]) continue;

                // Get uv on face at edge - need to check sense setting currently 0;
                status = EG_getEdgeUV(bodyFace[faceIndex], edge, 0, t, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                params[0] = data[0];
                params[1] = data[2];

                // Get derivative along face
                status = EG_evaluate(bodyFace[faceIndex], params, data);
                if (status != EGADS_SUCCESS) goto cleanup;

                // Get new face normal
                (void) cross_DoubleVal(&data[3], &data[6], normal2);

                mag = dot_DoubleVal(normal2, normal2);
                mag = sqrt(mag);

                normal2[0] = fabs(normal2[0])/mag;
                normal2[1] = fabs(normal2[1])/mag;
                normal2[2] = fabs(normal2[2])/mag;

                if (fabs(normal[0] - normal2[0]) <= residual &&
                    fabs(normal[1] - normal2[1]) <= residual &&
                    fabs(normal[2] - normal2[2]) <= residual) {

                    *coplanarFlag = true;
                    loopIndex = numChildrenLoop;
                    faceIndex = numBodyFace;
                    break;
                }
            }
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in check_edgeInCoplanarFace, status %d\n",
               status);

    if (bodyFace != NULL) EG_free(bodyFace);

    return status;
}

#ifdef DEFINED_BUT_NOT_USED
static int check_nodeNormalHist(double normal[3], int *numNormal,
                                double *normalHist[], int *normalExist)
{
    int status; // Function return status

    int i; // Indexing

    double mag;
    double residual = 1E-6; // Tolerence residual for normal comparison

    double scaledNormal[3];

    *normalExist = (int) false;

    mag = dot_DoubleVal(normal, normal);
    mag = sqrt(mag);

    scaledNormal[0] = fabs(normal[0])/mag;
    scaledNormal[1] = fabs(normal[1])/mag;
    scaledNormal[2] = fabs(normal[2])/mag;

    // Does the normal already exist
    for (i = 0; i < *numNormal; i++) {

        // Have a miss-step between numNormal and normalHist
        if (normalHist == NULL) return CAPS_NULLVALUE;

        if (fabs(scaledNormal[0] - (*normalHist)[3*i + 0]) <= residual &&
                fabs(scaledNormal[1] - (*normalHist)[3*i + 1]) <= residual &&
                fabs(scaledNormal[2] - (*normalHist)[3*i + 2]) <= residual) {

            *normalExist = (int) true;
            break;
        }
    }

    if (*normalExist == (int) true) {
        status = CAPS_SUCCESS;
        goto cleanup;
    }

    *numNormal += 1;

    *normalHist = (double *) EG_reall(*normalHist, (*numNormal)*3*sizeof(double));
    if (*normalHist == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    //printf("Normal = %f %f %f\n", normal[0], normal[1], normal[2]);

    (*normalHist)[3*(*numNormal-1) + 0] = scaledNormal[0];
    (*normalHist)[3*(*numNormal-1) + 1] = scaledNormal[1];
    (*normalHist)[3*(*numNormal-1) + 2] = scaledNormal[2];

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in check_nodeNormalHist, status %d\n",
               status);

    return status;
}
#endif

// Get the configuration sensitivity at a given point and write out the SNORMDT card
static int astros_getConfigurationSens(FILE *fp,
                                       void *aimInfo,
                                       int numDesignVariable,
                                       feaDesignVariableStruct *feaDesignVariable,
                                       feaFileFormatStruct *feaFileFormat,
                                       int numGeomIn,
                                       capsValue *geomInVal,
                                       ego tess, int topoType, int topoIndex,
                                       int pointIndex,
                                       double *pointNorm, int patchID,
                                       meshNodeStruct node)
{
    int status = 0; // Function return

    int i, j; // Indexing

    const char *geomInName;

    int numPoint;
    double *dxyz = NULL;

    double mag, snormDT[3];

    if (fp                == NULL) return CAPS_NULLVALUE;
    if (aimInfo           == NULL) return CAPS_NULLVALUE;
    if (feaDesignVariable == NULL) return CAPS_NULLVALUE;
    if (geomInVal         == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat     == NULL) return CAPS_NULLVALUE;
    if (pointNorm         == NULL) return CAPS_NULLVALUE;

    // Loop through design variables
    for (i = 0; i < numDesignVariable; i++) {

        for (j = 0; j < numGeomIn; j++) {

            status = aim_getName(aimInfo, j+1, GEOMETRYIN, &geomInName);
            if (status != CAPS_SUCCESS) goto cleanup;

            if (strcmp(feaDesignVariable[i].name, geomInName) == 0) break;
        }

        // If name is found in Geometry inputs skip design variables
        if (j >= numGeomIn) continue;

        if(aim_getGeomInType(aimInfo, j+1) != 0) {
            printf("Error: Geometric sensitivity not available for CFGPMTR = %s\n",
                   geomInName);
            status = CAPS_NOSENSITVTY;
            goto cleanup;
        }


        //printf("Geometric sensitivity name = %s\n", feaDesignVariable[j].name);

        if (dxyz) EG_free(dxyz);
        dxyz = NULL;

        if (geomInVal[j].length > 1) {
            printf("Warning: Can NOT write SNORMDT cards for multidimensional design variables!\n");
            continue;
        }

        status = aim_setSensitivity(aimInfo, geomInName, 1, 1);
        if (status != CAPS_SUCCESS) goto cleanup;


        status = aim_getSensitivity(aimInfo, tess, topoType, topoIndex,
                                    &numPoint, &dxyz);
        if (status != CAPS_SUCCESS) goto cleanup;
        if (dxyz == NULL) {
            status = CAPS_NULLVALUE;
            goto cleanup;
        }

//        printf("TopoType = %d, TopoIndex = %d\n", topoType, topoIndex);
//        printf("numPoint %d, pointIndex %d\n", numPoint, pointIndex);

        if (pointIndex > numPoint) {
            status = CAPS_BADINDEX;
            goto cleanup;
        }

        // Normalize in incoming normal vector - just in case
        mag = dot_DoubleVal(pointNorm, pointNorm);

        pointNorm[0] /= sqrt(mag);
        pointNorm[1] /= sqrt(mag);
        pointNorm[2] /= sqrt(mag);

        // Get the scalar projection of the configuration sensitivity on the normal
        mag = dot_DoubleVal(&dxyz[3*(pointIndex-1) + 0], pointNorm);

        // Get the vector project of the configuration on the normal
        snormDT[0] = mag*pointNorm[0];
        snormDT[1] = mag*pointNorm[1];
        snormDT[2] = mag*pointNorm[2];

        printf(">>> Writing SNORMDT cards\n");
        status = astros_writeSNORMDTCard(fp, feaDesignVariable[i].designVariableID,
                                         node, snormDT, patchID, feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;

    } // Design variables

cleanup:

    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in astros_getConfigurationSens, status %d\n",
             status);

    if (dxyz != NULL) EG_free(dxyz);

    return status;
}

// Write boundary normals for shape sensitivities - only valid for modifications made by Bob Canfield to Astros
static int astros_getBoundaryNormal(FILE *fp,
                                    void *aimInfo,
                                    int numDesignVariable,
                                    feaDesignVariableStruct *feaDesignVariable,
                                    int numGeomIn,
                                    capsValue *geomInVal,
                                    meshStruct *feaMesh,
                                    feaFileFormatStruct *feaFileFormat)
{
    int status; // Function return status

    int i, j, k, faceIndex, loopIndex, edgeIndex; // Indexing

    int numMesh;
    meshStruct *mesh = NULL;
    int nodeOffSet = 0; //Keep track of global node indexing offset due to combining multiple meshes

    ego body;

    int numBodyNode, numBodyEdge, numBodyFace;
    ego *bodyNode = NULL, *bodyEdge = NULL, *bodyFace = NULL;

    int numChildrenLoop, numChildrenEdge, numChildrenNode;
    ego *childrenLoop, *childrenEdge, *childrenNode;

    int numPoint;
    int tessState, pointLocalIndex, nodeOrEdge, pointTopoIndex, edgeTopoIndex;

    int len;
    const double *xyz, *t;

    ego geomRef;
    int oclass, mtypeFace, mtype, *senses, *dummySenses;
    double data[18], params[2];

    double normEdge[3], normFace[3], normBoundary[3], mag;

    int coplanarFlag;

    //int numNormalHist; //, normalExistFlag;
    //double *normalHist = NULL;

    if (fp == NULL || feaMesh == NULL || feaFileFormat == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    // Are we dealing with a single mesh or a combined mesh
    if (feaMesh->egadsTess != NULL) {
        numMesh = 1;
        mesh = feaMesh;

    } else {

        numMesh = feaMesh->numReferenceMesh;
    }

    if (numMesh == 0) {
        printf("No bodies with tessellations found!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    for (i = 0; i < numMesh; i++ ) {

        if (mesh == NULL) mesh = &feaMesh->referenceMesh[i];

        status = EG_statusTessBody(mesh->egadsTess, &body,
                                   &tessState, &numPoint);
        if (status != EGADS_SUCCESS) goto cleanup;

        if (bodyNode != NULL) EG_free(bodyNode);
        if (bodyEdge != NULL) EG_free(bodyEdge);
        if (bodyFace != NULL) EG_free(bodyFace);

        status = EG_getBodyTopos(body, NULL, NODE, &numBodyNode, &bodyNode);
        if (status < EGADS_SUCCESS) goto cleanup;
        if (bodyNode == NULL) {
            status = CAPS_NULLOBJ;
            goto cleanup;
        }

        status = EG_getBodyTopos(body, NULL, EDGE, &numBodyEdge, &bodyEdge);
        if (status != EGADS_SUCCESS) goto cleanup;
        if (bodyEdge == NULL) {
            status = CAPS_NULLOBJ;
            goto cleanup;
        }

        status = EG_getBodyTopos(body, NULL, FACE, &numBodyFace, &bodyFace);
        if (status != EGADS_SUCCESS) goto cleanup;
        if (bodyFace == NULL) {
            status = CAPS_NULLOBJ;
            goto cleanup;
        }

        for (j = 0; j < numPoint; j++) {

            status = EG_getGlobal(mesh->egadsTess, j+1,
                                  &pointLocalIndex, &pointTopoIndex, NULL);
            if (status != EGADS_SUCCESS) goto cleanup;

            //printf("J = %d, PointLocalIndex %d, TopoIndex %d\n", j+1, pointLocalIndex, pointTopoIndex);

            if (pointLocalIndex < 0) continue; // Don't care about face nodes; only want edge and node nodes

            nodeOrEdge = pointLocalIndex;
            edgeTopoIndex = pointTopoIndex;

            // Has this point been included in the mesh?
            for (k = 0; k < mesh->numNode; k++) {
                if (mesh->node[k].nodeID != j+1) continue;
                break;
            }

            if (k >= mesh->numNode) continue; // Point isn't in the mesh - it has been removed

            // Set our node normal history settings
            /*
            numNormalHist = 0;
            if (normalHist != NULL) EG_free(normalHist);
            normalHist = NULL;
            */

            // Loop through the faces and find what edges
            for (faceIndex = 0; faceIndex < numBodyFace; faceIndex++) {

                status = EG_getTopology(bodyFace[faceIndex],
                                        &geomRef, &oclass, &mtypeFace, data,
                                        &numChildrenLoop, &childrenLoop, &senses);
                if (status != EGADS_SUCCESS) goto cleanup;

                for (loopIndex = 0; loopIndex < numChildrenLoop; loopIndex++) {

                    status = EG_getTopology(childrenLoop[loopIndex],
                                            &geomRef, &oclass, &mtype, data,
                                            &numChildrenEdge, &childrenEdge, &senses);
                    if (status != EGADS_SUCCESS) goto cleanup;

                    for (edgeIndex = 0; edgeIndex < numChildrenEdge; edgeIndex++) {

                        if (nodeOrEdge > 0)  { // Edge point

                            if (bodyEdge[edgeTopoIndex -1] !=
                                childrenEdge[edgeIndex]) continue;

                        } else { // Node

                            status = EG_getTopology(childrenEdge[edgeIndex],
                                                    &geomRef, &oclass, &mtype,
                                                    data, &numChildrenNode,
                                                    &childrenNode, &dummySenses);
                            if (status != EGADS_SUCCESS) goto cleanup;

                            pointLocalIndex = -1;

                            status = EG_indexBodyTopo(body, childrenEdge[edgeIndex]);
                            if (status < EGADS_SUCCESS) goto cleanup;

                            edgeTopoIndex = status;

                            if (numChildrenNode == 1 || numChildrenNode == 2) {

                                status = EG_getTessEdge(mesh->egadsTess,
                                                        edgeTopoIndex, &len, &xyz, &t);
                                if (status != EGADS_SUCCESS) goto cleanup;

                                if (bodyNode[pointTopoIndex -1] == childrenNode[0])
                                    pointLocalIndex = 1;

                                if (numChildrenNode > 1) {
                                    if (bodyNode[pointTopoIndex -1] == childrenNode[1])
                                        pointLocalIndex = len;
                                }

                            } else {
                                printf("Warning: Number of nodes = %d  for edge index %d\n",
                                       numChildrenNode, edgeIndex);
                                continue;
                            }

                            if (pointLocalIndex < 0) continue;
                        }

//                        printf("J = %d, egdeTopoIndex %d, pointLocalIndex %d\n", j+1, edgeTopoIndex, pointLocalIndex);

                        if (pointLocalIndex < 0) {
                            printf("Unable to determine pointLocalIndex\n");
                            status = CAPS_NOTFOUND;
                            goto cleanup;
                        }

                        // Get t - along edge
                        status = EG_getTessEdge(mesh->egadsTess,
                                                edgeTopoIndex, &len, &xyz, &t);
                        if (status != EGADS_SUCCESS) goto cleanup;

                        params[0] = t[pointLocalIndex-1];
                        //status = EG_invEvaluate(bodyEdge[pointLocalIndex-1], mesh->node[k].xyz,
                        //      params, data);
                        //if (status != EGADS_SUCCESS) goto cleanup;

                        // Check to see if edge is part of a co-planar face
                        status = check_edgeInCoplanarFace(bodyEdge[edgeTopoIndex-1],
                                                          params[0], faceIndex, body,
                                                          &coplanarFlag);
                        if (status != CAPS_SUCCESS) goto cleanup;

                        if (coplanarFlag == (int) true) continue;

                        // Get derivative along edge
                        status = EG_evaluate(bodyEdge[edgeTopoIndex-1], params,
                                             data);
                        if (status != EGADS_SUCCESS) goto cleanup;

                        normEdge[0] = data[3];
                        normEdge[1] = data[4];
                        normEdge[2] = data[5];

                        if (senses[edgeIndex]  < 0) {
                            normEdge[0] *= -1;
                            normEdge[1] *= -1;
                            normEdge[2] *= -1;
                        }

                        /*
                                                printf("Edge normal = %f %f %f, MTypeFace %d, %d %d\n", normEdge[0], normEdge[1], normEdge[2],
                                                                                                                                                                mtypeFace, SFORWARD, SREVERSE);
                         */

                        // Get uv on face at edge - need to check sense setting currently 0; in general the co-planar check
                        //  should catch this I think so it shouldn't be an issue
                        status = EG_getEdgeUV(bodyFace[faceIndex],
                                              bodyEdge[edgeTopoIndex-1], 0,
                                              params[0], data);
                        if (status != EGADS_SUCCESS) goto cleanup;

                        params[0] = data[0];
                        params[1] = data[2];

                        // Get derivative along face
                        status = EG_evaluate(bodyFace[faceIndex], params, data);
                        if (status != EGADS_SUCCESS) goto cleanup;

                        // Get face normal
                        (void) cross_DoubleVal(&data[3], &data[6], normFace);

                        // Test face orientation
                        if (mtypeFace == SREVERSE) {
                            normFace[0] *= -1;
                            normFace[1] *= -1;
                            normFace[2] *= -1;
                        }

                        // Get normal boundary
                        (void) cross_DoubleVal(normEdge, normFace, normBoundary);

                        mag = dot_DoubleVal(normBoundary, normBoundary);
                        mag = sqrt(mag);

                        normBoundary[0] = normBoundary[0]/mag;
                        normBoundary[1] = normBoundary[1]/mag;
                        normBoundary[2] = normBoundary[2]/mag;

                        /* Want to repeat normals - for now
                        if (nodeOrEdge == 0) {
                            status = check_nodeNormalHist(normBoundary, &numNormalHist, &normalHist, &normalExistFlag);
                            if (status != CAPS_SUCCESS) goto cleanup;

                            if (normalExistFlag == (int) true) continue;
                        }
                         */

                        printf(">>> Writing SNORM card\n");
                        status = astros_writeSNORMCard(fp, feaMesh->node[k+nodeOffSet],
                                                       normBoundary, edgeIndex+1, 1,
                                                       feaFileFormat);
                        if (status != CAPS_SUCCESS) goto cleanup;

                        printf(">>> Getting SNORMDT information\n");
                        status = astros_getConfigurationSens(fp,
                                                             aimInfo,
                                                             numDesignVariable,
                                                             feaDesignVariable,
                                                             feaFileFormat,
                                                             numGeomIn,
                                                             geomInVal,
                                                             mesh->egadsTess,
                                                             -1, edgeTopoIndex,
                                                             pointLocalIndex,
                                                             normBoundary,
                                                             faceIndex+1,
                                                             feaMesh->node[k+nodeOffSet]);
                        printf(">>> Done with SNORMDT information\n");
                        if (status != CAPS_SUCCESS) goto cleanup;


                    } // Children edge loop

                } // Children loop loop

            } // Face loop

        } // Point loop

        nodeOffSet += mesh->numNode;

        mesh = NULL;

    } // Mesh loop

    status = CAPS_SUCCESS;
    goto cleanup;

cleanup:
    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in astros_getBoundaryNormal, status %d\n",
             status);

    if (bodyNode != NULL) EG_free(bodyNode);
    if (bodyEdge != NULL) EG_free(bodyEdge);
    if (bodyFace != NULL) EG_free(bodyFace);

    //if (normalHist != NULL) EG_free(normalHist);

    return status;
}

// Write geometric parametrization - only valid for modifications made by Bob Canfield to Astros
int astros_writeGeomParametrization(FILE *fp,
                                    void *aimInfo,
                                    int numDesignVariable,
                                    feaDesignVariableStruct *feaDesignVariable,
                                    int numGeomIn,
                                    capsValue *geomInVal,
                                    meshStruct *feaMesh,
                                    feaFileFormatStruct *feaFileFormat)
{
    int status; // Function return status

    int i, j, k, m; // row, col; // Indexing

    int numMesh;
    meshStruct *mesh = NULL;

    int nodeOffSet = 0; //Keep track of global node indexing offset due to combining multiple meshes

    const char *geomInName;
    int numPoint;
    double *xyz = NULL;

    if (fp                == NULL) return CAPS_NULLVALUE;
    if (aimInfo           == NULL) return CAPS_NULLVALUE;
    if (feaDesignVariable == NULL) return CAPS_NULLVALUE;
    if (geomInVal         == NULL) return CAPS_NULLVALUE;
    if (feaMesh           == NULL) return CAPS_NULLVALUE;
    if (feaFileFormat     == NULL) return CAPS_NULLVALUE;

    // Are we dealing with a single mesh or a combined mesh
    if (feaMesh->egadsTess != NULL) {
        numMesh = 1;
        mesh = feaMesh;

    } else {
        //printf("Using reference meshes\n");
        numMesh = feaMesh->numReferenceMesh;
    }

    if (numMesh == 0) {
        printf("No bodies with tessellations found!\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Write sensitivity files for each  mesh
    for (i = 0; i < numMesh; i++ ) {
        printf(">>> Parametrization on mesh %d\n", i);
        if (mesh == NULL) mesh = &feaMesh->referenceMesh[i];

        for (j = 0; j < numDesignVariable; j++) {

            for (k = 0; k < numGeomIn; k++) {

                status = aim_getName(aimInfo, k+1, GEOMETRYIN, &geomInName);
                if (status != CAPS_SUCCESS) goto cleanup;

                if (strcmp(feaDesignVariable[j].name, geomInName) == 0) break;
            }

            // If name isn't found in Geometry inputs skip design variables
            if (k >= numGeomIn) continue;

            if(aim_getGeomInType(aimInfo, k+1) != 0) {
                printf("Error: Geometric sensitivity not available for CFGPMTR = %s\n",
                       geomInName);
                status = CAPS_NOSENSITVTY;
                goto cleanup;
            }

            printf("Geometric sensitivity name = %s\n", geomInName);
            //printf("Geometric sensitivity name = %s\n", feaDesignVariable[j].name);

            if (xyz != NULL) EG_free(xyz);
            xyz = NULL;

            if (geomInVal[k].length == 1) {
                printf(">>> Getting sensitivity\n");
                status = aim_tessSensitivity(aimInfo,
                                             geomInName,
                                             1, 1,
                                             mesh->egadsTess,
                                             &numPoint, &xyz);
                printf(">>> Back from getting sensitivity\n");
                AIM_STATUS(aimInfo, status, "Sensitivity for: %s\n", geomInName);

                if (numPoint != mesh->numNode) {
                    printf("Error: the number of nodes returned by aim_senitivity does NOT match the surface mesh!\n");
                    status = CAPS_MISMATCH;
                    goto cleanup;
                }

                if (xyz != NULL)
                    for (m = 0; m < mesh->numNode; m++) {

                        if (mesh->node[m].nodeID != m+1) {
                            printf("Error: Node Id %d is out of order (%d). No current fix!\n",
                                   mesh->node[m].nodeID, m+1);
                            status = CAPS_MISMATCH;
                            goto cleanup;
                        }
                        printf(">>> Write DVGRID cards\n");
                        status = astros_writeDVGRIDCard(fp,
                                                        feaDesignVariable[j].designVariableID,
                                                        feaMesh->node[m+nodeOffSet],
                                                        1.0,
                                                        &xyz[3*m + 0], feaFileFormat);
                        if (status != CAPS_SUCCESS) goto cleanup;
                    }

            } else {

                printf("Warning: Can NOT write DVGRID cards for multidimensional design variables!\n");
                continue;

                /*
                                for (row = 0; row < geomInVal[k].nrow; row++) {
                                        for (col = 0; col < geomInVal[k].ncol; col++) {

                                                if (xyz != NULL) EG_free(xyz);
                                                xyz = NULL;

                                                status = aim_tessSensitivity(aimInfo,
                                                                                                 geomInName,
                                                                                                 row+1, col+1, // row, col
                                                                                                 mesh->egadsTess,
                                                                                                 &numPoint, &xyz);
                                                AIM_STATUS(aimInfo, status, "Sensitivity for: %s\n", geomInName);


                                                if (numPoint != mesh->numNode) {
                                                        printf("Error: the number of nodes returned by aim_senitivity does NOT match the surface mesh!\n");
                                                        status = CAPS_MISMATCH;
                                                        goto cleanup;
                                                }

                                                for (m = 0; m < mesh->numNode; m++) {

                                                        if (mesh->node[m].nodeID != m+1) {
                                                                printf("Error: Node Id %d is out of order (%d). No current fix!\n", mesh->node[m].nodeID, m+1);
                                                                status = CAPS_MISMATCH;
                                                                goto cleanup;
                                                        }

                                                        status = astros_writeDVGRIDCard(fp,
                                                                                                                        feaDesignVariable[j].designVariableID,
                                                                                                                        feaMesh->node[m+nodeOffSet],
                                                                                                                        1.0,
                                                                                                                        &xyz[3*m + 0], feaFileFormat);
                                                        if (status != CAPS_SUCCESS) goto cleanup;
                                                }
                                        } // Col
                                } // Row
                 */
            } // Size of geometry
        } // Design variables

        nodeOffSet += mesh->numNode;

        mesh = NULL;

    } // Mesh
    printf(">>> Done with DVGRID cards\n");

    // Write out SNorm card
    printf(">>> Getting SNORM data\n");
    status = astros_getBoundaryNormal(fp,
                                      aimInfo,
                                      numDesignVariable,
                                      feaDesignVariable,
                                      numGeomIn,
                                      geomInVal,
                                      feaMesh, feaFileFormat);
    printf(">>> Done with SNORM data\n");
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("Error: Premature exit in astros_writeGeomParametrization, status %d\n", status);

        if (xyz != NULL) EG_free(xyz);

        return status;

}

// Write a mesh contained in the mesh structure in Astros format (*.bdf)
int astros_writeMesh(void *aimInfo,
                     char *fname,
                     int asciiFlag, // 0 for binary, anything else for ascii
                     meshStruct *mesh,
                     feaFileTypeEnum gridFileType,
                     int numDesignVariable,
                     feaDesignVariableStruct feaDesignVariable[],
                     double scaleFactor) // Scale factor for coordinates
{
    int status; // Function return status
    FILE *fp = NULL;
    int i, j;
    int coordID, propertyID;
    double xyz[3];
    char *filename = NULL;
    char fileExt[] = ".bdf";
    feaFileTypeEnum formatType;
    int *cp = NULL;
    void *tm = NULL;
    double *tmax = NULL;

    feaMeshDataStruct *feaData;

    // Design variables
    int foundDesignVar, designIndex;
    double maxDesignVar = 0.0;

    if (mesh == NULL) return CAPS_NULLVALUE;

    printf("\nWriting Astros grid and connectivity file (in large field format) ....\n");

    if (asciiFlag == 0) {
        printf("\tBinary output is not currently supported for working with Astros\n");
        printf("\t..... switching to ASCII!\n");
        //asciiFlag = 1;
    }

    if (scaleFactor <= 0) {
        printf("\tScale factor for mesh must be > 0! Defaulting to 1!\n");
        scaleFactor = 1;
    }

    filename = (char *) EG_alloc((strlen(fname) + 1 + strlen(fileExt)) *sizeof(char));
    if (filename == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    sprintf(filename,"%s%s", fname, fileExt);

    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        printf("\tUnable to open file: %s\n", filename);

        status = CAPS_IOERR;
        goto cleanup;
    }

    if (gridFileType == LargeField) {
        fprintf(fp,"$---1A--|-------2-------|-------3-------|-------4-------|-------5-------|-10A--|\n");
        fprintf(fp,"$---1B--|-------6-------|-------7-------|-------8-------|-------9-------|-10B--|\n");
    }
    else if (gridFileType == FreeField) {
        fprintf(fp,"$-------------------------------------------------------------------------------\n");
    }
    else {
        fprintf(fp,"$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|\n");
    }

    // Write nodal coordinates
    for (i = 0; i < mesh->numNode; i++) {

        coordID = 0; // blank by default

        if (mesh->node[i].analysisType == MeshStructure) {
            feaData = (feaMeshDataStruct *) mesh->node[i].analysisData;
            coordID = feaData->coordID;
        }

        if (coordID != 0) {
            cp = &coordID;
        }
        else {
            cp = NULL;
        }

        xyz[0] = mesh->node[i].xyz[0]*scaleFactor;
        xyz[1] = mesh->node[i].xyz[1]*scaleFactor;
        xyz[2] = mesh->node[i].xyz[2]*scaleFactor;

        // printf("Writing GRID %d\n", mesh->node[i].nodeID);
        status = astrosCard_grid(fp,
                                 &mesh->node[i].nodeID, // id
                                 cp, // cp
                                 xyz, // xi
                                 NULL, // cd, blank
                                 NULL, // ps, blank
                                 gridFileType);
        if (status != CAPS_SUCCESS) goto cleanup;

    }

    if (gridFileType == LargeField) {
        formatType = SmallField;
    }
    else {
        formatType = gridFileType;
    }

    if (formatType == LargeField) {
        fprintf(fp,"$---1A--|-------2-------|-------3-------|-------4-------|-------5-------|-10A--|\n");
        fprintf(fp,"$---1B--|-------6-------|-------7-------|-------8-------|-------9-------|-10B--|\n");
    }
    else if (formatType == FreeField) {
        fprintf(fp,"$-------------------------------------------------------------------------------\n");
    }
    else {
        fprintf(fp,"$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|\n");
    }

    for (i = 0; i < mesh->numElement; i++) {

        // If we have a volume mesh skip the surface elements
        if (mesh->meshType == VolumeMesh) {
            if (mesh->element[i].elementType != Tetrahedral &&
                mesh->element[i].elementType != Pyramid     &&
                mesh->element[i].elementType != Prism       &&
                mesh->element[i].elementType != Hexahedral) continue;
        }

        // Grab Structure specific related data if available
        if (mesh->element[i].analysisType == MeshStructure) {
            feaData = (feaMeshDataStruct *) mesh->element[i].analysisData;
            propertyID = feaData->propertyID;
            coordID = feaData->coordID;
        } else {
            propertyID = mesh->element[i].markerID;
            coordID = 0;
        }

        // If got structure specific related data, set tm = coordID
        if (coordID != 0) {
            tm = &coordID;
        }
        else {
            tm = NULL;
        }

        // Check for design minimum area
        foundDesignVar = (int) false;
        for (designIndex = 0; designIndex < numDesignVariable; designIndex++) {
            for (j = 0; j < feaDesignVariable[designIndex].numPropertyID; j++) {

                if (feaDesignVariable[designIndex].propertySetID[j] == propertyID) {
                    foundDesignVar = (int) true;

                    maxDesignVar = feaDesignVariable[designIndex].upperBound;

                    // If 0.0 don't do anything
                    if (maxDesignVar == 0.0) foundDesignVar = (int) false;

                    break;
                }
            }

            if (foundDesignVar == (int) true) break;
        }

        // If found design variable, set tmax = upper bound
        if (foundDesignVar) {
            tmax = &maxDesignVar;
        }
        else {
            tmax = NULL;
        }

        if ( mesh->element[i].elementType == Line) { // Need to add subType for bar and beam .....

            // printf("Writing CROD %d\n", mesh->element[i].elementID);
            status = astrosCard_crod(fp,
                                     &mesh->element[i].elementID,
                                     &propertyID,
                                     mesh->element[i].connectivity,
                                     tmax,
                                     formatType);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        if ( mesh->element[i].elementType == Triangle) {

            // printf("Writing CTRIA3 %d\n", mesh->element[i].elementID);
            status = astrosCard_ctria3(fp,
                                       &mesh->element[i].elementID, // eid
                                       &propertyID, // pid
                                       mesh->element[i].connectivity, // Gi
                                       Integer, tm, // tm,
                                       NULL, // zoff
                                       tmax, // tmax
                                       NULL, // Ti
                                       formatType);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        if ( mesh->element[i].elementType == Quadrilateral) {

            // printf("Writing CQUAD4 %d\n", mesh->element[i].elementID);
            status = astrosCard_cquad4(fp,
                                       &mesh->element[i].elementID, // eid
                                       &propertyID, // pid
                                       mesh->element[i].connectivity, // Gi
                                       Integer, tm, // tm,
                                       NULL, // zoff
                                       tmax, // tmax
                                       NULL, // Ti
                                       formatType);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        if ( mesh->element[i].elementType == Tetrahedral)   {

            printf("\tWarning: Astros doesn't support tetrahedral elements - skipping element %d\n",
                   mesh->element[i].elementID);
        }

        if ( mesh->element[i].elementType == Pyramid) {

            printf("\tWarning: Astros doesn't support pyramid elements - skipping element %d\n",
                   mesh->element[i].elementID);
        }

        if ( mesh->element[i].elementType == Prism)   {

            printf("\tWarning: Astros doesn't support prism elements - skipping element %d\n",
                   mesh->element[i].elementID);
        }

        if ( mesh->element[i].elementType == Hexahedral) {

            status = astrosCard_cihex1(fp,
                                       &mesh->element[i].elementID, // eid
                                       &propertyID, // pid
                                       mesh->element[i].connectivity, // Gi
                                       formatType);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    // fprintf(fp,"$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|\n");

    printf("Finished writing Astros grid file\n\n");

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit in astros_writeMesh, status = %d\n", status);

    if (filename != NULL) EG_free(filename);

    if (fp != NULL) fclose(fp);

    return status;
}
