/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "abaqusUtils.h"

#include <string.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"    // Meshing utilities
#include "arrayUtils.h"   // Array utilities

#ifdef WIN32
#define strcasecmp  stricmp
#endif


#ifdef HAVE_PYTHON

#include "Python.h" // Bring in Python API

#include "abaqusFILReader.h" // Bring in Cython generated header file

#ifndef CYTHON_PEP489_MULTI_PHASE_INIT
#define CYTHON_PEP489_MULTI_PHASE_INIT (PY_VERSION_HEX >= 0x03050000)
#endif

#if CYTHON_PEP489_MULTI_PHASE_INIT
static int abaqusFILReader_Initialized = (int)false;
#endif
#endif


// Write a Abaqus Analysis card from a feaAnalysis structure
int abaqus_writeAnalysisCard(void *aimInfo,
                             FILE *fp,
                             const int numLoad,
                             const feaLoadStruct feaLoad[],
                             const feaAnalysisStruct *feaAnalysis,
                             const meshStruct *mesh) {

    int status;

    int i, j;

    // TODO: add ability to set nonlinear geometry and numIncrements ==100
    fprintf(fp,"*STEP, NAME=%s, NLGEOM=NO, INC=100\n", feaAnalysis->name);

    // Eigenvalue
    if (feaAnalysis->analysisType == Modal) {
        fprintf(fp, "*FREQUENCY, ");
        if (feaAnalysis->extractionMethod != NULL)  fprintf(fp, "EIGENSOLVER=%s,", feaAnalysis->extractionMethod);

        if (feaAnalysis->eigenNormalization != NULL) fprintf(fp, "NORMALIZATION=%s", feaAnalysis->eigenNormalization);

        fprintf(fp,"\n");

        if (feaAnalysis->numDesiredEigenvalue != 0) fprintf(fp, "%d,",feaAnalysis->numDesiredEigenvalue);
        else fprintf(fp, ",");

        fprintf(fp, "%f,", feaAnalysis->frequencyRange[0]);

        if (feaAnalysis->frequencyRange[1] != 0) fprintf(fp, "%f", feaAnalysis->frequencyRange[1]);

        fprintf(fp, "\n");
    }

    // Static
    if (feaAnalysis->analysisType == Static) {
        fprintf(fp, "*STATIC\n");

        fprintf(fp, "%f, ", 1.0); // Initial Time increment
        fprintf(fp, "%f, ", 1.0); // Time period
        fprintf(fp, "%f, ", 1E-5); // Minimum Time Increment
        fprintf(fp, "%f", 1.0); // Maximum Time Increment


        fprintf(fp, "\n");
    }

    // Write out loads for the case
    for (i = 0; i < feaAnalysis->numLoad; i++) {
        AIM_NOTNULL(feaLoad, aimInfo, status);
//        int numLoad;     // Number of loads in the analysis
//            int *loadSetID; // List of the load IDSs

        for (j = 0; j < numLoad; j++) {
            if (feaAnalysis->loadSetID[i] == feaLoad[j].loadID) break;
        }

        if (j >= numLoad) {
            printf("Unable to find load ID %d for analysis %s\n", feaAnalysis->loadSetID[i], feaAnalysis->name);
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        status = abaqus_writeLoadCard(fp, &feaLoad[j], mesh);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    //fprintf(fp, "*RESTART, WRITE, FREQUENCY=1\n");
    //fprintf(fp, "*CONTACT FILE\n");
    //fprintf(fp, "*ENERGY FILE\n");
    //fprintf(fp, "*MODAL FILE\n");
    //fprintf(fp, "*EL FILE\n");
    fprintf(fp, "*NODE PRINT\n");
    fprintf(fp, "U,\n");
    //fprintf(fp, "RF,\n");
    fprintf(fp, "*EL PRINT\n");
    fprintf(fp, "MISES,\n");
    //fprintf(fp, "*SECTION FILE\n");
    fprintf(fp, "*END STEP\n");

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("Error: Status %d during abaqus_writeAnalysisCard\n", status);

        return status;
}


// Write Abaqus constraint card from a feaConstraint structure
int abaqus_writeConstraintCard(FILE *fp, const feaConstraintStruct *feaConstraint) {

    int i; // Indexing

    if (feaConstraint->constraintType == Displacement) {

        fprintf(fp,"*BOUNDARY, TYPE=DISPLACEMENT\n");
        for (i = 0; i < feaConstraint->numGridID; i++ ) {

            if (feaConstraint->dofConstraint <= 10)     fprintf(fp,"%d, %d, %f\n",feaConstraint->gridIDSet[i],
                                                                                  feaConstraint->dofConstraint,
                                                                                  feaConstraint->gridDisplacement);

            if (feaConstraint->dofConstraint == 123456) fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 1,6,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 12345)  fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 1,5,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 1234)   fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 1,4,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 123)    fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 1,3,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 12)     fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 1,2,feaConstraint->gridDisplacement);

            if (feaConstraint->dofConstraint == 23456) fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 2,6,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 2345)  fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 2,5,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 234)   fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 2,4,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 23)    fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 2,3,feaConstraint->gridDisplacement);

            if (feaConstraint->dofConstraint == 3456) fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 3,6,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 345)  fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 3,5,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 34)   fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 3,4,feaConstraint->gridDisplacement);

            if (feaConstraint->dofConstraint == 456) fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 4,6,feaConstraint->gridDisplacement);
            if (feaConstraint->dofConstraint == 45)  fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 4,5,feaConstraint->gridDisplacement);

            if (feaConstraint->dofConstraint == 56) fprintf(fp,"%d, %d, %d, %f\n",feaConstraint->gridIDSet[i], 5,6,feaConstraint->gridDisplacement);
        }
    }
	if (feaConstraint->constraintType == ZeroDisplacement) {

	    fprintf(fp,"*BOUNDARY\n");

	    for (i = 0; i < feaConstraint->numGridID; i++ ) {

            if (feaConstraint->dofConstraint <= 10)     fprintf(fp,"%d, %d\n",feaConstraint->gridIDSet[i],
                                                                              feaConstraint->dofConstraint);

            if (feaConstraint->dofConstraint == 123456) fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 1,6);
            if (feaConstraint->dofConstraint == 12345)  fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 1,5);
            if (feaConstraint->dofConstraint == 1234)   fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 1,4);
            if (feaConstraint->dofConstraint == 123)    fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 1,3);
            if (feaConstraint->dofConstraint == 12)     fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 1,2);

            if (feaConstraint->dofConstraint == 23456) fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 2,6);
            if (feaConstraint->dofConstraint == 2345)  fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 2,5);
            if (feaConstraint->dofConstraint == 234)   fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 2,4);
            if (feaConstraint->dofConstraint == 23)    fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 2,3);

            if (feaConstraint->dofConstraint == 3456) fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 3,6);
            if (feaConstraint->dofConstraint == 345)  fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 3,5);
            if (feaConstraint->dofConstraint == 34)   fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 3,4);

            if (feaConstraint->dofConstraint == 456) fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 4,6);
            if (feaConstraint->dofConstraint == 45)  fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 4,5);

            if (feaConstraint->dofConstraint == 56) fprintf(fp,"%d, %d, %d\n",feaConstraint->gridIDSet[i], 5,6);
	    }
	}

	return CAPS_SUCCESS;
}

// Write a Abaqus Material card from a feaMaterial structure
int abaqus_writeMaterialCard(FILE *fp, const feaMaterialStruct *feaMaterial) {

    fprintf(fp,"*MATERIAL, NAME=%s\n", feaMaterial->name);

    if (feaMaterial->density != 0) fprintf(fp,"*DENSITY\n%f\n", feaMaterial->density);

	// Isotropic
	if (feaMaterial->materialType == Isotropic) {

	    if (feaMaterial->youngModulus != 0 &&
	        feaMaterial->poissonRatio != 0) fprintf(fp,"*ELASTIC, TYPE=Isotropic\n%f,%f\n",feaMaterial->youngModulus,
	                                                                                       feaMaterial->poissonRatio);
	    else {

	        if (feaMaterial->shearModulus != 0) fprintf(fp,"*ELASTIC, TYPE=Shear\n%f\n",feaMaterial->shearModulus);
	    }
	}

	// Orthotropic
	if (feaMaterial->materialType == Orthotropic) {

	    fprintf(fp,"*ELASTIC, TYPE=ENGINEERING CONSTANTS\n");

	    fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f,\n", feaMaterial->youngModulus,
	                                                        feaMaterial->youngModulusLateral,
	                                                        feaMaterial->youngModulusLateral,
	                                                        feaMaterial->poissonRatio,
	                                                        feaMaterial->poissonRatio,
	                                                        feaMaterial->poissonRatio,
	                                                        feaMaterial->shearModulus,
	                                                        feaMaterial->shearModulusTrans1Z);
	    fprintf(fp, "%f\n", feaMaterial->shearModulusTrans2Z);
	}

	return CAPS_SUCCESS;
}

// Write a Abaqus Property card from a feaProperty structure
int abaqus_writePropertyCard(FILE *fp, const feaPropertyStruct *feaProperty) {

	//          1D Elements        //

	// Rod
	if (feaProperty->propertyType ==  Rod) {

	    /*
        fprintf(fp,"%-8s", "PROD");

        tempString = convert_integerToString(feaProperty->propertyID, 7, 1);

        tempString = convert_integerToString(feaProperty->materialID, 7, 1);

        tempString = convert_doubleToString( feaProperty->crossSecArea, 7, 1);

        tempString = convert_doubleToString( feaProperty->torsionalConst, 7, 1);

        tempString = convert_doubleToString( feaProperty->torsionalStressReCoeff, 7, 1);

        tempString = convert_doubleToString( feaProperty->massPerLength, 7, 1);
        */

	}

	// Bar
	if (feaProperty->propertyType ==  Bar) {
	    /*
	    fprintf(fp,"%-8s", "PBAR");

		tempString = convert_integerToString(feaProperty->propertyID, 7, 1);

        tempString = convert_integerToString(feaProperty->materialID, 7, 1);

        tempString = convert_doubleToString( feaProperty->crossSecArea, 7, 1);

        tempString = convert_doubleToString( feaProperty->zAxisInertia, 7, 1);

        tempString = convert_doubleToString( feaProperty->yAxisInertia, 7, 1);

        tempString = convert_doubleToString( feaProperty->torsionalConst, 7, 1);
	    */

	}

	//          2D Elements        //

	// Shell
	if (feaProperty->propertyType == Shell) {

	    fprintf(fp,"*SHELL SECTION, ELSET=%s, MATERIAL=%s\n",feaProperty->name, feaProperty->materialName);

	    fprintf(fp,"%f\n", feaProperty->membraneThickness);

	    /*
        if (feaProperty->materialBendingID != 0) {
            tempString = convert_integerToString(feaProperty->materialBendingID, 7, 1);
        }

        tempString = convert_doubleToString(feaProperty->bendingInertiaRatio, 7, 1);

        if (feaProperty->materialShearID != 0) {
            tempString = convert_integerToString(feaProperty->materialShearID, 7, 1);
        } else {// Print a blank
        }

        tempString = convert_doubleToString(feaProperty->shearMembraneRatio, 7, 1);

        tempString = convert_doubleToString(feaProperty->massPerArea, 7, 1);
        */
	}


    //          3D Elements        //

    // Solid
    if (feaProperty->propertyType == Solid) {

        /*
        fprintf(fp,"%-8s", "PSOLID");

        tempString = convert_integerToString(feaProperty->propertyID, 7, 1);

        tempString = convert_integerToString(feaProperty->materialID, 7, 1);
        */

    }

	return CAPS_SUCCESS;
}

// Write a Abaqus Load card from a feaLoad structure
int abaqus_writeLoadCard(FILE *fp, const feaLoadStruct *feaLoad, const meshStruct *mesh)
{

    int i, j; // Indexing
    int elemIDindex; // used in PLOAD element ID check logic
    double avgPressure;

    if (fp == NULL) return CAPS_IOERR;
    if (feaLoad == NULL) return CAPS_NULLVALUE;

    // Concentrated force at a grid point
    if (feaLoad->loadType ==  GridForce) {

        for (i = 0; i < feaLoad->numGridID; i++) {

            fprintf(fp,"%s\n", "*CLOAD");
            for (j = 0; j < 3; j++)
              fprintf(fp, "%d, %d, %3.16e\n", feaLoad->gridIDSet[i], j+1, feaLoad->directionVector[j]*feaLoad->forceScaleFactor);
        }
    }

    // Concentrated moment at a grid point
    if (feaLoad->loadType ==  GridMoment) {
        printf("LoadType GridMoment isn't supported\n");
        return CAPS_BADVALUE;
//        for (i = 0; i < feaLoad->numGridID; i++) {
//            fprintf(fp,"%-8s", "MOMENT");
//
//            tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_integerToString(feaLoad->gridIDSet[i], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_integerToString(feaLoad->coordSystemID, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString(feaLoad->momentScaleFactor, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString(feaLoad->directionVector[0], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString(feaLoad->directionVector[1], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString(feaLoad->directionVector[2], fieldWidth, 1);
//            fprintf(fp, "%s%s\n",delimiter, tempString);
//            EG_free(tempString);
//        }
    }

    // Gravitational load
    if (feaLoad->loadType ==  Gravity) {

        fprintf(fp,"%s\n", "*DLOAD");

        fprintf(fp,",GRAV, %3.16e, %3.16e, %3.16e, %3.16e\n", feaLoad->gravityAcceleration,
                                                              feaLoad->directionVector[0],
                                                              feaLoad->directionVector[1],
                                                              feaLoad->directionVector[2]);
    }

    // Pressure load
    if (feaLoad->loadType ==  Pressure) {

        fprintf(fp,"%s\n", "*DLOAD");

        fprintf(fp,"%s, P, %3.16e\n", feaLoad->name, feaLoad->pressureForce);
    }

    // Pressure load at element Nodes
    if (feaLoad->loadType ==  PressureDistribute) {

        printf("LoadType PressureDistribute isn't supported\n");
        return CAPS_BADVALUE;
//        for (i = 0; i < feaLoad->numElementID; i++) {
//            fprintf(fp,"%-8s", "PLOAD4");
//
//            tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_integerToString(feaLoad->elementIDSet[i], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            for (j = 0; j < 4; j++) {
//                tempString = convert_doubleToString(feaLoad->pressureDistributeForce[j], fieldWidth, 1);
//                fprintf(fp, "%s%s",delimiter, tempString);
//                EG_free(tempString);
//            }
//
//            fprintf(fp,"\n");
//        }
    }

    // Pressure load at element Nodes - different distribution per element
    if (feaLoad->loadType ==  PressureExternal) {

        //printf("LoadType PressureExternal isn't supported\n");
        //return CAPS_BADVALUE;

        fprintf(fp,"%s\n", "*DLOAD");
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

            avgPressure = 0;
            for (j = 0; j < mesh_numMeshConnectivity(mesh->element[elemIDindex].elementType); j++) {
                avgPressure += feaLoad->pressureMultiDistributeForce[4*i + j]*feaLoad->loadScaleFactor;
            }
            avgPressure /= mesh_numMeshConnectivity(mesh->element[elemIDindex].elementType);

            //fprintf(fp,"%s, P, %.6f\n", feaLoad->name, avgPressure);
            fprintf(fp,"%d, P, %3.16e\n", mesh->element[elemIDindex].elementID, avgPressure);


            if (mesh->element[elemIDindex].elementType == Quadrilateral ||
                    mesh->element[elemIDindex].elementType == Triangle) {

                for (j = 0; j <  mesh_numMeshConnectivity(mesh->element[elemIDindex].elementType); j++) {

//                    tempString = convert_integerToString(mesh->element[elemIDindex].connectivity[j], fieldWidth, 1);
//                    fprintf(fp, "%s%s",delimiter, tempString);
//                    EG_free(tempString);
                }

            } else {
                printf("Unsupported element type for PLOAD in astrosAIM!\n");
                return CAPS_BADVALUE;
            }
        }
    }

    // Rotational velocity
    if (feaLoad->loadType ==  Rotational) {
        printf("LoadType Rotational isn't supported\n");
        return CAPS_BADVALUE;
//        for (i = 0; i < feaLoad->numGridID; i++) {
//            fprintf(fp,"%-8s", "RFORCE");
//
//            tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_integerToString(feaLoad->gridIDSet[i], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_integerToString(feaLoad->coordSystemID, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString( feaLoad->angularVelScaleFactor, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString( feaLoad->directionVector[0], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString( feaLoad->directionVector[1], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString( feaLoad->directionVector[2], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            fprintf(fp," %7s%-7s\n", "", "+RF");
//            fprintf(fp,"%-8s", "+RF");
//
//            tempString = convert_doubleToString( feaLoad->angularAccScaleFactor, fieldWidth, 1);
//            fprintf(fp, "%s%s\n",delimiter, tempString);
//            EG_free(tempString);
//        }
    }
    // Thermal load at a grid point
    if (feaLoad->loadType ==  Thermal) {
        printf("LoadType Thermal isn't supported\n");
        return CAPS_BADVALUE;
//        fprintf(fp,"%-8s", "TEMPD");
//
//        tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
//        fprintf(fp, "%s%s",delimiter, tempString);
//        EG_free(tempString);
//
//        tempString = convert_doubleToString(feaLoad->temperatureDefault, fieldWidth, 1);
//        fprintf(fp, "%s%s\n",delimiter, tempString);
//        EG_free(tempString);
//
//        for (i = 0; i < feaLoad->numGridID; i++) {
//            fprintf(fp,"%-8s", "TEMP");
//
//            tempString = convert_integerToString(feaLoad->loadID, fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_integerToString(feaLoad->gridIDSet[i], fieldWidth, 1);
//            fprintf(fp, "%s%s",delimiter, tempString);
//            EG_free(tempString);
//
//            tempString = convert_doubleToString(feaLoad->temperature, fieldWidth, 1);
//            fprintf(fp, "%s%s\n",delimiter, tempString);
//            EG_free(tempString);
//        }
    }

    return CAPS_SUCCESS;
}


// Write a Abaqus Element Set card
int abaqus_writeElementSet(void *aimInfo, FILE *fp, char *setName,
                           int numElement, int elementSet[],
                           int numName, char *name[])
{
    int status=CAPS_SUCCESS;

    int i=0;

    int comma;

    // Nothing to do
    if (numElement == 0 && numName == 0) return CAPS_SUCCESS;

    AIM_NOTNULL(setName, aimInfo, status);
    AIM_NOTNULL(elementSet, aimInfo, status);

    fprintf(fp, "*ELSET, ELSET=%s\n", setName);

    comma = (int) false;
    for (i = 0; i < numElement; i++) {

        if (i % 16 == 0 && i>0) {
            fprintf(fp, "\n");
            comma = (int) false;
        }

        if (comma) fprintf(fp, ", ");

        fprintf(fp, "%d", elementSet[i]);

        comma = (int) true;
    }

    if (numName != 0 && name != NULL) {
        for (; i < numElement+numName; i++) {

            if (i % 16 == 0 && i>0) {
                fprintf(fp, "\n");
                comma = (int) false;
            }

            if (comma) fprintf(fp, ", ");

            fprintf(fp, "%s", name[i-numElement]);

            comma = (int) true;
        }
    }

    fprintf(fp, "\n");

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Status %d during abaqus_writeElementSet\n", status);

    return status;
}

// Write a Abaqus Node Set card
int abaqus_writeNodeSet(void *aimInfo, FILE *fp, char *setName,
                        int numNode, int nodeSet[],
                        int numName, char *name[])
{
    int status=CAPS_SUCCESS;

    int i=0;

    int comma;

    // Nothing to do
    if (numNode == 0 && numName == 0) return CAPS_SUCCESS;

    AIM_NOTNULL(setName, aimInfo, status);
    AIM_NOTNULL(nodeSet, aimInfo, status);

    fprintf(fp, "*NSET, NSET=%s\n", setName);

    comma = (int) false;
    for (i = 0; i < numNode; i++) {

        if (i % 16 == 0 && i>0) {
            fprintf(fp, "\n");
            comma = (int) false;
        }

        if (comma) fprintf(fp, ", ");

        fprintf(fp, "%d", nodeSet[i]);

        comma = (int) true;
    }

    if (numName != 0 && name != NULL) {
        for (; i < numNode+numName; i++) {

            if (i % 16 == 0 && i>0) {
                fprintf(fp, "\n");
                comma = (int) false;
            }

            if (comma) fprintf(fp, ", ");

            fprintf(fp, "%s", name[i-numNode]);

            comma = (int) true;
        }
    }

    fprintf(fp, "\n");

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Status %d during abaqus_writeNodeSet\n", status);

    return status;
}

// Write element and node sets for various components of the problem
int abaqus_writeAllSets(void *aimInfo, FILE *fp, const feaProblemStruct *feaProblem) {
    int status;
    int i, j;

    const meshStruct *mesh;
    const feaPropertyStruct *property;
    const feaLoadStruct *load;
    const feaMeshDataStruct *feaData;

    int propertyID;

    int numElement;
    int *elementSet=NULL;

    if (feaProblem == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    // TODO: Add other stuff?
    //      // Constraint information
    //      int numConstraint;
    //      feaConstraintStruct *feaConstraint; // size = [numConstraint]
    //
    //      // Support information
    //      int numSupport;
    //      feaSupportStruct *feaSupport; // size = [numSupport]
    //

    fprintf(fp, "**\n**Node and Element Sets\n**\n");

    mesh = &feaProblem->feaMesh;

    status = array_allocIntegerVector(mesh->numElement, 0, &elementSet);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(elementSet, aimInfo, status);

    for (j = 0; j < feaProblem->numProperty; j++) {

        numElement = 0;

        property = &feaProblem->feaProperty[j];

        for (i = 0; i < mesh->numElement; i++){
            // Grab Structure specific related data if available
            if (mesh->element[i].analysisType == MeshStructure) {
                feaData = (feaMeshDataStruct *) mesh->element[i].analysisData;
                propertyID = feaData->propertyID;
            } else {
                propertyID = mesh->element[i].markerID;
            }

            if (propertyID != property->propertyID) continue;

            elementSet[numElement] = mesh->element[i].elementID;

            numElement += 1;
        }

        if (numElement != 0) {

            status = abaqus_writeElementSet(aimInfo, fp, property->name, numElement, elementSet,  0, NULL);
            if (status != CAPS_SUCCESS) goto cleanup;

        } else{
            printf("Warning: No elements found for property %s\n", property->name);
        }
    }

    if (elementSet != NULL) EG_free(elementSet);
    elementSet = NULL;

    for (i = 0; i < feaProblem->numLoad; i++) {

        load = &feaProblem->feaLoad[i];

        status = abaqus_writeNodeSet(aimInfo, fp, load->name, load->numGridID, load->gridIDSet,  0, NULL);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = abaqus_writeElementSet(aimInfo, fp, load->name, load->numElementID, load->elementIDSet,  0, NULL);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("Error: Status %d during abaqus_writeAllSets\n", status);

        if (elementSet != NULL) EG_free(elementSet);

        return status;
}


// Read data from a Abaqus DAT file and load it into a dataMatrix[numGridPoint][6]
// where variables are U1, U2, U3, UR1, UR2, UR3
int abaqus_readDATDisplacement(void *aimInfo, const char *filename,
                               int numGridPoint, int **nodeID, DOUBLE_6 **dataMatrix)
{
    int status; // Function return

    int i, j; // Indexing

    FILE *fp = NULL;

    size_t linecap = 0;
    char *line = NULL; // Temporary line holder

    const char *header = "       NODE FOOT-  U1             U2             U3             UR1            UR2            UR3";

    const int numVariable = 6; // U1, U2, U3, UR1, UR2, UR3
//    int numDataRead = 0;

    printf("Reading Abaqus DAT file - extracting Displacements!\n");

    fp = aim_fopen(aimInfo, filename, "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open Output file: %s!", filename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Allocate dataMatrix array
    AIM_FREE(*nodeID);
    AIM_FREE(*dataMatrix);
    AIM_ALLOC(*nodeID, numGridPoint, int, aimInfo, status);
    AIM_ALLOC(*dataMatrix, numGridPoint, DOUBLE_6, aimInfo, status);
    for (i = 0; i < numGridPoint; i++) {
      (*nodeID)[i] = -1;
      for (j = 0; j < numVariable; j++)
        (*dataMatrix)[i][j] = 0;
    }

    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {
        if (line == NULL) continue;

        // Look for start of Eigen-Vector
        if (strncmp(header, line, strlen(header)) == 0) {

            // Fast forward 2 lines
            for (i = 0; i < 2; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through the file and fill up the data matrix
            for (i = 0; i < numGridPoint; i++) {

                if (fscanf(fp, "%d", &(*nodeID)[i]) != 1) break;
                for (j = 0; j < numVariable; j++) {

                    if (fscanf(fp, "%lf", &(*dataMatrix)[i][j]) != 1) break;
                    //numDataRead++;
                }
            }

            break;
        }
    }

//    if (numDataRead/numVariable != numGridPoint) {
//        AIM_ERROR(aimInfo, "Failed to read %d grid points. Only found %d.",
//                  numGridPoint, numDataRead/numVariable);
//        status = CAPS_IOERR;
//        goto cleanup;
//    }

    status = CAPS_SUCCESS;
cleanup:
    if (fp != NULL) fclose(fp);
    if (line != NULL) free(line);

    return status;
}



// Read data from a Abaqus DAT file and load it into a dataMatrix[numGridPoint][2]
// where variables are Grid Id, Mises
int abaqus_readDATMises(void *aimInfo, const char *filename,
                        int numElement, int **elemID, double **elemData)
{
    int status; // Function return

    int i; // Indexing

    FILE *fp = NULL;

    size_t linecap = 0;
    char *line = NULL; // Temporary line holder
    int elem = 1;
    int prev = -1;
    int navg = 0;
    int pt, sec;
    double avg = 0, mises;

    const char *header = "    ELEMENT  PT SEC FOOT-   MISES";

    int numDataRead = 0;

    printf("Reading Abaqus DAT file - extracting Mises!\n");

    fp = aim_fopen(aimInfo, filename, "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open Output file: %s!", filename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Allocate arrays
    AIM_FREE(*elemID);
    AIM_FREE(*elemData);
    AIM_ALLOC(*elemID, numElement, int, aimInfo, status);
    AIM_ALLOC(*elemData, numElement, double, aimInfo, status);
    for (i = 0; i < numElement; i++) {
      (*elemID)[i] = -1;
      (*elemData)[i] = 0;
    }

    // Loop through the file again and pull out data
    while (getline(&line, &linecap, fp) >= 0) {
        if (line == NULL) continue;

        // Look for start of Mises
        if (strncmp(header, line, strlen(header)) == 0) {

            // Fast forward 2 lines
            for (i = 0; i < 2; i++) {
                status = getline(&line, &linecap, fp);
                if (status < 0) break;
            }

            // Loop through the file and fill up the data matrix
            elem = -1;
            avg = navg = 0;
            while (getline(&line, &linecap, fp) >= 0) {

                prev = elem;
                if (sscanf(line, "%d %d %d %lf", &elem, &pt, &sec, &mises) != 4) break;
                if (elem != prev && prev > 0) {
                  (*elemID)[numDataRead] = prev;
                  (*elemData)[numDataRead] = avg/navg;
                  navg = 0;
                  avg = 0;
                  numDataRead++;
                }
                navg++;
                avg += mises;
            }
        }
    }

//    if (numDataRead != numElement) {
//        AIM_ERROR(aimInfo, "Failed to read %d elements. Only found %d.\n",
//                  numElement, numDataRead);
//        status = CAPS_IOERR;
//        goto cleanup;
//    }

    status = CAPS_SUCCESS;
cleanup:
    if (fp != NULL) fclose(fp);
    if (line != NULL) free(line);

    return status;
}


int abaqus_readFIL(void *aimInfo, char *filename, int field, int *numData,  double ***dataMatrix) {

    int status;

#ifdef HAVE_PYTHON
    PyObject* mobj = NULL;
#if CYTHON_PEP489_MULTI_PHASE_INIT
    PyModuleDef *mdef = NULL;
    PyObject *modname = NULL;
#endif
#endif
    AIM_NOTNULL(filename, aimInfo, status);
    *numData = 0;
    *dataMatrix = NULL;

#ifdef HAVE_PYTHON
    // Flag to see if Python was already initialized - i.e. the AIM was called from within Python
    int initPy = (int) false;

    // Python thread state variable for new interpreter if needed
    //PyThreadState *newThread = NULL;

    printf("\nUsing Python to read Abaqus FIL file\n");

    // Initialize python
    if (Py_IsInitialized() == 0) {
        printf("\tInitializing Python for Abaqus AIM\n\n");
        Py_InitializeEx(0);
        initPy = (int) true;
    } else {

        /*
      // If python is already initialized - create a new interpreter to run the namelist creator
         printf("\tInitializing new Python thread for FUN3D AIM\n\n");
         PyEv    al_InitThreads();
         newThread = Py_NewInterpreter();
         (void) PyThreadState_Swap(newThread);
        */
        initPy = (int) false;
    }

    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();

     // Taken from "main" by running cython with --embed
    #if PY_MAJOR_VERSION < 3
        initabaqusFILReader();
    #elif CYTHON_PEP489_MULTI_PHASE_INIT
           if (abaqusFILReader_Initialized == (int)false || initPy == (int)true) {
               abaqusFILReader_Initialized = (int)true;
               mobj = PyInit_abaqusFILReader();

               if (!PyModule_Check(mobj)) {
                mdef = (PyModuleDef *) mobj;
                modname = PyUnicode_FromString("abaqusFILReader");
                mobj = NULL;
                if (modname) {
                    mobj = PyModule_NewObject(modname);
                    Py_CLEAR(modname);
                    if (mobj) PyModule_ExecDef(mobj, mdef);
                  }
            }
        }
    #else
        mobj = PyInit_abaqusFILReader();
    #endif

    if (PyErr_Occurred()) {
        PyErr_Print();
        #if PY_MAJOR_VERSION < 3
        if (Py_FlushLine()) PyErr_Clear();
        #endif
        /* Release the thread. No Python API allowed beyond this point. */
        PyGILState_Release(gstate);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (field == 0)
      status = abaqus_parseFILDisplacement(filename, numData, dataMatrix);
    else if (field == 1)
      status = abaqus_parseFILvonMises(filename, numData, dataMatrix);
    else {
      AIM_ERROR(aimInfo, "Unkonwn Field ID %d", field);
      status = CAPS_NOTIMPLEMENT;
    }

    if (status == -1) {
        printf("\tWarning: Python error occurred while reading FIL file: %s\n", filename);
    } else {
        printf("\tDone reading FIL file with Python\n");
    }

    Py_CLEAR(mobj);

    if (PyErr_Occurred()) {
        PyErr_Print();
        #if PY_MAJOR_VERSION < 3
        if (Py_FlushLine()) PyErr_Clear();
        #endif
        /* Release the thread. No Python API allowed beyond this point. */
        PyGILState_Release(gstate);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    /* Release the thread. No Python API allowed beyond this point. */
    PyGILState_Release(gstate);

    // Close down python
    if (initPy == (int) false) {
        /*
        printf("\n");
        printf("\tClosing new Python thread\n");
        Py_EndInterpreter(newThread);
        (void) PyThreadState_Swap(NULL); // This function call is probably not necessary;
         */
    } else {
        printf("\tClosing Python\n");
        Py_Finalize(); // Do not finalize if the AIM was called from Python
    }

#else
    if (field != 0) status = CAPS_NOTIMPLEMENT; // suppress lint
    AIM_ERROR(aimInfo, "Abaqus AIM must be compiled with Python to read displacements");
    status = CAPS_NOTIMPLEMENT;
#endif


cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Status %d during abaqus_readFIL\n", status);

    //if (elementSet != NULL) EG_free(elementSet);

    return status;
}
