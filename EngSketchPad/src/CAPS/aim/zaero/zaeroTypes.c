#include "capsTypes.h"
#include "zaeroTypes.h"

#include "aimUtil.h"
#include "feaUtils.h"

/* "constructors" / "destructors" */

int initiate_zaeroLinearFlutterStruct(zaeroLinearFlutterStruct *flutter) {

    if (flutter == NULL) return CAPS_NULLVALUE;

    flutter->function = UnknownFlutterFunction;

    flutter->boundaryCondition = NULL;

    flutter->density = 0.0;

    flutter->refVelocity = 0.0;

    flutter->numVelocities = 0;
    flutter->velocities = NULL;

    flutter->numModes = 0;

    flutter->numDampingFreq = 0;
    flutter->dampingFreq = NULL;

    flutter->numDampingValues = 0;
    flutter->dampingValues = NULL;
    flutter->dampingUnits = NULL;

    flutter->printFlag = 0;

    return CAPS_SUCCESS;
}

int destroy_zaeroLinearFlutterStruct(zaeroLinearFlutterStruct *flutter) {

    if (flutter == NULL) return CAPS_NULLVALUE;

    AIM_FREE(flutter->boundaryCondition);
    AIM_FREE(flutter->velocities);
    AIM_FREE(flutter->dampingFreq);
    AIM_FREE(flutter->dampingValues);
    AIM_FREE(flutter->dampingUnits);

    initiate_zaeroLinearFlutterStruct(flutter);

    return CAPS_SUCCESS;
}

int initiate_zaeroTrimVarStruct(zaeroTrimVarStruct *variable) {

    if (variable == NULL) return CAPS_NULLVALUE;

    variable->name = NULL;
    variable->label = NULL;
    variable->value = NULL;
    variable->lowerLimit = 0.0;
    variable->upperLimit = 0.0;
    variable->dmi = NULL;
    variable->initial = 0.0;

    return CAPS_SUCCESS;
}

int destroy_zaeroTrimVarStruct(zaeroTrimVarStruct *variable) {

    if (variable == NULL) return CAPS_NULLVALUE;

    AIM_FREE(variable->name);
    AIM_FREE(variable->label);
    AIM_FREE(variable->value);
    AIM_FREE(variable->dmi);

    initiate_zaeroTrimVarStruct(variable);

    return CAPS_SUCCESS;
}

int initiate_zaeroTrimStruct(zaeroTrimStruct *trim) {

   if (trim == NULL) return CAPS_NULLVALUE;

   trim->dynamicPressure = 0.0;
   trim->vectorToCG[0] = 0.0;
   trim->vectorToCG[1] = 0.0;
   trim->vectorToCG[2] = 0.0;
   trim->gravityAcceleration = 0.0;
   trim->weight = 0.0;
   trim->weightMomentOfInertia[0] = 0.0;
   trim->weightMomentOfInertia[1] = 0.0;
   trim->weightMomentOfInertia[2] = 0.0;
   trim->weightMomentOfInertia[3] = 0.0;
   trim->weightMomentOfInertia[4] = 0.0;
   trim->weightMomentOfInertia[5] = 0.0;
   trim->accUnit = NULL;
   trim->numVariables = 0;
   trim->variables = NULL;

   return CAPS_SUCCESS;
}

int destroy_zaeroTrimStruct(zaeroTrimStruct *trim) {

    int i;

    if (trim == NULL) return CAPS_NULLVALUE;

    AIM_FREE(trim->accUnit);

    if (trim->variables != NULL) {
        for (i = 0; i < trim->numVariables; i++) {
            if (trim->variables[i] != NULL) {
              AIM_FREE(trim->variables[i]);
            }
        }
        AIM_FREE(trim->variables);
    }

    initiate_zaeroTrimStruct(trim);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) a zaero subcase structure
int initiate_zaeroSubcaseStruct(zaeroSubcaseStruct *zaeroSubcase) {

    if (zaeroSubcase == NULL) return CAPS_NULLVALUE;

    zaeroSubcase->name = NULL;
    zaeroSubcase->subtitle = NULL;
    zaeroSubcase->subcaseID = 0;
    zaeroSubcase->label = NULL;
    zaeroSubcase->uaicName = NULL;

    zaeroSubcase->analysisID = 0;
    zaeroSubcase->analysisType = UnknownAnalysis;

    // discipline
    zaeroSubcase->disciplineType = UnknownDiscipline;
    zaeroSubcase->discipline = NULL;

    // initiate_zaeroUAICStruct(&zaeroSubcase->uaic);

    return CAPS_SUCCESS;
}

// Destroy the subcase discipline according to discipline type
int destroy_zaeroDiscipline(zaeroSubcaseStruct *zaeroSubcase) {

    int status = CAPS_SUCCESS;

    if (zaeroSubcase->discipline != NULL) {
        switch (zaeroSubcase->disciplineType) {
          case LinearFlutter:
            status = destroy_zaeroLinearFlutterStruct(zaeroSubcase->discipline);
            break;
          case ParamFlutter:
            status = CAPS_NOTIMPLEMENT;
            break;
          case Aeroservoelastic:
            status = CAPS_NOTIMPLEMENT;
            break;
          case StaticAeroelastic:
            status = destroy_zaeroTrimStruct(zaeroSubcase->discipline);
            break;
          case ManeuverLoads:
            status = CAPS_NOTIMPLEMENT;
            break;
          case EjectionLoads:
            status = CAPS_NOTIMPLEMENT;
            break;
          case GustLoads:
            status = CAPS_NOTIMPLEMENT;
            break;
          case MFTGustLoads:
            status = CAPS_NOTIMPLEMENT;
            break;
          case NonLinearFlutter:
            status = CAPS_NOTIMPLEMENT;
            break;
          case UnknownDiscipline:
          default:
            printf("Unknown disciplineType %d and discipline is not NULL\n",
                   zaeroSubcase->disciplineType);
            status = CAPS_BADVALUE;
        }
        if (status != CAPS_SUCCESS) {
          printf("[WARNING] Status %d during destroy_zaeroDiscipline\n",
                 status);
        }
    }

    AIM_FREE(zaeroSubcase->discipline);

    return status;
}

// Destroy (0 out all values and NULL all pointers) zaero subcase structure
int destroy_zaeroSubcaseStruct(zaeroSubcaseStruct *zaeroSubcase) {

    int status;

    if (zaeroSubcase == NULL) return CAPS_NULLVALUE;

    AIM_FREE(zaeroSubcase->name);
    AIM_FREE(zaeroSubcase->subtitle);
    AIM_FREE(zaeroSubcase->label);
    AIM_FREE(zaeroSubcase->uaicName);

    status = destroy_zaeroDiscipline(zaeroSubcase);
    if (status != CAPS_SUCCESS) return status;

    // reset all variables
    initiate_zaeroSubcaseStruct(zaeroSubcase);

    return CAPS_SUCCESS;
}

int initiate_zaeroHFGStruct(zaeroHFGStruct *hfg) {

    if (hfg == NULL) return CAPS_NULLVALUE;

    hfg->XZSymmetry = NULL;
    hfg->flip = NULL;

    hfg->refChord = 0.0;
    hfg->refSpan = 0.0;
    hfg->refArea = 0.0;

    hfg->refCenter[0] = 0.0;
    hfg->refCenter[1] = 0.0;
    hfg->refCenter[2] = 0.0;

    return CAPS_SUCCESS;
}

int destroy_zaeroHFGStruct(zaeroHFGStruct *hfg) {

    if (hfg == NULL) return CAPS_NULLVALUE;

    AIM_FREE(hfg->XZSymmetry);
    AIM_FREE(hfg->flip);

    // reset all variables
    initiate_zaeroHFGStruct(hfg);

    return  CAPS_SUCCESS;
}

int initiate_zaeroFEMStruct(zaeroFEMStruct *fem) {

    if (fem == NULL) return CAPS_NULLVALUE;

    fem->filename = NULL;
    fem->form = NULL;
    fem->boundary = NULL;
    fem->suport[0] = 0;
    fem->suport[1] = 0;
    fem->printFlag = 0;
    fem->asetFlag = 0;

    return CAPS_SUCCESS;
}

int destroy_zaeroFEMStruct(zaeroFEMStruct *fem) {

    if (fem == NULL) return CAPS_NULLVALUE;

    AIM_FREE(fem->filename);
    AIM_FREE(fem->form);
    AIM_FREE(fem->boundary);

    initiate_zaeroFEMStruct(fem);

    return CAPS_SUCCESS;
}

int initiate_zaeroUAICStruct(zaeroUAICStruct *uaic) {

    if (uaic == NULL) return CAPS_NULLVALUE;

    uaic->name = NULL;

    uaic->machNumber = 0.0;
    uaic->methodFlag = 0;
    uaic->aicFilename[0] = '\0';

    uaic->numReducedFreq = 0;
    uaic->reducedFreq = NULL;

    uaic->saveFlag = 0;
    uaic->printFlag = 0;

    return CAPS_SUCCESS;
}

int destroy_zaeroUAICStruct(zaeroUAICStruct *uaic) {

    if (uaic == NULL) return CAPS_NULLVALUE;

    AIM_FREE(uaic->name);
    AIM_FREE(uaic->reducedFreq);

    // reset all variables
    initiate_zaeroUAICStruct(uaic);

    return CAPS_SUCCESS;
}

int initiate_zaeroSplineStruct(zaeroSplineStruct *spline) {

    if (spline == NULL) return CAPS_NULLVALUE;

    spline->method = 1;
    spline->eps = 1.0e-6;
    spline->attachFlex = 0.0;

    return CAPS_SUCCESS;
}

int destroy_zaeroSplineStruct(zaeroSplineStruct *spline) {

    if (spline == NULL) return CAPS_NULLVALUE;

    // reset all variables
    initiate_zaeroSplineStruct(spline);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) a zaero problem structure
int initiate_zaeroProblemStruct(zaeroProblemStruct *zaeroProblem) {

    if (zaeroProblem == NULL) return CAPS_NULLVALUE;

    zaeroProblem->title = NULL;
    zaeroProblem->echo = NULL;
    zaeroProblem->numSubcases = 0;
    zaeroProblem->subcases = NULL;

    // FEM module data
    zaeroProblem->numFEMs = 0;
    initiate_zaeroFEMStruct(&zaeroProblem->FEMs[0]);
    initiate_zaeroFEMStruct(&zaeroProblem->FEMs[1]);

    // Aerodynamics
    zaeroProblem->numAero = 0;
    zaeroProblem->feaAero = NULL;

    zaeroProblem->numTrimVariables = 0;
    zaeroProblem->trimVariables = NULL;

    // initialize modules
    initiate_zaeroHFGStruct(&zaeroProblem->hfg);
    initiate_zaeroSplineStruct(&zaeroProblem->spline);
    zaeroProblem->numUAICs = 0;
    zaeroProblem->UAICs = NULL;

    // // artifacts
    // initiate_zaeroArtifactsStruct(&zaeroProblem->artifacts);

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) zaero problem structure
int destroy_zaeroProblemStruct(zaeroProblemStruct *zaeroProblem) {

    int i; // Indexing

    int status; // Status return

    if (zaeroProblem == NULL) return CAPS_NULLVALUE;

    // Title
    AIM_FREE(zaeroProblem->title);

    // Echo
    AIM_FREE(zaeroProblem->echo);

    // Modules
    for (i = 0; i < zaeroProblem->numFEMs; i++) {
      status = destroy_zaeroFEMStruct(&zaeroProblem->FEMs[i]);
      if (status != CAPS_SUCCESS) {
          printf("[WARNING] Status %d during destroy_zaeroFEMStruct\n",
                 status);
      }
    }

    status = destroy_zaeroHFGStruct(&zaeroProblem->hfg);
    if (status != CAPS_SUCCESS) {
        printf("[WARNING] Status %d during destroy_zaeroHFGStruct\n",
               status);
    }

    status = destroy_zaeroSplineStruct(&zaeroProblem->spline);
    if (status != CAPS_SUCCESS) {
        printf("[WARNING] Status %d during destroy_zaeroSplineStruct\n",
               status);
    }

    // UAICs
    if (zaeroProblem->UAICs != NULL) {
        for (i = 0; i < zaeroProblem->numUAICs; i++) {
            status = destroy_zaeroUAICStruct(&zaeroProblem->UAICs[i]);
            if (status != CAPS_SUCCESS) {
                printf("[WARNING] Status %d during destroy_zaeroUAICStruct\n",
                       status);
            }
        }
        AIM_FREE(zaeroProblem->UAICs);
    }

    // Subcases
    if (zaeroProblem->subcases != NULL) {
        for (i = 0; i < zaeroProblem->numSubcases; i++) {
            status = destroy_zaeroSubcaseStruct(&zaeroProblem->subcases[i]);
            if (status != CAPS_SUCCESS) {
                printf("[WARNING] Status %d during destroy_zaeroSubcaseStruct\n",
                       status);
            }
        }
        AIM_FREE(zaeroProblem->subcases);
    }

    // Aerodynamics
    if (zaeroProblem->feaAero != NULL) {

        for (i = 0; i < zaeroProblem->numAero; i++) {
            status = destroy_feaAeroStruct(&zaeroProblem->feaAero[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaAeroStruct\n", status);
        }
    }
    AIM_FREE(zaeroProblem->feaAero);
    zaeroProblem->numAero = 0;

    // Trim variables
    if (zaeroProblem->trimVariables != NULL) {
        for (i = 0; i < zaeroProblem->numTrimVariables; i++) {
            status = destroy_zaeroTrimVarStruct(&zaeroProblem->trimVariables[i]);
            if (status != CAPS_SUCCESS) {
                printf("[WARNING] Status %d during destroy_zaeroTrimVarStruct\n",
                       status);
            }
        }
        AIM_FREE(zaeroProblem->trimVariables);
    }

    // // Artifacts
    // status = destroy_zaeroArtifactsStruct(&zaeroProblem->artifacts);
    // if (status != CAPS_SUCCESS)  {
    //   printf("[WARNING] Status %d during destroy_zaeroArtifactsStruct!\n", status);
    // }

    // reset all variables
    initiate_zaeroProblemStruct(zaeroProblem);

    return CAPS_SUCCESS;
}

int initiate_zaeroArtifactsStruct(zaeroArtifactsStruct *artifacts) {

    if (artifacts == NULL) return CAPS_NULLVALUE;

    artifacts->input = NULL;
    artifacts->output = NULL;
    artifacts->numAIC = 0;
    artifacts->aic = NULL;

    return CAPS_SUCCESS;
}

int destroy_zaeroArtifactsStruct(zaeroArtifactsStruct *artifacts) {

    int i;

    if (artifacts == NULL) return CAPS_NULLVALUE;

    AIM_FREE(artifacts->input);
    AIM_FREE(artifacts->output);

    if (artifacts->aic != NULL) {
        for (i = 0; i < artifacts->numAIC; i++) {
            if (artifacts->aic[i].name != NULL) {
              AIM_FREE(artifacts->aic[i].name);
            }
            if (artifacts->aic[i].value != NULL) {
              AIM_FREE(artifacts->aic[i].value);
            }
        }
        AIM_FREE(artifacts->aic);
    }

    initiate_zaeroArtifactsStruct(artifacts);

    return CAPS_SUCCESS;
}
