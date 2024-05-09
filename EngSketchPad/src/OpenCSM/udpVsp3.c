/*
 ************************************************************************
 *                                                                      *
 * udpVsp3 -- read a .vsp3 file (or a .stp file exported from vsp)      *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2024  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

#define NUMUDPARGS 3
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME( IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define KEEPTEMPS(IUDP)  ((int    *) (udps[IUDP].arg[1].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"filename", "keeptemps", "recycle",   };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,   ATTRINT,     ATTRRECYCLE, };
static int    argIdefs[NUMUDPARGS] = {0,          0,           0,           };
static double argDdefs[NUMUDPARGS] = {0.,         0.,          0.,          };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#ifdef GRAFIC
   #include "grafic.h"
#endif

#ifdef WIN32
    #define  SLASH       '\\'
    #define  SLEEP(msec)  Sleep(msec)
#else
    #include <unistd.h>
    #define  SLASH       '/'
    #define  SLEEP(msec)  usleep(1000*msec)
#endif

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           EPS20           1.0e-20
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

/* message (to be shared by all functions) */
static char message[1024];

static int processStepFile(ego context, char filename[], ego *ebody);


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     lenname, ipmtr, ptype, nrow, ncol, i;
    double  value, dot;
    char    filename[MAX_LINE_LEN+1], *vsp3_root, command[1024], pname[MAX_NAME_LEN];
    FILE    *fp_vspscript;
    void    *modl;
    modl_T  *MODL;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename( 0) = %s\n", FILENAME( 0));
    printf("keeptemps(0) = %d\n", KEEPTEMPS(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size == 0) {
        snprintf(message, 1023, "\"filename\" must be given");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 1023, "\"keeptemps\" must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename( %d) = %s\n", numUdp, FILENAME( numUdp));
    printf("keeptemps(%d) = %d\n", numUdp, KEEPTEMPS(numUdp));
#endif

    /* process based upon type file filename given */
    strncpy(filename, FILENAME(numUdp), MAX_LINE_LEN);
    lenname = strlen(filename);

    /* since OpenCSM converts all forward slashes in filenames to backslashes
       (before we enter the UDP), we need to convert them back since anglescript
       (in vspscript) only allows forward slashes in filenames */
    for (i = 0; i < lenname; i++) {
        if (filename[i] == SLASH) filename[i] = '/';
    }

    /* filename is a .stp file */
    if        (lenname > 4 && strcmp(&filename[lenname-4], ".stp" ) == 0) {
        status = processStepFile(context, filename, ebody);
        CHECK_STATUS(processStepFile);

    /* filename is a .vsp3 file */
    } else if (lenname > 5 && strcmp(&filename[lenname-5], ".vsp3") == 0) {

        /* get pointer to OpenCSM MODL */
        status = EG_getUserPointer(context, (void**)(&(modl)));
        CHECK_STATUS(EG_getUserPointer);

        MODL = (modl_T *)modl;
        if (MODL->perturb != NULL) {
            MODL = MODL->perturb;
        }

        /* create the TeMpVsP3.vspscript file */
        fp_vspscript = fopen("TeMpVsP3.vspscript", "w");
        if (fp_vspscript == NULL) {
            snprintf(message, 1023, "could not create \"TeMpVsP3.vspscript\"");
            status = EGADS_NOTFOUND;
            goto cleanup;
        }

        /* .vspscript prolog */
        fprintf(fp_vspscript, "void main()\n");
        fprintf(fp_vspscript, "{\n");
        fprintf(fp_vspscript, "    string vspname = \"%s\";\n", filename);
        fprintf(fp_vspscript, "    ReadVSPFile(vspname);\n\n");

        /* update vsp UserParms:ESP_Group:* from DESPMTRs */
        fprintf(fp_vspscript, "    string user_ctr = FindContainer(\"UserParms\", 0);\n");
        fprintf(fp_vspscript, "    SilenceErrors();\n");
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            status = ocsmGetPmtr(MODL, ipmtr, &ptype, &nrow, &ncol, pname);
            CHECK_STATUS(ocsmGetPmtr);

            if (ptype == OCSM_DESPMTR && nrow == 1 && ncol == 1) {
                status = ocsmGetValu(MODL, ipmtr, 1, 1, &value, &dot);
                CHECK_STATUS(ocsmGetValu);

                for (i = 0; i < STRLEN(pname); i++) {
                    if (pname[i] == ':') pname[i] = '.';    // OpenCSM uses ':' but VPS uses '.'
                }

                fprintf(fp_vspscript, "    SetParmVal(FindParm(user_ctr, \"%s\", \"ESP_Group\"), %20.14e);\n", pname, value);
            }
        }
        fprintf(fp_vspscript, "    PrintOnErrors();\n\n");

        /* .vspscript epilog */
        fprintf(fp_vspscript, "    string veh_id = GetVehicleID();\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"SplitSurfs\", \"STEPSettings\"), 1);\n\n");

        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelID\", \"STEPSettings\"), 1);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelName\", \"STEPSettings\"), 1);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelSurfNo\", \"STEPSettings\"), 1);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelDelim\", \"STEPSettings\"), DELIM_COMMA);\n\n");

        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"SplitSubSurfs\", \"STEPSettings\"), 0);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"MergePoints\", \"STEPSettings\"), 0);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"ToCubic\", \"STEPSettings\"), 0);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"TrimTE\", \"STEPSettings\"), 0);\n");
        fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"ExportPropMainSurf\", \"STEPSettings\"), 0);\n\n");

        fprintf(fp_vspscript, "    string stpname = \"TeMpVsP3.stp\";\n");
        fprintf(fp_vspscript, "    ExportFile(stpname, SET_ALL, EXPORT_STEP);\n");
        fprintf(fp_vspscript, "}\n");

        fclose(fp_vspscript);

        /* exeucute vspscript */
        vsp3_root = getenv("VSP3_ROOT");
        if (vsp3_root != NULL) {
            snprintf(command, 1023, "%s%cvspscript -script TeMpVsP3.vspscript", vsp3_root, SLASH);
        } else {
            snprintf(command, 1023, "vspscript -script TeMpVsP3.vspscript");
        }

        printf("\n====================\nRunning: %s\n", command);
        system(command);
        SLEEP(1000);
        printf("vspscript has completed\n====================\n\n");

        /* process the .stp file */
        status = processStepFile(context, "TeMpVsP3.stp", ebody);
        CHECK_STATUS(processStepFile);

        /* clean up temporary files */
        if (KEEPTEMPS(0) == 0) {
            remove("TeMpVsP3.stp");
            remove("TeMpVsP3.vspscript");
        }

    /* unknown filename type */
    } else {
        snprintf(message, 1023, "\"%s\" is not a .stp or .vsp3 file", filename);
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (strlen(message) > 0) {
        MALLOC(*string, char, 1024);
        strncpy(*string, message, 1023);
    } else if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                      *
 ************************************************************************
 */

int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
   /*@unused@*/int    npnt,             /* (in)  number of points */
   /*@unused@*/int    entType,          /* (in)  OCSM entity type */
   /*@unused@*/int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
   /*@unused@*/double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        status = EGADS_NOTMODEL;
        goto cleanup;
    }

    /* this routine is not written yet */
    status = EGADS_NOLOAD;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   processStepFile - read .stp file nd create Model with the Bodys    *
 *                                                                      *
 ************************************************************************
 */

static int
processStepFile(ego    context,         /* (in)  EGADS context */
                char   filename[],      /* (in)  filename */
                ego    *ebody)          /* (out) Model */
{
    int     status = EGADS_SUCCESS;

    int     attrType, attrLen, oclass, mtype, nchild, ichild, *senses;
    int     i, j, igeomID, igeomName, isurfNum, nbody=0, nface=0, nnn, *sss;
    int     nchild2, nface2, periodic, *header;
    CINT    *tempIlist;
    double  data[18], uvrange[4], bbox[6], *rdata, mprops[14];
    double  xmin=+HUGEQ, xmax=-HUGEQ, ymin=+HUGEQ, ymax=-HUGEQ, zmin=+HUGEQ, zmax=-HUGEQ, eps06=EPS06;
    CDOUBLE *tempRlist;
    char    oldGeomID[80], geomID[80], oldGeomName[80], geomName[80], oldSurfNum[80], surfNum[80], bodyName[80];
    CCHAR   *tempClist;
    ego     emodel, eref, *echilds, *efaces=NULL, *newBodys=NULL, newModel, *eee, eshell;
    ego     *efaces2, esurf, *echilds2;

    ROUTINE(processStepFile);

    /* --------------------------------------------------------------- */

    /* read the .stp file (created by exporting from OpenVSP) */
    status = EG_loadModel(context, 0, filename, &emodel);
    CHECK_STATUS(EG_loadModel);

    /* get the Bodys contained in emodel */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    printf("There are %d children in emodel\n", nchild);

    MALLOC(efaces,   ego, nchild);
    MALLOC(newBodys, ego, nchild);

    nbody          = 0;
    oldGeomID[0]   = '\0';
    oldGeomName[0] = '\0';
    oldSurfNum[0]  = '\0';

    /* get the Name of each Body */
    for (ichild = 0; ichild < nchild; ichild++) {

        /* if the Body does not have a Name, skip it */
        status = EG_attributeRet(echilds[ichild], "Name", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status != EGADS_SUCCESS || attrType != ATTRSTRING) {
            printf("Skipping   Child %3d (does not have Name)\n", ichild+1);
            continue;
        } else {
            /* if the Body has no area, skip it */
            status = EG_getMassProperties(echilds[ichild], mprops);
            CHECK_STATUS(EG_getMassProperties);

            if (mprops[1] < EPS12) {
                printf("Skipping   Child %3d (area=%12.6ef)\n", ichild+1, mprops[1]);
                continue;
            }

            printf("Processing Child %3d (%s), nface=%d\n", ichild+1, tempClist, nface);
        }


        /* extract the geomID from the Name */
        igeomID   = 0;
        geomID[0] = '\0';

        for (i = 0; i < attrLen; i++) {
            if (tempClist[i] == ',') {
                i++;
                break;
            }

            geomID[igeomID  ] = tempClist[i];
            geomID[igeomID+1] = '\0';
            igeomID++;
        }

        /* extract the geomName from the Name */
        igeomName   = 0;
        geomName[0] = '\0';

        for ( ; i < attrLen; i++) {
            if (tempClist[i] == ',') {
                i++;
                break;
            }

            if (tempClist[i] == ' ') {
                if (igeomName == 0) {
                    continue;
                } else {
                    geomName[igeomName  ] = '_';
                    geomName[igeomName+1] = '\0';
                    igeomName++;
                }
            } else {
                geomName[igeomName  ] = tempClist[i];
                geomName[igeomName+1] = '\0';
                igeomName++;
            }
        }

        /* extract the surfNum from the Name */
        isurfNum   = 0;
        surfNum[0] = '\0';

        for ( ; i < attrLen; i++) {
            if (tempClist[i] == ',') {
                i++;
                break;
            }

            if (tempClist[i] == ' ') {
                if (isurfNum == 0) {
                    continue;
                } else {
                    surfNum[isurfNum  ] = '_';
                    surfNum[isurfNum+1] = '\0';
                    isurfNum++;
                }
            } else {
                surfNum[isurfNum  ] = tempClist[i];
                surfNum[isurfNum+1] = '\0';
                isurfNum++;
            }
        }

        /* if this has a different component name than the previous Body,
           make the SolidBody (or SheetBody if it is open) from the FaceBodys
           that have been previously processed (if any) */
        if (strcmp(geomID, oldGeomID) != 0 || strcmp(surfNum, oldSurfNum) != 0) {
            if (nface > 0) {
                eps06 = EPS06 * MAX(MAX(xmax-xmin,ymax-ymin),MAX(zmax-zmin,100));

                if (nface > 1) {
                    status = EG_sewFaces(nface, efaces, eps06, 1, &newModel);
                    CHECK_STATUS(EG_sewFaces);

                    status = EG_getTopology(newModel, &eref, &oclass, &mtype, data,
                                            &nnn, &eee, &sss);
                    CHECK_STATUS(EG_getTopology);

                    if (nnn != 1) {
                        printf("EG_sewFaces(nface=%d) generated %d Bodys.  Only using first Body\n", nface, nnn);
                        ocsmPrintEgo(newModel);
                    }

                    status = EG_copyObject(eee[0], NULL, &(newBodys[nbody++]));
                    CHECK_STATUS(EG_copyObject);

                    status = EG_deleteObject(newModel);
                    CHECK_STATUS(EG_deleteObject);
                } else {
                    status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL,
                                             1, efaces, NULL, &eshell);
                    CHECK_STATUS(EG_makeTopology);

                    status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL,
                                             1, &eshell, NULL, &(newBodys[nbody++]));
                    CHECK_STATUS(EG_makeTopology);
                }

                snprintf(bodyName, 79, "%s.%s:%d", oldGeomName, oldSurfNum, nbody);

                status = EG_attributeAdd(newBodys[nbody-1], "_name", ATTRSTRING, STRLEN(bodyName),
                                         NULL, NULL, bodyName);
                CHECK_STATUS(EG_attributeAdd);

                status = EG_attributeAdd(newBodys[nbody-1], "_vspBody", ATTRINT, 1, &nbody, NULL, NULL);
                CHECK_STATUS(EG_attributeAdd);

                printf("   Made Body %3d (%s) with %d Faces\n", nbody-1, bodyName, nface);

                nface = 0;
                xmin  = +HUGEQ;
                xmax  = -HUGEQ;
                ymin  = +HUGEQ;
                ymax  = -HUGEQ;
                zmin  = +HUGEQ;
                zmax  = -HUGEQ;
            }
            strcpy(oldGeomID,   geomID  );
            strcpy(oldGeomName, geomName);
            strcpy(oldSurfNum,  surfNum );
        }

        /* since older versions of OpenVSP can create knot vectors with
           jumps near the end, we will extract the Surface from the Face,
           adjust the knot vector, and remake the Surface and Face */
        status = EG_getBodyTopos(echilds[ichild], NULL, FACE, &nface2, &efaces2);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_getTopology(efaces2[0], &esurf, &oclass, &mtype, data,
                                &nchild2, &echilds2, &senses);
        CHECK_STATUS(EG_getTopology);

        EG_free(efaces2);

        status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &header, &rdata);
        CHECK_STATUS(EG_getGeometry);

        /* fix jumps in knot vectors (fixes error in OpenVSP) */
        for (i = 1; i < header[3]; i++) {
            if (rdata[i]-rdata[i-1] > 1.01) {
                for (j = i; j < header[3]; j++) {
                    rdata[j] -= 1;
                }
                i--;
            }
        }

        for (i = 1; i < header[6]; i++) {
            if (rdata[header[3]+i]-rdata[header[3]+i-1] > 1.01) {
                for (j = i; j < header[6]; j++) {
                    rdata[header[3]+j] -= 1;
                }
                i--;
            }
        }

        status = EG_makeGeometry(context, oclass, mtype, NULL, header, rdata, &esurf);
        CHECK_STATUS(EG_makeGeometry);

        EG_free(header);
        EG_free(rdata );

        status = EG_getRange(esurf, uvrange, &periodic);
        CHECK_STATUS(EG_getRange);

        status = EG_makeFace(esurf, SFORWARD, uvrange, &(efaces[nface]));
        CHECK_STATUS(EG_makeFace);

        status = EG_attributeAdd(efaces[nface], "_vspID", ATTRSTRING, strlen(oldGeomID), NULL, NULL, oldGeomID);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_attributeAdd(efaces[nface], "_vspFace", ATTRINT, 1, &nface, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_getBoundingBox(efaces[nface], bbox);
        CHECK_STATUS(EG_getBoundingBox);

        xmin = MIN(xmin, bbox[0]);
        ymin = MIN(ymin, bbox[1]);
        zmin = MIN(zmin, bbox[2]);
        xmax = MAX(xmax, bbox[3]);
        ymax = MAX(ymax, bbox[4]);
        zmax = MAX(zmax, bbox[5]);

        nface++;
    } /* next Body from the .stp file */

    /* make Sheet/SolidBodys from last FaceBody(s) */
    if (nface > 0) {
        eps06 = EPS06 * MAX(MAX(xmax-xmin,ymax-ymin),MAX(zmax-zmin,100));

        if (nface > 1) {
            status = EG_sewFaces(nface, efaces, eps06, 1, &newModel);
            CHECK_STATUS(EG_sewFaces);

            status = EG_getTopology(newModel, &eref, &oclass, &mtype, data,
                                    &nnn, &eee, &sss);
            CHECK_STATUS(EG_getTopology);

            if (nnn != 1) {
                printf("EG_sewFaces(nface=%d) generated %d Bodys.  Only using first Body\n", nface, nnn);
                ocsmPrintEgo(newModel);
            }

            status = EG_copyObject(eee[0], NULL, &(newBodys[nbody++]));
            CHECK_STATUS(EG_copyObject);

            status = EG_deleteObject(newModel);
            CHECK_STATUS(EG_deleteObject);

        } else {
            status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL,
                                     1, efaces, NULL, &eshell);
            CHECK_STATUS(EG_makeTopology);

            status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL,
                                     1, &eshell, NULL, &(newBodys[nbody++]));
            CHECK_STATUS(EG_makeTopology);
        }

        snprintf(bodyName, 79, "%s.%s:%d", oldGeomName, oldSurfNum, nbody);

        status = EG_attributeAdd(newBodys[nbody-1], "_name", ATTRSTRING, STRLEN(bodyName),
                                 NULL, NULL, bodyName);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_attributeAdd(newBodys[nbody-1], "_vspBody", ATTRINT, 1, &nbody, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        printf("   Made Body %d (%s) with %d Faces\n", nbody-1, bodyName, nface);
    }

    /* make a Model to return to OpenCSM */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, nbody, newBodys, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

cleanup:
    FREE(efaces  );
    FREE(newBodys);

    return status;
}
