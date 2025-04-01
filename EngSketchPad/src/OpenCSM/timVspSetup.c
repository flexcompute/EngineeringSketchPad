/*
 ************************************************************************
 *                                                                      *
 * timVspSetup -- Tool Integration Module for setting up vsp3           *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2025  John F. Dannenhoffer, III (Syracuse University)
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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "OpenCSM.h"
#include "tim.h"
#include "wsserver.h"

#define CINT     const int
#define CDOUBLE  const double
#define CCHAR    const char

static int       outLevel   = 1;

#ifdef WIN32
    #define  SLASH       '\\'
#else
    #define  SLASH       '/'
#endif


/***********************************************************************/
/*                                                                     */
/*   timLoad - open a tim instance                                     */
/*                                                                     */
/***********************************************************************/

int
timLoad(esp_T *ESP,                     /* (in)  pointer to ESP structure */
        void  *text)                    /* (in)  calling text string */
{
    int    status=0;                    /* (out) return status */

    int    i, j, same;
    char   buffer[512], *vsp3_root, command[1024];
    char   *vsp3name=NULL, *udcname=NULL, *vspShortName=NULL;
    FILE   *fp_vsp3=NULL, *fp_vspscript=NULL;

    ROUTINE(timLoad(vspSetup));

    /* --------------------------------------------------------------- */

    SPLINT_CHECK_FOR_NULL(text);

    outLevel = ocsmSetOutLevel(-1);

    if (ESP == NULL) {
        printf("ERROR:: cannot run timVspSetup without serveESP\n");
        status = EGADS_SEQUERR;
        goto cleanup;
    }

    /* create the vspSetup_T structure */
    if (ESP->nudata >= MAX_TIM_NESTING) {
        printf("ERROR:: cannot nest more than %d TIMs\n", MAX_TIM_NESTING);
        exit(0);
    }

    ESP->nudata++;
    MALLOC(ESP->udata[ESP->nudata-1], char, 1);

    strcpy(ESP->timName[ESP->nudata-1], "vspSetup");

    /* initialize the structure */
    GetToken(text, 0, '%', &vsp3name);
    GetToken(text, 1, '%', &udcname );

    SPLINT_CHECK_FOR_NULL(vsp3name);
    SPLINT_CHECK_FOR_NULL(udcname );

    /* make sure vsp3 file exists */
    fp_vsp3 = fopen(vsp3name, "r");
    if (fp_vsp3 == NULL) {
        snprintf(buffer, 511, "timQuit|vspSetup|ERROR:: \"%s\" does not exist", vsp3name);
        tim_bcst("vspSetup", buffer);
//$$$        status = EGADS_NOTFOUND;
        goto cleanup;
    } else {
        fclose(fp_vsp3);
        fp_vsp3 = NULL;
    }

    /* the vspShortName is either the same as vsp3name or a relative filename */
    MALLOC(vspShortName, char, STRLEN(vsp3name)+4);

    strcpy(vspShortName, vsp3name);

    /* find the parts of vsp3name and udcname that are common */
    for (i = 0; i < MIN(STRLEN(vsp3name), STRLEN(udcname)); i++) {
        if (vsp3name[i] == udcname[i]) continue;

        /* now that they disagree (at character i), make sure that
           there are no SLASHes in what is left */
        same = 1;

        for (j = i; j < STRLEN(vsp3name); j++) {
            if (vsp3name[j] == SLASH) {
                same = 0;
                break;
            }
        }

        for (j = i; j < STRLEN(udcname); j++) {
            if (udcname[j] == SLASH) {
                same = 0;
                break;
            }
        }

        /* decrease i until we get back to a SLASH or get to the beginning */
        while (i >= 0) {
            if (vsp3name[i] == SLASH) {
                i++;
                break;
            }
            i--;
        }
        if (i < 0) i = 0;

        /* if same==1, then we can remove the common characters
           in vspShortName */
        if (same == 1) {
            strcpy(vspShortName, "$/");
            strcat(vspShortName, &vsp3name[i]);
        }

        break;
    }

    /* write the .vspscript that will create the .udc file */
    fp_vspscript = fopen("TeMpVsP3.vspscript", "w");
    if (fp_vspscript == NULL) {
        snprintf(buffer, 511, "timQuit|vspSetup|ERROR:: could not create TeMpVsP3.vspscript");
        tim_bcst("vspSetup", buffer);
//$$$        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    fprintf(fp_vspscript, "void main()\n");
    fprintf(fp_vspscript, "{\n");
    fprintf(fp_vspscript, "    string vspName      = \"%s\";\n", vsp3name    );
    fprintf(fp_vspscript, "    string vspShortName = \"%s\";\n", vspShortName);
    fprintf(fp_vspscript, "    string udcName      = \"%s\";\n", udcname     );
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "    // return if unable to open the .udc file\n");
    fprintf(fp_vspscript, "    file udcFile;\n");
    fprintf(fp_vspscript, "    if (udcFile.open(udcName, \"w\") < 0) return;\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "    // prolog\n");
    fprintf(fp_vspscript, "    udcFile.writeString(\"# generated by VspSetup from \" + vspName + \"\\n\\n\");\n");
    fprintf(fp_vspscript, "    udcFile.writeString(\"INTERFACE . ALL\\n\\n\");\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "    // DESPMTR/LBOUND/UOUND statements\n");
    fprintf(fp_vspscript, "    ReadVSPFile(vspName);\n");
    fprintf(fp_vspscript, "    array<string> @id_arr = GetAllUserParms();\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "    for (int i = 0; i < int(id_arr.size()); i++) {\n");
    fprintf(fp_vspscript, "        string parmName   = GetParmName(      id_arr[i]);\n");
    fprintf(fp_vspscript, "        string groupName  = GetParmGroupName( id_arr[i]);\n");
    fprintf(fp_vspscript, "        double parmValue  = GetParmVal(       id_arr[i]);\n");
    fprintf(fp_vspscript, "        double parmUlimit = GetParmUpperLimit(id_arr[i]);\n");
    fprintf(fp_vspscript, "        double parmLlimit = GetParmLowerLimit(id_arr[i]);\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "        // only process User Parms if in ESP_Group\n");
    fprintf(fp_vspscript, "        if (groupName != \"ESP_Group\") continue;\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "        // convert \".\" to \":\" in parmName\n");
    fprintf(fp_vspscript, "        string newParmName = \"\";\n");
    fprintf(fp_vspscript, "        for (uint j = 0; j < parmName.length(); j++) {\n");
    fprintf(fp_vspscript, "            if (parmName.substr(j,1) == \".\") {\n");
    fprintf(fp_vspscript, "                newParmName += \":\";\n");
    fprintf(fp_vspscript, "            } else {\n");
    fprintf(fp_vspscript, "                newParmName += parmName.substr(j,1);\n");
    fprintf(fp_vspscript, "            }\n");
    fprintf(fp_vspscript, "        }\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "        udcFile.writeString(\"   DESPMTR   \" + newParmName + \"    \" + parmValue  + \"\\n\");\n");
    fprintf(fp_vspscript, "        udcFile.writeString(\"   LBOUND    \" + newParmName + \"    \" + parmLlimit + \"\\n\");\n");
    fprintf(fp_vspscript, "        udcFile.writeString(\"   UBOUND    \" + newParmName + \"    \" + parmUlimit + \"\\n\\n\");\n");
    fprintf(fp_vspscript, "    }\n");
    fprintf(fp_vspscript, "\n");
    fprintf(fp_vspscript, "    // epilog\n");
    fprintf(fp_vspscript, "    udcFile.writeString(\"UDPRIM    vsp3    filename $\" + vspShortName + \"\\n\\n\");\n");
    fprintf(fp_vspscript, "    udcFile.writeString(\"END\\n\");\n");
    fprintf(fp_vspscript, "    \n");
    fprintf(fp_vspscript, "    udcFile.close();\n");
    fprintf(fp_vspscript, "}\n");

    fclose(fp_vspscript);
    fp_vspscript = NULL;

    /* execute the vspscript */
    /* exeucute vspscript */
    vsp3_root = getenv("VSP3_ROOT");
    if (vsp3_root != NULL) {
        snprintf(command, 1023, "%s%cvspscript -script TeMpVsP3.vspscript", vsp3_root, SLASH);
    } else {
        snprintf(command, 1023, "vspscript -script TeMpVsP3.vspscript");
    }

    printf("\n====================\nRunning: %s\n", command);
    system(command);
    printf("vspscript has completed\n====================\n\n");

    snprintf(buffer, 511, "timQuit|vspSetup|");
    tim_bcst("vspSetup", buffer);

cleanup:
    FREE(vsp3name    );
    FREE(udcname     );
    FREE(vspShortName);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timMesg - get command, process, and return response               */
/*                                                                     */
/***********************************************************************/

int
timMesg(esp_T *ESP,                     /* (in)  pointer to ESP structure */
/*@unused@*/char  command[])            /* (in)  command */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    i;

    ROUTINE(timMesg(vspSetup));

    /* --------------------------------------------------------------- */

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "vspSetup") != 0) {
        printf("WARNING:: TIM on top of stack is not \"vspSetup\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timSave - save tim data and close tim instance                    */
/*                                                                     */
/***********************************************************************/

int
timSave(esp_T *ESP)                     /* (in)  pointer to ESP structure */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    i;

    ROUTINE(timSave(vspSetup));

    /* --------------------------------------------------------------- */

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "vspSetup") != 0) {
        printf("WARNING:: TIM on top of stack is not \"vspSetup\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    }

    /* cleanup */
    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

    tim_bcst("vspSetup", "timSave|vspSetup|");

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   timQuit - close tim instance without saving                       */
/*                                                                     */
/***********************************************************************/

int
timQuit(esp_T *ESP,                     /* (in)  pointer to ESP structure */
/*@unused@*/int   unload)               /* (in)  flag to unload */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int     i;

    ROUTINE(timQuit(vspSetup));

    /* --------------------------------------------------------------- */

    if (ESP->nudata <= 0) {
        goto cleanup;
    } else if (strcmp(ESP->timName[ESP->nudata-1], "vspSetup") != 0) {
        printf("WARNING:: TIM on top of stack is not \"vspSetup\"\n");
        for (i = 0; i < ESP->nudata; i++) {
            printf("   timName[%d]=%s\n", i, ESP->timName[i]);
        }
        goto cleanup;
    }

    /* cleanup */
    FREE(ESP->udata[ESP->nudata-1]);
    ESP->timName[   ESP->nudata-1][0] = '\0';
    ESP->nudata--;

    tim_bcst("vspSetup", "timQuit|vspSetup|");

cleanup:
    return status;
}
