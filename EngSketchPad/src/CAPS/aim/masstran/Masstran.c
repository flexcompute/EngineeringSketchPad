/*
 ************************************************************************
 *                                                                      *
 * Masstran -- read .bdf file and compute mass properties               *
 *                                                                      *
 *              Initial version by Adam Steward @ Syracuse University   *
 *              Modified by        Dannenhoffer @ Syracuse University   *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2012/2021  John F. Dannenhoffer, III (Syracuse University)
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
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define NINT(A)   (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))
#define MIN(A,B)  (((A) < (B)) ? (A) : (B))
#define MAX(A,B)  (((A) < (B)) ? (B) : (A))

/*
1-------2-------3-------4-------5-------6-------7-------8-------9-------10------
GRID    GID     CP      X1      X2      X3      CD      PS      SEID

GRID*   GID             CP              X1              X2              +
*       X3              CD              PS              SEID

CROD    EID     PID     G1      G2

CTRIA3  EID     PID     G1      G2      G3

CQUAD4  EID     PID     G1      G2      G3      G4

CSHEAR  EID     PID     G1      G2      G3      G4

CTETRA  EID     PID     G1      G2      G3      G4

CPYRAM  EID     PID     G1      G2      G3      G4      G5

CPENTA  EID     PID     G1      G2      G3      G4      G5      G6

CHEXA   EID     PID     G1      G2      G3      G4      G5      G6      +CH
+CH     G7      G8

PSHELL  PID     MID1    T       MID2    12/T..3 MID3    TS/T    NSM
        Z1      Z2      MID4

MAT1    MID     E       G       NU      RHO     A       TREF    GE
        ST      SC      SS      MCSID

*/

typedef struct {
    int    cp;        /* coordinate ID (bias-1) */ 
    double x1;        /* location of point in cp */
    double x2;
    double x3;
    int    cd;        /* coordinate system for displacements */
    int    ps;        /* single-point constraint */
    int    seid;      /* super-element ID (bias-1) */
} grid_T;

typedef struct {
    char   type[8];   /* CTRIA3, CQUAD4, CSHEAR, CTETRA, CPYRAM, CPENTA, or CHEXA */
    int    pid;       /* property ID (bias-1) */
    int    g1;        /* grid     ID (bias-1) */
    int    g2;
    int    g3;
    int    g4;
    int    g5;
    int    g6;
    int    g7;
    int    g8;
} elem_T;

typedef struct {
    int    mid;       /* material ID (bias-1) */
 
   double T;         /* thickness */
} prop_T;

typedef struct {
    double E;         /* Young's mdoulus */
    double G;         /* shear modulus */
    double NU;        /* Poisson ratio */
    double RHO;       /* mass density */
} matl_T;


/*
 ************************************************************************
 *                                                                      *
 *   nextCard - read next card                                          *
 *                                                                      *
 ************************************************************************
 */

void nextCard(FILE *fp, char name[], double fields[])
{
    int  ichar, jchar, ifield;
    char card[161], field[17];

    /* initialize the outputs */
    name[0] = '\0';
    for (ifield = 0; ifield < 11; ifield++) {
        fields[ifield] = 0;
    }

    /* read the next card (but skip comments and blank cards) */
    while (!feof(fp)) {
        (void) fgets(card, 160, fp);
        if (feof(fp)) return;
        if (card[0] != '$' && card[0] != ' ' && card[0] != 10) break;
    }

    /* if it contains commas, then we have free format */
    if (strstr(card, ",") != NULL) {

        /* add a comma to the end */
        ichar = strlen(card);
        card[ichar  ] = ',';
        card[ichar+1] = '\0';

        /* start reading frombeginning */
        ichar = 0;

        /* pull out name */
        jchar = 0;
        while (ichar < strlen(card)) {
            if (card[ichar] != ',') {
                if (card[ichar] != ' ') {
                    name[  jchar] = card[ichar++];
                    name[++jchar] ='\0';
                } else {
                    ichar++;
                }
            } else {
                ichar++;
                break;
            }
        }

        /* pull out the fields */
        for (ifield = 2; ifield < 11; ifield++) {
        
            jchar = 0;
            while (ichar < strlen(card)) {
                if (card[ichar] != ',') {
                    field[  jchar] = card[ichar++];
                    field[++jchar] = '\0';
                } else {
                    ichar++;
                    fields[ifield] = strtod(field, NULL);
                    break;
                }
            }
        }

    /* if title contains a *, then we are long format */
    } else if (strstr(card, "*") != NULL) {
        /* pull out name*/
        name[0] = '\0';
        for (ichar = 0; ichar < MIN(8, strlen(card)); ichar++) {
            if (card[ichar] != ' ') {
                name[ichar  ] = card[ichar];
                name[ichar+1] = '\0';
            } else {
                break;
            }
        }

        /* pull out fields on first card */
        for (ifield = 2; ifield <= 5; ifield++) {
            strncpy(field, &(card[ifield*16-24]), 16);
            field[16] = '\0';

            fields[ifield] = strtod(field, NULL);
        }

        /* pull out fields on second card */
        (void) fgets(card, 160, fp);
        for (ifield = 6; ifield <= 9; ifield++) {
            strncpy(field, &(card[ifield*16-88]), 16);
            field[16] = '\0';

            fields[ifield] = strtod(field, NULL);
        }

    /* otherwise we have short format */
    } else {
        /* pull out name*/
        name[0] = '\0';
        for (ichar = 0; ichar < MIN(8, strlen(card)); ichar++) {
            if (card[ichar] != ' ') {
                name[ichar  ] = card[ichar];
                name[ichar+1] = '\0';
            } else {
                break;
            }
        }

        /* pull out fields */
        for (ifield = 2; ifield <= 10; ifield++) {
            strncpy(field, &(card[ifield*8-8]), 8);
            field[8] = '\0';

            fields[ifield] = strtod(field, NULL);
        }
    }
}


/*
 ************************************************************************
 *                                                                      *
 *   magCross - compute magnitude of cross product                      *
 *                                                                      *
 ************************************************************************
 */

double mag_cross(double a[], double b[])
{
    double c[3];

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];

    return sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
}


/*
 ************************************************************************
 *                                                                      *
 *   get_areaTri - compute area of CTRIA3                               *
 *                                                                      *
 ************************************************************************
 */

double get_areaTri(double x[], double y[], double z[])
{
    double vec1[3], vec2[3];

    vec1[0] = x[1] - x[0];
    vec1[1] = y[1] - y[0];
    vec1[2] = z[1] - z[0];

    vec2[0] = x[2] - x[0];
    vec2[1] = y[2] - y[0];
    vec2[2] = z[2] - z[0];

    return mag_cross(vec1, vec2) / 2.0;
}


/*
 ************************************************************************
 *                                                                      *
 *   get_areaQuad - compute area of CQUAD4                              *
 *                                                                      *
 ************************************************************************
 */

double get_areaQuad(double x[], double y[], double z[])
{
    double vec1[3], vec2[3];

    vec1[0] = x[2] - x[0];
    vec1[1] = y[2] - y[0];
    vec1[2] = z[2] - z[0];

    vec2[0] = x[3] - x[1];
    vec2[1] = y[3] - y[1];
    vec2[2] = z[3] - z[1];

    return mag_cross(vec1, vec2) / 2.0;
}


/*
 ************************************************************************
 *                                                                      *
 *   main - min program                                                 *
 *                                                                      *
 ************************************************************************
 */

int main(int argc, char *argv[])
{
    int    status = EXIT_SUCCESS;

    char   name[9];
    double fields[11];

    double area   = 0;
    double cgxmom = 0;
    double cgymom = 0;
    double cgzmom = 0;
    double cxmom  = 0;
    double cymom  = 0;
    double czmom  = 0;
    double weight = 0;
    double Ixx    = 0;
    double Ixy    = 0;
    double Ixz    = 0;
    double Iyy    = 0;
    double Iyz    = 0;
    double Izz    = 0;
    double Cx, Cy, Cz, CGx, CGy, CGz;
    double xelem[4], yelem[4], zelem[4], myArea, thick, density;
    double xcent, ycent, zcent;
    
    FILE *fp;

    int igid, ngid=0;
    int ieid, neid=0;
    int ipid, npid=0;
    int imid, nmid=0;

    grid_T *grid=NULL;
    elem_T *elem=NULL;
    prop_T *prop=NULL;
    matl_T *matl=NULL;

    if (argc < 2) {
        printf("Proper usage: Masstran filename\n");
        status = EXIT_FAILURE;
        goto cleanup;
    }

    fp = fopen(argv[1], "r");
    if (fp == NULL) {
        printf("Error opening file \"%s\"\n", argv[1]);
        status = EXIT_FAILURE;
        goto cleanup;
    }

    /* find the size of required arrays by reading cards until while file is read */
    while (!feof(fp)) {
        nextCard(fp, name, fields);

        if (strcmp(name, "GRID" ) == 0 ||
            strcmp(name, "GRID*") == 0   ) {
            igid = NINT(fields[2]);
            if (igid > ngid) ngid = igid;
        } else if (strcmp(name, "CROD") == 0) {

        } else if (strcmp(name, "CTRIA3" ) == 0) {
            ieid = NINT(fields[2]);
            if (ieid > neid) neid = ieid;
            ipid = NINT(fields[3]);
            if (ipid > npid) npid = ipid;
        } else if (strcmp(name, "CQUAD4" ) == 0) {
            ieid = NINT(fields[2]);
            if (ieid > neid) neid = ieid;
            ipid = NINT(fields[3]);
            if (ipid > npid) npid = ipid;
        } else if (strcmp(name, "PSHELL") == 0) {
            ipid = NINT(fields[2]);
            if (ipid > npid) npid = ipid;
            imid = NINT(fields[3]);
            if (imid > nmid) nmid = imid;
        } else if (strcmp(name, "MAT1") == 0) {
            imid = NINT(fields[2]);
            if (imid > nmid) nmid= imid;
        } else if (name[0] == '$') {

        } else if (name[0] == '\0') {

        } else {
            printf("SKIPPING \"%s\"\n", name);
        }
    }    

    printf("\nProblem size:\n");
    printf("     ngid=%5d\n", ngid);
    printf("     neid=%5d\n", neid);
    printf("     npid=%5d\n", npid);
    printf("     nmid=%5d\n", nmid);

    /* allocate the needed arrays */
    grid = (grid_T *) malloc((ngid+1)*sizeof(grid_T));
    elem = (elem_T *) malloc((neid+1)*sizeof(elem_T));
    prop = (prop_T *) malloc((npid+1)*sizeof(prop_T));
    matl = (matl_T *) malloc((nmid+1)*sizeof(matl_T));

    if (grid == NULL || elem == NULL || prop == NULL || matl == NULL) {
        printf("malloc problem\n");
        status = EXIT_FAILURE;
        goto cleanup;
    }

    /* initialize the structures */
    for (igid = 0; igid <= ngid; igid++) {
        grid[igid].cp = -1;
    }
    for (ieid = 0; ieid <= neid; ieid++) {
        elem[ieid].type[0] = '\0';
    }
    for (ipid = 0; ipid <= npid; ipid++) {
        prop[ipid].mid = 0;
    }
    for (imid = 0; imid <= nmid; imid++) {
        matl[imid].E = -1;
    }

    /* re-read input file and fill structures (which are accessed bias-1) */
    rewind(fp);

    while (!feof(fp)) {
        nextCard(fp, name, fields);

        if (strcmp(name, "GRID" ) == 0 ||
            strcmp(name, "GRID*") == 0   ) {
            igid            = NINT(fields[2]);
            grid[igid].cp   = NINT(fields[3]);
            grid[igid].x1   =      fields[4];
            grid[igid].x2   =      fields[5];
            grid[igid].x3   =      fields[6];
            grid[igid].cd   = NINT(fields[7]);
            grid[igid].ps   = NINT(fields[8]);
            grid[igid].seid = NINT(fields[9]);
        } else if (strcmp(name, "CTRIA3" ) == 0) {
            ieid            = NINT(fields[2]);
            elem[ieid].pid  = NINT(fields[3]);
            elem[ieid].g1   = NINT(fields[4]);
            elem[ieid].g2   = NINT(fields[5]);
            elem[ieid].g3   = NINT(fields[6]);

            strcpy(elem[ieid].type, "CTRIA3");
        } else if (strcmp(name, "CQUAD4" ) == 0) {
            ieid            = NINT(fields[2]);
            elem[ieid].pid  = NINT(fields[3]);
            elem[ieid].g1   = NINT(fields[4]);
            elem[ieid].g2   = NINT(fields[5]);
            elem[ieid].g3   = NINT(fields[6]);
            elem[ieid].g4   = NINT(fields[7]);

            strcpy(elem[ieid].type, "CQUAD4");
        } else if (strcmp(name, "PSHELL") == 0) {
            ipid            = NINT(fields[2]);
            prop[ipid].mid  = NINT(fields[3]);
            prop[ipid].T    =      fields[4];
        } else if (strcmp(name, "MAT1") == 0) {
            imid            = NINT(fields[2]);
            matl[imid].E    =      fields[3];
            matl[imid].G    =      fields[4];
            matl[imid].NU   =      fields[5];
            matl[imid].RHO  =      fields[6];
        }
    }

    fclose(fp);

    /* check for data structure coherency */
    for (ieid = 1; ieid <= neid; ieid++) {
        if        (strcmp(elem[ieid].type, "CTRIA3") == 0) {
            ipid = elem[ieid].pid;
            if        (ipid < 0 || ipid > npid) {
                printf("illegal elem[ieid=%d].pid=%d\n", ieid, ipid);
                status = EXIT_FAILURE;
            } else if (ipid == 0) {
            } else if (prop[ipid].mid < 0 || prop[ipid].mid > nmid) {
                printf("illegal elem[ieid=%d].pid=%d\n", ieid, elem[ieid].pid);
                status = EXIT_FAILURE;
            }

            imid = prop[ipid].mid;
            if (imid < 0 || imid > nmid) {
                printf("illegal prop[ipid=%d].mid=%d\n", ipid, imid);
                status = EXIT_FAILURE;
            } else if (imid == 0) {
            } else if (matl[imid].E < 0) {
                printf("illegal prop[ipid=%d].mid=%d\n", ipid, imid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g1;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g1=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g1=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g2;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g2=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g2=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g3;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g3=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g3=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }
        } else if (strcmp(elem[ieid].type, "CQUAD4") == 0) {
            ipid = elem[ieid].pid;
            if        (ipid < 0 || ipid > npid) {
                printf("illegal elem[ieid=%d].pid=%d\n", ieid, ipid);
                status = EXIT_FAILURE;
            } else if (ipid == 0) {
            } else if (prop[ipid].mid < 0 || prop[ipid].mid > nmid) {
                printf("illegal elem[ieid=%d].pid=%d\n", ieid, elem[ieid].pid);
                status = EXIT_FAILURE;
            }

            imid = prop[ipid].mid;
            if (imid < 0 || imid > nmid) {
                printf("illegal prop[ipid=%d].mid=%d\n", ipid, imid);
                status = EXIT_FAILURE;
            } else if (imid == 0) {
            } else if (matl[imid].E < 0) {
                printf("illegal prop[ipid=%d].mid=%d\n", ipid, imid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g1;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g1=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g1=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g2;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g2=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g2=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g3;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g3=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g3=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }

            igid = elem[ieid].g4;
            if        (igid < 0 || igid > ngid) {
                printf("illegal elem[ieid=%d].g4=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            } else if (grid[igid].cp < 0) {                
                printf("illegal elem[ieid=%d].g4=%d\n", ieid, igid);
                status = EXIT_FAILURE;
            }
        }
    }

    if (status == EXIT_FAILURE) goto cleanup;

    /* compute sums */
    for (ieid = 1; ieid <= neid; ieid++) {
        if (strcmp(elem[ieid].type, "CTRIA3") == 0) {
            xelem[0] = grid[elem[ieid].g1].x1;
            yelem[0] = grid[elem[ieid].g1].x2;
            zelem[0] = grid[elem[ieid].g1].x3;
            xelem[1] = grid[elem[ieid].g2].x1;
            yelem[1] = grid[elem[ieid].g2].x2;
            zelem[1] = grid[elem[ieid].g2].x3;
            xelem[2] = grid[elem[ieid].g3].x1;
            yelem[2] = grid[elem[ieid].g3].x2;
            zelem[2] = grid[elem[ieid].g3].x3;

            xcent  = (xelem[0] + xelem[1] + xelem[2]) / 3;
            ycent  = (yelem[0] + yelem[1] + yelem[2]) / 3;
            zcent  = (zelem[0] + zelem[1] + zelem[2]) / 3;
            
            myArea = get_areaTri(xelem, yelem, zelem);

            thick   = 1;
            density = 1;
            if (elem[ieid].pid > 0) {
                thick = prop[elem[ieid].pid].T;
                if (prop[elem[ieid].pid].mid > 0) {
                    density = matl[prop[elem[ieid].pid].mid].RHO;
                }
            }

            area   += myArea;
            weight += myArea * density * thick;

            cxmom  += xcent * myArea;
            cymom  += ycent * myArea;
            czmom  += zcent * myArea;
            
            cgxmom += xcent * myArea * density * thick;
            cgymom += ycent * myArea * density * thick;
            cgzmom += zcent * myArea * density * thick;

            Ixx    += (ycent * ycent + zcent * zcent) * myArea;
            Ixy    -= (xcent         * ycent        ) * myArea;
            Ixz    -= (xcent         * zcent        ) * myArea;
            Iyy    += (xcent * xcent + zcent * zcent) * myArea;
            Iyz    -= (ycent         * zcent        ) * myArea;
            Izz    += (xcent * xcent + ycent * ycent) * myArea;
        } else if (strcmp(elem[ieid].type,"CQUAD4") == 0) {
            xelem[0] = grid[elem[ieid].g1].x1;
            yelem[0] = grid[elem[ieid].g1].x2;
            zelem[0] = grid[elem[ieid].g1].x3;
            xelem[1] = grid[elem[ieid].g2].x1;
            yelem[1] = grid[elem[ieid].g2].x2;
            zelem[1] = grid[elem[ieid].g2].x3;
            xelem[2] = grid[elem[ieid].g3].x1;
            yelem[2] = grid[elem[ieid].g3].x2;
            zelem[2] = grid[elem[ieid].g3].x3;
            xelem[3] = grid[elem[ieid].g4].x1;
            yelem[3] = grid[elem[ieid].g4].x2;
            zelem[3] = grid[elem[ieid].g4].x3;

            xcent  = (xelem[0] + xelem[1] + xelem[2] + xelem[3]) / 4;
            ycent  = (yelem[0] + yelem[1] + yelem[2] + yelem[3]) / 4;
            zcent  = (zelem[0] + zelem[1] + zelem[2] + zelem[3]) / 4;
            
            myArea = get_areaQuad(xelem, yelem, zelem);

            thick   = 1;
            density = 1;
            if (elem[ieid].pid > 0) {
                thick = prop[elem[ieid].pid].T;
                if (prop[elem[ieid].pid].mid > 0) {
                    density = matl[prop[elem[ieid].pid].mid].RHO;
                }
            }

            area   += myArea;
            weight += myArea * density * thick;

            cxmom  += xcent * myArea;
            cymom  += ycent * myArea;
            czmom  += zcent * myArea;
            
            cgxmom += xcent * myArea * density * thick;
            cgymom += ycent * myArea * density * thick;
            cgzmom += zcent * myArea * density * thick;

            Ixx    += (ycent * ycent + zcent * zcent) * myArea;
            Ixy    -= (xcent         * ycent        ) * myArea;
            Ixz    -= (xcent         * zcent        ) * myArea;
            Iyy    += (xcent * xcent + zcent * zcent) * myArea;
            Iyz    -= (ycent         * zcent        ) * myArea;
            Izz    += (xcent * xcent + ycent * ycent) * myArea;
        }
    }

    /* compute statistics for whole Body */
    Cx  = cxmom / area;
    Cy  = cymom / area;
    Cz  = czmom / area;

    if (weight > 0) {
        CGx = cgxmom / weight;
        CGy = cgymom / weight;
        CGz = cgzmom / weight;
    } else {
        CGx = Cx;
        CGy = Cy;
        CGz = Cz;
    }

    Ixx -= area * (Cy * Cy + Cz * Cz);
    Ixy += area *  Cx * Cy;
    Ixz += area *  Cx * Cz;
    Iyy -= area * (Cx * Cx + Cz * Cz);
    Iyz += area *  Cy * Cz;
    Izz -= area * (Cx * Cx + Cy * Cy);

    /* report statistics */
    printf("\nStatistics:\n");
    printf("Area:     %lf\n", area  );
    printf("Weight:   %lf\n", weight);
    printf("Centroid: %lf, %lf, %lf\n", Cx,  Cy,  Cz );
    printf("CG:       %lf, %lf, %lf\n", CGx, CGy, CGz);
    printf("Ixx:      %lf\n", Ixx);
    printf("Iyy:      %lf\n", Iyy);
    printf("Izz:      %lf\n", Izz);
    printf("Ixy:      %lf\n", Ixy);
    printf("Ixz:      %lf\n", Ixz);
    printf("Iyz:      %lf\n", Iyz);

cleanup:
    if (grid != NULL) free(grid);
    if (elem != NULL) free(elem);
    if (prop != NULL) free(prop);
    if (matl != NULL) free(matl);

    exit(status);
}

