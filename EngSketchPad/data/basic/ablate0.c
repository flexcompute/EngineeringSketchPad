/*
 *******************************************************************
 *                                                                 *
 * ablate0 [noise] -- create ablate0.cloud file for Plugs          *
 *                                                                 *
 *              Written by John Dannenhoffer @ Syracuse University *
 *                                                                 *
 *******************************************************************
 */

/*
 * Copyright (C) 2013/2023  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later
 *     version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free
 *    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *    Boston, MA  02110-1301  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  MAX_BUMP  10


/******************************************************************/
/*                                                                */
/*   main - main program                                          */
/*                                                                */
/******************************************************************/

int
main(int       argc,               /* (in)  number of arguments */
     char      *argv[])            /* (in)  array of arguments */
{

    int    ipnt, npnt=10000;
    double x, y, z, r, noise, temp;

    int    nbump=0, ibump, btype[MAX_BUMP];
    double xcent[MAX_BUMP], ycent[MAX_BUMP], rad[MAX_BUMP], depth[MAX_BUMP];
    
    FILE   *fp;

    /* get inputs from user */
    printf("enter npnt: ");
    scanf("%d", &npnt);

    while (nbump < MAX_BUMP-1) {
        printf("enter 1 for step, 2 for gaussian: ");
        scanf("%d", &btype[nbump]);
        if (btype[nbump] != 1 && btype[nbump] != 2) break;

        printf("enter xcent (0-4): ");
        scanf("%lf", &xcent[nbump]);

        printf("enter ycent (0-3): ");
        scanf("%lf", &ycent[nbump]);

        printf("enter rad: ");
        scanf("%lf", &rad[nbump]);

        printf("enter depth: ");
        scanf("%lf", &depth[nbump]);

        nbump++;
    }

    printf("enter noise: ");
    scanf("%lf", &noise);

    srand(12345);

    fp = fopen("ablate0.cloud", "w");
    fprintf(fp, "%5d%5d ablate0.cloud\n", npnt, 0);

    for (ipnt = 0; ipnt < npnt; ipnt++) {

        /* point on Face 6 */
        x = 0.02 + 3.96 * (double)(rand()) / (double)(RAND_MAX);
        y = 0.02 + 2.96 * (double)(rand()) / (double)(RAND_MAX);
        z = 2.00;

        /* add bumps */
        for (ibump = 0; ibump < nbump; ibump++) {
            r = sqrt((x-xcent[ibump])*(x-xcent[ibump]) + (y-ycent[ibump])*(y-ycent[ibump]));

            if (r < rad[ibump]) {
                if (btype[ibump] == 1) {
                    z += depth[ibump];
                } else {
                    z += depth[ibump] * exp(-rad[ibump]/(rad[ibump]*rad[ibump]-r*r));
                }
            }
        }

        /* add noise */
        temp = (double)(rand()) / (double)(RAND_MAX);
        z += noise * (temp - 0.5);

        fprintf(fp, "%12.6f %12.6f %12.6f\n", x, y, z);
    }

    fprintf(fp, "%5d%5d end\n", 0, 0);
    fclose(fp);

    return EXIT_SUCCESS;
}
        
    
