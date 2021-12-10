/*
 ************************************************************************
 *                                                                      *
 * Demo1_server -- demonstration of EGADS and distributed EGADS-Lite    *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2017  John F. Dannenhoffer, III (Syracuse University)
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
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mpi.h"

#include "egads.h"
#include "OpenCSM.h"

static int massProps(ego etess, int iface, double props[]);


/***********************************************************************/
/*                                                                     */
/*   main program - server                                             */
/*                                                                     */
/***********************************************************************/

int
main(int    argc,                       /* (in)  number of arguments */
     char   *argv[])                    /* (in)  array  of arguments */
{
    int          status, myRank, numRanks, source, dest, tag;
    int          ibody, atype, alen, nnode, nedge, nface;
    int          senses, buildTo, builtTo, nbody, myNbytes, iface;
    int          Nnode, Npnt;
    const int    *tempIlist;
    double       rootProps[5], myProps[5], totProps[5];
    const double *tempRlist;
    char         message[100], *stream;
    const char   *tempClist;
    MPI_Status   MPIstat;
    MPI_Comm     myCommWorld;
    void         *modl;
    modl_T       *MODL;
    ego          ebodys=NULL, emodel;
    size_t       nbytes;

    /* --------------------------------------------------------------- */

    /* make sure that we get a .csm filename */
    assert (argc >= 2);

    /* start MPI */
    status = MPI_Init(&argc, &argv);
    assert (status == 0);

    /* make a copy of MPI_COMM_WORLD so that MPI operations
       in this application do not collide with those in a
       library, etc. */
    status = MPI_Comm_dup(MPI_COMM_WORLD, &myCommWorld);
    assert (status == 0);

    /* determine my rank */
    status = MPI_Comm_rank(myCommWorld, &myRank);
    assert (status == 0);

    /* determine the number of ranks */
    status = MPI_Comm_size(myCommWorld, &numRanks);
    assert (status == 0);

    /* this (server) process should only be rank 0 */
    assert (myRank == 0);

    /* load .csm file */
    status = ocsmLoad(argv[1], &modl);
    assert (status == 0);

    MODL = (modl_T*) modl;

    /* build BRep */
    nbody   = 0;
    buildTo = 0;
    status = ocsmBuild(modl, buildTo, &builtTo, &nbody, NULL);
    assert (status == 0);

    /* use the last Body on the stack */
    for (ibody = MODL->nbody; ibody > 0; ibody--) {
        if (MODL->body[ibody].onstack == 1) {
            ebodys = MODL->body[ibody].ebody;
            senses = SFORWARD;

            break;
        }
    }
    assert (ebodys != NULL);

    /* generate a stream to send to the clients */
    status = EG_makeTopology(MODL->context, NULL, MODEL, 0,
                             NULL, 1, &ebodys, &senses, &emodel);
    assert (status == 0);

    status = EG_exportModel(emodel, &nbytes, &stream);
    assert (status == 0);

    printf("Broadcasting to clients...\n");

    /* broadcast size of stream so that clients can malloc arrays */
    myNbytes = (int)nbytes;
    status = MPI_Bcast((void*)(&myNbytes), 1, MPI_INT, 0, myCommWorld);
    assert (status == 0);

    /* broadcast the stream to the clients */
    status = MPI_Bcast((void*)stream, (int)nbytes, MPI_CHAR, 0, myCommWorld);
    assert (status == 0);

    EG_free(stream);

    /* compute the mass properties */
    rootProps[0] = 0;
    rootProps[1] = 0;
    rootProps[2] = 0;
    rootProps[3] = 0;
    rootProps[4] = 0;
    
    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
        status = massProps(MODL->body[ibody].etess, iface, rootProps);
        assert (status == 0);
    }
    rootProps[2] /= rootProps[0];
    rootProps[3] /= rootProps[0];
    rootProps[4] /= rootProps[0];

    /* accumulate the volume integrals that were computed by the clients */
    myProps[0] = 0;
    myProps[1] = 0;
    myProps[2] = 0;
    myProps[3] = 0;
    myProps[4] = 0;

    status = MPI_Reduce(myProps, totProps, 5, MPI_DOUBLE, MPI_SUM, 0, myCommWorld);
    assert (status == 0);

    totProps[2] /= totProps[0];
    totProps[3] /= totProps[0];
    totProps[4] /= totProps[0];

    /* compare the accumulated volume integrals with the one computed locally */
    printf("Comparing results...\n");
    printf("volume %10.5f %10.5f (%12.4e)\n", rootProps[0], totProps[0], rootProps[0]-totProps[0]);
    printf("area   %10.5f %10.5f (%12.4e)\n", rootProps[1], totProps[1], rootProps[1]-totProps[1]);
    printf("xcg    %10.5f %10.5f (%12.4e)\n", rootProps[2], totProps[2], rootProps[2]-totProps[2]);
    printf("ycg    %10.5f %10.5f (%12.4e)\n", rootProps[3], totProps[3], rootProps[3]-totProps[3]);
    printf("zcg    %10.5f %10.5f (%12.4e)\n", rootProps[4], totProps[4], rootProps[4]-totProps[4]);

    /* remove the copy of my world communicator */
    status = MPI_Comm_free(&myCommWorld);
    assert (status == 0);

    /* shut down MPI */
    status = MPI_Finalize();
    assert (status == 0);
}


/***********************************************************************/
/*                                                                     */
/*   massProps - compute mass properties for an ego                    */
/*                                                                     */
/***********************************************************************/

static int
massProps(ego    etess,                 /* (in)  tessellation object */
          int    iface,                 /* (in)  Face index (bias-1) */
          double props[])               /* (both) properties
                                                 [0] = volume
                                                 [1] = surface area/length
                                                 [2] = xint
                                                 [3] = yint
                                                 [4] = zint */
{
    int       status = SUCCESS;         /* return status */

    int           npnt, ntri, itri, ip0, ip1, ip2;
    const int     *ptype, *pindx, *tris, *tric;
    double        volume, area, xint, yint, zint;
    double        xarea, yarea, zarea, xcent, ycent, zcent;
    const double  *xyz, *uv;

    /* --------------------------------------------------------------- */

    /* initialize the summations */
    volume = 0;
    area   = 0;
    xint   = 0;
    yint   = 0;
    zint   = 0;

    /* loop through Traiangle on given Face */
    status = EG_getTessFace(etess, iface,
                            &npnt, &xyz, &uv, &ptype, &pindx,
                            &ntri, &tris, &tric);
    if (status != 0) goto cleanup;

    for (itri = 0; itri < ntri; itri++) {
        ip0 = tris[3*itri  ] - 1;
        ip1 = tris[3*itri+1] - 1;
        ip2 = tris[3*itri+2] - 1;

        xarea = ((xyz[3*ip1+1]-xyz[3*ip0+1]) * (xyz[3*ip2+2]-xyz[3*ip0+2])
               - (xyz[3*ip1+2]-xyz[3*ip0+2]) * (xyz[3*ip2+1]-xyz[3*ip0+1]));
        yarea = ((xyz[3*ip1+2]-xyz[3*ip0+2]) * (xyz[3*ip2  ]-xyz[3*ip0  ])
               - (xyz[3*ip1  ]-xyz[3*ip0  ]) * (xyz[3*ip2+2]-xyz[3*ip0+2]));
        zarea = ((xyz[3*ip1  ]-xyz[3*ip0  ]) * (xyz[3*ip2+1]-xyz[3*ip0+1])
               - (xyz[3*ip1+1]-xyz[3*ip0+1]) * (xyz[3*ip2  ]-xyz[3*ip0  ]));

        xcent = (xyz[3*ip0  ] + xyz[3*ip1  ] + xyz[3*ip2  ]);
        ycent = (xyz[3*ip0+1] + xyz[3*ip1+1] + xyz[3*ip2+1]);
        zcent = (xyz[3*ip0+2] + xyz[3*ip1+2] + xyz[3*ip2+2]);

        area   += sqrt(xarea*xarea         + yarea*yarea         + zarea*zarea        );
        volume +=     (xarea*xcent         + yarea*ycent         + zarea*zcent        );
        xint   +=     (xarea*xcent*xcent/2 + yarea*ycent*xcent   + zarea*zcent*xcent  );
        yint   +=     (xarea*xcent*ycent   + yarea*ycent*ycent/2 + zarea*zcent*ycent  );
        zint   +=     (xarea*xcent*zcent   + yarea*ycent*zcent   + zarea*zcent*zcent/2);
    }

    /* add the mass properties into the array */
    props[ 0] += volume / 18;
    props[ 1] += area   /  2;
    props[ 2] += xint   / 54;
    props[ 3] += yint   / 54;
    props[ 4] += zint   / 54;

cleanup:
    return status;
}
