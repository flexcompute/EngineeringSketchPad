/*
 ************************************************************************
 *                                                                      *
 * simple_server -- demonstration of EGADS and distributed EGADS-Lite   *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *                       Bob Haimes        @ MIT                        *
 *                                                                      *
 ************************************************************************
 */


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mpi.h"

#include "egads.h"



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
  int           status = EGADS_SUCCESS;       /* return status */
  
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


/***********************************************************************/
/*                                                                     */
/*   main program - server                                             */
/*                                                                     */
/***********************************************************************/

int
main(int    argc,                       /* (in)  number of arguments */
     char   *argv[])                    /* (in)  array  of arguments */
{
  int          status, myRank, numRanks, *senses;
  int          atype, alen, nnode, nface, oclass, mtype;
  int          nbody, myNbytes, iface;
  int          Nnode, Npnt;
  double       rootProps[5], myProps[5], totProps[5], data[4];
  char         *stream;
  const int    *tempIlist;
  const double *tempRlist;
  const char   *tempClist;
  MPI_Status   MPIstat;
  MPI_Comm     myCommWorld;
  ego          context, emodel, eref, etess, *ebodys;
  size_t       nbytes;
  
  /* --------------------------------------------------------------- */
  
  /* make sure that we get a filename */
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
  
  /* open EGADS */
  status = EG_open(&context);
  assert (status == 0);
  
  /* load file */
  status = EG_loadModel(context, 0, argv[1], &emodel);
  assert (status == 0);
  
  status = EG_exportModel(emodel, &nbytes, &stream);
  assert (status == 0);
/*
  {
    FILE *fp;
    
    fp = fopen("stream.lite", "wb");
    fwrite(stream, sizeof(char), nbytes, fp);
    fclose(fp);
  }
*/
  printf("Broadcasting to clients...\n");
  
  /* broadcast size of stream so that clients can malloc arrays */
  myNbytes = (int) nbytes;
  status = MPI_Bcast((void *)(&myNbytes), 1, MPI_INT, 0, myCommWorld);
  assert (status == 0);
  
  /* broadcast the stream to the clients */
  status = MPI_Bcast((void *)stream, (int)nbytes, MPI_CHAR, 0, myCommWorld);
  assert (status == 0);
  
  EG_free(stream);
  
  /* extract the Body from the Model */
  status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                          data, &nbody, &ebodys, &senses);
  assert (status == 0);
  
  /* make the tessellation */
  status = EG_attributeRet(ebodys[0], "_tParams", &atype, &alen,
                           &tempIlist, &tempRlist, &tempClist);
  assert (status == 0);
  status  = EG_makeTessBody(ebodys[0], (double *) tempRlist, &etess);
  assert (status == 0);
  status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, NULL);
  assert (status == 0);
  
  /* compute the mass properties */
  rootProps[0] = 0;
  rootProps[1] = 0;
  rootProps[2] = 0;
  rootProps[3] = 0;
  rootProps[4] = 0;
  
  for (iface = 1; iface <= nface; iface++) {
    status = massProps(etess, iface, rootProps);
    assert (status == 0);
  }
  rootProps[2] /= rootProps[0];
  rootProps[3] /= rootProps[0];
  rootProps[4] /= rootProps[0];
  
  if (numRanks == 2) {
    /* receive the Nodes from the server */
    double *xyz_client, xyz_server[18];
    int    inode, nchild;
    ego    *enodes, *echilds;
    
    printf("Comparing Nodes...\n");
    
    status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, &enodes);
    assert (status == 0);
    
    xyz_client = (double *) malloc(3*nnode*sizeof(double));
    assert (xyz_client != NULL);
    
    status = MPI_Recv((void *)(&Nnode),   1,       MPI_INT,    MPI_ANY_SOURCE,
                      100, myCommWorld, &MPIstat);
    assert (status == 0);
    
    assert (Nnode == nnode);
    
    status = MPI_Recv((void *)xyz_client, 3*nnode, MPI_DOUBLE, MPI_ANY_SOURCE,
                      200, myCommWorld, &MPIstat);
    assert (status == 0);
    
    for (inode = 0; inode < nnode; inode++) {
      status = EG_getTopology(enodes[inode], &eref, &oclass, &mtype,
                              xyz_server, &nchild, &echilds, &senses);
      assert (status == 0);
      
      if (xyz_server[0] != xyz_client[3*inode  ] ||
          xyz_server[1] != xyz_client[3*inode+1] ||
          xyz_server[2] != xyz_client[3*inode+2]   ) {
        printf("node %3d  %12.5f %12.5f (%12.4e)\n", inode,
               xyz_server[0], xyz_client[3*inode  ],
               xyz_server[0]- xyz_client[3*inode  ]);
        printf("          %12.5f %12.5f (%12.4e)\n",
               xyz_server[1], xyz_client[3*inode+1],
               xyz_server[1]- xyz_client[3*inode+1]);
        printf("          %12.5f %12.5f (%12.4e)\n",
               xyz_server[2], xyz_client[3*inode+2],
               xyz_server[2]- xyz_client[3*inode+2]);
        
        assert (fabs(xyz_server[0]-xyz_client[0]) < 1e-15);
        assert (fabs(xyz_server[1]-xyz_client[1]) < 1e-15);
        assert (fabs(xyz_server[2]-xyz_client[2]) < 1e-15);
      }
    }
    
    EG_free(enodes);
    free(xyz_client);
  }
  
  if (numRanks == 2) {
    /* receive the Edges from the server */
    int          nedge, iedge, npnt, ipnt;
    double       *xyz_client, *t_client;
    const double *xyz_server, *t_server;
    
    printf("Comparing Edges...\n");
    
    status = EG_getBodyTopos(ebodys[0], NULL, EDGE, &nedge, NULL);
    assert (status == 0);
    
    for (iedge = 1; iedge <= nedge; iedge++) {
      status = EG_getTessEdge(etess, iedge, &npnt, &xyz_server, &t_server);
      assert (status == 0);
      
      xyz_client = (double *) malloc(3*npnt*sizeof(double));
      assert (xyz_client != NULL);
      t_client = (double *) malloc(npnt*sizeof(double));
      assert (t_client != NULL);
      
      status = MPI_Recv((void *)(&Npnt),    1,      MPI_INT,    MPI_ANY_SOURCE,
                        300+iedge, myCommWorld, &MPIstat);
      assert (status == 0);
      
      assert (Npnt == npnt);
      
      status = MPI_Recv((void *)xyz_client, 3*npnt, MPI_DOUBLE, MPI_ANY_SOURCE,
                        400+iedge, myCommWorld, &MPIstat);
      assert (status == 0);
      
      status = MPI_Recv((void *)t_client, npnt, MPI_DOUBLE, MPI_ANY_SOURCE,
                        700+iedge, myCommWorld, &MPIstat);
      assert (status == 0);
      
      for (ipnt = 0; ipnt < npnt; ipnt++) {
        if (xyz_server[3*ipnt  ] != xyz_client[3*ipnt  ] ||
            xyz_server[3*ipnt+1] != xyz_client[3*ipnt+1] ||
            xyz_server[3*ipnt+2] != xyz_client[3*ipnt+2]   ) {
          printf("edge %3d %5d  %12.5f %12.5f (%12.4e)\n", iedge, ipnt,
                 xyz_server[3*ipnt  ], xyz_client[3*ipnt  ],
                 xyz_server[3*ipnt  ]- xyz_client[3*ipnt  ]);
          printf("                %12.5f %12.5f (%12.4e)\n",
                 xyz_server[3*ipnt+1], xyz_client[3*ipnt+1],
                 xyz_server[3*ipnt+1]- xyz_client[3*ipnt+1]);
          printf("                %12.5f %12.5f (%12.4e)\n",
                 xyz_server[3*ipnt+2], xyz_client[3*ipnt+2],
                 xyz_server[3*ipnt+2]- xyz_client[3*ipnt+2]);
          
          assert (fabs(xyz_server[3*ipnt  ]-xyz_client[3*ipnt  ]) < 1e-15);
          assert (fabs(xyz_server[3*ipnt+1]-xyz_client[3*ipnt+1]) < 1e-15);
          assert (fabs(xyz_server[3*ipnt+2]-xyz_client[3*ipnt+2]) < 1e-15);
        }
        if (t_server[ipnt] !=  t_client[ipnt]) {
          printf("edge %3d %5d  %12.5f %12.5f (%12.4e)\n", iedge, ipnt,
                 t_server[ipnt], t_client[ipnt], t_server[ipnt]-t_client[ipnt]);
          assert (fabs(t_server[ipnt]-t_client[ipnt]) < 1e-15);
        }
      }
      
      free(t_client);
      free(xyz_client);
    }
  }
  
  if (numRanks == 2) {
    /* receive the Faces from the server */
    int          npnt, ipnt, ntri;
    const int    *ptype, *pindx, *tris, *tric;
    double       *xyz_client;
    const double *xyz_server, *uv_server;
    ego          *faces;
    
    printf("Comparing Faces...\n");
    
    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &faces);
    assert (status == 0);
    
    for (iface = 1; iface <= nface; iface++) {
      status = EG_getTessFace(etess, iface, &npnt, &xyz_server, &uv_server,
                              &ptype, &pindx, &ntri, &tris, &tric);
      assert (status == 0);
      
      status = MPI_Recv((void *)(&Npnt),    1,      MPI_INT,    MPI_ANY_SOURCE,
                        500+iface, myCommWorld, &MPIstat);
      assert (status == 0);
      
      xyz_client = (double *) malloc(3*Npnt*sizeof(double));
      assert (xyz_client != NULL);
      
      if (Npnt != npnt) printf("face %3d:  Npnts = %d %d\n", iface, Npnt, npnt);
//    assert (Npnt == npnt);
      
      status = MPI_Recv((void *)xyz_client, 3*Npnt, MPI_DOUBLE, MPI_ANY_SOURCE,
                        600+iface, myCommWorld, &MPIstat);
      assert (status == 0);
      
/*
      EG_getTopology(faces[iface-1], &eref, &oclass, &mtype,
                     uvbox, &nchild, &children, &senses);
      printf("face %3d: surface mtype = %d\n", iface, eref->mtype);
      
      for (ipnt = 0; ipnt < npnt; ipnt++) {
        if (xyz_server[3*ipnt  ] != xyz_client[3*ipnt  ] ||
            xyz_server[3*ipnt+1] != xyz_client[3*ipnt+1] ||
            xyz_server[3*ipnt+2] != xyz_client[3*ipnt+2]   ) {
          printf("face %3d %5d  %12.5f %12.5f (%12.4e)\n", iface, ipnt,
                 xyz_server[3*ipnt  ], xyz_client[3*ipnt  ],
                 xyz_server[3*ipnt  ]- xyz_client[3*ipnt  ]);
          printf("                %12.5f %12.5f (%12.4e)\n",
                 xyz_server[3*ipnt+1], xyz_client[3*ipnt+1],
                 xyz_server[3*ipnt+1]- xyz_client[3*ipnt+1]);
          printf("                %12.5f %12.5f (%12.4e)",
                 xyz_server[3*ipnt+2], xyz_client[3*ipnt+2],
                 xyz_server[3*ipnt+2]- xyz_client[3*ipnt+2]);
          printf("   %d %d\n", ptype[ipnt], pindx[ipnt]);
          
          assert (fabs(xyz_server[3*ipnt  ]-xyz_client[3*ipnt  ]) < 1e-14);
          assert (fabs(xyz_server[3*ipnt+1]-xyz_client[3*ipnt+1]) < 1e-14);
          assert (fabs(xyz_server[3*ipnt+2]-xyz_client[3*ipnt+2]) < 1e-14);
        }
      }
*/
      free(xyz_client);
    }
    EG_free(faces);
  }
  
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
  printf("volume %10.5f %10.5f (%12.4e)\n", rootProps[0], totProps[0],
                                            rootProps[0]- totProps[0]);
  printf("area   %10.5f %10.5f (%12.4e)\n", rootProps[1], totProps[1],
                                            rootProps[1]- totProps[1]);
  printf("xcg    %10.5f %10.5f (%12.4e)\n", rootProps[2], totProps[2],
                                            rootProps[2]- totProps[2]);
  printf("ycg    %10.5f %10.5f (%12.4e)\n", rootProps[3], totProps[3],
                                            rootProps[3]- totProps[3]);
  printf("zcg    %10.5f %10.5f (%12.4e)\n", rootProps[4], totProps[4],
                                            rootProps[4]- totProps[4]);
  
  /* remove the copy of my world communicator */
  status = MPI_Comm_free(&myCommWorld);
  assert (status == 0);
  
  /* shut down MPI */
  status = MPI_Finalize();
  assert (status == 0);
  
  return 0;
}
