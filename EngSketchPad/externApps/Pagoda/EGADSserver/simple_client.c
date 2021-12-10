/*
 ************************************************************************
 *                                                                      *
 * simple_client -- demonstration of EGADS and distributed EGADS-Lite   *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *                       Bob Haimes @ MIT                               *
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
/*   main program - client                                             */
/*                                                                     */
/***********************************************************************/

int
main(int    argc,                       /* (in)  number of arguments */
     char   *argv[])                    /* (in)  array  of arguments */
{
  int          status, myRank, numRanks;
  int          nbytes, oclass, mtype, nchild ,*senses, nface, iface, atype, alen;
  const int    *tempIlist;
  double       data[18], myProps[5], totProps[5];
  const double *tempRlist;
  char         *stream;
  const char   *tempClist;
  ego          context, emodel, eref, *ebodys, etess;
  MPI_Comm     myCommWorld;
  
  /* --------------------------------------------------------------- */
  
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
  
  /* this (client) process should not be rank 0 */
  assert (myRank != 0);
  
  /* determine the number of ranks */
  status = MPI_Comm_size(myCommWorld, &numRanks);
  assert (status == 0);
  
  /* receive (via broadcast) the size of the stream */
  status = MPI_Bcast((void *)(&nbytes), 1, MPI_INT, 0, myCommWorld);
  assert (status == 0);
  
  /* allocate a buffer to receive the stream */
  stream = (char *) EG_alloc(nbytes+1);
  assert (stream != NULL);
  
  /* receive (via broadcast) the stream */
  status = MPI_Bcast((void *)stream, nbytes, MPI_CHAR, 0, myCommWorld);
  assert (status == 0);
  
  /* import the stream into EGADSlite */
  status = EG_open(&context);
  assert (status == 0);
  
  status = EG_importModel(context, nbytes, stream, &emodel);
  assert (status == 0);
  EG_free(stream);
  
  /* extract the Body from the Model */
  status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                          data, &nchild, &ebodys, &senses);
  assert (status == 0);
  
  /* tessellate the Body with the tessellation parameters from _tParams */
  status = EG_attributeRet(ebodys[0], "_tParams", &atype, &alen,
                           &tempIlist, &tempRlist, &tempClist);
  assert (status == 0);
  
  status = EG_makeTessBody(ebodys[0], (double *) tempRlist, &etess);
  assert (status == 0);
  
  if (numRanks == 2) {
    /* make a list of Nodes and send to server */
    double *xyz;
    int    nnode, inode;
    ego    *enodes, *echilds;
    
    status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, &enodes);
    assert (status == 0);
    
    xyz = (double *) malloc(3*nnode*sizeof(double));
    assert (xyz != NULL);
    
    for (inode = 0; inode < nnode; inode++) {
      status = EG_getTopology(enodes[inode], &eref, &oclass, &mtype,
                              data, &nchild, &echilds, &senses);
      assert (status == 0);
      
      xyz[3*inode  ] = data[0];
      xyz[3*inode+1] = data[1];
      xyz[3*inode+2] = data[2];
    }
    
    status = MPI_Send((void *)(&nnode), 1,       MPI_INT,    0,
                      100, myCommWorld);
    assert (status == 0);
    
    status = MPI_Send((void *)xyz,      3*nnode, MPI_DOUBLE, 0,
                      200, myCommWorld);
    assert (status == 0);
    
    EG_free(enodes);
    free(xyz);
  }
  
  if (numRanks == 2) {
    /* send each of the Edges to the server */
    int          nedge, iedge, npnt;
    const double *xyz, *t;
    
    status = EG_getBodyTopos(ebodys[0], NULL, EDGE, &nedge, NULL);
    assert (status == 0);
    
    for (iedge = 1; iedge <= nedge; iedge++) {
      status = EG_getTessEdge(etess, iedge,
                              &npnt, &xyz, &t);
      assert (status == 0);
      
      status = MPI_Send((void *)(&npnt), 1,      MPI_INT,    0,
                        300+iedge, myCommWorld);
      assert (status == 0);
      
      status = MPI_Send((void *)xyz,     3*npnt, MPI_DOUBLE, 0,
                        400+iedge, myCommWorld);
      assert (status == 0);
      status = MPI_Send((void *)t,         npnt, MPI_DOUBLE, 0,
                        700+iedge, myCommWorld);
      assert (status == 0);
    }
  }
  
  if (numRanks == 2) {
    /* send each of the Faces to the server */
    int          npnt, ntri;
    const int    *ptype, *pindx, *tris, *tric;
    const double *xyz, *uv;
    
    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, NULL);
    assert (status == 0);
    
    for (iface = 1; iface <= nface; iface++) {
      status = EG_getTessFace(etess, iface,
                              &npnt, &xyz, &uv, &ptype, &pindx,
                              &ntri, &tris, &tric);
      assert (status == 0);
      
      status = MPI_Send((void *)(&npnt), 1,      MPI_INT,    0,
                        500+iface, myCommWorld);
      assert (status == 0);
      
      status = MPI_Send((void *)xyz,     3*npnt, MPI_DOUBLE, 0,
                        600+iface, myCommWorld);
      assert (status == 0);
    }
  }
  
  /* compute the volume intergals on my portion of the configuration */
  myProps[0] = 0;
  myProps[1] = 0;
  myProps[2] = 0;
  myProps[3] = 0;
  myProps[4] = 0;
  
  status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, NULL);
  assert (status == 0);
  
  for (iface = myRank; iface <= nface; iface += (numRanks-1)) {
    status = massProps(etess, iface, myProps);
    assert (status == 0);
  }
  
  /* send my portion of the integrals back to the server */
  status = MPI_Reduce(myProps, totProps, 5, MPI_DOUBLE,
                      MPI_SUM, 0, myCommWorld);
  assert (status == 0);
  
  /* remove the copy of my world communicator */
  status = MPI_Comm_free(&myCommWorld);
  assert (status == 0);
  
  /* shut down MPI */
  status = MPI_Finalize();
  assert (status == 0);
  
  return 0;
}
