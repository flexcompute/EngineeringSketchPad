/*
 * OpenCSM/EGADS User Defined Primitive include
 *
 */

/*
 * Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
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

#define OCSM_NODE 600
#define OCSM_EDGE 601
#define OCSM_FACE 602

/* initialize & get info about the list of arguments */
extern int
udp_initialize(char   primName[],
               int    *nArgs,
               char   **name[],
               int    *type[],
               int    *idefault[],
               double *ddefault[]);

/* get number of Bodys expected in call to udp_executePrim */
extern int
udp_numBodys(char primName[]);

/* get list of Bodys associated with primitive */
extern int
udp_bodyList(char primName[], ego body, const int *bodyList[]);

/* set the argument list back to default */
extern int
udp_clrArguments(char primName[]);

/* clean the udp cache */
extern int
udp_clean(char primName[],
          ego  body);

/* set an argument -- characters are converted based on type */
extern int
udp_setArgument(char primName[],
                char name[],
                void *value,
                int  nvalue,
      /*@null@*/char message[]);

/* execute */
extern int
udp_executePrim(char primName[],
                ego  context,
                ego  *body,             // destroy with EG_deleteObject
                int  *nMesh,            // =0 for no 3D mesh
                char *string[]);        // freeable

/* get an output parameter -- characters are converted based on type */
extern int
udp_getOutput(char primName[],
              ego  body,
              char name[],
              void *value,
              char message[]);

/* return the OverSet Mesh */
extern int
udp_getMesh(char   primName[],
            ego    body,
            int    iMesh,               // bias-1
            int    *imax,
            int    *jmax,
            int    *kmax,
            double *mesh[]);             // freeable

/* set a design velocity -- characters are converted based on type */
extern int
udp_setVelocity(char   primName[],
                ego    body,
                char   name[],
                double value[],
                int    nvalue);

/* return sensitivity derivatives for the "real/sens" argument */
extern int
udp_sensitivity(char   primName[],
                ego    body,
                int    npts,
                int    entType,
                int    entIndex,
      /*@null@*/double uvs[],
                double vels[]);

/* unload and cleanup for all */
extern void
udp_cleanupAll();

/* additional attribute type */
#define ATTRREALSEN 4
#define ATTRFILE    5

/* define needed for WIN32 */
#ifdef WIN32
    #define snprintf _snprintf
#endif

/* prototypes for routines to be provided with a udp */
int udpInitialize(int    *nArgs,
                  char   **namex[],
                  int    *typex[],
                  int    *idefault[],
                  double *ddefault[]);

int udpReset(int flag);

int udpClean(ego ebody);

int udpSet(char name[],
           void *values,
           int  nvalue,
           char message[]);

int udpGet(ego  ebody,
           char name[],
           void *value,
           char message[]);

int udpVel(ego    ebody,
           char   name[],
           double value[],
           int    nvalue);

int udpExecute(ego  context,
               ego  *ebody,
               int  *nMesh,
               char *string[]);         // freeable

int udpSensitivity(ego    ebody,
                   int    npts,
                   int    entType,
                   int    entIndex,
                   double uvs[],
                   double vels[]);

int udpMesh(ego body,
            int  imesh,
            int  *imax,
            int  *jmax,
            int  *kmax,
            double *mesh[]);            // freeable
