/*
 ************************************************************************
 *                                                                      *
 * udpDumpCATIA -- udp file to dump STEP/IGES for CATIA                 *
 *                                                                      *
 *            Written by Marshall Galbraith @ MIT                       *
 *            Patterned after code written by                           *
 *            John Dannenhoffer @ Syracuse University and               *
 *            Bob Haimes  @ MIT                                         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2025  Marshall Galbraith @ (MIT)
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

#include "egads_dot.h"

extern "C" {
#define NUMUDPARGS 1
#define NUMUDPINPUTBODYS -999
#include "udpUtilities.h"
}

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)  ((char   *) (udps[IUDP].arg[0].val))

/* data about possible arguments */
static const char*  argNames[NUMUDPARGS] = {"filename", };
static       int    argTypes[NUMUDPARGS] = {ATTRSTRING, };
static       int    argIdefs[NUMUDPARGS] = {0,          };
static       double argDdefs[NUMUDPARGS] = {0.,         };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
 udpGet, udpVel, udpClean, udpMesh */
extern "C" {
#include "udpUtilities.c"
}

#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include <algorithm>

/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */
extern "C"
int
udpExecute(ego  emodel,                 /* (in)  Model containing Bodies */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, nbody, nface, nloop, nedge, nnode, *senses, *loop_senses, ndump = 0;
    double  data[4];
    int     atype, alen, wsense[1] = {SFORWARD};
    const int    *pint;
    const double *preal;
    const char   *pstr;
    char    *message=NULL;
    ego     context=NULL, eref, eshell, eloop, *eloop_edges=NULL, *eloops=NULL, *ebodys=NULL, *ebodys_copy=NULL, *edumpBodys=NULL;
    ego     *efaces=NULL, *eedges=NULL, *enodes=NULL, edumpModel=NULL;
    udp_T   *udps = *Udps;

    std::vector<std::pair<std::string,std::vector<ego>>> loops;
    std::vector<std::pair<std::string,std::vector<ego>>> shells;

    std::map<std::string,std::vector<ego>> named_faces;
    std::map<std::string,std::vector<ego>> named_edges;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 1024);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size <= 1) {
      snprintf(message, 1024, "FILENAME must be given");
      status = EGADS_RANGERR;
      goto cleanup;
    }

    if (emodel->oclass != MODEL) {
      snprintf(message, 1024, "Expecting a model!");
      status = EGADS_RANGERR;
      goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nbody, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    /* group topology based on the 'Name' attribute */
    for (int ibody = 0; ibody < nbody; ibody++) {

      if (ebodys[ibody]->mtype == WIREBODY) {
        status = EG_getTopology(ebodys[ibody], &eref, &oclass, &mtype, data,
                                &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getTopology(eloops[0], &eref, &oclass, &mtype, data,
                                &nedge, &eloop_edges, &senses);
        CHECK_STATUS(EG_getTopology);

        for (int iedge = 0; iedge < nedge; iedge++) {

          status = EG_attributeRet(eloop_edges[iedge], "Name", &atype, &alen, &pint, &preal, &pstr);
          if (status == EGADS_NOTFOUND) {
            if (!loops.empty() && loops.back().first == "") {
              loops.back().second.push_back(eloop_edges[iedge]);
            }
            else {
              loops.push_back( std::make_pair<std::string,std::vector<ego>>("", {eloop_edges[iedge]}) );
            }
            continue;
          }
          CHECK_STATUS(EG_attributeRet);

          if (atype != ATTRSTRING) {
            snprintf(message, 1024, "Body %d Edge %d 'Name' attribute is not a string", ibody+1, iedge+1);
            status = EGADS_ATTRERR;
            goto cleanup;
          }

          if (!loops.empty() && loops.back().first == pstr) {
            loops.back().second.push_back(eloop_edges[iedge]);
          } else {
            loops.push_back( std::make_pair<std::string,std::vector<ego>>(pstr, {eloop_edges[iedge]}) );
          }
        }

      } else {

        //===================================================================

        status = EG_getBodyTopos(ebodys[ibody], NULL, EDGE, &nedge, &eedges);
        CHECK_STATUS(EG_getBodyTopos);

        for (int iedge = 0; iedge < nedge; iedge++) {

          status = EG_attributeRet(eedges[iedge], "Name", &atype, &alen, &pint, &preal, &pstr);
          if (status == EGADS_NOTFOUND) {
            //named_edges[""].push_back(eedges[iedge]);
            continue;
          }
          CHECK_STATUS(EG_attributeRet);

          if (atype != ATTRSTRING) {
            snprintf(message, 1024, "Body %d Edge %d 'Name' attribute is not a string", ibody+1, iedge+1);
            status = EGADS_ATTRERR;
            FREE(eedges);
            goto cleanup;
          }

          named_edges[pstr].push_back(eedges[iedge]);
        }

        FREE(eedges);

        for (std::pair<std::string,std::vector<ego>> name_edge : named_edges) {

          std::list<std::set<ego>> edge_neighbors;

          for (ego edge : name_edge.second) {
            if (edge->mtype == DEGENERATE) continue;

            edge_neighbors.resize(edge_neighbors.size()+1);
            edge_neighbors.back().insert(edge);

            status = EG_getTopology(edge, &eref, &oclass, &mtype, data,
                                    &nnode, &enodes, &senses);
            CHECK_STATUS(EG_getTopology);

            for (int inode = 0; inode < nnode; inode++) {

              status = EG_getBodyTopos(ebodys[ibody], enodes[inode], EDGE, &nedge, &eedges);
              CHECK_STATUS(EG_getBodyTopos);

              for (int iedge = 0; iedge < nedge; iedge++) {
                if (eedges[iedge]->mtype == DEGENERATE) continue;
                if (eedges[iedge] == edge) continue;
                if (std::find(name_edge.second.begin(), name_edge.second.end(), eedges[iedge]) == name_edge.second.end()) continue;

                edge_neighbors.back().insert(eedges[iedge]);
              }

              FREE(eedges);
            }
          }

          for (std::list<std::set<ego>>::iterator it = edge_neighbors.begin(); it != edge_neighbors.end(); ) {

            bool found = false;
            for (std::list<std::set<ego>>::iterator suburb = it; suburb != edge_neighbors.end(); suburb++) {
              if (it == suburb) continue;

              for (ego edge : *suburb) {
                if (it->find(edge) != it->end()) {
                  it->insert(suburb->begin(), suburb->end());
                  edge_neighbors.erase(suburb);
                  suburb = it;
                  found = true;
                  break;
                }
              }
            }
            if (!found) it++;
          }

          for (const std::set<ego>& suburb : edge_neighbors ) {

            // custom sorting function object for sorting based on nodes
            struct
            {
              bool operator()(const ego& L, const ego& R) const
              {
                int nnode, oclass, mtype, *senses;
                double data[3];
                ego eref, *enodesL, *enodesR;
                EG_getTopology(L, &eref, &oclass, &mtype, data, &nnode, &enodesL, &senses);
                EG_getTopology(R, &eref, &oclass, &mtype, data, &nnode, &enodesR, &senses);

                return enodesL[0] == enodesR[0] ||
                       enodesL[0] == enodesR[1] ||
                       enodesL[1] == enodesR[0] ||
                       enodesL[1] == enodesR[1];
              }
            } nodeMatch;

            std::vector<ego> edges(suburb.begin(), suburb.end());
            std::vector<ego> sorted_edges;
            sorted_edges.push_back(edges.back());
            edges.pop_back();
            while (!edges.empty()) {
              for (std::size_t i = 0; i < sorted_edges.size(); i++)
                if (nodeMatch(edges.back(), sorted_edges[i])) {
                  sorted_edges.push_back(edges.back());
                  edges.pop_back();
                  if (edges.empty()) break;
                  i = 0;
                }
              if (!edges.empty()) std::rotate(edges.begin(), edges.begin()+1, edges.end());
            }

            std::pair<std::string,std::vector<ego>> pair(name_edge.first, sorted_edges);
            loops.emplace_back( pair );
          }
        }

        named_edges.clear();

        //===================================================================

        // Gather shells with common "Name"
        status = EG_getBodyTopos(ebodys[ibody], NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        for (int iface = 0; iface < nface; iface++) {

          status = EG_attributeRet(efaces[iface], "Name", &atype, &alen, &pint, &preal, &pstr);
          if (status == EGADS_NOTFOUND) {
            named_faces[""].push_back(efaces[iface]);
            continue;
          }
          CHECK_STATUS(EG_attributeRet);

          if (atype != ATTRSTRING) {
            snprintf(message, 1024, "Body %d Face %d 'Name' attribute is not a string", ibody+1, iface+1);
            status = EGADS_ATTRERR;
            FREE(efaces);
            goto cleanup;
          }

          named_faces[pstr].push_back(efaces[iface]);
        }

        FREE(efaces);

        for (std::pair<std::string,std::vector<ego>> name_face : named_faces) {

          std::list<std::set<ego>> face_neighbors;

          for (ego face : name_face.second) {

            face_neighbors.resize(face_neighbors.size()+1);
            face_neighbors.back().insert(face);

            status = EG_getTopology(face, &eref, &oclass, &mtype, data,
                                    &nloop, &eloops, &loop_senses);
            CHECK_STATUS(EG_getTopology);

            for (int iloop = 0; iloop < nloop; iloop++) {
              status = EG_getTopology(eloops[iloop], &eref, &oclass, &mtype, data,
                                      &nedge, &eloop_edges, &senses);
              CHECK_STATUS(EG_getTopology);

              for (int iedge = 0; iedge < nedge; iedge++) {

                status = EG_getBodyTopos(ebodys[ibody], eloop_edges[iedge], FACE, &nface, &efaces);
                CHECK_STATUS(EG_getBodyTopos);

                for (int iface = 0; iface < nface; iface++) {
                  if (efaces[iface] == face) continue;
                  if (std::find(name_face.second.begin(), name_face.second.end(), efaces[iface]) == name_face.second.end()) continue;

                  face_neighbors.back().insert(efaces[iface]);
                }

                FREE(efaces);
              }
            }
          }

          for (std::list<std::set<ego>>::iterator it = face_neighbors.begin(); it != face_neighbors.end(); ) {

            bool found = false;
            for (std::list<std::set<ego>>::iterator suburb = it; suburb != face_neighbors.end(); suburb++) {
              if (it == suburb) continue;

              for (ego face : *suburb) {
                if (it->find(face) != it->end()) {
                  it->insert(suburb->begin(), suburb->end());
                  face_neighbors.erase(suburb);
                  suburb = it;
                  found = true;
                  break;
                }
              }
            }
            if (!found) it++;
          }

          for (std::set<ego> suburb : face_neighbors ) {
            std::vector<ego> faces(suburb.begin(), suburb.end());
            std::pair<std::string,std::vector<ego>> pair(name_face.first, faces);
            shells.emplace_back( pair );
          }
        }

        named_faces.clear();
      }
    }


    ndump = shells.size() + loops.size();
    MALLOC(edumpBodys, ego, ndump);
    for (int ibody = 0; ibody < ndump; ibody++) edumpBodys[ibody] = NULL;
    ndump = 0;


    for (std::size_t iwire = 0; iwire < loops.size(); iwire++) {
      std::pair<std::string,std::vector<ego>>& loop = loops[iwire];

      if (loop.second.size() == 1) {
        // Single Edge, use makeTopology
        status = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, loop.second.size(), loop.second.data(), wsense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeAdd(eloop, "Name", ATTRSTRING, 1, NULL, NULL, loop.first.c_str());
        CHECK_STATUS(EG_attributeAdd);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, &edumpBodys[ndump++]);
        EG_deleteObject(eloop);
        CHECK_STATUS(EG_makeTopology);

      } else {
        // Multiple Edge, use makeNmWireBody
        status = EG_makeNmWireBody(loop.second.size(), loop.second.data(), 0, &edumpBodys[ndump]);
        CHECK_STATUS(EG_makeNmWireBody);

        status = EG_getTopology(edumpBodys[ndump++], &eref, &oclass, &mtype, data,
                                &nloop, &eloops, &loop_senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_attributeAdd(eloops[0], "Name", ATTRSTRING, 1, NULL, NULL, loop.first.c_str());
        CHECK_STATUS(EG_attributeAdd);
      }
    }

    for (std::pair<std::string,std::vector<ego>> shell : shells) {
      status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL, shell.second.size(), shell.second.data(), NULL, &eshell);
      CHECK_STATUS(EG_makeTopology);

      if (shell.first != "") {
        status = EG_attributeAdd(eshell, "Name", ATTRSTRING, 1, NULL, NULL, shell.first.c_str());
        CHECK_STATUS(EG_attributeAdd);
      }

      status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1, &eshell, NULL, &edumpBodys[ndump++]);
      EG_deleteObject(eshell);
      CHECK_STATUS(EG_makeTopology);
    }

    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, ndump, edumpBodys, NULL, &edumpModel);
    CHECK_STATUS(EG_makeTopology);

    remove(FILENAME(numUdp));
    status = EG_saveModel(edumpModel, FILENAME(numUdp));
    CHECK_STATUS(EG_saveModel);

    udps[numUdp].ebody = *ebody = edumpModel;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:

  if (strlen(message) > 0) {
    *string = message;
    printf("%s\n", message);
  } else if (status != EGADS_SUCCESS) {
    FREE(message);
    *string = udpErrorStr(status);
  } else {
    FREE(message);
  }

  /* cleanup temporary ego's */
  EG_deleteObject(context);

  FREE(ebodys_copy);
  FREE(efaces);
  FREE(eedges);
  FREE(edumpBodys);

  return status;
}


/*
 *****************************************************************************
 *                                                                           *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                           *
 *****************************************************************************
 */

extern "C"
int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
  int    status = EGADS_SUCCESS;

  int    iudp, judp;

#ifdef DEBUG
  printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
         (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
#endif

  /* check that ebody matches one of the ebodys */
  iudp = 0;
  for (judp = 1; judp <= numUdp; judp++) {
    if (ebody == udps[judp].ebody) {
      iudp = judp;
      break;
    }
  }
  if (iudp <= 0) {
    return EGADS_NOTMODEL;
  }

  return status;
}

