/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             vtk 3D Mesh Writer Code
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include "aimUtil.h"
#include "aimMesh.h"

#include "vtkWriter.h"

static int surfaceMeshWrite(void *aimInfo, aimMesh *mesh);
static int generalMeshWrite(void *aimInfo, aimMesh *mesh);

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

const char *meshExtension()
{
/*@-observertrans@*/
  return MESHEXTENSION;
/*@+observertrans@*/
}

/* =================================================================== */

int meshWrite(void *aimInfo, aimMesh *mesh)
{
    int  status; // Function return status
    int  nLine=0, nTri=0, nQuad=0;
    int  nTet=0, nPyramid=0, nPrism=0, nHex=0;
    int  igroup;
    aimMeshData *meshData = NULL;

    /* ---------------------------------------- */

    printf("\nWriting vtk file ....\n");

    if (mesh == NULL) return CAPS_NULLVALUE;
    if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
    if (mesh->meshData == NULL) return CAPS_NULLVALUE;

    meshData = mesh->meshData;

    if (meshData->dim != 2 && meshData->dim != 3) {
        AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", meshData->dim);
        return CAPS_BADVALUE;
    }

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
        if (meshData->elemGroups[igroup].order != 1) {
            AIM_ERROR(aimInfo, "vtk files only supports linear mesh elements! group %d order = %d",
                      igroup, meshData->elemGroups[igroup].order);
            status = CAPS_IOERR;
            goto cleanup;
        }

        // count the number of element types
        if        (meshData->elemGroups[igroup].elementTopo == aimLine   ) {
            nLine    += meshData->elemGroups[igroup].nElems;
        } else if (meshData->elemGroups[igroup].elementTopo == aimTri    ) {
            nTri     += meshData->elemGroups[igroup].nElems;
        } else if (meshData->elemGroups[igroup].elementTopo == aimQuad   ) {
            nQuad    += meshData->elemGroups[igroup].nElems;
        } else if (meshData->elemGroups[igroup].elementTopo == aimTet    ) {
            nTet     += meshData->elemGroups[igroup].nElems;
        } else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) {
            nPyramid += meshData->elemGroups[igroup].nElems;
        } else if (meshData->elemGroups[igroup].elementTopo == aimPrism  ) {
            nPrism   += meshData->elemGroups[igroup].nElems;
        } else if (meshData->elemGroups[igroup].elementTopo == aimHex    ) {
            nHex     += meshData->elemGroups[igroup].nElems;
        } else {
            AIM_ERROR(aimInfo, "Unknown group %d element topology: %d", igroup+1, meshData->elemGroups[igroup].elementTopo);
            status = CAPS_MISMATCH;
            goto cleanup;
        }
    }

    /* create surface vtk file if there are only Triangles */
    if (nLine == 0 && nTri > 0 && nQuad == 0 && nTet == 0 && nPyramid == 0 && nPrism == 0 && nHex == 0) {
        status = surfaceMeshWrite(aimInfo, mesh);

    /* otherwise create general vtk file */
    } else {
        status = generalMeshWrite(aimInfo, mesh);
    }

cleanup:
    return status;
}

/* =================================================================== */

static int surfaceMeshWrite(void *aimInfo, aimMesh *mesh)
{
    int  status; // Function return status

    int         nTri=0;
    int         ivrtx, igroup, jgroup;

    int         ntess, itess, ngroup, ntemp, nedge, iedge, nface, iface, jface;
    int         npnt, npnt2, ntri, ntri2, itri, jtri;
    int         stat, attrType, attrLen, found;
    int         it0, it1, it2, ig0, ig1, ig2, IG0, IG1, IG2, tri[3];
    int         *edgeLeft=NULL, *edgeRite=NULL, *faceOffset=NULL;
    CINT        *ptype, *ptype2, *pindx, *pindx2, *tris, *tris2, *tric, *tric2;
    CINT        *tempIlist;
    CDOUBLE     *xyz, *xyz2, *uv, *uv2, *tempRlist;
    char        filename[PATH_MAX];
    char        **groups=NULL;
    CCHAR       *tempClist;
    ego         etess, ebody, *eedges=NULL, *efaces=NULL, *etemps=NULL;
    FILE        *fp=NULL;
    aimMeshData *meshData = NULL;

    /* ---------------------------------------- */

    printf("\nWriting VTK surface file ....\n");

    if (mesh == NULL) return CAPS_NULLVALUE;
    if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
    if (mesh->meshData == NULL) return CAPS_NULLVALUE;

    meshData = mesh->meshData;
    ntess = mesh->meshRef->nmap;

    /* determine the number of capsGroups found */
    ngroup = 0;

    for (itess = 0; itess < ntess; itess++) {
        etess = mesh->meshRef->maps[itess].tess;

        status = EG_statusTessBody(etess, &ebody, &stat, &npnt);
        AIM_STATUS(aimInfo, status, "EG_statusTessBody");

        if (efaces != NULL) {printf("line %d\n", __LINE__); exit(EXIT_FAILURE);}
        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        AIM_STATUS(aimInfo, status, "EG_getBodyTopos");

        for (iface = 0; iface < nface; iface++) {
            status = EG_attributeRet(efaces[iface], "capsGroup",
                                     &attrType, &attrLen, &tempIlist, &tempRlist, &tempClist);
            AIM_STATUS(aimInfo, status, "EG_attributeRet");

            if (status == EGADS_SUCCESS && attrType == ATTRSTRING) {
                found = 0;
                for (igroup = 0; igroup < ngroup; igroup++) {
                    if (strcmp(tempClist, groups[igroup]) == 0) {
                        found = 1;
                        break;
                    }
                }

                if (found == 0) {
                    ngroup++;
                    AIM_REALL(groups, ngroup, char*, aimInfo, status);
                    groups[ngroup-1] = NULL;

                    AIM_ALLOC(groups[ngroup-1], strlen(tempClist)+1, char, aimInfo, status);

                    strcpy(groups[ngroup-1], tempClist);
                }
            }
        }

        AIM_FREE(efaces);
    }

    printf("\ncapsGroups found:\n");
    for (igroup = 0; igroup < ngroup; igroup++) {
        printf("%3d: %s\n", igroup+1, groups[igroup]);
    }

    if (meshData->dim != 2 && meshData->dim != 3) {
        AIM_ERROR(aimInfo, "meshData dim = %d must be 2 or 3!!!", meshData->dim);
        return CAPS_BADVALUE;
    }

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
        if (meshData->elemGroups[igroup].order != 1) {
            AIM_ERROR(aimInfo, "vtk files only supports linear mesh elements! group %d order = %d",
                      igroup, meshData->elemGroups[igroup].order);
            status = CAPS_IOERR;
            goto cleanup;
        }

        // count the number of element types
        nTri += meshData->elemGroups[igroup].nElems;
    }

    /* open the file */
    (void) snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

    fp = fopen(filename, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
        return CAPS_IOERR;
    }

    /* write the header */
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    for (igroup = 0; igroup < ngroup; igroup++) {
        if (igroup > 0) {
            fprintf(fp, ", ");
        }
        fprintf(fp, "%d=\"%s\"", igroup+1, groups[igroup]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");

    /* write all of the verticies */
    fprintf(fp, "POINTS %d double\n", meshData->nVertex);
    if (meshData->dim == 3) {
        for (ivrtx = 0; ivrtx < meshData->nVertex; ivrtx++) {
            fprintf(fp, "%16.8f %16.8f %16.8f\n",
                    meshData->verts[ivrtx][0],
                    meshData->verts[ivrtx][1],
                    meshData->verts[ivrtx][2]);
        }
    } else {
        for (ivrtx = 0; ivrtx < meshData->nVertex; ivrtx++) {
            fprintf(fp, "%16.8f %16.8f %16.8f\n",
                    meshData->verts[ivrtx][0],
                    meshData->verts[ivrtx][1],
                    0.0);
        }
    }

    /* write out each of the polygons */
    fprintf(fp, "POLYGONS  %d  %d\n", nTri, 4*nTri);

    for (itess = 0; itess < ntess; itess++) {
        etess = mesh->meshRef->maps[itess].tess;

        /* find the number of Edges and Faces */
        status = EG_statusTessBody(etess, &ebody, &stat, &npnt);
        AIM_STATUS(aimInfo, status, "EG_statusTessBody");

        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, NULL);
        AIM_STATUS(aimInfo, status, "EG_getBodyTopos");

        for (iface = 1; iface <= nface; iface++) {
            status = EG_getTessFace(etess, iface,
                                    &npnt, &xyz, &uv, &ptype, &pindx,
                                    &ntri, &tris, &tric);
            AIM_STATUS(aimInfo, status, "EG_getTessFace");

            for (itri = 0; itri < ntri; itri++) {

                for (ivrtx = 0; ivrtx < 3; ivrtx++) {
                    status = EG_localToGlobal(etess, iface, tris[3*itri+ivrtx], &tri[ivrtx]);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");
                }

                fprintf(fp, "%d %8d %8d %8d\n", 3, tri[0]-1, tri[1]-1, tri[2]-1);
            }
        }
    }

    /* write out the body number and connectivity for each triangle */
    fprintf(fp, "CELL_DATA  %d\n", nTri);
    fprintf(fp, "SCALARS  group/connectivity  int  4\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    ntess = mesh->meshRef->nmap;
    for (itess = 0; itess < ntess; itess++) {
        etess = mesh->meshRef->maps[itess].tess;

        /* find the number of Edges and Faces */
        status = EG_statusTessBody(etess, &ebody, &stat, &npnt);
        AIM_STATUS(aimInfo, status, "EG_statusTessBody");

        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
        AIM_STATUS(aimInfo, status, "EG_getBodyTopos");
        AIM_NOTNULL(eedges, aimInfo, status);

        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        AIM_STATUS(aimInfo, status, "EG_getBodyTopos");
        AIM_NOTNULL(efaces, aimInfo, status);

        AIM_ALLOC(edgeLeft, nedge+1, int, aimInfo, status);
        AIM_ALLOC(edgeRite, nedge+1, int, aimInfo, status);

        AIM_ALLOC(faceOffset, nface+2, int, aimInfo, status);

        /* set up Edge and Face info */
        for (iedge = 1; iedge <= nedge; iedge++) {
            edgeLeft[iedge] = 0;
            edgeRite[iedge] = 0;

            status = EG_getBodyTopos(ebody, eedges[iedge-1], FACE, &ntemp, &etemps);
            AIM_STATUS(aimInfo, status, "EG_getBodyTopos");
            AIM_NOTNULL(etemps, aimInfo, status);

            if (ntemp == 2) {
                status = edgeLeft[iedge] = EG_indexBodyTopo(ebody, etemps[0]);
                if (status > CAPS_SUCCESS) status = CAPS_SUCCESS;
                AIM_STATUS(aimInfo, status, "EG_indexBodyTopo");

                status = edgeRite[iedge] = EG_indexBodyTopo(ebody, etemps[1]);
                if (status > CAPS_SUCCESS) status = CAPS_SUCCESS;
                AIM_STATUS(aimInfo, status, "EG_indexBodyTopo");
            }

            AIM_FREE(etemps);
        }

        faceOffset[1] = 0;
        for (iface = 1; iface <= nface; iface++) {
            status = EG_getTessFace(etess, iface,
                                    &npnt, &xyz, &uv, &ptype, &pindx,
                                    &ntri, &tris, &tric);
            AIM_STATUS(aimInfo, status, "EG_getTessFace");

            faceOffset[iface+1] = faceOffset[iface] + ntri;
        }

        /* process all Triangles in all Faces */
        for (iface = 1; iface <= nface; iface++) {

            /* get the group (or 0 if no capsGroup attribute) */
            igroup = 0;

            status = EG_attributeRet(efaces[iface-1], "capsGroup",
                                     &attrType, &attrLen, &tempIlist, &tempRlist, &tempClist);

            if (status == EGADS_SUCCESS && attrType == ATTRSTRING) {
                for (jgroup = 0; jgroup < ngroup; jgroup++) {
                    if (strcmp(tempClist, groups[jgroup]) == 0) {
                        igroup = jgroup + 1;
                        break;
                    }
                }
            }

            /* get the tessellation */
            status = EG_getTessFace(etess, iface,
                                    &npnt, &xyz, &uv, &ptype, &pindx,
                                    &ntri, &tris, &tric);
            AIM_STATUS(aimInfo, status, "EG_getTessFace");

            for (itri = 0; itri < ntri; itri++) {

                /* find neighbor it0 */
                if (tric[3*itri  ] > 0) {
                    it0 = tric[3*itri  ] + faceOffset[iface];
                } else if (tric[3*itri  ] < 0 && edgeLeft[-tric[3*itri  ]] == 0 &&
                                                 edgeRite[-tric[3*itri  ]] == 0   ) {
                    it0 = 0;            // non-manifold
                } else {
                    it0 = -1;

                    status = EG_localToGlobal(etess, iface, tris[3*itri+1], &ig1);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                    status = EG_localToGlobal(etess, iface, tris[3*itri+2], &ig2);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                    iedge = -tric[3*itri  ];
                    if        ( edgeLeft[iedge] == iface) {
                        jface = edgeRite[iedge];
                    } else if ( edgeRite[iedge] == iface) {
                        jface = edgeLeft[iedge];
                    } else {
                        printf("WE SHOULD NOT GET HERE %d: iface=%d, itri=%d, iedge=%d\n", __LINE__, iface, itri, iedge);
                        printf("edgeLeft=%d, edgeRite=%d\n", edgeLeft[iedge], edgeRite[iedge]);
                        printf("tris=%d %d %d, tric=%d %d %d\n",
                               tris[3*itri  ], tris[3*itri  ], tris[3*itri  ],
                               tric[3*itri  ], tric[3*itri  ], tric[3*itri  ]);
                        exit(0);
                    }

                    if (jface != 0) {
                        status = EG_getTessFace(etess, jface,
                                                &npnt2, &xyz2, &uv2, &ptype2, &pindx2,
                                                &ntri2, &tris2, &tric2);
                        AIM_STATUS(aimInfo, status, "EG_getTessFace");

                        for (jtri = 0; jtri < ntri2; jtri++) {
                            if (tric2[3*jtri] == -iedge || tric2[3*jtri+1] == -iedge || tric2[3*jtri+2] == -iedge) {
                                status = EG_localToGlobal(etess, jface, tris2[3*jtri  ], &IG0);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                status = EG_localToGlobal(etess, jface, tris2[3*jtri+1], &IG1);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                status = EG_localToGlobal(etess, jface, tris2[3*jtri+2], &IG2);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                if ((ig1 == IG2 && ig2 == IG1) ||
                                    (ig1 == IG1 && ig2 == IG0) ||
                                    (ig1 == IG0 && ig2 == IG2)   ) {
                                    it0 = (jtri+1) + faceOffset[jface];
                                    break;
                                }
                            }
                        }
                    } else {
                        it0 = 0;
                    }
                }

                /* find neighbor it1 */
                if (tric[3*itri+1] > 0) {
                    it1 = tric[3*itri+1] + faceOffset[iface];
                } else if (tric[3*itri+1] < 0 && edgeLeft[-tric[3*itri+1]] == 0 &&
                                                 edgeRite[-tric[3*itri+1]] == 0   ) {
                    it1 = 0;            // non-manifold
                } else {
                    it1 = -1;

                    status = EG_localToGlobal(etess, iface, tris[3*itri+2], &ig2);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                    status = EG_localToGlobal(etess, iface, tris[3*itri  ], &ig0);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                    iedge = -tric[3*itri+1];
                    if        ( edgeLeft[iedge] == iface) {
                        jface = edgeRite[iedge];
                    } else if ( edgeRite[iedge] == iface) {
                        jface = edgeLeft[iedge];
                    } else {
                        printf("WE SHOULD NOT GET HERE %d: iface=%d, itri=%d, iedge=%d\n", __LINE__, iface, itri, iedge);
                        printf("edgeLeft=%d, edgeRite=%d\n", edgeLeft[iedge], edgeRite[iedge]);
                        printf("tris=%d %d %d, tric=%d %d %d\n",
                               tris[3*itri+1], tris[3*itri+1], tris[3*itri+1],
                               tric[3*itri+1], tric[3*itri+1], tric[3*itri+1]);
                        exit(0);
                    }

                    if (jface != 0) {
                        status = EG_getTessFace(etess, jface,
                                                &npnt2, &xyz2, &uv2, &ptype2, &pindx2,
                                                &ntri2, &tris2, &tric2);
                        AIM_STATUS(aimInfo, status, "EG_getTessFace");

                        for (jtri = 0; jtri < ntri2; jtri++) {
                            if (tric2[3*jtri] == -iedge || tric2[3*jtri+1] == -iedge || tric2[3*jtri+2] == -iedge) {
                                status = EG_localToGlobal(etess, jface, tris2[3*jtri  ], &IG0);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                status = EG_localToGlobal(etess, jface, tris2[3*jtri+1], &IG1);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                status = EG_localToGlobal(etess, jface, tris2[3*jtri+2], &IG2);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                if ((ig2 == IG2 && ig0 == IG1) ||
                                    (ig2 == IG1 && ig0 == IG0) ||
                                    (ig2 == IG0 && ig0 == IG2)   ) {
                                    it1 = (jtri+1) + faceOffset[jface];
                                    break;
                                }
                            }
                        }
                    } else {
                        it1 = 0;
                    }
                }

                /* find neighbor it2 */
                if (tric[3*itri+2] > 0) {
                    it2 = tric[3*itri+2] + faceOffset[iface];
                } else if (tric[3*itri+2] < 0 && edgeLeft[-tric[3*itri+2]] == 0 &&
                                                 edgeRite[-tric[3*itri+2]] == 0   ) {
                    it2 = 0;            // non-manifold
                } else {
                    it2 = -1;

                    status = EG_localToGlobal(etess, iface, tris[3*itri  ], &ig0);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                    status = EG_localToGlobal(etess, iface, tris[3*itri+1], &ig1);
                    AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                    iedge = -tric[3*itri+2];
                    if        ( edgeLeft[iedge] == iface) {
                        jface = edgeRite[iedge];
                    } else if ( edgeRite[iedge] == iface) {
                        jface = edgeLeft[iedge];
                    } else {
                        printf("WE SHOULD NOT GET HERE %d: iface=%d, itri=%d, iedge=%d\n", __LINE__, iface, itri, iedge);
                        printf("edgeLeft=%d, edgeRite=%d\n", edgeLeft[iedge], edgeRite[iedge]);
                        printf("tris=%d %d %d, tric=%d %d %d\n",
                               tris[3*itri+2], tris[3*itri+2], tris[3*itri+2],
                               tric[3*itri+2], tric[3*itri+2], tric[3*itri+2]);
                        exit(0);
                    }

                    if (jface != 0) {
                        status = EG_getTessFace(etess, jface,
                                                &npnt2, &xyz2, &uv2, &ptype2, &pindx2,
                                                &ntri2, &tris2, &tric2);
                        AIM_STATUS(aimInfo, status, "EG_getTessFace");

                        for (jtri = 0; jtri < ntri2; jtri++) {
                            if (tric2[3*jtri] == -iedge || tric2[3*jtri+1] == -iedge || tric2[3*jtri+2] == -iedge) {
                                status = EG_localToGlobal(etess, jface, tris2[3*jtri  ], &IG0);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                status = EG_localToGlobal(etess, jface, tris2[3*jtri+1], &IG1);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                status = EG_localToGlobal(etess, jface, tris2[3*jtri+2], &IG2);
                                AIM_STATUS(aimInfo, status, "EG_localToGlobal");

                                if ((ig0 == IG2 && ig1 == IG1) ||
                                    (ig0 == IG1 && ig1 == IG0) ||
                                    (ig0 == IG0 && ig1 == IG2)   ) {
                                    it2 = (jtri+1) + faceOffset[jface];
                                    break;
                                }
                            }
                        }
                    } else {
                        it2 = 0;
                    }
                }

                /* add to file */
                fprintf(fp, "%8d %8d %8d %8d\n", igroup, it0-1, it1-1, it2-1);
            }
        }

        AIM_FREE(eedges);
        AIM_FREE(efaces);

        AIM_FREE(edgeLeft);
        AIM_FREE(edgeRite);
        AIM_FREE(faceOffset);
    }

    printf("Finished writing vtk file\n\n");

    status = CAPS_SUCCESS;

cleanup:
    for (igroup = 0; igroup < ngroup; igroup++) {
        AIM_FREE(groups[igroup]);
    }
    AIM_FREE(groups);

    AIM_FREE(eedges);
    AIM_FREE(efaces);
    AIM_FREE(etemps);

    AIM_FREE(edgeLeft);
    AIM_FREE(edgeRite);
    AIM_FREE(faceOffset);

    if (fp != NULL) fclose(fp);

    return status;
}

/* =================================================================== */

static int generalMeshWrite(void *aimInfo, aimMesh *mesh)
{
    int  status; // Function return status

    int         i, nPoint, nElems, elemType;
    int         nPoly=0, nTable=0;
    int         ivrtx, ipnt, igroup, ielem;
    char        filename[PATH_MAX];
    FILE        *fp=NULL;
    aimMeshData *meshData = NULL;

    /* ---------------------------------------- */

    printf("\nWriting VTK general file ....\n");

    if (mesh == NULL) return CAPS_NULLVALUE;
    if (mesh->meshRef  == NULL) return CAPS_NULLVALUE;
    if (mesh->meshData == NULL) return CAPS_NULLVALUE;

    meshData = mesh->meshData;

    for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
        if (meshData->elemGroups[igroup].order != 1) {
            AIM_ERROR(aimInfo, "VTK files only supports linear mesh elements! group %d order = %d",
                      igroup, meshData->elemGroups[igroup].order);
            status = CAPS_IOERR;
            goto cleanup;
        }

        // count the number of elements and connectivity table size
        nPoint = meshData->elemGroups[igroup].nPoint;
        nElems = meshData->elemGroups[igroup].nElems;

        nPoly += nElems;
        nTable += nElems*(nPoint+1);
    }

    /* open the file */
    (void) snprintf(filename, PATH_MAX, "%s%s", mesh->meshRef->fileName, MESHEXTENSION);

    fp = fopen(filename, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
        return CAPS_IOERR;
    }

    /* write the header */
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "written by CAPS\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    /* write all of the verticies */
    fprintf(fp, "POINTS %d double\n", meshData->nVertex);
    for (ivrtx = 0; ivrtx < meshData->nVertex; ivrtx++) {
        fprintf(fp, "%16.8f %16.8f %16.8f\n",
                meshData->verts[ivrtx][0],
                meshData->verts[ivrtx][1],
                meshData->verts[ivrtx][2]);
    }

    /* write out each of the cells */
    fprintf(fp, "CELLS  %d  %d\n", nPoly, nTable);
    if (meshData->elemMap != NULL) {

      /* Write elements in order from the mesh generator */
      for (i = 0; i < meshData->nTotalElems; i++) {

        igroup = meshData->elemMap[i][0];
        ielem  = meshData->elemMap[i][1];

        nPoint = meshData->elemGroups[igroup].nPoint;

        /* number of points in the element */
        fprintf(fp, "%d", nPoint);

        /* VTK is 0-based indexed */
        for (ipnt = 0; ipnt < nPoint; ipnt++) {
          fprintf(fp, " %d", meshData->elemGroups[igroup].elements[nPoint*ielem+ipnt]-1);
        }
        fprintf(fp, "\n");
      }

    } else {

      /* Write elements by grouping */
      for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

        nPoint = meshData->elemGroups[igroup].nPoint;
        nElems = meshData->elemGroups[igroup].nElems;

        for (ielem = 0; ielem < nElems; ielem++) {

          /* number of points in the element */
          fprintf(fp, "%d", nPoint);

          /* VTK is 0-based indexed */
          for (ipnt = 0; ipnt < nPoint; ipnt++) {
            fprintf(fp, " %d", meshData->elemGroups[igroup].elements[nPoint*ielem+ipnt]-1);
          }
          fprintf(fp, "\n");
        }
      }

    }

    /* write out each of the cell types */
    fprintf(fp, "CELL_TYPES  %d\n", nPoly);

    if (meshData->elemMap != NULL) {

        /* Write elements in order from the mesh generator */
        for (i = 0; i < meshData->nTotalElems; i++) {

            igroup = meshData->elemMap[i][0];

                 if (meshData->elemGroups[igroup].elementTopo == aimLine)    elemType = 3;
            else if (meshData->elemGroups[igroup].elementTopo == aimTri)     elemType = 5;
            else if (meshData->elemGroups[igroup].elementTopo == aimQuad)    elemType = 9;
            else if (meshData->elemGroups[igroup].elementTopo == aimTet)     elemType = 10;
            else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) elemType = 14;
            else if (meshData->elemGroups[igroup].elementTopo == aimPrism)   elemType = 13;
            else if (meshData->elemGroups[igroup].elementTopo == aimHex)     elemType = 12;
            else { AIM_ERROR(aimInfo, "Unknown element type"); status = CAPS_NOTIMPLEMENT; goto cleanup; }

            /* write the element type */
            fprintf(fp, "%d\n", elemType);
        }

    } else {

        /* Write elements by grouping */
        for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {

                 if (meshData->elemGroups[igroup].elementTopo == aimLine)    elemType = 3;
            else if (meshData->elemGroups[igroup].elementTopo == aimTri)     elemType = 5;
            else if (meshData->elemGroups[igroup].elementTopo == aimQuad)    elemType = 9;
            else if (meshData->elemGroups[igroup].elementTopo == aimTet)     elemType = 10;
            else if (meshData->elemGroups[igroup].elementTopo == aimPyramid) elemType = 14;
            else if (meshData->elemGroups[igroup].elementTopo == aimPrism)   elemType = 13;
            else if (meshData->elemGroups[igroup].elementTopo == aimHex)     elemType = 12;
            else { AIM_ERROR(aimInfo, "Unknown element type"); status = CAPS_NOTIMPLEMENT; goto cleanup; }

            nElems = meshData->elemGroups[igroup].nElems;

            /* write the element type */
            for (ielem = 0; ielem < nElems; ielem++) {
              fprintf(fp, "%d\n", elemType);
            }
        }

    }

    /* write out the capsGroup ID */
    fprintf(fp, "CELL_DATA  %d\n", nPoly);
    fprintf(fp, "SCALARS ID int  1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    if (meshData->elemMap != NULL) {

        /* Write element IDs in order from the mesh generator */
        for (i = 0; i < meshData->nTotalElems; i++) {

            igroup = meshData->elemMap[i][0];
            fprintf(fp, "%d\n", meshData->elemGroups[igroup].ID);
        }

    } else {

        /* Write element IDs by grouping */
        for (igroup = 0; igroup < meshData->nElemGroup; igroup++) {
            for (ielem = 0; ielem < meshData->elemGroups[igroup].nElems; ielem++) {
                fprintf(fp, "%d\n", meshData->elemGroups[igroup].ID);
            }
        }

    }

    printf("Finished writing VTK file\n\n");

    status = CAPS_SUCCESS;

cleanup:
    if (fp != NULL) fclose(fp);

    return status;
}
