#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "caps.h"


#include "meshUtils.h"    // Meshing utilities
#include "feaUtils.h"     // FEA utilities
#include "miscUtils.h"    // Miscellaneous utilities
#include "tecplotUtils.h" // Tecplot utilities

#include "hsmTypes.h" // Bring in hsm structures


// Initiate hsmMemory structure
int initiate_hsmMemoryStruct(hsmMemoryStruct *mem)
{

    /*
    mem->numNode = 0;
    mem->numElement = 0;
    mem->numBCEdge = 0;
    mem->numBCNode = 0;
    mem->numJoint = 0;
    */

    mem->numBCNode = 0;

    mem->parg     = NULL; // Global parameters

    mem->vars     = NULL; // Node primary variables

    mem->deps      = NULL; // Node dependent variables

    mem->pars     = NULL; // Node parameters

    mem->kelem    = NULL; // Element connectivity

    mem->pare     = NULL; // Edge-BC parameters

    mem->parp     = NULL; // Node-BC parameters

    mem->kbcedge  = NULL; // Node indices of BC edge

    mem->kbcnode  = NULL; // Node of BC node

    mem->lbcnode  = NULL; // Node BC type indicator

    mem->kjoint   = NULL; // Node Indices of the Joint

    return CAPS_SUCCESS;
}


// Destroy hsmMemory structure
int destroy_hsmMemoryStruct(hsmMemoryStruct *mem)
{

    /*
    mem->numNode = 0;
    mem->numElement = 0;
    mem->numBCEdge = 0;
    mem->numBCNode = 0;
    mem->numJoint = 0;
    */

    mem->numBCNode = 0;

    EG_free(mem->parg);     mem->parg    = NULL;

    EG_free(mem->vars);     mem->vars    = NULL;

    EG_free(mem->deps);    mem->deps    = NULL;

    EG_free(mem->pars);     mem->pars    = NULL;

    EG_free(mem->kelem);    mem->kelem   = NULL;

    EG_free(mem->pare);     mem->pare    = NULL;

    EG_free(mem->parp);     mem->parp    = NULL;

    EG_free(mem->kbcedge);  mem->kbcedge = NULL;

    EG_free(mem->kbcnode);  mem->kbcnode = NULL;

    EG_free(mem->lbcnode);  mem->lbcnode = NULL;

    EG_free(mem->kjoint);   mem->kjoint  = NULL;

    return CAPS_SUCCESS;
}


// Allocate hsmMemory structure
int allocate_hsmMemoryStruct(int numNode, int numElement, int maxDim,
                             hsmMemoryStruct *mem)
{
//int numBCEdge, int numBCNode, int numJoint)

    // Global parameters
    mem->parg = (double *) EG_alloc(LGTOT*sizeof(double));
    if (mem->parg == NULL) goto bail;

    // Set to zero
    memset(mem->parg, 0, LGTOT*sizeof(double));

    if (numNode > 0) {
        mem->vars      = (double *) EG_alloc(IVTOT*numNode*sizeof(double));
        if (mem->vars     == NULL) goto bail;
        memset(mem->vars, 0,                 IVTOT*numNode*sizeof(double));

        mem->deps      = (double *) EG_alloc(JVTOT*numNode*sizeof(double));
        if (mem->deps      == NULL) goto bail;
        memset(mem->deps , 0,                JVTOT*numNode*sizeof(double));

        mem->pars      = (double *) EG_alloc(LVTOT*numNode*sizeof(double));
        if (mem->pars     == NULL) goto bail;
        memset(mem->pars, 0,                 LVTOT*numNode*sizeof(double));
    }

    if (numElement > 0) {
        mem->kelem  = (int *) EG_alloc(4*numElement*sizeof(int));
        if (mem->kelem == NULL) goto bail;
        memset(mem->kelem, 0,          4*numElement*sizeof(int));
     }

    if (maxDim > 0) { // numBCEdge
        mem->pare     = (double *) EG_alloc(LETOT*maxDim*sizeof(double));
        if (mem->pare     == NULL) goto bail;
        memset(mem->pare, 0,                LETOT*maxDim*sizeof(double));

        mem->kbcedge  = (int *)    EG_alloc(    2*maxDim*sizeof(int));
        if (mem->kbcedge  == NULL) goto bail;
        memset(mem->kbcedge, 0,                 2*maxDim*sizeof(int));
    }

    if (maxDim > 0) { // numBCNode
        mem->parp     = (double *) EG_alloc(LPTOT*maxDim*sizeof(double));
        if (mem->parp     == NULL) goto bail;
        memset(mem->parp, 0,                LPTOT*maxDim*sizeof(double));

        mem->kbcnode  = (int *)    EG_alloc(      maxDim*sizeof(int));
        if (mem->kbcnode  == NULL) goto bail;
        memset(mem->kbcnode, 0,                   maxDim*sizeof(int));

        mem->lbcnode  = (int *)    EG_alloc(      maxDim*sizeof(int));
        if (mem->lbcnode  == NULL) goto bail;
        memset(mem->lbcnode, 0,                   maxDim*sizeof(int));
    }

    if (maxDim > 0) { // numJoint
        mem->kjoint   = (int *)    EG_alloc(    2*maxDim*sizeof(int));
        if (mem->kjoint   == NULL) goto bail;
        memset(mem->kjoint, 0,                  2*maxDim*sizeof(int));
    }

    return CAPS_SUCCESS;

    bail:
        return EGADS_MALLOC;
}


// Initiate hsmTempMemory structure
int initiate_hsmTempMemoryStruct(hsmTempMemoryStruct *mem)
{

    // Don't know what these are variables
    mem->bf = NULL;
    mem->bf_dj = NULL;

    mem->bm = NULL;
    mem->bm_dj = NULL;

    mem->resc = NULL;
    mem->resc_vars = NULL;

    mem->resp = NULL;
    mem->resp_vars = NULL;
    mem->resp_dvp = NULL;

    mem->ares = NULL;

    mem->dvars = NULL;

    mem->res = NULL;

    mem->rest = NULL;
    mem->rest_t = NULL;

    mem->resv = NULL;
    mem->resv_v = NULL;

    mem->ibx  = NULL;
    mem->kdvp = NULL;
    mem->ndvp = NULL;
    mem->frst = NULL;
    mem->idt  = NULL;
    mem->frstt= NULL;
    mem->kdv  = NULL;
    mem->ndv  = NULL;
    mem->frstv= NULL;

    mem->amat  = NULL;
    mem->amatt = NULL;
    mem->amatv = NULL;
    mem->ipp   = NULL;

    return CAPS_SUCCESS;

}


// Destroy hsmTempMemory structure
int destroy_hsmTempMemoryStruct(hsmTempMemoryStruct *mem)
{

    // Don't know what these are variables
    EG_free(mem->bf);        mem->bf = NULL;
    EG_free(mem->bf_dj);     mem->bf_dj = NULL;

    EG_free(mem->bm);        mem->bm = NULL;
    EG_free(mem->bm_dj);     mem->bm_dj = NULL;

    EG_free(mem->resc);      mem->resc = NULL;
    EG_free(mem->resc_vars); mem->resc_vars = NULL;

    EG_free(mem->resp);      mem->resp = NULL;
    EG_free(mem->resp_vars); mem->resp_vars = NULL;
    EG_free(mem->resp_dvp);  mem->resp_dvp = NULL;

    EG_free(mem->ares);      mem->ares = NULL;

    EG_free(mem->dvars);     mem->dvars = NULL;

    EG_free(mem->res);       mem->res = NULL;

    EG_free(mem->rest);      mem->rest = NULL;
    EG_free(mem->rest_t);    mem->rest_t = NULL;

    EG_free(mem->resv);      mem->resv= NULL;
    EG_free(mem->resv_v);    mem->resv_v = NULL;


    EG_free(mem->ibx);       mem->ibx = NULL;
    EG_free(mem->kdvp);      mem->kdvp = NULL;
    EG_free(mem->ndvp);      mem->ndvp = NULL;
    EG_free(mem->frst);      mem->frst = NULL;
    EG_free(mem->idt);       mem->idt = NULL;
    EG_free(mem->frstt);     mem->frstt = NULL;
    EG_free(mem->kdv);       mem->kdv = NULL;
    EG_free(mem->ndv);       mem->ndv = NULL;
    EG_free(mem->frstv);     mem->frstv = NULL;

    EG_free(mem->amat);      mem->amat = NULL;
    EG_free(mem->amatt);     mem->amatt = NULL;
    EG_free(mem->amatv);     mem->amatv = NULL;
    EG_free(mem->ipp);       mem->ipp = NULL;

    return CAPS_SUCCESS;
}


// Allocate hsmTempMemory structure
int allocate_hsmTempMemoryStruct(int numNode, int maxValence, int maxDim,
                                 hsmTempMemoryStruct *mem)
{

    mem->bf       = (double *) EG_alloc(  3*3*numNode*sizeof(double));
    if (mem->bf       == NULL) goto bail;
    memset(mem->bf, 0,                    3*3*numNode*sizeof(double));

    mem->bf_dj    = (double *) EG_alloc(3*3*3*numNode*sizeof(double));
    if (mem->bf_dj    == NULL) goto bail;
    memset(mem->bf_dj, 0,               3*3*3*numNode*sizeof(double));

    mem->bm       = (double *) EG_alloc(  3*3*numNode*sizeof(double));
    if (mem->bm       == NULL) goto bail;
    memset(mem->bm, 0,                    3*3*numNode*sizeof(double));

    mem->bm_dj    = (double *) EG_alloc(3*3*3*numNode*sizeof(double));
    if (mem->bm_dj    == NULL) goto bail;
    memset(mem->bm_dj, 0,               3*3*3*numNode*sizeof(double));

    mem->resc     = (double *) EG_alloc(IVTOT*numNode*sizeof(double));
    if (mem->resc     == NULL) goto bail;
    memset(mem->resc, 0,                IVTOT*numNode*sizeof(double));

    mem->resc_vars = (double *) EG_alloc(IVTOT*IVTOT*numNode*maxValence*sizeof(double));
    if (mem->resc_vars == NULL) goto bail;
    memset(mem->resc_vars, 0,            IVTOT*IVTOT*numNode*maxValence*sizeof(double));

    mem->resp     = (double *) EG_alloc(IRTOT*numNode*sizeof(double));
    if (mem->resp     == NULL) goto bail;
    memset(mem->resp, 0,                IRTOT*numNode*sizeof(double));

    mem->resp_vars = (double *) EG_alloc(IRTOT*IVTOT*maxValence*sizeof(double));
    if (mem->resp_vars == NULL) goto bail;
    memset(mem->resp_vars, 0,            IRTOT*IVTOT*maxValence*sizeof(double));

    mem->resp_dvp = (double *) EG_alloc(IRTOT*IRTOT*maxValence*numNode*sizeof(double));
    if (mem->resp_dvp == NULL) goto bail;
    memset(mem->resp_dvp, 0,            IRTOT*IRTOT*maxValence*numNode*sizeof(double));

    mem->ares     = (double *) EG_alloc(numNode*sizeof(double));
    if (mem->ares     == NULL) goto bail;
    memset(mem->ares, 0,                numNode*sizeof(double));

    mem->dvars     = (double *) EG_alloc(IVTOT*numNode*sizeof(double));
    if (mem->dvars     == NULL) goto bail;
    memset(mem->dvars, 0,                IVTOT*numNode*sizeof(double));

    mem->res      = (double *) EG_alloc(6*numNode*sizeof(double));
    if (mem->res      == NULL) goto bail;
    memset(mem->res, 0,                 6*numNode*sizeof(double));

    mem->rest     = (double *) EG_alloc(3*4*numNode*sizeof(double));
    if (mem->rest     == NULL) goto bail;
    memset(mem->rest, 0,                3*4*numNode*sizeof(double));

    mem->rest_t   = (double *) EG_alloc(3*3*maxValence*numNode*sizeof(double));
    if (mem->rest_t   == NULL) goto bail;
    memset(mem->rest_t, 0,              3*3*maxValence*numNode*sizeof(double));

    mem->resv     = (double *) EG_alloc(2*2*numNode*sizeof(double));
    if (mem->resv     == NULL) goto bail;
    memset(mem->resv, 0,                2*2*numNode*sizeof(double));

    mem->resv_v   = (double *) EG_alloc(2*2*maxValence*numNode*sizeof(double));
    if (mem->resv_v   == NULL) goto bail;
    memset(mem->resv_v, 0,              2*2*maxValence*numNode*sizeof(double));

    mem->ibx   = (int *) EG_alloc(6*maxDim*sizeof(int));
    if (mem->ibx   == NULL) goto bail;
    memset(mem->ibx, 0,           6*maxDim*sizeof(int));

    mem->kdvp  = (int *) EG_alloc(maxValence*numNode*sizeof(int));
    if (mem->kdvp  == NULL) goto bail;
    memset(mem->kdvp, 0,          maxValence*numNode*sizeof(int));

    mem->ndvp  = (int *) EG_alloc(numNode*sizeof(int));
    if (mem->ndvp  == NULL) goto bail;
    memset(mem->kdvp, 0,          numNode*sizeof(int));

    mem->frst  = (int *) EG_alloc((3*numNode+1)*sizeof(int));
    if (mem->frst  == NULL) goto bail;
    memset(mem->frst, 0,          (3*numNode+1)*sizeof(int));

    mem->idt   = (int *) EG_alloc((maxValence+1)*numNode*sizeof(int));
    if (mem->idt   == NULL) goto bail;
    memset(mem->idt, 0,           (maxValence+1)*numNode*sizeof(int));

    mem->frstt = (int *) EG_alloc((3*numNode+1)*sizeof(int));
    if (mem->frstt == NULL) goto bail;
    memset(mem->frstt, 0,         (3*numNode+1)*sizeof(int));

    mem->kdv   = (int *) EG_alloc(maxValence*numNode*sizeof(int));
    if (mem->kdv   == NULL) goto bail;
    memset(mem->kdv, 0,           maxValence*numNode*sizeof(int));

    mem->ndv   = (int *) EG_alloc(numNode*sizeof(int));
    if (mem->ndv   == NULL) goto bail;
    memset(mem->ndv, 0,           numNode*sizeof(int));

    mem->frstv = (int *) EG_alloc((3*numNode+1)*sizeof(int));
    if (mem->frstv == NULL) goto bail;
    memset(mem->frstv, 0,         (3*numNode+1)*sizeof(int));

    return CAPS_SUCCESS;

    bail:
        return EGADS_MALLOC;
}


// Convert an EGADS body to a boundary element model - disjointed at edges
int hsm_bodyToBEM(void *aimInfo,
                  ego    ebody,                        // (in)  EGADS Body
                  double paramTess[3],                 // (in)  Tessellation parameters
                  int    edgePointMin,                 // (in)  minimum points along any Edge
                  int    edgePointMax,                 // (in)  maximum points along any Edge
                  int    quadMesh,                     // (in)  only do tris-for faces
                  mapAttrToIndexStruct *attrMap,       // (in)  map from CAPSGroup names to indexes
                  mapAttrToIndexStruct *coordSystemMap,// (in)  map from CoordSystem names to indexes
                  mapAttrToIndexStruct *constraintMap, // (in)  map from CAPSConstraint names to indexes
                  mapAttrToIndexStruct *loadMap,       // (in)  map from CAPSLoad names to indexes
                  mapAttrToIndexStruct *transferMap,   // (in)  map from CAPSTransfer names to indexes
       /*@null@*/ mapAttrToIndexStruct *connectMap,    // (in)  map from CAPSConnect names to indexes
                  meshStruct *feaMesh)                 // (out) FEA mesh structure
{
    int status = 0; // Function return status

    int i, j, face, edge, patch; // Indexing

    // Body entities
    int numNode = 0, numEdge = 0, numFace = 0;
    ego *enodes=NULL, *eedges=NULL, *efaces=NULL;

    // Meshing
    int numPoint = 0, numTri = 0;
    const int  *triConn = NULL, *triNeighbor = NULL; // Triangle

    int numPatch = 0; // Patching
    int n1, n2;

    const double *xyz, *uv;

    const int *pointType = NULL, *pointTopoIndex = NULL, *pvindex = NULL, *pbounds = NULL;

    int       periodic, nchange, oclass, mtype, nchild, *senses;

    // Edge point distributions
    int    *points=NULL, *isouth=NULL, *ieast=NULL, *inorth=NULL, *iwest=NULL;
    double params[3], bbox[6], size, range[2], arclen, data[4];
    double *rpos=NULL;
    ego    eref, *echilds, eloop;

    int          n, cnt, iloop, iedge, eindex, last, numEdgePoints, nloop;
    int          nedge, *sen, *qints=NULL;
    double       uvbox[4];
    const double *xyzs = NULL, *ts = NULL;
    ego          *loops = NULL, *edges = NULL, *nodes = NULL;

    // Attributues
    const char *attrName;
    int         attrIndex, coordSystemIndex, loadIndex;

    int numElement = 0; // Number of elements

    feaMeshDataStruct *feaData;

    // Geometry calculations
    double result[18];

    // ---------------------------------------------------------------

    printf("Creating HSM BEM\n");

    // Get number of Nodes, Edges, and Faces in ebody
    status = EG_getBodyTopos(ebody, NULL, NODE, &numNode, &enodes);
    if (status < EGADS_SUCCESS) goto cleanup;

    status = EG_getBodyTopos(ebody, NULL, EDGE, &numEdge, &eedges);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (eedges == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    status = EG_getBodyTopos(ebody, NULL, FACE, &numFace, &efaces);
    if (status < EGADS_SUCCESS) goto cleanup;
    if (efaces == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    // Determine the nominal number of points along each Edge
    points = (int    *) EG_alloc((numEdge+1)     *sizeof(int   ));
    rpos   = (double *) EG_alloc((edgePointMax)*sizeof(double));
    if (points == NULL || rpos == NULL) {
        status = EGADS_MALLOC;
        printf("\tError in hsm_bodyToBEM: EG_alloc\n");
        goto cleanup;
    }

    status = EG_getBoundingBox(ebody, bbox);
    if (status < EGADS_SUCCESS) {
        printf("\tError in hsm_bodyToBEM: EG_getBoundingBox\n");
        goto cleanup;
    }

    size = sqrt( (bbox[3] - bbox[0]) * (bbox[3] - bbox[0])
                +(bbox[4] - bbox[1]) * (bbox[4] - bbox[1])
                +(bbox[5] - bbox[2]) * (bbox[5] - bbox[2]));

    params[0] = paramTess[0] * size;
    params[1] = paramTess[1] * size;
    params[2] = paramTess[2];

    status = EG_attributeAdd(ebody, ".tParam", ATTRREAL, 3,
                             NULL, params, NULL);
    if (status < EGADS_SUCCESS) {
        printf("\tError in hsm_bodyToBEM: EG_attributeAdd\n");
        goto cleanup;
    }

    for (i = 1; i <= numEdge; i++) {
        status = EG_getRange(eedges[i-1], range, &periodic);
        if (status < EGADS_SUCCESS) {
            printf("\tError in hsm_bodyToBEM: EG_getRange\n");
            goto cleanup;
        }

        status = EG_arcLength(eedges[i-1], range[0], range[1], &arclen);
        if (status < EGADS_SUCCESS) {
            printf("\tError in hsm_bodyToBEM: EG_arcLength\n");
            goto cleanup;
        }

      //points[i] = MIN(MAX(MAX(edgePointMin,2), (int)(1+arclen/params[0])),
      //                edgePointMax);
        points[i] = (int) min_DoubleVal(
                              max_DoubleVal(
                                  max_DoubleVal( (double) edgePointMin, 2.0),
                                                 (double) (1+arclen/params[0])),
                                 (double) edgePointMax);
    }

    // make arrays for "opposite" sides of four-sided Faces (with only one loop)
    isouth = (int *) EG_alloc((numFace+1)*sizeof(int));
    ieast  = (int *) EG_alloc((numFace+1)*sizeof(int));
    inorth = (int *) EG_alloc((numFace+1)*sizeof(int));
    iwest  = (int *) EG_alloc((numFace+1)*sizeof(int));

    if (isouth == NULL ||
        ieast == NULL  ||
        inorth == NULL ||
        iwest == NULL   ) {

        status = EGADS_MALLOC;
        printf("\tError in hsm_bodyToBEM: EG_alloc\n");
        goto cleanup;
    }

    for (i = 1; i <= numFace; i++) {
        isouth[i] = 0;
        ieast [i] = 0;
        inorth[i] = 0;
        iwest [i] = 0;

        status = EG_getTopology(efaces[i-1],
                                &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        if (status < EGADS_SUCCESS) goto cleanup;

        if (nchild != 1) continue;

        eloop = echilds[0];
        status = EG_getTopology(eloop, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        if (status < EGADS_SUCCESS) goto cleanup;

        if (nchild != 4) continue;

        isouth[i] = status = EG_indexBodyTopo(ebody, echilds[0]);
        if (status < EGADS_SUCCESS) goto cleanup;

        ieast[i]  = status = EG_indexBodyTopo(ebody, echilds[1]);
        if (status < EGADS_SUCCESS) goto cleanup;

        inorth[i] = status = EG_indexBodyTopo(ebody, echilds[2]);
        if (status < EGADS_SUCCESS) goto cleanup;

        iwest[i]  = status = EG_indexBodyTopo(ebody, echilds[3]);
        if (status < EGADS_SUCCESS) goto cleanup;
    }

    // make "opposite" sides of four-sided Faces (with only one loop) match
    nchange = 1;
    for (i = 0; i < 20; i++) {
        nchange = 0;

        for (face = 1; face <= numFace; face++) {
            if (isouth[face] <= 0 || ieast[face] <= 0 ||
                inorth[face] <= 0 || iwest[face] <= 0   ) continue;

            if        (points[iwest[face]] < points[ieast[face]]) {
                points[iwest[face]] = points[ieast[face]];
                nchange++;

            } else if (points[ieast[face]] < points[iwest[face]]) {
                points[ieast[face]] = points[iwest[face]];
                nchange++;
            }

            if        (points[isouth[face]] < points[inorth[face]]) {
                points[isouth[face]] = points[inorth[face]];
                nchange++;
            } else if (points[inorth[face]] < points[isouth[face]]) {
                points[inorth[face]] = points[isouth[face]];
                nchange++;
            }
        }
        if (nchange == 0) break;
    }
    if (nchange > 0) {
        status = -999;
        goto cleanup;
    }

    // mark the Edges with points[iedge] evenly-spaced points
    for (edge = 1; edge <= numEdge; edge++) {
        for (i = 1; i < points[edge]-1; i++) {

            rpos[i-1] = (double)(i) / (double)(points[edge]-1);
        }

        if (points[edge] == 2) {
            i = 0;
            status = EG_attributeAdd(eedges[edge-1], ".rPos", ATTRINT,
                                     1, &i, NULL, NULL);
            if (status < EGADS_SUCCESS) goto cleanup;
        } else {
            status = EG_attributeAdd(eedges[edge-1], ".rPos", ATTRREAL,
                                     points[edge]-2, NULL, rpos, NULL);
            if (status < EGADS_SUCCESS) goto cleanup;
        }
    }

    // Make tessellation
    status = EG_makeTessBody(ebody, params, &feaMesh->egadsTess);
    if (status != EGADS_SUCCESS) {
        printf("\tError in hsm_bodyToBEM: EG_makeTessBody\n");
        goto cleanup;
    }

    // Make Quads on each four-sided Face
    params[0] = 0;
    params[1] = 0;
    params[2] = 0;

    // If making quads on faces lets setup an array to keep track of which faces have been quaded.
    if (quadMesh == (int) true) {
        if (numFace > 0) {
            AIM_ALLOC(qints, numFace, int, aimInfo, status);

            // Set default to 0
            for (face = 0; face < numFace; face++)
                qints[face] = 0;
        }
    }

    if (quadMesh == (int) true) {
        for (face = 1; face <= numFace; face++) {
            if (iwest[face] <= 0) continue;

            status = EG_makeQuads(feaMesh->egadsTess, params, face);
            if (status < EGADS_SUCCESS) goto cleanup;
        }
    }

    // Set the mesh type information
    feaMesh->meshType = SurfaceMesh;
    feaMesh->analysisType = MeshStructure;

    feaMesh->numNode = 0; // This should be redundant
    feaMesh->numElement = 0;

    // Fill element information

    if (quadMesh == (int) true && numFace > 0) {
        printf("\tGetting quads for BEM!\n");

        // Turn off meshQuick guide if you are getting quads
        feaMesh->meshQuickRef.useStartIndex = (int) false;
    } else {
        feaMesh->meshQuickRef.useStartIndex = (int) true;
        feaMesh->meshQuickRef.startIndexTriangle = numElement;
    }
    feaMesh->meshQuickRef.useStartIndex = (int) false;

    // Get Tris and Quads from faces
    for (face = 0; face < numFace; face++) {

        status = retrieve_CAPSIgnoreAttr(efaces[face], &attrName);
        if (status == CAPS_SUCCESS) {
            printf("\tcapsIgnore attribute found for face - %d!! - NOT currently allowed\n",
                   face+1);
            status = CAPS_BADVALUE;
            goto cleanup; // This isn't currently allowed - see element->connectivity in permutation array
            //continue;
        }

        status = retrieve_CAPSGroupAttr(efaces[face], &attrName);
        if (status != CAPS_SUCCESS) {
            printf("Error: no capsGroup attribute found for face - %d!!\n",
                   face+1);
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(attrMap, attrName, &attrIndex);
        if (status != CAPS_SUCCESS) {
            printf("Error: capsGroup name %s not found in attribute to index map\n",
                   attrName);
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(coordSystemMap, attrName,
                                         &coordSystemIndex);
        if (status != CAPS_SUCCESS) coordSystemIndex = 0;

        loadIndex = CAPSMAGIC;
        status = retrieve_CAPSLoadAttr(efaces[face], &attrName);
        if (status == CAPS_SUCCESS) {

            status = get_mapAttrToIndexIndex(loadMap, attrName, &loadIndex);

            if (status != CAPS_SUCCESS) {
                printf("Error: capsLoad name %s not found in attribute to index map\n",
                       attrName);
                goto cleanup;
            }
        }

        if (quadMesh == (int) true) {
            status = EG_getQuads(feaMesh->egadsTess, face+1,
                                 &numPoint, &xyz, &uv, &pointType,
                                 &pointTopoIndex, &numPatch);
            if (status < EGADS_SUCCESS) goto cleanup;

        } else numPatch = -1;

        if ((numPatch > 0) && (qints != NULL)) {

            qints[face] = 0;
            for (patch = 1; patch <= numPatch; patch++) {

                status = EG_getPatch(feaMesh->egadsTess, face+1,
                                     patch, &n1, &n2, &pvindex, &pbounds);
                if (status < EGADS_SUCCESS) goto cleanup;
                if (pvindex == NULL) {
                    status = CAPS_NULLVALUE;
                    goto cleanup;
                }

                for (j = 1; j < n2; j++) {
                    for (i = 1; i < n1; i++) {
                        numElement += 1;

                        feaMesh->meshQuickRef.numQuadrilateral += 1;
                        feaMesh->numElement = numElement;

                        feaMesh->element = (meshElementStruct *)
                                           EG_reall(feaMesh->element,
                                                    feaMesh->numElement*
                                                    sizeof(meshElementStruct));

                        if (feaMesh->element == NULL) {
                            status = EGADS_MALLOC;
                            goto cleanup;
                        }

                        status = initiate_meshElementStruct(&feaMesh->element[numElement-1],
                                                            feaMesh->analysisType);
                        if (status != CAPS_SUCCESS) goto cleanup;

                        qints[face] += 1;

                        feaMesh->element[numElement-1].elementType = Quadrilateral;

                        feaMesh->element[numElement-1].elementID = numElement;

                        status = mesh_allocMeshElementConnectivity(&feaMesh->element[numElement-1]);
                        if (status != CAPS_SUCCESS) goto cleanup;

                        feaMesh->element[numElement-1].connectivity[0] =
                                    pvindex[(i-1)+n1*(j-1)] + feaMesh->numNode;
                        feaMesh->element[numElement-1].connectivity[1] =
                                    pvindex[(i  )+n1*(j-1)] + feaMesh->numNode;
                        feaMesh->element[numElement-1].connectivity[2] =
                                    pvindex[(i  )+n1*(j  )] + feaMesh->numNode;
                        feaMesh->element[numElement-1].connectivity[3] =
                                    pvindex[(i-1)+n1*(j  )] + feaMesh->numNode;

                        feaMesh->element[numElement-1].markerID = attrIndex;

                        feaData = (feaMeshDataStruct *)
                                  feaMesh->element[numElement-1].analysisData;

                        feaData->propertyID = attrIndex;
                        feaData->coordID = coordSystemIndex;
                        feaData->loadIndex = loadIndex;
                    }
                }
            }
        } else {
            status = EG_getTessFace(feaMesh->egadsTess, face+1,
                                    &numPoint, &xyz, &uv, &pointType,
                                    &pointTopoIndex,
                                    &numTri, &triConn, &triNeighbor);
            if (status < EGADS_SUCCESS) goto cleanup;
            if (triConn   == NULL) {
                status = CAPS_NULLVALUE;
                goto cleanup;
            }

            feaMesh->element = (meshElementStruct *)
                               EG_reall(feaMesh->element,
                                        (feaMesh->numElement+numTri)*
                                        sizeof(meshElementStruct));

            if (feaMesh->element == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (i= 0; i < numTri; i++) {

                numElement += 1;

                feaMesh->meshQuickRef.numTriangle += 1;
                feaMesh->numElement = numElement;

                status = initiate_meshElementStruct(&feaMesh->element[numElement-1],
                                                    feaMesh->analysisType);
                if (status != CAPS_SUCCESS) goto cleanup;

                feaMesh->element[numElement-1].elementType = Triangle;

                feaMesh->element[numElement-1].elementID = numElement;

                status = mesh_allocMeshElementConnectivity(&feaMesh->element[numElement-1]);
                if (status != CAPS_SUCCESS) goto cleanup;

                feaMesh->element[numElement-1].connectivity[0] =
                                           triConn[3*i + 0] + feaMesh->numNode;
                feaMesh->element[numElement-1].connectivity[1] =
                                           triConn[3*i + 1] + feaMesh->numNode;
                feaMesh->element[numElement-1].connectivity[2] =
                                           triConn[3*i + 2] + feaMesh->numNode;

                feaMesh->element[numElement-1].markerID = attrIndex;

                feaData = (feaMeshDataStruct *)
                           feaMesh->element[numElement-1].analysisData;

                feaData->propertyID = attrIndex;
                feaData->coordID = coordSystemIndex;
                feaData->loadIndex = loadIndex;

            }
        }

        // also extract all edge elements from the face
        status = EG_getTopology(efaces[face], &eref, &oclass, &mtype, uvbox,
                                &nloop, &loops, &senses);
        if (status != EGADS_SUCCESS) goto cleanup;
        if (loops == NULL) {
            status = CAPS_NULLVALUE;
            goto cleanup;
        }

        for (cnt = iloop = 0; iloop < nloop; iloop++) {
            status = EG_getTopology(loops[iloop], &eref, &oclass, &mtype,
                                    NULL, &nedge, &edges, &senses);
            if (status != EGADS_SUCCESS) goto cleanup;
            if (edges == NULL) {
                status = CAPS_NULLVALUE;
                goto cleanup;
            }

            last = cnt;
            for (iedge = 0; iedge < nedge; iedge++) {
                status = EG_getTopology(edges[iedge], &eref, &oclass,
                                        &mtype, range, &n, &nodes, &sen);
                if (status != EGADS_SUCCESS) goto cleanup;
                if (mtype == DEGENERATE) continue;

                // get the load information on the edge
                loadIndex = CAPSMAGIC;
                status    = retrieve_CAPSLoadAttr(edges[iedge], &attrName);
                if (status == CAPS_SUCCESS) {

                    status = get_mapAttrToIndexIndex(loadMap, attrName,
                                                     &loadIndex);

                    if (status != CAPS_SUCCESS) {
                        printf("Error: capsLoad name %s not found in attribute to index map\n",
                               attrName);
                        goto cleanup;
                    }
                }

                eindex = EG_indexBodyTopo(ebody, edges[iedge]);
                if (eindex < EGADS_SUCCESS) goto cleanup;

                status = EG_getTessEdge(feaMesh->egadsTess, eindex,
                                        &numEdgePoints, &xyzs, &ts);
                if (status != EGADS_SUCCESS) goto cleanup;

                feaMesh->element = (meshElementStruct *)
                                   EG_reall(feaMesh->element,
                                            (feaMesh->numElement+numEdgePoints-1)*
                                            sizeof(meshElementStruct));

                if (feaMesh->element == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                for (i = 0; i < numEdgePoints-1; i++, cnt++) {
                    status = initiate_meshElementStruct(&feaMesh->element[numElement],
                                                        feaMesh->analysisType);
                    if (status != CAPS_SUCCESS) goto cleanup;

                    feaMesh->element[numElement].elementType = Line;

                    feaMesh->element[numElement].elementID = numElement+1;

                    status = mesh_allocMeshElementConnectivity(&feaMesh->element[numElement]);
                    if (status != CAPS_SUCCESS) goto cleanup;

                    feaMesh->element[numElement].connectivity[0] = cnt+1 + feaMesh->numNode;
                    feaMesh->element[numElement].connectivity[1] = cnt+2 + feaMesh->numNode;

                    feaMesh->element[numElement].markerID = attrIndex;

                    feaData = (feaMeshDataStruct *) feaMesh->element[numElement].analysisData;

                    feaData->propertyID = attrIndex;
                    feaData->coordID = coordSystemIndex;
                    feaData->loadIndex = loadIndex;

                    numElement += 1;
                }
                feaMesh->numElement = numElement;
            }
            feaMesh->element[numElement-1].connectivity[1] = last+1 + feaMesh->numNode;
        }

        // Get node information
        feaMesh->node = (meshNodeStruct *)
                        EG_reall(feaMesh->node,
                                 (feaMesh->numNode+numPoint)*
                                 sizeof(meshNodeStruct));
        if (feaMesh->node == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }
        if ((pointType == NULL) || (pointTopoIndex == NULL)) {
            status = CAPS_NULLVALUE;
            goto cleanup;
        }

        for (i = 0; i < numPoint; i++) {

            j = i + feaMesh->numNode;

            status = initiate_meshNodeStruct(&feaMesh->node[j],
                                             feaMesh->analysisType);
            if (status != CAPS_SUCCESS) goto cleanup;

            feaMesh->node[j].nodeID = j+1;
            feaMesh->node[j].xyz[0] = xyz[3*i+0];
            feaMesh->node[j].xyz[1] = xyz[3*i+1];
            feaMesh->node[j].xyz[2] = xyz[3*i+2];

            // Get geometry data for node
            feaMesh->node[j].geomData = (meshGeomDataStruct *)
                                        EG_alloc(sizeof(meshGeomDataStruct));

            if (feaMesh->node[j].geomData == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            status = initiate_meshGeomDataStruct(feaMesh->node[j].geomData);
            if (status != CAPS_SUCCESS) goto cleanup;

            feaMesh->node[j].geomData->type = pointType[i];
            feaMesh->node[j].geomData->topoIndex = pointTopoIndex[i];

            // Want the face index to be set for topoIndex
            if (feaMesh->node[j].geomData->topoIndex < 0)
                feaMesh->node[j].geomData->topoIndex = face+1;

            feaMesh->node[j].geomData->uv[0] = uv[2*i + 0];
            feaMesh->node[j].geomData->uv[1] = uv[2*i + 1];

            status = EG_evaluate(efaces[face], (const double *)
                                 &feaMesh->node[j].geomData->uv, result);
            if (status != EGADS_SUCCESS) goto cleanup;

            // U
            feaMesh->node[j].geomData->firstDerivative[0] = result[3];
            feaMesh->node[j].geomData->firstDerivative[1] = result[4];
            feaMesh->node[j].geomData->firstDerivative[2] = result[5];

            // V
            feaMesh->node[j].geomData->firstDerivative[3] = result[6];
            feaMesh->node[j].geomData->firstDerivative[4] = result[7];
            feaMesh->node[j].geomData->firstDerivative[5] = result[8];

            // Get attributes
            feaData = (feaMeshDataStruct *) feaMesh->node[j].analysisData;
/*@-nullpass@*/
            status = fea_setFEADataPoint(efaces, eedges, enodes,
                                         attrMap,
                                         coordSystemMap,
                                         constraintMap,
                                         loadMap,
                                         transferMap,
                                         connectMap,
                                         NULL,
                                         feaMesh->node[j].geomData->type,
                                         feaMesh->node[j].geomData->topoIndex,
                                         feaData);
/*@+nullpass@*/
            if (status != CAPS_SUCCESS) goto cleanup;
            feaData->propertyID = attrIndex;
        }
        feaMesh->numNode += numPoint;
    }

    if (qints != NULL) {
        status = EG_attributeAdd(feaMesh->egadsTess, ".mixed", ATTRINT, numFace, qints, NULL, NULL);
        AIM_STATUS(aimInfo, status);
    }


cleanup:

    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in hsm_bodyToBEM, status %d\n", status);

    AIM_FREE(iwest );
    AIM_FREE(inorth);
    AIM_FREE(ieast );
    AIM_FREE(isouth);
    AIM_FREE(rpos  );
    AIM_FREE(points );

    AIM_FREE(enodes);
    AIM_FREE(eedges);
    AIM_FREE(efaces);
    AIM_FREE(qints);

    return status;
}


// Write hsm data to a Tecplot file
int hsm_writeTecplot(void *aimInfo, char *projectName, meshStruct feaMesh,
                     hsmMemoryStruct *hsmMemory, int permutation[])
{

    int status; // Function return status

    int i, j, k, m, elem; // Indexing

    meshElementStruct *element; // Not freeable

    int stringLength;
    char *filename = NULL;
    char fileExt[] = ".dat";
    int numElement = 0;

    // Tecplot output
    const char * variableName[] = {"x" , "y" , "z",
                                   "x'", "y'", "z'",
                                   "x Displacement", "y Displacement", "z Displacement",
                                   "x'<sub>Material normal</sub>",
                                   "y'<sub>Material normal</sub>",
                                   "z'<sub>Material normal</sub>",
                                   "drilling rotation DOF",
                                   "e<sub>1,x</sub>", "e<sub>1,y</sub>", "e<sub>1,z</sub>",
                                   "e<sub>2,x</sub>", "e<sub>2,y</sub>", "e<sub>2,z</sub>",
                                   "n<sub>x</sub>", "n<sub>y</sub>", "n<sub>z</sub>",
                                   "Strain, <greek>e</greek><sub>11</sub>",
                                   "Strain, <greek>e</greek><sub>22</sub>",
                                   "Strain, <greek>e</greek><sub>12</sub>",
                                   "Curv. Change, <greek>k</greek><sub>11</sub>",
                                   "Curv. Change, <greek>k</greek><sub>22</sub>",
                                   "Curv. Change, <greek>k</greek><sub>12</sub>",
                                   "Stress, f<sub>11</sub>", "Stress, f<sub>22</sub>",
                                   "Stress, f<sub>12</sub>",
                                   "Stress Mom., m<sub>11</sub>",
                                   "Stress Mom., m<sub>22</sub>",
                                   "Stress Mom., m<sub>12</sub>",
                                   "Shear Stress, f<sub>1n</sub>",
                                   "Shear Stress, f<sub>2n</sub>",
                                   "Tilt Angle, <greek>g</greek><sub>1</sub>",
                                   "Tilt Angle, <greek>g</greek><sub>2</sub>",
                                    };
    // Needs to be consistent with variableName
    int numOutVariable = sizeof(variableName)/sizeof(const char *);

    double **dataMatrix = NULL;
    int *connectMatrix =  NULL;

    if (permutation == NULL) return CAPS_NULLVALUE;

    stringLength = strlen(projectName) + strlen(fileExt) + 1;

    // Write out Tecplot files
    filename = (char *) EG_alloc(stringLength*sizeof(char));
    if (filename == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    strcpy(filename, projectName);
    strcat(filename, fileExt);

    // Lets allocate our Tecplot matrix
    dataMatrix    = (double **) EG_alloc(numOutVariable*sizeof(double));

    numElement = feaMesh.meshQuickRef.numTriangle +
                 feaMesh.meshQuickRef.numQuadrilateral;
    connectMatrix = (int *) EG_alloc(4*numElement*sizeof(int));

    if (dataMatrix == NULL  || connectMatrix == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < numOutVariable; i++) dataMatrix[i] = NULL;

    for (i = 0; i < numOutVariable; i++) {

        dataMatrix[i] = (double *) EG_alloc(feaMesh.numNode*sizeof(double));

        if (dataMatrix[i] == NULL) { // If allocation failed ....

            status =  EGADS_MALLOC;
            goto cleanup;
        }
    }

    // Set the data
    for (i = 0; i < feaMesh.numNode; i++) {

      k = permutation[i]-1;
      m = 0;

      // XYZ
      dataMatrix[m++][i] = hsmMemory->pars[k*LVTOT + lvr0x];
      dataMatrix[m++][i] = hsmMemory->pars[k*LVTOT + lvr0y];
      dataMatrix[m++][i] = hsmMemory->pars[k*LVTOT + lvr0z];

      // Deformed XYZ
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivrx];
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivry];
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivrz];

      // Displacement
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivrx] -
                           hsmMemory->pars[k*LVTOT + lvr0x];
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivry] -
                           hsmMemory->pars[k*LVTOT + lvr0y];
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivrz] -
                           hsmMemory->pars[k*LVTOT + lvr0z];

      // Unit material-normal vector of deformed geometry
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivdx];
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivdy];
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivdz];

      // drilling rotation DOF
      dataMatrix[m++][i] = hsmMemory->vars[k*IVTOT + ivps];


      // e - local basis unit tangential vector 1
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve1x];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve1y];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve1z];

      // e - local basis unit tangential vector 2
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve2x];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve2y];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve2z];

      // n - local basis unit normal vector
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvnx];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvny];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvnz];

      // eps - strain
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve11];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve22];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jve12];

      // kap - curvature change
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvk11];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvk22];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvk12];

      // f - stress resultant
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvf11];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvf22];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvf12];

      // m - stress-moment resultant
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvm11];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvm22];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvm12];

      // fn - transverse shear stress resultant
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvf1n];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvf2n];

      // gam - n tilt angle in e* direction
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvga1];
      dataMatrix[m++][i] = hsmMemory->deps[k*JVTOT + jvga2];
    }

    // Set the data and connectivity
    elem = 0;
    for (i = 0; i < feaMesh.numElement; i++) {

        element = &feaMesh.element[i];

        if (element->elementType == Line) continue;

        if (element->elementType != Triangle && element->elementType !=
            Quadrilateral) {
            printf("Unsupported element type\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        if (element->connectivity == NULL) {
            status = CAPS_NULLVALUE;
            goto cleanup;
        }

        for (j = 0; j < 4; j++) {

            if (element->elementType == Triangle && j == 3) {

                m = element->connectivity[j-1];

            } else {

                m = element->connectivity[j];
            }

            connectMatrix[4*elem+j] = m;
        }
        elem++;
    }
/*@-nullpass@*/
    status = tecplot_writeFEPOINT(aimInfo,
                                  filename,
                                  "HSM solution to Tecplot",
                                  "HSM solution",
                                  numOutVariable,
                                  (char **) variableName,
                                  feaMesh.numNode,
                                  dataMatrix,
                                  NULL,
                                  numElement,
                                  connectMatrix,
                                  NULL);
/*@+nullpass@*/
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in hsm_writeTecplot status = %d\n",
                status);

    if (filename != NULL) EG_free(filename);

    if (dataMatrix != NULL) {
        for (i = 0; i < numOutVariable; i++) {
            if (dataMatrix[i] != NULL)  EG_free(dataMatrix[i]);
        }
        EG_free(dataMatrix);
    }

    if (connectMatrix !=  NULL) EG_free(connectMatrix);

    return status;
}


// Set global parameters in hsmMemory structure
int hsm_setGlobalParameter(feaProblemStruct feaProblem,
                           hsmMemoryStruct *hsmMemory)
{
    int loadIndex; // Indexing

    double normalize;

    feaLoadStruct *feaLoad;

    printf("NEED TO ADD MORE CODE TO hsm_setGlobalParameter\n");

    for (loadIndex = 0; loadIndex < feaProblem.numLoad; loadIndex++) {

        feaLoad = &feaProblem.feaLoad[loadIndex];

        normalize = dot_DoubleVal(feaLoad->directionVector, feaLoad->directionVector);
        normalize = sqrt(normalize);

        if (feaLoad->loadType == Gravity) {
            hsmMemory->parg[lggeex] =
              feaLoad->gravityAcceleration*feaLoad->directionVector[0]/normalize;
            hsmMemory->parg[lggeey] =
              feaLoad->gravityAcceleration*feaLoad->directionVector[1]/normalize;
            hsmMemory->parg[lggeez] =
              feaLoad->gravityAcceleration*feaLoad->directionVector[2]/normalize;
        }

        // Linear - velocity + acceleration  -- NO INPUTS YET
//        hsmMemory->parg[lgvelx] = 0.0;
//        hsmMemory->parg[lgvely] = 0.0;
//        hsmMemory->parg[lgvelz] = 0.0;
//
//        hsmMemory->parg[lgvacx] = 0.0;
//        hsmMemory->parg[lgvacy] = 0.0;
//        hsmMemory->parg[lgvacz] = 0.0;

        if (feaLoad->loadType == Rotational) {

            hsmMemory->parg[lgrotx] = feaLoad->angularVelScaleFactor;
            hsmMemory->parg[lgroty] = feaLoad->angularVelScaleFactor;
            hsmMemory->parg[lgrotz] = feaLoad->angularVelScaleFactor;

            hsmMemory->parg[lgracx] = feaLoad->angularAccScaleFactor;
            hsmMemory->parg[lgracy] = feaLoad->angularAccScaleFactor;
            hsmMemory->parg[lgracz] = feaLoad->angularAccScaleFactor;
        }
    }

// NOT SET YET
//             &   lgposx = 16,  !  R_X   position of xyz origin in XYZ (earth) axes
//             &   lgposy = 17,  !  R_Y
//             &   lgposz = 18,  !_ R_Z
//             &   lgephi = 19,  !  Phi   Euler angles of xyz frame
//             &   lgethe = 20,  !  Theta
//             &   lgepsi = 21,  !_ Psi
//             &   lggabx = 22,  !  g_X   gravity in XYZ (earth) axes, typically (0,0,-g)
//             &   lggaby = 23,  !  g_Y
//             &   lggabz = 24,  !_ g_Z)

    return CAPS_SUCCESS;
}


// Set surface on nodes in hsmMemory->pars structure
int hsm_setSurfaceParameter(feaProblemStruct feaProblem, int permutation[],
                            hsmMemoryStruct *hsmMemory)
{

    int status; // Function return status
    int i, j, k, m, loadIndex; // Indexing

    int numConnect;
    feaLoadStruct *feaLoad;

    for (loadIndex = 0; loadIndex < feaProblem.numLoad; loadIndex++) {

        feaLoad = &feaProblem.feaLoad[loadIndex];

        // Grid Force and Moment loads
        if (feaLoad->loadType == GridForce ||
            feaLoad->loadType == GridMoment) {
            continue; // processessed in hsm_setNodeBCParameter
        // Edge Force and Moment loads
        } else if (feaLoad->loadType == LineForce ||
                   feaLoad->loadType == LineMoment) {
            continue; // processessed in hsm_setEdgeBCParameter
        // Pressure
        } else if ((feaLoad->loadType == Pressure) ||
                   (feaLoad->loadType == PressureDistribute) ||
                   (feaLoad->loadType == PressureExternal)) {

            for (i = 0; i < feaLoad->numElementID; i++) {

                if (feaLoad->elementIDSet == NULL) {
                    printf("Error: NULL gridIDSet!\n");
                    status = CAPS_NULLVALUE;
                    goto cleanup;
                }

                for (j = 0; j < feaProblem.feaMesh.numElement; j++) {
                    if (feaProblem.feaMesh.element[j].elementID ==
                        feaLoad->elementIDSet[i]) break;
                }

                numConnect = mesh_numMeshConnectivity(feaProblem.feaMesh.element[j].elementType);
                for (m = 0; m < numConnect; m++) {

                    // Get index in hsmMemory
                    k =  permutation[feaProblem.feaMesh.element[j].connectivity[m]-1] -1;

                    //Loading - shell-following normal load/area
                    if (feaLoad->loadType == Pressure) {

                        hsmMemory->pars[k*LVTOT+lvqn] = feaLoad->pressureForce;
                    }

                    // Check for element size
                    if ((feaLoad->loadType == PressureDistribute) ||
                        (feaLoad->loadType == PressureExternal) ) {

                        if (numConnect > 4) {
                            printf("Error: Unsupported element type (connectivity length = %d) for load type PressureDistribute or PressureExternal\n",
                                   numConnect);
                            status = CAPS_NOTIMPLEMENT;
                            goto cleanup;
                        }
                    }

                    if (feaLoad->loadType == PressureDistribute) {

                        hsmMemory->pars[k*LVTOT+lvqn] = feaLoad->pressureDistributeForce[m];
                    }

                    if (feaLoad->loadType == PressureExternal) {
                        hsmMemory->pars[k*LVTOT+lvqn] = feaLoad->pressureMultiDistributeForce[4*i + m];
                    }
                }
            }
        } else {
            printf("Error: Unsupported load type - %d!\n", feaLoad->loadType);
            status = CAPS_NOTIMPLEMENT;
            goto cleanup;
        }

    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in hsm_setSurfaceParameter status = %d\n",
               status);

    return status;
}


// Set BC values on hsm-Edges (mesh segments) parameters in hsmMemory->pare structure
int hsm_setEdgeBCParameter(feaProblemStruct feaProblem, int permutation[],
                           hsmMemoryStruct *hsmMemory)
{

    int status; // Function return status
    int i, j, k, m, loadIndex; // Indexing

    int numConnect, numBCEdge = 0;
    double vector[3];
    feaLoadStruct *feaLoad;

    for (loadIndex = 0; loadIndex < feaProblem.numLoad; loadIndex++) {

        feaLoad = &feaProblem.feaLoad[loadIndex];

        // Edge Force and Moment loads
        if (feaLoad->loadType == LineForce ||
            feaLoad->loadType == LineMoment) {

            if (feaLoad->elementIDSet == NULL) {
                printf("Error: NULL elementIDSet!\n");
                status = CAPS_NULLVALUE;
                goto cleanup;
            }

            vector[0] = feaLoad->directionVector[0];
            vector[1] = feaLoad->directionVector[1];
            vector[2] = feaLoad->directionVector[2];

            for (i = 0; i < feaLoad->numElementID; i++) {

                 for (j = 0; j < feaProblem.feaMesh.numElement; j++) {
                     if (feaProblem.feaMesh.element[j].elementID == feaLoad->elementIDSet[i] ) break;
                 }

                 // self consistency check
                 if (feaProblem.feaMesh.element[j].elementType != Line) {
                     printf("Error: Edge Force/Moment applied to a non-Line element!\n");
                     status = CAPS_BADVALUE;
                     goto cleanup;
                 }

                 k = numBCEdge;
                 numBCEdge++;

                 numConnect = 2;
                 for (m = 0; m < numConnect; m++) {

                     // Get index in hsmMemory
                     hsmMemory->kbcedge[2*k + m] =
                         permutation[feaProblem.feaMesh.element[j].connectivity[m]-1];

                     if (feaLoad->loadType == LineForce) {

                         // Loading - f1  force/length  vector in xyz axes
                         hsmMemory->pare[k*LETOT+lef1x] = feaLoad->forceScaleFactor*vector[0]; // x
                         hsmMemory->pare[k*LETOT+lef1y] = feaLoad->forceScaleFactor*vector[1]; // y
                         hsmMemory->pare[k*LETOT+lef1z] = feaLoad->forceScaleFactor*vector[2]; // z

                         //           f2  force/length  vector in xyz axes
                         hsmMemory->pare[k*LETOT+lef2x] = feaLoad->forceScaleFactor*vector[0]; // x
                         hsmMemory->pare[k*LETOT+lef2y] = feaLoad->forceScaleFactor*vector[1]; // y
                         hsmMemory->pare[k*LETOT+lef2z] = feaLoad->forceScaleFactor*vector[2]; // z
                     }

                     if (feaLoad->loadType == LineMoment) {

                         // Load - m1  moment/length vector in xyz axes
                         hsmMemory->pare[k*LETOT+lem1x] = feaLoad->momentScaleFactor*vector[0]; // x
                         hsmMemory->pare[k*LETOT+lem1y] = feaLoad->momentScaleFactor*vector[1]; // y
                         hsmMemory->pare[k*LETOT+lem1z] = feaLoad->momentScaleFactor*vector[2]; // z

                         //        m2  moment/length vector in xyz axes
                         hsmMemory->pare[k*LETOT+lem2x] = feaLoad->momentScaleFactor*vector[0]; // x
                         hsmMemory->pare[k*LETOT+lem2y] = feaLoad->momentScaleFactor*vector[1]; // y
                         hsmMemory->pare[k*LETOT+lem2z] = feaLoad->momentScaleFactor*vector[2]; // z
                     }

                 } // End gridIDSet
            }
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in hsm_setEdgeBCParameter status = %d\n",
               status);

    return status;
}


// Set BC values on hsm-nodes (vertexes) in hsmMemory->parp structure
int hsm_setNodeBCParameter(feaProblemStruct feaProblem, int permutation[],
                           hsmMemoryStruct *hsmMemory)
{

    int status; // Function return status

    int i, j, k, m, constraintIndex, loadIndex; // Indexing

    int incrementFlag = (int) false;

    char numString[10];
    double vector[3];
    feaLoadStruct *feaLoad;
    feaConstraintStruct *feaConstraint;

    printf("NEED TO ADD MORE CODE TO hsm_setNodeBCParameter\n");

    hsmMemory->numBCNode = 0;

    // parp node forces and moment constraints
    for (loadIndex = 0; loadIndex < feaProblem.numLoad; loadIndex++) {

        feaLoad = &feaProblem.feaLoad[loadIndex];

        // Grid Force and Moment loads
        if (feaLoad->loadType == GridForce ||
            feaLoad->loadType == GridMoment) {

          if (feaLoad->gridIDSet == NULL) {
              printf("Error: NULL gridIDSet!\n");
              status = CAPS_NULLVALUE;
              goto cleanup;
          }

          vector[0] = feaLoad->directionVector[0];
          vector[1] = feaLoad->directionVector[1];
          vector[2] = feaLoad->directionVector[2];

          for (i = 0; i < feaLoad->numGridID; i++) {

              // Get index in hsmMemory
              k =  permutation[feaLoad->gridIDSet[i]-1];

              j = -1;
              if (hsmMemory->kbcnode[hsmMemory->numBCNode] == 0) {
                  j = hsmMemory->numBCNode;
              } else {
                  for (m = 0 ; m < hsmMemory->numBCNode; m++) {
                      if (hsmMemory->kbcnode[m] == k) {
                          j = m;
                          break;
                      }
                  }
              }

              // Index
              hsmMemory->kbcnode[j] = k;

              if (feaLoad->loadType == GridForce) {
                  hsmMemory->lbcnode[j] += lbcf;

                  // Loading - fixed-direction load
                  hsmMemory->parp[j*LPTOT+lpfx] = feaLoad->forceScaleFactor*vector[0]; // x
                  hsmMemory->parp[j*LPTOT+lpfy] = feaLoad->forceScaleFactor*vector[1]; // y
                  hsmMemory->parp[j*LPTOT+lpfz] = feaLoad->forceScaleFactor*vector[2]; // z
              }

              if (feaLoad->loadType == GridMoment) {
                  hsmMemory->lbcnode[j] += lbcm;

                  // Load - fixed-direction moment
                  hsmMemory->parp[j*LPTOT+lpmx] = feaLoad->momentScaleFactor*vector[0]; // x
                  hsmMemory->parp[j*LPTOT+lpmy] = feaLoad->momentScaleFactor*vector[1]; // y
                  hsmMemory->parp[j*LPTOT+lpmz] = feaLoad->momentScaleFactor*vector[2]; // z
              }
          }
       }
    }


    // parp node displacement constraints
    for (constraintIndex = 0; constraintIndex < feaProblem.numConstraint; constraintIndex++) {

        feaConstraint = &feaProblem.feaConstraint[constraintIndex];
        // Sort types { Displacement, ZeroDisplacement}

        // Zero -Displacement and Displacement
        if (feaConstraint->constraintType == ZeroDisplacement) {// ||
            //) feaConstraint->constraintType == Displacement) {

            for (i = 0; i < feaConstraint->numGridID; i++) {

                incrementFlag = (int) false;

                // Get index in hsmMemory
                k = permutation[feaConstraint->gridIDSet[i]-1];

                j = -1;
                if (hsmMemory->kbcnode[hsmMemory->numBCNode] == 0) {
                    j = hsmMemory->numBCNode;
                    incrementFlag = (int) true;
                } else {
                    for (m = 0 ; m < hsmMemory->numBCNode; m++) {
                        if (hsmMemory->kbcnode[m] == k) {
                            j = m;
                            break;
                        }
                    }
                }

                sprintf(numString, "%d", feaConstraint->dofConstraint);

                // Index
                hsmMemory->kbcnode[j] = k;
                k--; // convert to zero-based index

                // Boundary value 123456
                if (strstr(numString, "1") != NULL &&
                    strstr(numString, "2") != NULL &&
                    strstr(numString, "3") != NULL &&
                    strstr(numString, "4") != NULL &&
                    strstr(numString, "5") != NULL &&
                    strstr(numString, "6") != NULL ) {

                    hsmMemory->lbcnode[j] += lbcr3 + lbcd3;

                } else if (strstr(numString, "1") != NULL &&
                           strstr(numString, "2") != NULL &&
                           strstr(numString, "3") != NULL) {

                    hsmMemory->lbcnode[j] += lbcr3;

                } else {
                    printf("Error: DOF constraint %d, not supported yet\n",
                           feaConstraint->dofConstraint);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                if (feaConstraint->constraintType == ZeroDisplacement) {

                    // Fix position to undeformed coordinates
                    hsmMemory->parp[j*LPTOT+lprx] = hsmMemory->pars[k*LVTOT+lvr0x];
                    hsmMemory->parp[j*LPTOT+lpry] = hsmMemory->pars[k*LVTOT+lvr0y];
                    hsmMemory->parp[j*LPTOT+lprz] = hsmMemory->pars[k*LVTOT+lvr0z];

                    hsmMemory->parp[j*LPTOT+lpt1x] = hsmMemory->pars[k*LVTOT+lve01x];
                    hsmMemory->parp[j*LPTOT+lpt1y] = hsmMemory->pars[k*LVTOT+lve01y];
                    hsmMemory->parp[j*LPTOT+lpt1z] = hsmMemory->pars[k*LVTOT+lve01z];

                    hsmMemory->parp[j*LPTOT+lpt2x] = hsmMemory->pars[k*LVTOT+lve02x];
                    hsmMemory->parp[j*LPTOT+lpt2y] = hsmMemory->pars[k*LVTOT+lve02y];
                    hsmMemory->parp[j*LPTOT+lpt2z] = hsmMemory->pars[k*LVTOT+lve02z];

                } else {
                    printf("Error: Unsupported constraint type - %d!\n",
                           feaConstraint->constraintType);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                if (incrementFlag == (int) true) hsmMemory->numBCNode += 1;
            }

        } else {
            printf("Error: Unsupported constraint type - %d!\n",
                   feaConstraint->constraintType);
            status = CAPS_NOTIMPLEMENT;
            goto cleanup;
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in hsm_setNodeBCParameter status = %d\n",
               status);

    return status;
}
