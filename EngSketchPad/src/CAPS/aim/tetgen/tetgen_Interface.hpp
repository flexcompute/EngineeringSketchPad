#ifdef __cplusplus
extern "C" {
#endif

int tetgen_VolumeMesh(void *aimInfo,
                      meshInputStruct meshInput,
                      const mapAttrToIndexStruct *groupMap,
                      const char *fileName,
                      const int numSurfMesh,
                      meshStruct *surfaceMesh,
                      meshStruct *volumeMesh);

#ifdef __cplusplus
}
#endif
