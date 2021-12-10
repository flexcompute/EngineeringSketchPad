#ifdef __cplusplus
extern "C" {
#endif

int tetgen_VolumeMesh(void *aimInfo,
                      meshInputStruct meshInput,
                      const char *fileName,
                      meshStruct *surfaceMesh,
                      meshStruct *volumeMesh);

#ifdef __cplusplus
}
#endif
