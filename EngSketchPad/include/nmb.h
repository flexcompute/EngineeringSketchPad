#ifndef NMB_H
#define NMB_H

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif

__ProtoExt__ int  NMB_write( const ego model,      const char *name,
                             const int asciiOut,   const int  verbose,
                             const float modelSize );

#ifdef __cplusplus
}
#endif

#endif
