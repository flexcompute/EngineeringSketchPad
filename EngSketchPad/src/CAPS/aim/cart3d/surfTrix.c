#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* libCart3d defaults to WORD_BIT 32 unless told otherwise
 * this makes that assumption explicit and suppresses the warning
 */

#define WORD_BIT 32

#include <c3dio_lib.h>

int writeSurfTrix(const p_tsTriangulation p_config, const int nComps,
                  const char *const p_fileName, const int options)
{
  return io_writeSurfTrix(p_config, nComps, p_fileName, options);
}


int readSurfTrix(const char *const p_fileName, p_tsTriangulation *pp_config,
                 int *const p_nComps, const char *const p_compName,
                 const char *const p_vertDataNames,
                 const char *const p_triDataNames,
                 const int options)
{
  return io_readSurfTrix(p_fileName, pp_config, p_nComps, p_compName,
                         p_vertDataNames, p_triDataNames,
                         options);
}
