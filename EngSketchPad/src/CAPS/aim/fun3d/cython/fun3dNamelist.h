/* Generated by Cython 3.0.11 */

#ifndef __PYX_HAVE__fun3dNamelist
#define __PYX_HAVE__fun3dNamelist

#include "Python.h"
#include "egadsTypes.h"
#include "egads.h"
#include "egadsErrors.h"
#include "prm.h"
#include "capsTypes.h"
#include "caps.h"
#include "capsErrors.h"
#include "cfdTypes.h"
#include "aimUtil.h"

    #if defined(WIN32)
      #define PATH_MAX _MAX_PATH
    #else
      #include <limits.h>
    #endif

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#ifndef __PYX_HAVE_API__fun3dNamelist

#ifdef CYTHON_EXTERN_C
    #undef __PYX_EXTERN_C
    #define __PYX_EXTERN_C CYTHON_EXTERN_C
#elif defined(__PYX_EXTERN_C)
    #ifdef _MSC_VER
    #pragma message ("Please do not define the '__PYX_EXTERN_C' macro externally. Use 'CYTHON_EXTERN_C' instead.")
    #else
    #warning Please do not define the '__PYX_EXTERN_C' macro externally. Use 'CYTHON_EXTERN_C' instead.
    #endif
#else
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(_T) _T
#endif

__PYX_EXTERN_C int fun3d_writeNMLPython(void *, capsValue *, cfdBoundaryConditionStruct);

#endif /* !__PYX_HAVE_API__fun3dNamelist */

/* WARNING: the interface of the module init function changed in CPython 3.5. */
/* It now returns a PyModuleDef instance instead of a PyModule instance. */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initfun3dNamelist(void);
#else
/* WARNING: Use PyImport_AppendInittab("fun3dNamelist", PyInit_fun3dNamelist) instead of calling PyInit_fun3dNamelist directly from Python 3.5 */
PyMODINIT_FUNC PyInit_fun3dNamelist(void);

#if 0 && PY_VERSION_HEX >= 0x03050000 && (defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER) || (defined(__cplusplus) && __cplusplus >= 201402L))
#if defined(__cplusplus) && __cplusplus >= 201402L
[[deprecated("Use PyImport_AppendInittab(\"fun3dNamelist\", PyInit_fun3dNamelist) instead of calling PyInit_fun3dNamelist directly.")]] inline
#elif defined(__GNUC__) || defined(__clang__)
__attribute__ ((__deprecated__("Use PyImport_AppendInittab(\"fun3dNamelist\", PyInit_fun3dNamelist) instead of calling PyInit_fun3dNamelist directly."), __unused__)) __inline__
#elif defined(_MSC_VER)
__declspec(deprecated("Use PyImport_AppendInittab(\"fun3dNamelist\", PyInit_fun3dNamelist) instead of calling PyInit_fun3dNamelist directly.")) __inline
#endif
static PyObject* __PYX_WARN_IF_PyInit_fun3dNamelist_INIT_CALLED(PyObject* res) {
  return res;
}
#define PyInit_fun3dNamelist() __PYX_WARN_IF_PyInit_fun3dNamelist_INIT_CALLED(PyInit_fun3dNamelist())
#endif
#endif

#endif /* !__PYX_HAVE__fun3dNamelist */
