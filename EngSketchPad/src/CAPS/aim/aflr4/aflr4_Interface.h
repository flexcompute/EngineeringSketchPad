#include "meshTypes.h"  // Bring in mesh structurs
#include "capsTypes.h"  // Bring in CAPS types
#include "cfdTypes.h"   // Bring in cfd structures

// AFLR4 interface function - Modified from functions provided with
//	AFLR4_LIB source (aflr4_cad_geom_setup.c) written by David L. Marcum


enum aimInputs
{
  Proj_Name = 1,               /* index is 1-based */
  Mesh_Quiet_Flag,
  Mesh_Format,
  Mesh_ASCII_Flag,
  Mesh_Gen_Input_String,
  Ff_cdfr,
  Min_ncell,
  Mer_all,
  No_prox,
  Abs_min_scale,
  BL_Thickness,
  RE_l,
  Curv_factor,
  Erw_all,
  Max_scale,
  Min_scale,
  Mesh_Length_Factor,
  Mesh_Sizing,
  Multiple_Mesh,
  EGADS_Quad,
  NUMINPUT = EGADS_Quad        /* Total number of inputs */
};

#define AFLR4TESSFILE "aflr4_%d.eto"

extern int aflr4_Surface_Mesh(void *aimInfo,
                              int quiet,
                              int numBody, ego *bodies,
                              capsValue *aimInputs,
                              meshInputStruct meshInput);
