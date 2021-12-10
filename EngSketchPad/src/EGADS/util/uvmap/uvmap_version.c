#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_version.c,v 1.37 2021/08/28 04:58:37 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

void uvmap_version (
  char Compile_Date[],
  char Compile_OS[],
  char Version_Date[],
  char Version_Number[])
{
  // Put compile date, compile OS, version date, and version number in text
  // string.

strncpy (Compile_Date, "", 40);
strncpy (Compile_OS, "", 40);
strncpy (Version_Date, "08/27/21 @ 11:58PM", 40);
strncpy (Version_Number, "1.7.1", 40);

  strcpy (&(Compile_Date[40]), "\0");
  strcpy (&(Compile_OS[40]), "\0");
  strcpy (&(Version_Date[40]), "\0");
  strcpy (&(Version_Number[40]), "\0");

  return;

}
