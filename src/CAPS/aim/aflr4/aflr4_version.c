#include <ug/UG_LIB.h>
#include <aflr4/aflr4_version.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#define strtok_r   strtok_s
#endif

// the minimum version of AFLR4 API that aflr4AIM is designed for
const int AFLR4_MIN_VERSION[3] = { 10,17,5 };

int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  CHAR_133 Compile_Date;
  CHAR_133 Compile_OS;
  CHAR_133 Version_Date;
  CHAR_133 Version_Number;

  char *rest = NULL, *token = NULL;

  char *AFLR4 = NULL;

  int i = 0, status = 0;

  // get the version
  aflr4_version(Compile_Date, Compile_OS, Version_Date, Version_Number);

  rest = Version_Number;
  while ((token = strtok_r(rest, ".", &rest))) {
    if (i > 2) {
      aflr4_version(Compile_Date, Compile_OS, Version_Date, Version_Number);
      printf("error: AFLR4 version number %s has more than 3 integers. Please fix aflr4_version.c\n", Version_Number );
      return 1;
    }
    if ( atoi(token) < AFLR4_MIN_VERSION[i] ) {
      status = 1;
      break;
    } else if ( atoi(token) > AFLR4_MIN_VERSION[i] ) {
      break;
    }
    i++;
  }

  if (status != 0) {
    // get the version again as strtok corrupts it
    aflr4_version(Compile_Date, Compile_OS, Version_Date, Version_Number);
    AFLR4 = getenv("AFLR4");

    printf("\n");
    if ( AFLR4 != NULL )
      printf("Using AFLR4: %s\n", AFLR4);
    printf("error: AFLR4 version number %s is less than %d.%d.%d\n",
           Version_Number,
           AFLR4_MIN_VERSION[0],
           AFLR4_MIN_VERSION[1],
           AFLR4_MIN_VERSION[2]);
    printf("\n");
    return 1;
  }

  return status;
}
