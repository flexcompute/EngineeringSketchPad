#include <ug/UG_LIB.h>
#include <aflr3/aflr3_version.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#define strtok_r   strtok_s
#endif

// the minimum version of AFLR3 API that aflr3AIM is designed for
const int AFLR3_MIN_VERSION[3] = { 16,31,5 };

int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  CHAR_133 Compile_Date;
  CHAR_133 Compile_OS;
  CHAR_133 Version_Date;
  CHAR_133 Version_Number;

  char *rest = NULL, *token = NULL;

  char *AFLR = NULL;
  int i = 0, status = 0;

  // get the version
  aflr3_version(Compile_Date, Compile_OS, Version_Date, Version_Number);

  printf("%s\n", Version_Number);

  rest = Version_Number;
  while ((token = strtok_r(rest, ".", &rest))) {
    if (i > 2) {
      aflr3_version(Compile_Date, Compile_OS, Version_Date, Version_Number);
      printf("error: AFLR3 version number %s has more than 3 integers. Please fix aflr3_version.c\n", Version_Number );
      return 1;
    }
    if ( atoi(token) < AFLR3_MIN_VERSION[i] ) {
      status = 1;
      break;
    } else if ( atoi(token) > AFLR3_MIN_VERSION[i] ) {
      break;
    }
    i++;
  }

  if (status != 0) {
    // get the version again as strtok corrupts it
    aflr3_version(Compile_Date, Compile_OS, Version_Date, Version_Number);
    AFLR = getenv("AFLR");

    printf("\n");
    if ( AFLR != NULL )
      printf("Using AFLR: %s\n", AFLR);
    printf("error: AFLR3 version number %s is less than %d.%d.%d\n",
           Version_Number,
           AFLR3_MIN_VERSION[0],
           AFLR3_MIN_VERSION[1],
           AFLR3_MIN_VERSION[2]);
    printf("\n");
  }

  return status;
}
