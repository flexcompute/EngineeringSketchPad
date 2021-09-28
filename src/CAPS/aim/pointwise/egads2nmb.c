#include "egads.h"
#include "nmb.h"
#include "string.h"


static void showUsage(const char *exe)
{
    const char *fmt =
        "Usage: %s [options] EGADSfileName [NMBfileName]\n"
        "\n"
        "  Converts the EGADS file EGADSfileName to an NMB file written\n"
        "  to NMBfileName.\n"
        "\n"
        "  If NMBfileName is not specified, a default name is constructed\n"
        "  from EGADSfileName with the extension changed to 'format.nmb'.\n"
        "  For example, file.ext is saved as file.a.nmb for ASCII and\n"
        "  file.b.nmb for BINARY.\n"
        "\n"
        "  options:\n"
        "    -b|--binary  Export nmb file in BINARY format (the default).\n"
        "    -a|--ascii   Export nmb file in ASCII format.\n"
        "    -v|--verbose Dump detailed conversion information.\n"
        "    --modelsize <value> Define GE ModelSize.\n"
        "    -h|--help    Display this help and stop.\n"
        "";
    printf(fmt, exe);
}
#define SHOWUSAGE showUsage(argv[0])



int main(int argc, char **argv)
{
  int        ii, maj, min, stat;
  char       *dot;
  const char *OCCrev;
  ego        context, model;
  const char *egadsFilename = 0;
  const char *nmbFilename   = 0;
  int        asciiOut       = 0;
  int        verbose        = 0;
  int        help           = 0;
  float      modelSize      = 1000.0;

  for (ii = 1; ii < argc; ++ii) {
    if (0 == strcmp("-a", argv[ii]) || 0 == strcmp("--ascii", argv[ii])) {
      asciiOut = 1;
    }
    else if (0 == strcmp("-b", argv[ii]) || 0 == strcmp("--binary", argv[ii])) {
      asciiOut = 0;
    }
    else if (0 == strcmp("-v", argv[ii]) || 0 == strcmp("--verbose", argv[ii])) {
      verbose  = 1;
    }
    else if (0 == strcmp("--modelsize", argv[ii])) {
      ++ii;
      if (ii >= argc) {
        printf("Error: --modelsize switch requires value\n");
        SHOWUSAGE;
        return 1;
      }
      if (1 != sscanf(argv[ii], "%f", &modelSize)) {
        printf("Error: --modelsize switch requires value\n");
        SHOWUSAGE;
        return 1;
      }
    }
    else if (0 == strcmp("-h", argv[ii]) || 0 == strcmp("--help", argv[ii])) {
      help = 1;
    }
    else if ('-' == argv[ii][0]) {
      printf("Error: Unknown switch %s\n", argv[ii]);
      SHOWUSAGE;
      return 1;
    }
    else if (0 == egadsFilename) {
      // first non-switch is EGADS input file name
      egadsFilename = argv[ii];
    }
    else if (0 == nmbFilename) {
      // second non-switch is NMB output file name
      nmbFilename = argv[ii];
    }
    else {
      printf("Error: Unexpected filename %s\n", argv[ii]);
      SHOWUSAGE;
      return 1;
    }
  }

  if (help) {
      SHOWUSAGE;
      return 0;
  }

  if (0 == egadsFilename) {
      printf("Error: EGADSfileName required\n");
      SHOWUSAGE;
      return 1;
  }

  if (0 == nmbFilename) {
      // The nmb filename not provided. Use a default name:
      //    filename.ext --> filename.nmb
      #define sz  2048
      static char defaultNmbFilename[sz] = "";
      strncpy(defaultNmbFilename, egadsFilename, sz - 8);
      dot = strrchr(defaultNmbFilename, '.');
      if (0 != dot) {
          *dot = '\0';
      }
      strcat(defaultNmbFilename, (asciiOut ? ".a" : ".b"));
      strcat(defaultNmbFilename, ".nmb");
      nmbFilename = defaultNmbFilename;
  }
  if (verbose) {
      printf("Input : '%s'\n", egadsFilename);
      printf("Output: '%s' (%s)\n", nmbFilename,
                (asciiOut ? "ASCII" : "Binary"));
  }

  /* look at EGADS revision, open and load Model */
  EG_revision(&maj, &min, &OCCrev);
  if (verbose) {
      printf("\n Using EGADS %2d.%02d %s\n\n", maj, min, OCCrev);
  }

  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open = %d\n", stat);
    return 1;
  }

  stat = EG_loadModel(context, 0, egadsFilename, &model);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_loadModel = %d\n", stat);
    EG_close(context);
    return 1;
  }
  
  stat = NMB_write(model, nmbFilename, asciiOut, verbose, modelSize);
  if (stat != EGADS_SUCCESS) {
      printf(" NMB_write = %d\n", stat);
  }

  EG_deleteObject(model);
  EG_close(context);
  
  return 0;
}
