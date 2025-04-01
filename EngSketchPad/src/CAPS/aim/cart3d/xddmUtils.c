
#include "xddm.h"
#include "xddmUtils.h"

#include <stdio.h>

/* deprecated functions from old xddm library */
int
xddm_addAeroFunForce(p_tsXddmAFun p_afun,
                     const char *const name,
                     const int force,
                     const int frame,
                     const int J,
                     const int N,
                     const double target,
                     const double weight,
                     const int bnd,
                     const char *const comp)
{
  char tmp[512];
  const char *header = "\n"
                       "#         Name    Force   Frame    J      N    Target   Weight  Bound  GMP Comp\n"
                       "#        (String) (0,1,2) (0,1) (0,1,2) (int)  (dble)   (dble)   (0)\n"
                       "#------------------------------------------------------------------------------\n";
  const char *cmp = comp;

  if (p_afun->p_text == NULL) {
    /* add the header information */
    p_afun->p_text = malloc((strlen(header)+1)*sizeof(char));
    if (p_afun->p_text == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
    strcpy(p_afun->p_text, header);
  }

  if (comp == NULL) {
    cmp = "entire";
  }

  snprintf(tmp, 512, "optForce   %7s   %d      %d      %d      %d      %6g  %6g   %d    %s\n",
                     name, force, frame, J, N, target, weight, bnd, cmp);

  p_afun->p_text = realloc(p_afun->p_text, (strlen(p_afun->p_text) + strlen(tmp) + 1)*sizeof(char));
  if (p_afun->p_text == NULL) {
    ERR("malloc failed\n");
    return -1;
  }

  strcat(p_afun->p_text, tmp);


  return 0;
}

int
xddm_addAeroFunMoment_Point(p_tsXddmAFun p_afun,
                            const char *const name,
                            const int index,
                            const int moment,
                            const int frame,
                            const int J,
                            const int N,
                            const double target,
                            const double weight,
                            const int bnd,
                            const char *const comp)
{
  char tmp[512];
  const char *header =  "\n"
                        "#                  Name   Index  Moment  Frame   J     N   Target  Weight  Bound  GMP_Comp\n"
                        "#                (String) (int) (0,1,2)  (0,1) (0,1) (int) (dble)  (dble)  (0)\n"
                        "#---------------------------------------------------------------------------------------\n";

  const char *cmp = comp;

  if (p_afun->p_text == NULL) {
    /* add the header information */
    p_afun->p_text = malloc((strlen(header)+1)*sizeof(char));
    if (p_afun->p_text == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
    strcpy(p_afun->p_text, header);
  }

  if (comp == NULL) {
    cmp = "entire";
  }

  snprintf(tmp, 512, "optMoment_Point %7s    %d      %d        %d      %d    %d  %6g %6g     %d    %s\n",
                     name, index, moment, frame, J, N, target, weight, bnd, cmp);

  p_afun->p_text = realloc(p_afun->p_text, (strlen(p_afun->p_text) + strlen(tmp) + 1)*sizeof(char));
  if (p_afun->p_text == NULL) {
    ERR("malloc failed\n");
    return -1;
  }

  strcat(p_afun->p_text, tmp);

  return 0;
}


int
xddm_addAeroFunLoD(p_tsXddmAFun p_afun,
                   const char *const name,
                   const int frame,
                   const int J,
                   const int N,
                   const double A,
                   const double bias,
                   const double target,
                   const double weight,
                   const int bnd,
                   const char *const comp)
{
  char tmp[512];
  const char *header = "\n"
                       "# L/D -> SIGN(CL)*ABS(CL)^A/(CD+Bias) in Aero Frame\n"
                       "#     -> SIGN(CN)*ABS(CN)^A/(CA+Bias) in Body Frame\n"
                       "# Format:\n"
                       "#      Name   Frame   J     N     A     Bias  Target  Weight  Bound  GMP_Comp\n"
                       "#    (String) (0,1) (0,1) (int) (dble) (dble) (dble)  (dble)   (0)\n"
                       "#----------------------------------------------------------------------------\n";
  const char *cmp = comp;

  if (p_afun->p_text == NULL) {
    /* add the header information */
    p_afun->p_text = malloc((strlen(header)+1)*sizeof(char));
    if (p_afun->p_text == NULL) {
      ERR("malloc failed\n");
      return -1;
    }
    strcpy(p_afun->p_text, header);
  }

  if (comp == NULL) {
    cmp = "entire";
  }

  snprintf(tmp, 512, "optLD  %7s   %d      %d   %d   %6g  %6g   %6g    %6g   %d    %s\n",
                     name, frame, J, N, A, bias, target, weight, bnd, cmp);

  p_afun->p_text = realloc(p_afun->p_text, (strlen(p_afun->p_text) + strlen(tmp) + 1)*sizeof(char));
  if (p_afun->p_text == NULL) {
    ERR("malloc failed\n");
    return -1;
  }

  strcat(p_afun->p_text, tmp);

  return 0;
}


