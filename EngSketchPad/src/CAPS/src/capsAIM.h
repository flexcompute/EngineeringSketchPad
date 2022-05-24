/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Function Interface Include
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */



/* initialize & get info about the list of arguments */
extern int
aim_Initialize(aimContext *cntxt,
               const char *analysisName,
               int        *qeFlag,       /* query/execute flag (input/output) */
    /*@null@*/ const char *unitSys,      /* Unit System requested */
    /*@null@*/ void       *aimStruc,     /* the AIM context */
               int        *major,        /* the returned major version */
               int        *minor,        /* the returned minor version */
               int        *nIn,          /* returned number of inputs */
               int        *nOut,         /* returned number of outputs */
               int        *nField,       /* returned number of DataSet fields */
               char       ***fnames,     /* returned pointer to field strings */
               int        **franks,      /* returned pointer to field ranks */
               int        **fInOut,      /* returned pointer to field in/out */
               void       **instStore);  /* returned AIM storage */

/* get analysis index */
extern int
aim_Index(aimContext cntxt, const char *analysisName);

/* fill in the discrete data for a Bound */
extern int
aim_Discr(aimContext cntxt,
          const char *analysisName,
          char       *bname,            /* Bound name */
          capsDiscr  *discr);           /* the structure to fill */
  
/* frees up data in a discrete structure */
extern int
aim_FreeDiscr(aimContext cntxt,
              const char *analysisName,
              capsDiscr  *discr);       /* the structure to free up */
  
/* locate an element in the mesh */
extern int
aim_LocateElement(aimContext cntxt,
                  const char *analysisName,
                  capsDiscr  *discr,    /* the input discrete structure */
                  double     *params,   /* the input global parametric space */
                  double     *param,    /* the requested position */
                  int        *bIndex,   /* the returned body index */
                  int        *eIndex,   /* the returned element index */
                  double     *bary);    /* the barycentric coordinates */

extern int
aim_LocateElIndex(aimContext cntxt,
                  int        index,
                  capsDiscr  *discr,    /* the input discrete structure */
                  double     *params,   /* the input global parametric space */
                  double     *param,    /* the requested position */
                  int        *bIndex,   /* the returned body index */
                  int        *eIndex,   /* the returned element index */
                  double     *bary);    /* the barycentric coordinates */

/* input information for the AIM */
extern int
aim_Inputs(aimContext cntxt,
           const char *analysisName,
/*@null@*/ void       *instStore,       /* instance storage */
/*@null@*/ void       *aimStruc,        /* the AIM context */
           int        index,            /* the input index [1-nIn] */
           char       **ainame,         /* pointer to the returned name */
           capsValue  *defaultVal);     /* pointer to default value (filled) */

/* parse input data & generate input file(s) */
extern int
aim_PreAnalysis(aimContext cntxt,
                const char *analysisName,
     /*@null@*/ void       *instStore,  /* instance storage */
                void       *aimStruc,   /* the AIM context */
     /*@null@*/ capsValue  *inputs);    /* complete suite of analysis inputs */

/* execute the analysis */
extern int
aim_Execute(aimContext cntxt,
            const char *analysisName,
 /*@null@*/ void       *instStore,      /* instance storage */
            void       *aimStruc,       /* the AIM context */
            int        *state);         /* the state of the execution */

#ifdef ASYNCEXEC
/* check the analysis execution */
extern int
aim_Check(aimContext cntxt,
          const char *analysisName,
/*@null@*/void       *instStore,        /* instance storage */
          void       *aimStruc,         /* the AIM context */
          int        *state);           /* the state of the execution */
#endif

/* output information for the AIM */
extern int
aim_Outputs(aimContext cntxt,
            const char *analysisName,
 /*@null@*/ void       *instStore,      /* instance storage */
 /*@null@*/ void       *aimStruc,       /* the AIM context */
            int        index,           /* the output index [1-nOut] */
            char       **aoname,        /* pointer to the returned name */
            capsValue  *formVal);       /* pointer to form/units (filled) */

/* parse output data / file(s) */
extern int
aim_PostAnalysis(aimContext cntxt,
                 const char *analysisName,
      /*@null@*/ void       *instStore,  /* instance storage */
                 void       *aimStruc,   /* the AIM context */
                 int        restart,     /* 0 - normal, 1 - restart */
      /*@null@*/ capsValue  *inputs);    /* complete suite of analysis inputs */

/* calculate output information */
extern int
aim_CalcOutput(aimContext cntxt,
               const char *analysisName,
    /*@null@*/ void       *instStore,   /* instance storage */
               void       *aimStruc,    /* the AIM context */
               int        index,        /* the output index [1-nOut] */
               capsValue  *value);      /* pointer to value struct to fill */

/* data transfer using the discrete structure */
extern int
aim_Transfer(aimContext cntxt,
             const char *analysisName,
             capsDiscr  *discr,         /* the input discrete structure */
             const char *name,          /* the field name */
             int        npts,           /* number of points to be filled */
             int        rank,           /* the rank of the data */
             double     *data,          /* pointer to the memory to be filled */
             char       **units);       /* the units string returned */

/* interpolation */
extern int
aim_Interpolation(aimContext cntxt,
                  const char *analysisName,
                  capsDiscr  *discr,    /* the input discrete structure */
                  const char *name,     /* the dataset name */
                  int        bIndex,    /* the input discretized body (1-bias) */
                  int        eIndex,    /* the input element (1-bias) */
                  double     *bary,     /* the barycentric coordinates */
                  int        rank,      /* the rank of the data */
                  double     *data,     /* global discrete support */
                  double     *result);  /* the result (rank in length) */

extern int
aim_InterpolIndex(aimContext cntxt,
                  int        index,
                  capsDiscr  *discr,    /* the input discrete structure */
                  const char *name,     /* the dataset name */
                  int        bIndex,    /* the input discretized body (1-bias) */
                  int        eIndex,    /* the input element (1-bias) */
                  double     *bary,     /* the barycentric coordinates */
                  int        rank,      /* the rank of the data */
                  double     *data,     /* global discrete support */
                  double     *result);  /* the result (rank in length) */

/* reverse differentiated interpolation */
extern int
aim_InterpolateBar(aimContext cntxt,
                   const char *analysisName,
                   capsDiscr  *discr,   /* the input discrete structure */
                   const char *name,    /* the dataset name */
                   int        bIndex,   /* the input discretized body (1-bias) */
                   int        eIndex,   /* the input element (1-bias) */
                   double     *bary,    /* the barycentric coordinates */
                   int        rank,     /* the rank of the data */
                   double     *r_bar,   /* input d(objective)/d(result) */
                   double     *d_bar);  /* returned d(objective)/d(data) */

extern int
aim_InterpolIndBar(aimContext cntxt,
                   int       index,
                   capsDiscr  *discr,   /* the input discrete structure */
                   const char *name,    /* the dataset name */
                   int        bIndex,   /* the input discretized body (1-bias) */
                   int        eIndex,   /* the input element (1-bias) */
                   double     *bary,    /* the barycentric coordinates */
                   int        rank,     /* the rank of the data */
                   double     *r_bar,   /* input d(objective)/d(result) */
                   double     *d_bar);  /* returned d(objective)/d(data) */

/* integration */
extern int
aim_Integration(aimContext cntxt,
                const char *analysisName,
                capsDiscr  *discr,      /* the input discrete structure */
                const char *name,       /* the dataset name */
                int        bIndex,      /* the input discretized body (1-bias) */
                int        eIndex,      /* the input element (1-bias) */
                int        rank,        /* the rank of the data */
     /*@null@*/ double     *data,       /* global discrete support */
                double     *result);    /* the result (rank in length) */

extern int
aim_IntegrIndex(aimContext cntxt,
                int        index,
                capsDiscr  *discr,      /* the input discrete structure */
                const char *name,       /* the dataset name */
                int        bIndex,      /* the input discretized body (1-bias) */
                int        eIndex,      /* the input element (1-bias) */
                int        rank,        /* the rank of the data */
     /*@null@*/ double     *data,       /* global discrete support */
                double     *result);    /* the result (rank in length) */

/* reverse differentiated integration */
extern int
aim_IntegrateBar(aimContext cntxt,
                 const char *analysisName,
                 capsDiscr  *discr,     /* the input discrete structure */
                 const char *name,      /* the dataset name */
                 int        bIndex,      /* the input discretized body (1-bias) */
                 int        eIndex,     /* the input element (1-bias) */
                 int        rank,       /* the rank of the data */
                 double     *r_bar,     /* input d(objective)/d(result) */
                 double     *d_bar);    /* returned d(objective)/d(data) */

extern int
aim_IntegrIndBar(aimContext cntxt,
                 int        index,
                 capsDiscr  *discr,     /* the input discrete structure */
                 const char *name,      /* the dataset name */
                 int        bIndex,      /* the input discretized body (1-bias) */
                 int        eIndex,     /* the input element (1-bias) */
                 int        rank,       /* the rank of the data */
                 double     *r_bar,     /* input d(objective)/d(result) */
                 double     *d_bar);    /* returned d(objective)/d(data) */

/* backdoor communication */
extern int
aim_Backdoor(aimContext cntxt,
             const char *analysisName,
  /*@null@*/ void       *instStore,     /* instance storage */
             void       *aimStruc,      /* the AIM context */
             const char *JSONin,        /* the input(s) */
             char       **JSONout);     /* the output(s) */

/* cleanup instance storage */
extern void
aim_cleanup(aimContext cntxt,
            const char *analysisName,
            void       *instStore);     /* instance storage */

/* unload the AIM */
extern void
aim_cleanupAll(aimContext *cntxt);
