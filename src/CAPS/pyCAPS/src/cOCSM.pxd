#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

cimport cEGADS

from libc.stdio cimport FILE

cdef extern from "OpenCSM.h":

    # Structures
    
    # "Node" is a 0-D topological entity in a Body 
    ctypedef struct node_T:
        int           nedge;  # number of indicent Edges 
        double        x;      # x-coordinate 
        double        y;      # y-coordinate 
        double        z;      # z-coordinate 
        int           ibody;  # Body index (1-nbody) 
        #grat_T        gratt;  # GRatt of the Node 
        double        *dxyz;  # tessellation velocity (or NULL) 
        cEGADS.ego           enode;  # EGADS node object 
   
    
    # "Edge" is a 1-D topological entity in a Body 
    ctypedef struct edge_T:
        int           itype;  # Edge type 
        int           ibeg;   # Node at beginning 
        int           iend;   # Node at end 
        int           ileft;  # Face on the left 
        int           irite;  # Face on the rite 
        int           nface;  # number of incident Faces 
        int           ibody;  # Body index (1-nbody) 
        int           iford;  # face-order 
        int           imark;  # value of "mark" Attribute (or -1) 
        #grat_T        gratt;  # GRatt of the Edge 
        double        *dxyz;  # tessellation velocity (or NULL) 
        double        *dt;    # parametric   velocity (or NULL) 
        int           globid; # global ID (bias-1) 
        cEGADS.ego    eedge;  # EGADS edge object 
   
    
    # "Face" is a 2-D topological entity in a Body 
    ctypedef struct face_T:
        int           ibody;  # Body index (1-nbody) 
        int           iford;  # face-order 
        int           imark;  # value of "mark" Attribute (or -1) 
        #grat_T        gratt;  # GRatt of the Face 
        void          *eggdata; # pointer to external grid generator data 
        double        *dxyz;  # tessellation velocity (or NULL) 
        double        *duv;   # parametric   velocity (or NULL) 
        int           globid; # global ID (bias-1) 
        cEGADS.ego     eface;  # EGADS face object 
    

    ctypedef struct body_T:
        int           ibrch;   # Branch associated with Body 
        int           brtype;  # Branch type (see below) 
        int           ileft;   # left parent Body (or 0) 
        int           irite;   # rite parent Body (or 0) 
        int           ichld;   # child Body (or 0 for root) 
        int           igroup;  # Group number 
       # varg_T        arg[10]; # array of evaluated arguments (actually use 1-9) 
    
        cEGADS.ego           ebody;   # EGADS Body         object(s) 
        cEGADS.ego           etess;   # EGADS Tessellation object(s) 
        int           npnts;   # total number of unique points 
        int           ntris;   # total number of triangles 
    
        int           onstack; # =1 if on stack (and returned); =0 otherwise 
        int           hasdots; # =1 if an argument has a dot; =2 if UDPARG is changed; =0 otherwise 
        int           botype;  # Body type (see below) 
        int           nnode;   # number of Nodes 
        node_T        *node;   # array  of Nodes 
        int           nedge;   # number of Edges 
        edge_T        *edge;   # array  of Edges 
        int           nface;   # number of Faces 
        face_T        *face;   # array  of Faces 
        void          *cache;  # (blind) structure for caching sensitivity info 
       # grat_T        gratt;   # GRatt of the Nodes 
  
    ctypedef struct modl_T:
    
        int           nattr # number of global Attributes
    
        int           nbrch # number of Branches
        brch_T        *brch # array  of Branches

        int           npmtr # number of Parameters
        pmtr_T        *pmtr # array  of Parameters

        
        int           nbody # number of Bodys
        body_T        *body # array  of Bodys
        
        modl_T *perturb  # model of perturbed body for sensitivty 
        modl_T *basemodl # base MODL while creating perturbation 
        double dtime     # time step in sensitivity 
                                #  0.001 = initial value 
                                # -2     = problem with previous attempt to create perturb 
                                
    ctypedef struct brch_T:                     
        char          *filename #filename where Branch is defined 
        int           linenum  # line number in file where Branch is defined
    
    ctypedef struct pmtr_T:
        char          *name # name of Parameter
        int           type  # Parameter type (see below)
        int           scope # associated scope (nominally 0)
        int           nrow  # number of rows    (=0 for string)
        int           ncol  # number of columns (=0 for string)
        double        *value# current value(s)
        double        *dot  # current velocity(s)
        double        *lbnd # lower Bound(s)
        double        *ubnd # upper Bound(s)
        char          *str  # string value

                                
    
    #************************************************************************
    #*                                                                      *
    #* Callable routines                                                    *
    #*                                                                      *
    #************************************************************************
     
    # return current version 
    cdef int ocsmVersion(int   *imajor,          # (out) major version number 
                         int   *iminor)          # (out) minor version number 
    
    # set output level 
    cdef int ocsmSetOutLevel(int    ilevel)      # (in)  output level: 
                                                 #       =0 warnings and errors only 
                                                 #       =1 nominal (default) 
                                                 #       =2 debug 
    
    # create a MODL by reading a .csm file 
    cdef int ocsmLoad(char   filename[],         # (in)  file to be read (with .csm) 
                      void   **modl)             # (out) pointer to MODL 
    
    # load dictionary from dictname 
    cdef int ocsmLoadDict(void   *modl,          # (in)  pointer to MODL
                          char   dictname[])     # (in)  file that contains dictionary

    # get a list of all .csm, .cpc. and .udc files
    cdef int ocsmGetFilelist(void   *modl,       # (in)  pointer to MODL 
                             char   *filelist[]) # (out) bar-sepatared list of files 
                                                 #      must be freed by user
                                        
    # save a MODL to a file 
    cdef int ocsmSave(void   *modl,              # (in)  pointer to MODL 
                      char   filename[])         # (in)  file to be written (with extension) 
                                                 #       .csm -> write outer .csm file 
                                                 #       .cpc -> write checkpointed .csm file 
                                                 #       .udc -> write a .udc file 
    
    # copy a MODL 
    cdef int ocsmCopy(void   *srcModl,           # (in)  pointer to source MODL 
                      void   **newModl)          # (out) pointer to new    MODL 
    
    # free up all storage associated with a MODL 
    cdef int ocsmFree(void   *modl)              # (in)  pointer to MODL (or NULL) 
    
    # get info about a MODL 
    cdef int ocsmInfo(void   *modl,              # (in)  pointer to MODL 
                      int    *nbrch,             # (out) number of Branches 
                      int    *npmtr,             # (out) number of Parameters 
                      int    *nbody)             # (out) number of Bodys 
    
    
    # check that Branches are properly ordered 
    cdef int ocsmCheck(void   *modl)            # (in)  pointer to MODL 
    
    # build Bodys by executing the MODL up to a given Branch 
    cdef int ocsmBuild(void   *modl,             # (in)  pointer to MODL 
                       int    buildTo,           # (in)  last Branch to execute (or 0 for all, or -1 for no recycling) 
                       int    *builtTo,          # (out) last Branch executed successfully 
                       int    *nbody,            # (in)  number of entries allocated in body[]  (out) number of Bodys on the stack 
                       int    body[])            # (out) array  of Bodys on the stack (LIFO) (at least nbody long) 
    
    # create a perturbed MODL 
    cdef int ocsmPerturb(void   *modl,           # (in)  pointer to MODL 
                         int    npmtrs,          # (in)  numner of perturbed Parameters (or 0 to remove) 
                         int    ipmtrs[],        # (in)  array of Parameter indices (1-npmtr) 
                         int    irows[],         # (in)  array of row       indices (1-nrow) 
                         int    icols[],         # (in)  array of column    indices (1-ncol) 
                         double values[])        # (in)  array of perturbed values 
    
    # create a new Branch 
    cdef int ocsmNewBrch(void   *modl,           # (in)  pointer to MODL 
                         int    iafter,          # (in)  Branch index (0-nbrch) after which to add 
                         int    type,            # (in)  Branch type (see below) 
                         char   filename[],      # (in)  filename where Branch is defined
                         #int    fileindx,        # (in)  index into filelist (bias-0) 
                         int    linenum,         # (in)  file number         (bias-1) 
                         char   arg1[],          # (in)  Argument 1 (or NULL) 
                         char   arg2[],          # (in)  Argument 2 (or NULL) 
                         char   arg3[],          # (in)  Argument 3 (or NULL) 
                         char   arg4[],          # (in)  Argument 4 (or NULL) 
                         char   arg5[],          # (in)  Argument 5 (or NULL) 
                         char   arg6[],          # (in)  Argument 6 (or NULL) 
                         char   arg7[],          # (in)  Argument 7 (or NULL) 
                         char   arg8[],          # (in)  Argument 8 (or NULL) 
                         char   arg9[])          # (in)  Argument 9 (or NULL) 
    
    # get info about a Branch 
    cdef int ocsmGetBrch(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         int    *type,           # (out) Branch type (see below) 
                         int    *bclass,         # (out) Branch class (see below) 
                         int    *actv,           # (out) Branch Activity (see below) 
                         int    *ichld,          # (out) ibrch of child (or 0 if root) 
                         int    *ileft,          # (out) ibrch of left parent (or 0) 
                         int    *irite,          # (out) ibrch of rite parent (or 0) 
                         int    *narg,           # (out) number of Arguments 
                         int    *nattr)          # (out) number of Attributes 
    
    # set activity for a Branch 
    cdef int ocsmSetBrch(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         int    actv)            # (in)  Branch activity (see below) 
    
    # delete a Branch (or whole Sketch if SKBEG) 
    cdef int ocsmDelBrch(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch)           # (in)  Branch index (1-nbrch) 
    
    # print Branches to file 
    cdef int ocsmPrintBrchs(void   *modl,        # (in)  pointer to MODL 
                            FILE   *fp)          # (in)  pointer to FILE 
    
    # get an Argument for a Branch 
    cdef int ocsmGetArg(void   *modl,            # (in)  pointer to MODL 
                        int    ibrch,            # (in)  Branch index (1-nbrch) 
                        int    iarg,             # (in)  Argument index (1-narg) 
                        char   defn[],           # (out) Argument definition (at least 129 long) 
                        double *value,           # (out) Argument value 
                        double *dot)             # (out) Argument velocity 
    
    # set an Argument for a Branch 
    cdef int ocsmSetArg(void   *modl,            # (in)  pointer to MODL 
                        int    ibrch,            # (in)  Branch index (1-nbrch) 
                        int    iarg,             # (in)  Argument index (1-narg) 
                        char   defn[])           # (in)  Argument definition 
    
    # return an Attribute for a Branch by index 
    cdef int ocsmRetAttr(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         int    iattr,           # (in)  Attribute index (1-nattr) 
                         char   aname[],         # (out) Attribute name  (at least 129 long) 
                         char   avalue[])        # (out) Attribute value (at least 129 long) 
    
    # get an Attribute for a Branch by name 
    cdef int ocsmGetAttr(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) or 0 for global 
                         char   aname[],         # (in)  Attribute name 
                         char   avalue[])        # (out) Attribute value (at least 129 long) 
    
    # set an Attribute for a Branch 
    cdef int ocsmSetAttr(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) or 0 for global 
                         char   aname[],         # (in)  Attribute name 
                         char   avalue[])        # (in)  Attribute value (or blank to delete) 
    
    # return a Csystem for a Branch by index 
    cdef int ocsmRetCsys(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         int    icsys,           # (in)  Csystem index (1-nattr) 
                         char   cname[],         # (out) Csystem name  (at least 129 long) 
                         char   cvalue[])        # (out) Csystem value (at least 129 long) 
    
    # get a Csystem for a Branch by name 
    cdef int ocsmGetCsys(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         char   cname[],         # (in)  Csystem name 
                         char   cvalue[])        # (out) Csystem value (at least 129 long) 
    
    # set a Csystem for a Branch 
    cdef int ocsmSetCsys(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch)  
                         char   cname[],         # (in)  Csystem name 
                         char   cvalue[])        # (in)  Csystem value (or blank to delete) 
    
    # print global Attributes to file 
    cdef int ocsmPrintAttrs(void   *modl,        # (in)  pointer to MODL 
                            FILE   *fp)          # (in)  pointer to FILE 
    
    # get the name of a Branch 
    cdef int ocsmGetName(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         char   name[])          # (out) Branch name (at least 129 long) 
    
    # set the name for a Branch 
    cdef int ocsmSetName(void   *modl,           # (in)  pointer to MODL 
                         int    ibrch,           # (in)  Branch index (1-nbrch) 
                         char   name[])          # (in)  Branch name 
    
    # get string data associated with a Sketch 
    cdef int ocsmGetSketch(void   *modl,         # (in)  pointer to MODL 
                           int    ibrch,         # (in)  Branch index (1-nbrch) within Sketch 
                           int    maxlen,        # (in)  length of begs, vars, cons, and segs 
                           char   begs[],        # (out) string with SKBEG info 
                                                 #       "xargxvalyargyvalzargzval" 
                           char   vars[],        # (out) string with Sketch variables 
                                                 #       "x1y1d1x2 ... dn" 
                           char   cons[],        # (out) string with Sketch constraints 
                                                 #       "type1index1_1index2_1value1 ... valuen" 
                                                 #       index1 and index2 are bias-1 
                           char   segs[])        # (out) string with Sketch segments 
                                                 #       "type1ibeg1iend1 ... iendn" 
                                                 #       ibeg and iend are bias-1 
    
    # solve for new Sketch variables 
    cdef int ocsmSolveSketch(void   *modl,       # (in)  pointer to MODL 
                             char   vars_in[],   # (in)  string with Sketch variables 
                             char   cons[],      # (in)  string with Sketch constraints 
                             char   vars_out[])  # (out) string (1024 long) with new Sketch variables 
    
    # overwrite Branches associated with a Sketch 
    cdef int ocsmSaveSketch(void   *modl,        # (in)  pointer to MODL 
                            int    ibrch,        # (in)  Branch index (1-nbrch) within Sketch 
                            char   vars[],       # (in)  string with Sketch variables 
                            char   cons[],       # (in)  string with Sketch constraints 
                            char   segs[])       # (in)  string with Sketch segments 
    
    # create a new Parameter 
    cdef int ocsmNewPmtr(void   *modl,           # (in)  pointer to MODL 
                         char   name[],          # (in)  Parameter name 
                         int    type,            # (in)  Parameter type 
                         int    nrow,            # (in)  number of rows 
                         int    ncol)            # (in)  number of columns 
    
    # delete a Parameter 
    cdef int ocsmDelPmtr(void   *modl,           # (in)  pointer to MODL 
                         int    ipmtr)           # (in)  Parameter index (1-npmtr) 
    
    # find (or create) a Parameter 
    cdef int ocsmFindPmtr(void   *modl,          # (in)  pointer to MODL 
                          char   name[],         # (in)  Parameter name 
                          int    type,           # (in)  Parameter type 
                          int    nrow,           # (in)  number of rows 
                          int    ncol,           # (in)  number of columns 
                          int    *ipmtr)         # (out) Parameter index (bias-1) 
    
    # get info about a Parameter 
    cdef int ocsmGetPmtr(void   *modl,           # (in)  pointer to MODL 
                         int    ipmtr,           # (in)  Parameter index (1-npmtr) 
                         int    *type,           # (out) Parameter type 
                         int    *nrow,           # (out) number of rows 
                         int    *ncol,           # (out) number of columns 
                         char   name[])          # (out) Parameter name (at least 33 long) 
    
    # print external Parameters to file 
    cdef int ocsmPrintPmtrs(void   *modl,        # (in)  pointer to MODL 
                            FILE   *fp)          # (in)  pointer to FILE 
    
    # get the Value of a Parameter 
    cdef int ocsmGetValu(void   *modl,           # (in)  pointer to MODL 
                         int    ipmtr,           # (in)  Parameter index (1-npmtr) 
                         int    irow,            # (in)  row    index (1-nrow) 
                         int    icol,            # (in)  column index (1-ncol) 
                         double *value,          # (out) Parameter value 
                         double *dot)            # (out) Parameter velocity 
    
    # get the Value of a string Parameter 
    cdef int ocsmGetValuS(void   *modl,          # (in)  pointer to MODL 
                          int    ipmtr,          # (in)  Parameter index (1-npmtr) 
                          char   str[])          # (out) Parameter (string) value 
    
    # set a Value for a Parameter 
    cdef int ocsmSetValu(void   *modl,           # (in)  pointer to MODL 
                         int    ipmtr,           # (in)  Parameter index (1-npmtr) 
                         int    irow,            # (in)  row    index (1-nrow) 
                         int    icol,            # (in)  column index (1-ncol) or 0 for index
                         char   defn[])          # (in)  definition of Value 
    
    # set a (double) Value for a Parameter 
    cdef int ocsmSetValuD(void   *modl,          # (in)  pointer to MODL 
                          int    ipmtr,          # (in)  Parameter index (1-npmtr) 
                          int    irow,           # (in)  row    index (1-nrow) 
                          int    icol,           # (in)  column index (1-ncol) or 0 for index
                          double value)          # (in)  value to set 
    
    # get the Bounds of a Parameter 
    cdef int ocsmGetBnds(void   *modl,           # (in)  pointer to MODL 
                         int    ipmtr,           # (in)  Parameter index (1-npmtr) 
                         int    irow,            # (in)  row    index (1-nrow) 
                         int    icol,            # (in)  column index (1-ncol) 
                         double *lbound,         # (out) lower Bound 
                         double *ubound)         # (out) upper Bound 
    
    # set sensitivity FD time step (or select analytic) 
    cdef int ocsmSetDtime(void   *modl,          # (in)  pointer to MODL 
                          double dtime)          # (in)  time step (or 0 to choose analytic) 
    
    # set the velocity for a Parameter 
    cdef int ocsmSetVel(void   *modl,            # (in)  pointer to MODL 
                        int    ipmtr,            # (in)  Parameter index (1-npmtr) or 0 for all 
                        int    irow,             # (in)  row    index (1-nrow)     or 0 for all 
                        int    icol,             # (in)  column index (1-ncol)     or 0 for index 
                        char   defn[])           # (in)  definition of Velocity 
    
    # set the (double) velocity for a Parameter 
    cdef int ocsmSetVelD(void   *modl,           # (in)  pointer to MODL 
                         int    ipmtr,           # (in)  Parameter index (1-npmtr) or 0 for all 
                         int    irow,            # (in)  row    index (1-nrow)     or 0 for all 
                         int    icol,            # (in)  column index (1-ncol)     or 0 for index 
                         double dot)             # (in)  velocity to set 
    
    # get the parametric coordinates on an Edge or Face 
    cdef int ocsmGetUV(void   *modl,             # (in)  pointer to MODL 
                       int    ibody,             # (in)  Body index (bias-1) 
                       int    seltype,           # (in)  OCSM_EDGE, or OCSM_FACE 
                       int    iselect,           # (in)  Edge or Face index (bias-1) 
                       int    npnt,              # (in)  number of points 
                       double xyz[],             # (in)  coordinates (NULL or 3*npnt in length) 
                       double uv[])              # (out) para coords (1*npnt or 2*npnt in length) 
    
    # get the coordinates on a Node, Edge, or Face 
    cdef int ocsmGetXYZ(void   *modl,            # (in)  pointer to MODL 
                        int    ibody,            # (in)  Body index (bias-1) 
                        int    seltype,          # (in)  OCSM_NODE, OCSM_EDGE, or OCSM_FACE 
                        int    iselect,          # (in)  Node, Edge, or Face index (bias-1) 
                        int    npnt,             # (in)  number of points 
                        double uv[],             # (in)  para coords (NULL, 1*npnt, or 2*npnt) 
                        double xyz[])           # (out) coordinates (3*npnt in length) 
    
    # get the unit normals for a Face 
    cdef int ocsmGetNorm(void   *modl,           # (in)  pointer to MODL 
                         int    ibody,           # (in)  Body index (bias-1) 
                         int    iface,           # (in)  Face index (bias-1) 
                         int    npnt,            # (in)  number of points 
                         double uv[],            # (in)  para coords (NULL or 2*npnt in length) 
                         double norm[])          # (out) normals (3*npnt in length) 
    
    # get the velocities of coordinates on a Node, Edge, or Face 
    cdef int ocsmGetVel(void   *modl,            # (in)  pointer to MODL 
                        int    ibody,            # (in)  Body index (bias-1) 
                        int    seltype,          # (in)  OCSM_NODE, OCSM_EDGE, or OCSM_FACE 
                        int    iselect,          # (in)  Node, Edge, or Face index (bias-1) 
                        int    npnt,             # (in)  number of points 
                        double uv[],             # (in)  para coords
                                                     #   NULL           for OCSM_NODE
                                                     #   NULL or 1*npnt for OCSM_EDGE
                                                     #   NULL or 2*npnt for OCSM_FACE 
                        double vel[])            # (out) velocities (in pre-allocated array)
                                                     #    3      for OCSM_NODE
                                                     #   3*npnt for OCEM_EDGE
                                                     #   3*npnt for OCSM_FACE 
    
    # set up alternative tessellation by an external grid generator 
    cdef int ocsmSetEgg(void   *modl,            # (in)  pointer to MODL 
                        char   *eggname)         # (in)  name of dynamically-loadable file 
    
    # get the tessellation velocities on a Node, Edge, or Face 
    cdef int ocsmGetTessVel(void   *modl,        # (in)  pointer to MODL 
                            int    ibody,        # (in)  Body index (bias-1) 
                            int    seltype,      # (in)  OCSM_NODE, OCSM_EDGE, or OCSM_FACE 
                            int    iselect,      # (in)  Node, Edge, or Face index (bias-1) 
                            const double *dxyz[])     # (out) pointer to storage containing velocities 
    
    # get info about a Body 
    cdef int ocsmGetBody(void   *modl,           # (in)  pointer to MODL 
                         int    ibody,           # (in)  Body index (1-nbody) 
                         int    *type,           # (out) Branch type (see below) 
                         int    *ichld,          # (out) ibody of child (or 0 if root) 
                         int    *ileft,          # (out) ibody of left parent (or 0) 
                         int    *irite,          # (out) ibody of rite parent (or 0) 
                         double vals[],          # (out) array  of Argument values (at least 10 long) 
                         int    *nnode,          # (out) number of Nodes 
                         int    *nedge,          # (out) number of Edges 
                         int    *nface)          # (out) number of Faces 
    
    # print all Bodys to file 
    cdef int ocsmPrintBodys(void   *modl,        # (in)  pointer to MODL 
                            FILE   *fp)          # (in)  pointer to FILE 
    
    # print the BRep associated with a specific Body 
    cdef int ocsmPrintBrep(void   *modl,         # (in)  pointer to MODL 
                           int    ibody,         # (in)  Body index (1-nbody) 
                           FILE   *fp)           # (in)  pointer to File 
    
    # evaluate an expression 
    cdef int ocsmEvalExpr(void   *modl,          # (in)  pointer to MODL 
                          char   expr[],         # (in)  expression 
                          double *value,         # (out) value 
                          double *dot,           # (out) velocity 
                          char   str[])          # (out) value if string-valued (w/o leading $) 
    
    # print the contents of an EGADS ego 
    cdef void ocsmPrintEgo( cEGADS.ego    obj)          # (in)  EGADS ego 
    
    # convert an OCSM code to text 
    #@observer@
    cdef char *ocsmGetText(int    icode)        # (in)  code to look up 
    
    # convert text to an OCSM code 
    cdef int ocsmGetCode(char   *text)          # (in)  text to look up 
    
    #Definitions
    cdef int OCSM_DESPMTR  "OCSM_DESPMTR"
    cdef int OCSM_CFGPMTR  "OCSM_CFGPMTR"
    cdef int OCSM_CONPMTR  "OCSM_CONPMTR"
    cdef int OCSM_LOCALVAR "OCSM_LOCALVAR"
    cdef int OCSM_OUTPMTR  "OCSM_OUTPMTR"
    cdef int OCSM_UNKNOWN  "OCSM_UNKNOWN"

    cdef int OCSM_ACTIVE     "OCSM_ACTIVE"   # Branch activities
    cdef int OCSM_SUPPRESSED "OCSM_SUPPRESSED"
    cdef int OCSM_INACTIVE   "OCSM_INACTIVE"
    cdef int OCSM_DEFERRED   "OCSM_DEFERRED"

    cdef int OCSM_NODE  "OCSM_NODE"   # Selector types
    cdef int OCSM_EDGE  "OCSM_EDGE"
    cdef int OCSM_FACE  "OCSM_FACE"
    cdef int OCSM_BODY  "OCSM_BODY"

    # Errors                             
    cdef int OCSM_FILE_NOT_FOUND 
    cdef int OCSM_ILLEGAL_STATEMENT 
    cdef int OCSM_NOT_ENOUGH_ARGS 
    cdef int OCSM_NAME_ALREADY_DEFINED 
    cdef int OCSM_NESTED_TOO_DEEPLY
    cdef int OCSM_IMPROPER_NESTING  
    cdef int OCSM_NESTING_NOT_CLOSED    
    cdef int OCSM_NOT_MODL_STRUCTURE
    cdef int OCSM_PROBLEM_CREATING_PERTURB      
    
    cdef int OCSM_MISSING_MARK                  
    cdef int OCSM_INSUFFICIENT_BODYS_ON_STACK   
    cdef int OCSM_WRONG_TYPES_ON_STACK           
    cdef int OCSM_DID_NOT_CREATE_BODY            
    cdef int OCSM_CREATED_TOO_MANY_BODYS          
    cdef int OCSM_TOO_MANY_BODYS_ON_STACK         
    cdef int OCSM_ERROR_IN_BODYS_ON_STACK         
    cdef int OCSM_MODL_NOT_CHECKED              
    cdef int OCSM_NEED_TESSELLATION              
    
    cdef int OCSM_BODY_NOT_FOUND                
    cdef int OCSM_FACE_NOT_FOUND               
    cdef int OCSM_EDGE_NOT_FOUND                  
    cdef int OCSM_NODE_NOT_FOUND
    cdef int OCSM_ILLEGAL_VALUE                   
    cdef int OCSM_ILLEGAL_ATTRIBUTE               
    cdef int OCSM_ILLEGAL_CSYSTEM                 
    
    cdef int OCSM_SKETCH_IS_OPEN                  
    cdef int OCSM_SKETCH_IS_NOT_OPEN              
    cdef int OCSM_COLINEAR_SKETCH_POINTS          
    cdef int OCSM_NON_COPLANAR_SKETCH_POINTS      
    cdef int OCSM_TOO_MANY_SKETCH_POINTS          
    cdef int OCSM_TOO_FEW_SPLINE_POINTS           
    cdef int OCSM_SKETCH_DOES_NOT_CLOSE           
    cdef int OCSM_SELF_INTERSECTING               
    cdef int OCSM_ASSERT_FAILED                   
    
    cdef int OCSM_ILLEGAL_CHAR_IN_EXPR            
    cdef int OCSM_CLOSE_BEFORE_OPEN               
    cdef int OCSM_MISSING_CLOSE                   
    cdef int OCSM_ILLEGAL_TOKEN_SEQUENCE          
    cdef int OCSM_ILLEGAL_NUMBER                  
    cdef int OCSM_ILLEGAL_PMTR_NAME               
    cdef int OCSM_ILLEGAL_FUNC_NAME               
    cdef int OCSM_ILLEGAL_TYPE                    
    cdef int OCSM_ILLEGAL_NARG                    
    
    cdef int OCSM_NAME_NOT_FOUND                  
    cdef int OCSM_NAME_NOT_UNIQUE                 
    cdef int OCSM_PMTR_IS_DESPMTR
    cdef int OCSM_PMTR_IS_CONPMTR
    cdef int OCSM_PMTR_IS_OUTPMTR
    cdef int OCSM_PMTR_IS_LOCALVAR
    cdef int OCSM_WRONG_PMTR_TYPE                 
    cdef int OCSM_FUNC_ARG_OUT_OF_BOUNDS          
    cdef int OCSM_VAL_STACK_UNDERFLOW              
    cdef int OCSM_VAL_STACK_OVERFLOW               
    
    cdef int OCSM_ILLEGAL_BRCH_INDEX              
    cdef int OCSM_ILLEGAL_PMTR_INDEX              
    cdef int OCSM_ILLEGAL_BODY_INDEX               
    cdef int OCSM_ILLEGAL_ARG_INDEX              
    cdef int OCSM_ILLEGAL_ACTIVITY                 
    cdef int OCSM_ILLEGAL_MACRO_INDEX              
    cdef int OCSM_ILLEGAL_ARGUMENT                
    cdef int OCSM_CANNOT_BE_SUPPRESSED            
    cdef int OCSM_STORAGE_ALREADY_USED            
    cdef int OCSM_NOTHING_PREVIOUSLY_STORED       
    
    cdef int OCSM_SOLVER_IS_OPEN                  
    cdef int OCSM_SOLVER_IS_NOT_OPEN              
    cdef int OCSM_TOO_MANY_SOLVER_VARS            
    cdef int OCSM_UNDERCONSTRAINED                
    cdef int OCSM_OVERCONSTRAINED                 
    cdef int OCSM_SINGULAR_MATRIX                 
    cdef int OCSM_NOT_CONVERGED                   
    
    cdef int OCSM_UDP_ERROR1                      
    cdef int OCSM_UDP_ERROR2                      
    cdef int OCSM_UDP_ERROR3                      
    cdef int OCSM_UDP_ERROR4                      
    cdef int OCSM_UDP_ERROR5                      
    cdef int OCSM_UDP_ERROR6                      
    cdef int OCSM_UDP_ERROR7                      
    cdef int OCSM_UDP_ERROR8                      
    cdef int OCSM_UDP_ERROR9                      
    
    cdef int OCSM_OP_STACK_UNDERFLOW              
    cdef int OCSM_OP_STACK_OVERFLOW               
    cdef int OCSM_RPN_STACK_UNDERFLOW             
    cdef int OCSM_RPN_STACK_OVERFLOW              
    cdef int OCSM_TOKEN_STACK_UNDERFLOW          
    cdef int OCSM_TOKEN_STACK_OVERFLOW           
    cdef int OCSM_UNSUPPORTED                    
    cdef int OCSM_INTERNAL_ERROR                 