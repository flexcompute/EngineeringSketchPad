
// These need to be consistent with "index.inc"

//---- pointers for primary variables  vars( . k )   (before projection)
enum hsm_vars {
  ivrx ,   //  r_x     position vector
  ivry ,   //  r_y
  ivrz ,   //_ r_z
  ivdx ,   //  d_x     material director vector
  ivdy ,   //  d_y
  ivdz ,   //_ d_z
  ivps ,   //  psi     drilling rotation DOF
  IVTOT    //  7
};

//---- pointers for equation residuals, projected variables
enum hsm_res {
  irf1 ,   //  f1   r1
  irf2 ,   //  f2   r2
  irfn ,   //_ fn   r3
  irm1 ,   //  m1   d1
  irm2 ,   //_ m2   d2
  irmn ,   //_ mn   psi
  IRTOT    //  6
};

//---- pointers for input parameters  pars( . k )
enum hsm_pars {
  lvr0x ,   //  r0    undeformed global position vector
  lvr0y ,   //
  lvr0z ,   //_
  lve01x,   //  e0_1  local basis unit tangential vector 1
  lve01y,   //
  lve01z,   //_
  lve02x,   //  e0_2  local basis unit tangential vector 2
  lve02y,   //
  lve02z,   //_
  lvn0x ,   //  n0    local basis unit normal vector
  lvn0y ,   //
  lvn0z ,   //_
  lva11 ,   //  extension/shear stiffness matrix A in e01,e02 axes
  lva22 ,   //
  lva12 ,   //
  lva16 ,   //
  lva26 ,   //
  lva66 ,   //_
  lvb11 ,   //  extension/bending stiffness matrix B in e01,e02 axes
  lvb22 ,   //
  lvb12 ,   //
  lvb16 ,   //
  lvb26 ,   //
  lvb66 ,   //_
  lvd11 ,   //  bending stiffness matrix D in e01,e02 axes
  lvd22 ,   //
  lvd12 ,   //
  lvd16 ,   //
  lvd26 ,   //
  lvd66 ,   //_
  lva55 ,   //  transverse shear stiffness, e01,n0 component
  lva44 ,   //_ transverse shear stiffness, e02,n0 component
  lvq1  ,   //  shell-following normal force/area  q_12n
  lvq2  ,   //
  lvqn  ,   //
  lvqx  ,   //  fixed-direction force/area  q_xyz
  lvqy  ,   //
  lvqz  ,   //_
  lvtx  ,   //  fixed-direction moment/area tau_xyz
  lvty  ,   //
  lvtz  ,   //_
  lvmu  ,   //_ mass/area  mu
  lvtsh ,   //_ shell thickness  h
  lvzrf ,   //_ reference-plane location  zeta_ref =  -1 .. +1
  LVTOT     //  44
};

//---- pointers for global parameters  parg( . )
enum hsm_parg {
  lggeex,  //  g_x
  lggeey,  //  g_y
  lggeez,  //_ g_z
  lgvelx,  //  U_x   velocity  U
  lgvely,  //  U_y
  lgvelz,  //_ U_z
  lgrotx,  //  W_x   rotation rate  Omega
  lgroty,  //  W_y
  lgrotz,  //_ W_z
  lgvacx,  //  Udot_x   linear acceleration  U_dot
  lgvacy,  //  Udot_y
  lgvacz,  //_ Udot_z
  lgracx,  //  Wdot_x   angular acceleration  Omega_dot
  lgracy,  //  Wdot_y
  lgracz,  //_ Wdot_z
  lgposx,  //  R_X   position of xyz origin in XYZ (earth) axes
  lgposy,  //  R_Y
  lgposz,  //_ R_Z
  lgephi,  //  Phi   Euler angles of xyz frame
  lgethe,  //  Theta
  lgepsi,  //_ Psi
  lggabx,  //  g_X   gravity in XYZ (earth) axes, typically (0,0,-g)
  lggaby,  //  g_Y
  lggabz,  //_ g_Z
  LGTOT    //  24
};

//---- pointers for edge-load BC parameters for edge nodes 1,2  pare( . l )
enum hsm_pare {
  lef1x,    //  f1  force/length  vector in xyz axes
  lef1y,    //
  lef1z,    //_
  lef2x,    //  f2  force/length  vector in xyz axes
  lef2y,    //
  lef2z,    //_
  lef1t,    //  f1  force/length  vector in tln axes
  lef1l,    //
  lef1n,    //_
  lef2t,    //  f2  force/length  vector in tln axes
  lef2l,    //
  lef2n,    //_
  lem1x,    //  m1  moment/length vector in xyz axes
  lem1y,    //
  lem1z,    //_
  lem2x,    //  m2  moment/length vector in xyz axes
  lem2y,    //
  lem2z,    //_
  lem1t,    //  m1  moment/length vector in tln axes
  lem1l,    //
  lem1n,    //_
  lem2t,    //  m2  moment/length vector in tln axes
  lem2l,    //
  lem2n,    //_
  LETOT     //  24
};

//---- pointers for node BC parameters  parp( . l )
//-    (only some of this data is used, depending on the BC flags for the node)
enum hsm_parp {
  lprx ,    //  rBC   specified boundary position
  lpry ,    //
  lprz ,    //_
  lpn1x,    //  nBC1  specified boundary normal vector
  lpn1y,    //
  lpn1z,    //_
  lpn2x,    //  nBC2  specified boundary normal vector
  lpn2y,    //
  lpn2z,    //_
  lpt1x,    //  tBC1  specified surface-tangent edge-normal vector
  lpt1y,    //
  lpt1z,    //_
  lpt2x,    //  tBC2  specified surface-tangent edge-normal vector
  lpt2y,    //
  lpt2z,    //_
  lpfx ,    //  fixed-direction point force,  F_xyz
  lpfy ,    //
  lpfz ,    //_
  lpmx ,    //  fixed-direction point moment, M_xyz
  lpmy ,    //
  lpmz ,    //_
  LPTOT     //  21
};

//---- pointers for secondary variables from post-processing  deps( . k )
enum hsm_deps {
  jve1x,   //  e_1   local basis unit tangential vector 1
  jve1y,   //
  jve1z,   //_
  jve2x,   //  e_2   local basis unit tangential vector 2
  jve2y,   //
  jve2z,   //_
  jvnx ,   //  n     local basis unit normal vector
  jvny ,   //
  jvnz ,   //_
  jve11,   //  e_11    eps11 strain
  jve22,   //  e_22    eps22 strain
  jve12,   //_ e_12    eps12 strain
  jvk11,   //  k_11    kap11 curvature change
  jvk22,   //  k_22    kap22 curvature change
  jvk12,   //_ k_12    kap12 curvature change
  jvf11,   //  f_11    f11 stress resultant
  jvf22,   //  f_22    f22 stress resultant
  jvf12,   //_ f_12    f12 stress resultant
  jvm11,   //  m_11    m11 stress-moment resultant
  jvm22,   //  m_22    m22 stress-moment resultant
  jvm12,   //_ m_12    m12 stress-moment resultant
  jvf1n,   //  f_1n    f1n transverse shear stress resultant
  jvf2n,   //_ f_2n    f2n transverse shear stress resultant
  jvga1,   //  gam_1   n tilt angle in e1 direction
  jvga2,   //_ gam_2   n tilt angle in e2 direction
  JVTOT    //  25
};

enum hsm_lbc {
  lbcr1 = 1,   //  restrict displacement normal to a surface
  lbcr2 = 2,   //  restrict displacement normal to two surfaces
  lbcr3 = 3,   //  restrict all displacement
  lbcd1 = 10,  //  restrict drilling DOF
  lbcd2 = 20,  //  restrict director lean in 1 direction, and drilling DOF
  lbcd3 = 30,  //  restrict director lean in 2 directions, and drilling DOF
  lbcf  = 100, //  force applied to node
  lbcm  = 1000 //  moment applied to node
};

typedef struct {
    /*
    int numNode;
    int numElement;
    int numBCEdge;
    int numBCNode;
    int numJoint;
    */

    int numBCNode;

    double *parg; // Global parameters [LGTOT];

    // Note although indicated as two-d arrays, arrays are packed one-dimensionally

    double *vars; // Node primary variables [IVTOT, numNode]
    double *deps; // Node dependent variables [JVTOT, numNode]

    int *kelem;   // Element connectivity [1:4, numElement]

    double *pars; // (par-surface) Surface parameters [LVTOT, numNode]
    double *pare; // (par-edge)    Edge-BC parameters [LETOT, numBCEdge]
    double *parp; // (par-point)   Node-BC parameters [LPTOT, numBCNode]

    int *kbcedge; // Node indices of BC edge [1:2, numBCEdge]

    int *kbcnode; // Node of BC node [numBCNode]
    int *lbcnode; // Node BC type indicator [numBCNode]
    int *kjoint;  // Node Indices of the Joint [1:2, numJoint]

} hsmMemoryStruct;

typedef struct {

    // Don't know what these are variables
    double *bf;
    double *bf_dj;

    double *bm;
    double *bm_dj;

    double *resc;
    double *resc_vars;

    double *resp;
    double *resp_vars;
    double *resp_dvp;

    double *ares;

    double *dvars;

    double *res;

    double *rest;
    double *rest_t;

    double *resv;
    double *resv_v;

    int *ibx;
    int *kdvp;
    int *ndvp;
    int *frst;
    int *idt;
    int *frstt;
    int *kdv;
    int *ndv;
    int *frstv;

    double *amat;
    double *amatt;
    double *amatv;
    int *ipp;

} hsmTempMemoryStruct;

