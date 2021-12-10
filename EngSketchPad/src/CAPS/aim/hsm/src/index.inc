
c---- pointers for primary variables  vars( . k )   (before projection)
      parameter (
     &   ivrx  =  1,   !  r_x     position vector
     &   ivry  =  2,   !  r_y
     &   ivrz  =  3,   !_ r_z
     &   ivdx  =  4,   !  d_x     material director vector
     &   ivdy  =  5,   !  d_y
     &   ivdz  =  6,   !_ d_z
     &   ivps  =  7,   !  psi     drilling rotation DOF
     &   ivtot =  7 )  ! 

c---- pointers for equation residuals, projected variables
      parameter (
     &   irf1  =  1,   !   f1   r1 
     &   irf2  =  2,   !   f2   r2
     &   irfn  =  3,   !_  fn   r3
     &   irm1  =  4,   !   m1   d1 
     &   irm2  =  5,   !_  m2   d2
     &   irmn  =  6,   !_  mn   psi
     &   irtot =  6 )  ! 

c---- pointers for input parameters  pars( . k )
      parameter (
     &   lvr0x =  1,   !  r0    undeformed global position vector
     &   lvr0y =  2,   ! 
     &   lvr0z =  3,   !_
     &   lve01x=  4,   !  e0_1  local basis unit tangential vector 1
     &   lve01y=  5,   ! 
     &   lve01z=  6,   !_
     &   lve02x=  7,   !  e0_2  local basis unit tangential vector 2
     &   lve02y=  8,   ! 
     &   lve02z=  9,   !_
     &   lvn0x = 10,   !  n0    local basis unit normal vector
     &   lvn0y = 11,   ! 
     &   lvn0z = 12,   !_
     &   lva11 = 13,   !  extension/shear stiffness matrix A in e01,e02 axes
     &   lva22 = 14,   !
     &   lva12 = 15,   !
     &   lva16 = 16,   !
     &   lva26 = 17,   !
     &   lva66 = 18,   !_
     &   lvb11 = 19,   !  extension/bending stiffness matrix B in e01,e02 axes
     &   lvb22 = 20,   !
     &   lvb12 = 21,   !
     &   lvb16 = 22,   !
     &   lvb26 = 23,   !
     &   lvb66 = 24,   !_
     &   lvd11 = 25,   !  bending stiffness matrix D in e01,e02 axes
     &   lvd22 = 26,   !
     &   lvd12 = 27,   !
     &   lvd16 = 28,   !
     &   lvd26 = 29,   !
     &   lvd66 = 30,   !_
     &   lva55 = 31,   !  transverse shear stiffness, e01,n0 component
     &   lva44 = 32,   !_ transverse shear stiffness, e02,n0 component
     &   lvq1  = 33,   !  shell-following normal force/area  q_12n
     &   lvq2  = 34,   !
     &   lvqn  = 35,   !
     &   lvqx  = 36,   !  fixed-direction force/area  q_xyz
     &   lvqy  = 37,   ! 
     &   lvqz  = 38,   !_
     &   lvtx  = 39,   !  fixed-direction moment/area tau_xyz
     &   lvty  = 40,   !
     &   lvtz  = 41,   !_
     &   lvmu  = 42,   !_ mass/area  mu
     &   lvtsh = 43,   !_ shell thickness  h
     &   lvzrf = 44,   !_ reference-plane location  zeta_ref =  -1 .. +1
     &   lvtot = 44 )  !

cc---- pointers for global variables  glv( . )
c      parameter (
c     &   igtot = 1  )

c---- pointers for global parameters  parg( . )
      parameter (
     &   lggeex =  1,  !  g_x
     &   lggeey =  2,  !  g_y
     &   lggeez =  3,  !_ g_z
     &   lgvelx =  4,  !  U_x   velocity  U
     &   lgvely =  5,  !  U_y
     &   lgvelz =  6,  !_ U_z
     &   lgrotx =  7,  !  W_x   rotation rate  Omega
     &   lgroty =  8,  !  W_y
     &   lgrotz =  9,  !_ W_z
     &   lgvacx = 10,  !  Udot_x   linear acceleration  U_dot
     &   lgvacy = 11,  !  Udot_y
     &   lgvacz = 12,  !_ Udot_z
     &   lgracx = 13,  !  Wdot_x   angular acceleration  Omega_dot
     &   lgracy = 14,  !  Wdot_y
     &   lgracz = 15,  !_ Wdot_z
     &   lgposx = 16,  !  R_X   position of xyz origin in XYZ (earth) axes
     &   lgposy = 17,  !  R_Y
     &   lgposz = 18,  !_ R_Z
     &   lgephi = 19,  !  Phi   Euler angles of xyz frame
     &   lgethe = 20,  !  Theta
     &   lgepsi = 21,  !_ Psi
     &   lggabx = 22,  !  g_X   gravity in XYZ (earth) axes, typically (0,0,-g)
     &   lggaby = 23,  !  g_Y
     &   lggabz = 24,  !_ g_Z
     &   lgtot  = 24 )

c---- pointers for edge-load BC parameters for edge nodes 1,2  pare( . l )
      parameter (
     &   lef1x =  1,    !  f1  force/length  vector in xyz axes
     &   lef1y =  2,    ! 
     &   lef1z =  3,    !_
     &   lef2x =  4,    !  f2  force/length  vector in xyz axes
     &   lef2y =  5,    ! 
     &   lef2z =  6,    !_
     &   lef1t =  7,    !  f1  force/length  vector in tln axes
     &   lef1l =  8,    !     
     &   lef1n =  9,    !_    
     &   lef2t = 10,    !  f2  force/length  vector in tln axes
     &   lef2l = 11,    !     
     &   lef2n = 12,    !_    
     &   lem1x = 13,    !  m1  moment/length vector in xyz axes
     &   lem1y = 14,    !
     &   lem1z = 15,    !_
     &   lem2x = 16,    !  m2  moment/length vector in xyz axes
     &   lem2y = 17,    !
     &   lem2z = 18,    !_
     &   lem1t = 19,    !  m1  moment/length vector in tln axes
     &   lem1l = 20,    !
     &   lem1n = 21,    !_    
     &   lem2t = 22,    !  m2  moment/length vector in tln axes
     &   lem2l = 23,    !
     &   lem2n = 24,    !_    
     &   letot = 24 )

c---- pointers for node BC parameters  parp( . l )
c-    (only some of this data is used, depending on the BC flags for the node)
      parameter (
     &   lprx  =  1,    !  rBC   specified boundary position
     &   lpry  =  2,    ! 
     &   lprz  =  3,    !_
     &   lpn1x =  4,    !  nBC1  specified boundary normal vector
     &   lpn1y =  5,    ! 
     &   lpn1z =  6,    !_
     &   lpn2x =  7,    !  nBC2  specified boundary normal vector
     &   lpn2y =  8,    ! 
     &   lpn2z =  9,    !_
     &   lpt1x = 10,    !  tBC1  specified surface-tangent edge-normal vector
     &   lpt1y = 11,    !
     &   lpt1z = 12,    !_
     &   lpt2x = 13,    !  tBC2  specified surface-tangent edge-normal vector
     &   lpt2y = 14,    !
     &   lpt2z = 15,    !_
     &   lpfx  = 16,    !  fixed-direction point force,  F_xyz
     &   lpfy  = 17,    ! 
     &   lpfz  = 18,    !_
     &   lpmx  = 19,    !  fixed-direction point moment, M_xyz
     &   lpmy  = 20,    ! 
     &   lpmz  = 21,    !_
     &   lptot = 21 )

c---- pointers for secondary variables from post-processing  deps( . k )
      parameter (
     &   jve1x =  1,   !  e_1   local basis unit tangential vector 1
     &   jve1y =  2,   !
     &   jve1z =  3,   !_
     &   jve2x =  4,   !  e_2   local basis unit tangential vector 2
     &   jve2y =  5,   !
     &   jve2z =  6,   !_
     &   jvnx  =  7,   !  n     local basis unit normal vector
     &   jvny  =  8,   !
     &   jvnz  =  9,   !_
     &   jve11 = 10,   !  e_11    eps11 strain
     &   jve22 = 11,   !  e_22    eps22 strain
     &   jve12 = 12,   !_ e_12    eps12 strain
     &   jvk11 = 13,   !  k_11    kap11 curvature change
     &   jvk22 = 14,   !  k_22    kap22 curvature change
     &   jvk12 = 15,   !_ k_12    kap12 curvature change
     &   jvf11 = 16,   !  f_11    f11 stress resultant
     &   jvf22 = 17,   !  f_22    f22 stress resultant
     &   jvf12 = 18,   !_ f_12    f12 stress resultant
     &   jvm11 = 19,   !  m_11    m11 stress-moment resultant
     &   jvm22 = 20,   !  m_22    m22 stress-moment resultant
     &   jvm12 = 21,   !_ m_12    m12 stress-moment resultant
     &   jvf1n = 22,   !  f_1n    f1n transverse shear stress resultant
     &   jvf2n = 23,   !_ f_2n    f2n transverse shear stress resultant
     &   jvga1 = 24,   !  gam_1   n tilt angle in e1 direction
     &   jvga2 = 25,   !_ gam_2   n tilt angle in e2 direction
     &   jvtot = 25 )  ! 

