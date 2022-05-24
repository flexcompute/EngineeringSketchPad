/*
 ************************************************************************
 *                                                                      *
 * OpenCSM.h -- header for OpenCSM.c                                    *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

#ifndef _OPENCSM_H_
#define _OPENCSM_H_

#define OCSM_MAJOR_VERSION  1
#define OCSM_MINOR_VERSION 21

#define MAX_NAME_LEN       64           /* maximum chars in name */
#define MAX_EXPR_LEN      512           /* maximum chars in expression */
#define MAX_FILENAME_LEN  256           /* maximum chars in filename */
#define MAX_LINE_LEN     2048           /* maximum chars in line in input file */
#define MAX_STR_LEN      4096           /* maximum chars in any string */
#define MAX_STRVAL_LEN    256           /* maximum chars in any string value */
#define MAX_STACK_SIZE   4096           /* maximum size of stack */
#define MAX_NESTING        20           /* maximum number of nested patterns, ifthens, macros, or UDCs */
#define MAX_SKETCH_SIZE  1024           /* maximum points in Sketch */
#define MAX_SOLVER_SIZE   256           /* maximum variable in solver */
#define MAX_NUM_SKETCHES  100           /* maximum number of Sketches in rule/loft/blend */
#define MAX_NUM_MACROS    100           /* maximum number of storage locations */

#include <time.h>
#include "egads.h"

/*
 ************************************************************************
 *                                                                      *
 * Definitions                                                          *
 *                                                                      *
 ************************************************************************
 */

/*
Node
    - a location in 3D space that serves as the terminus for one or
      more Edges

Edge
    - is associated with a 3D curve (if not degenerate)
    - has a range of parametric coordinates, where tmin <= t <= tmax
    - the positive orientation is going from tmin to tmax
    - has a Node for tmin and for tmax
    - if the Nodes at tmin and tmax are the same, the Edge forms a
      closed Loop (that is, is periodic) or is degenerate (if tmin
      equals tmax); otherwise it is open

Loop
    - free standing collection of one or more connected Edges with
      associated senses
    - if the Loop is closed, each of the corresponding Nodes is
      associated with exactly two Edges
    - if the Loop is open, the intermediate Nodes are each associated
      with two Edges and the Nodes at the ends each correspond to one
      Edge
    - the sense of the Loop is associated with the order of the Edges
      in the Loop and their associated senses

Face
    - a surface bounded by one or more Loops with associated senses
    - there may be only one outer Loop (sense = 1) and any number
      of inner Loops (sense = -1)
    - all associated Loops must be closed

Shell
    - a collection of one of more connected Faces
    - if all the Edges associated with a Shell are used by exactly two
      Faces in the Shell, the Shell is closed (manifold) and it
      segregates regionds of 3D space; otherwise the Shell is open

NODE  (NodeBody)
    - a single Node (which can be used as the terminus in operations
      such as RULE, BLEND, and LOFT)
    - formed by creating an empty Sketch

WIRE  (WireBody)
    - a single Loop

SHEET  (SheetBody)
    - a single Shell that can be either non-manifold (open) or manifold
      (closed)

SOLID  (SolidBody)
    - a manifold collection of one or more closed Shells with associated
      senses
    - there may be only one outer Shell (sense = 1) and any number of
      inner Shells (sense = -1)
 */

/*
 ************************************************************************
 *                                                                      *
 * CSM file format                                                      *
 *                                                                      *
 ************************************************************************
 */

/*
The .csm file contains a series of statements.

If a line contains a hash (#), all characters starting at the hash
   are ignored.

If a line contains a backslash, all characters starting at the
   backslash are ignored and the next line is appended; spaces at
   the beginning of the next line are treated normally.

All statements begin with a keyword (described below) and must
   contain at least the indicated number of arguments.

The keywords may either be all lowercase or all UPPERCASE.

Any CSM statement can be used except the INTERFACE statement.

Blocks of statements must be properly nested.  The Blocks are bounded
   by PATBEG/PATEND, IFTHEN/ELSEIF/ELSE/ENDIF, SOLBEG/SOLEND,
   and CATBEG/CATEND.

Extra arguments in a statement are discarded.  If one wants to add
   a comment, it is recommended to begin it with a hash (#) in case
   optional arguments are added in future releases.

Any statements after an END statement are ignored.

All arguments must not contain any spaces or must be enclosed
   in a pair of double quotes (for example, "a + b").

Parameters are evaluated in the order that they appear in the
   file, using MATLAB-like syntax (see 'Expression rules' below).

During the build process, OpenCSM maintains a LIFO 'Stack' that
   can contain Bodys and Sketches.

The csm statements are executed in a stack-like way, taking their
   inputs from the Stack and depositing their results onto the Stack.

The default name for each Branch is 'Brch_xxxxxx', where xxxxxx
   is a unique sequence number.

Special characters:
   #          introduces comment
   "          ignore spaces until following "
   \          ignore this and following characters and concatenate next line
   <space>    separates arguments in .csm file (except between " and ")

   0-9        digits used in numbers and in names
   A-Z a-z    letters used in names
   _ : @      characters used in names (see rule for names)
   ? % =      characters used in strings
   .          decimal separator (used in numbers), introduces dot-suffixes
                 (in names)
   ,          separates function arguments and row/column in subscripts
   ;          multi-value item separator
   ( )        groups expressions and function arguments
   [ ]        specifies subscripts in form [row,column] or [index]
   { } < >    characters used in strings
   + - * / ^  arithmetic operators
   $          as first character, introduces a string that is terminated
                 by end-of-line or un-escaped plus, comma, or
                 close-parenthesis
   @          as first character, introduces @-parameters (see below)
   '          used to escape comma, plus, or close-parenthesis within
                 strings
   !          if first character of implicit string, ignore $! and treat
                 as an expression

   |          cannot be used (reserved for OpenCSM internals)
   ~          cannot be used (reserved for OpenCSM internals)
   &          cannot be used (reserved for OpenCSM internals)
*/

/*
 ************************************************************************
 *                                                                      *
 * CPC file format                                                      *
 *                                                                      *
 ************************************************************************
 */

/*
A .cpc file follows the rules of a .csm file, EXCEPT:

UDPRIM statments that refer to a UDC do not revert to the .udc, but
    instead read the .udc contents from the .cpc file

END statements that are part of in included UDC create END branches
    and do not stop the reading process
*/

/*
 ************************************************************************
 *                                                                      *
 * UDC file format                                                      *
 *                                                                      *
 ************************************************************************
 */

/*
A .udc file follows the rules of a .csm file, EXCEPT:

Zero or more INTERFACE statements must preceed any other non-comment
    CSM statement.

Any CSM statement can be used except the CFGPMTR, CONPMTR, DESPMTR,
    OUTPMTR, LBOUND, and UBOUND statements.  Note that these statements
    (except CONPMTR) may be used in include-type UDCs at global scope.

SET statements define parameters that are visible only within the .udc
   file (that is, parameters have local scope), except that they have
   their parent's scope if in an include-type UDC.

Parameters defined outside the .udc file are not available, except those
   passed in via INTERFACE statements.

.udc files can be nested to a depth of 10 levels.

.udc files are executed via a UDPRIM statement.
*/

/*
 ************************************************************************
 *                                                                      *
 * Valid CSM statements                                                 *
 *                                                                      *
 ************************************************************************
 */

/*
APPLYCSYS $csysName ibody=0
          use:    transforms Group on top of stack so that their
                      origins/orientations coincide with given csys
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  if ibody>0, use csys associated with that Body
                  if ibody==0, then search for csys backward from
                     next-to-last Body on stack
                  if ibody==-1, transform Body on top of stack so
                     that its csys is moved to the origin
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $body_not_found
                     $insufficient_bodys_on_stack
                     $name_not_found

ARC       xend yend zend dist $plane=xy
          use:    create a new circular arc to the new point, with a
                     specified distance between the mid-chord and mid-arc
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  $plane must be xy, yz, or zx
                  if dist>0, sweep is counterclockwise
                  sensitivity computed w.r.t. xend, yend, zend, dist
                  signals that may be thrown/caught:

ASSERT    arg1 arg2 toler=0 verify=0
          use:    return error if arg1 and arg2 differ
          pops:   -
          pushes: -
          notes:  if toler==0, set toler=1e-6
                  if toler<0, set toler=abs(arg1*toler)
                  if (abs(arg1-arg2) > toler) return an error
                  only executed if verify<=MODL->verify
                  cannot be followed by ATTRIBUTE or CSYSTEM

ATTRIBUTE $attrName attrValue
          use:    sets an Attribute for the Group on top of Stack
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  if first char of attrValue is '$', then string Attribute
                  elseif attrValue is a Parameter name, all its elements
                     are stored in Attribute
                  otherwise attrValue is a semicolon-separated list of
                     scalar numbers/expressions
                  does not create a Branch
                  if before first Branch that creates a Body,
                     the Attribute is a string-valued global Attribute
                  if after BLEND, BOX, CHAMFER, COMBINE, CONE, CONNECT,
                        CYLINDER, EXTRUDE, FILLET, HOLLOW, IMPORT, LOFT,
                        RESTORE, REVOLVE, RULE, SPHERE, SWEEP, TORUS,
                        or UDPRIM
                     the Attribute is added to the Body and its Faces
                  else
                     the Attribute is only added to the Body
                  is applied to selected Nodes, Edges, or Faces if after a
                     SELECT statement

BEZIER    x y z
          use:    add a Bezier control point
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  sensitivity computed w.r.t. x, y, z
                  signals that may be thrown/caught:

BLEND     begList=0 endList=0 reorder=0 oneFace=0 periodic=0
          use:    create a Body by blending through Xsects since Mark
          pops:   Xsect1 ... Mark
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  all Xsects must have the same number of Edges
                  if all Xsects are NodeBodys
                     a WireBody is created
                  elseif all Xsects are WireBodys (or a NodeBody at one end)
                     a SheetBody is created
                  else
                     a SolidBody is created
                  if the first Xsect is a point
                      if begList is 0
                          pointed end is created
                      elseif begList contains 8 values
                          begList contains rad1;dx1;dy1;dz1;rad2;dx2;dy2;dz2
                          rounded end is created
                  elseif first Xsect is a WireBody
                      created SheetBody is open at the beginning
                  elseif first Xsect is a SheetBody
                      if begList is 0
                          created Body included SheetBody at its beginning
                      elseif begList contains 2 values and first is -1
                          begList contains -1;aspect
                          rounded end with approximately given aspect ratio
                  if first Xsect is a WireBody or SheetBody
                      if begList contains 3 values
                          begList describes inward tangency at beginning
                  if the last Xsect is a point
                      if endList is 0
                          pointed end is created
                      elseif endList contains 8 values
                          endList contains rad1;dx1;dy1;dz1;rad2;dx2;dy2;dz2
                          rounded end is created
                  elseif last Xsect is a WireBody
                      created SheetBody is open at the end
                  elseif last Xsect is a SheetBody
                      if endList is 0
                          created Body included SheetBody at its end
                      elseif endList contains 2 values and first is -1
                          endList contains -1;aspect
                          rounded end with approximately given aspect ratio
                  if last Xsect is a WireBody or SheetBody
                      if endList contains 3 values
                          endList describes inward tangency at end
                  if begList!=0 and endList!=0, there must be at least
                     three interior Xsects
                  interior Xsects can be repeated once for C1 continuity
                  interior Xsects can be repeated twice for C0 continuity
                  if reorder!=0 then Xsects are reordered to minimize Edge
                     lengths in the direction between Xsects
                  first Xsect is unaltered if reorder>0
                  last  Xsect is unaltered if reorder<0
                  if oneFace==1 then do not split at C0 (multiplicity=3)
                  if periodic=1 then connect the first and last Xsects
                  sensitivity computed w.r.t. begList, endList
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  Attributes on Xsects are maintained
                  face-order is: (base), (end), feat1:part1,
                     feat1:part2, ... feat2:part1, ...
                  signals that may be thrown/caught:
                     $error_in_bodys_on_stack
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

BOX       xbase ybase zbase dx dy dz
          use:    create a box SolidBody or planar SheetBody
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if one of dx, dy, or dz is zero, a SheetBody is created
                  if two of dx, dy, or dz is zero, a WireBody  is created
                  if dx, dy, dz      are all zero, a NodeBody  is creates
                  sensitivity computed w.r.t. xbase, ybase, zbase, dx, dy, dz
                  computes Face, Edge, and Node sensitivities analytically
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is: xmin, xmax, ymin, ymax, zmin, zmax
                  signals that may be thrown/caught:
                     $illegal_value

CATBEG    sigCode
          use:    execute Block of Branches if current signal matches
                     sigCode
          pops:   -
          pushes: -
          notes:  sigCode can be an integer or one of:
                     $all
                     $body_not_found
                     $colinear_sketch_points
                     $created_too_many_bodys
                     $did_not_create_body
                     $edge_not_found
                     $error_in_bodys_on_stack
                     $face_not_found
                     $file_not_found
                     $func_arg_out_of_bounds
                     $illegal_argument
                     $ilegal_attribute
                     $illegal_csystem
                     $illegal_pmtr_index
                     $illegal_pmtr_name
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $name_not_found
                     $node_not_found
                     $non_coplanar_sketch_points
                     $no_selection
                     $underconstrained
                     $overconstrained
                     $not_converged
                     $self_intersecting
                     $wrong_types_on_stack
                     $assert_failed
                  if sigCode does not match current signal, skip to matching
                     CATEND
                  Block contains all Branches up to matching CATEND
                  cannot be followed by ATTRIBUTE or CSYSTEM

CATEND
          use:    designates the end of a CATBEG Block
          pops:   -
          pushes: -
          notes:  inner-most Block must be a CATBEG Block
                  closes CATBEG Block
                  cannot be followed by ATTRIBUTE or CSYSTEM

CFGPMTR   $pmtrName value
          use:    define a configuration Paramater
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  statement may not be used in a .udc file
                  pmtrName must be in form 'name'
                  pmtrName must not start with '@'
                  pmtrName must not refer to an LOCALVAR/OUTPMTR/CONPMTR
                  pmtrName will be marked as CFGPMTR
                  pmtrName is used directly (without evaluation)
                  if value already exists, it is not overwritten
                  does not create a Branch
                  cannot be followed by ATTRIBUTE or CSYSTEM

CHAMFER   radius edgeList=0
          use:    apply a chamfer to a Body
          pops:   Body
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if previous operation is boolean, apply to all new Edges
                  edgeList=0 is the same as edgeList=[0;0]
                  edgeList is a multi-value Parameter or a semicolon-separated
                     list
                  pairs of edgeList entries are processed in order
                  pairs of edgeList entries are interpreted as follows:
                     col1  col2   meaning
                      =0    =0    add all Edges
                      >0    >0    add    Edges between iford=+icol1
                                                   and iford=+icol2
                      <0    <0    remove Edges between iford=-icol1
                                                   and iford=-icol2
                      >0    =0    add    Edges adjacent to iford=+icol1
                      <0    =0    remove Edges adjacent to iford=-icol1
                  sensitivity computed w.r.t. radius
                  sets up @-parameters
                  new Faces all receive the Branch's Attributes
                  face-order is based upon order that is returned from EGADS
                  signals that may be thrown/caught:
                     $illegal_argument
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

CIRARC    xon yon zon xend yend zend
          use:    create a new circular arc, using the previous point
                     as well as the two points specified
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  sensitivity computed w.r.t. xon, yon, zon, xend, yend, zend
                  signals that may be thrown/caught:

COMBINE   toler=0
          use:    combine Bodys since Mark into next higher type
          pops:   Body1 ... Mark
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  Mark must be set
                  if all Bodys since Mark are SheetBodys
                     create either a SolidBody from closed Shell or an
                     (open) SheetBody
                  elseif there is 1 planar WireBody that is closed
                     create SheetBody from Loop
                  elseif there is 1 planar WireBody that is open
                     create SheetBody from Loop after closing Loop first
                  elseif there are multiple planar WireBodys
                     create SheetBody from closed Loop
                  endif
                  if maxtol>0, then tolerance can be relaxed until successful
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $did_not_create_body
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

CONE      xvrtx yvrtx zvrtx xbase ybase zbase radius
          use:    create a cone Body
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. xvrtx, yvrtx, zvrtz, xbase, ybase,
                     zbase, radius
                  computes Face, Edge, and Node sensitivities analytically
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is: (empty), base, umin, umax
                     if x-aligned: umin=ymin, umax=ymax
                     if y-aligned: umin=xmax, umax=xmin
                     if z-aligned: umin=ymax, umax=ymin
                  signals that may be thrown/caught:
                     $illegal_value

CONNECT   faceList1 faceList2 edgeList1=0 edgeList2=0 toler=0
          use:    connects two Bodys with bridging Faces
          pops:   Body1 Body2
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  faceList1 and faceList2 must have the same length
                  edgeList1 and edgeList2 must have the same length
                  edgeList1[i] corresponds to edgeList2[i]
                  faceList1[i] corresponds to faceList2[i]
                  if Body1 is a Mark, left and rite Bodys are the same
                  if edgeLists are given
                      Body1 is either WireBody, SheetBody, or SolidBody
                      Body2 is same type as Body1
                      Body  is same type as Body1
                      Face in faceLists are removed
                      bridging Faces are made between edgeList pairs
                      a zero in an edgelist creates a degenerate Face
                  else
                      Body1 and Body2 must both be SolidBodys
                      Faces within each faceList must be contiguous
                      bridging Faces between exposed Edges are created
                  new Faces all receive the Branch's Attributes
                  sets up @-parameters
                  if edgeLists are given
                      face-order is same as edgeList
                  else
                      face-order is arbitrary
                  signals that may be thrown/caught:
                     $illegal_argument
                     $illegal_value
                     $insufficient_bodys_on_stack

CONPMTR   $pmtrName values
          use:    define a constant Parameter
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  statement may not be used in a .udc file
                  pmtrName must be in form 'name'
                  pmtrName must not start with '@'
                  pmtrName must not refer to an LOCALVAR/OUTPMTR/DESPMTR
                  pmtrName will be marked as CONPMTR
                  pmtrName is used directly (without evaluation)
                  pmtrName is available within .csm and .udc files
                  value(s) must be numbers
                  does not create a Branch
                  cannot be followed by ATTRIBUTE or CSYSTEM

CSYSTEM   $csysName csysList
          use:    attach a Csystem to Body on top of stack
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  if     csysList contains 9 entries:
                     {x0, y0, z0, dx1, dy1, dz1, dx2, dy2, dz2}
                     origin is at (x0,y0,z0)
                     dirn1  is in (dx1,dy1,dz1) direction
                     dirn2  is part of (dx2,dy2,dz2) that is orthog. to dirn1
                  elseif csysList contains 5 entries and first is positive
                     {+iface, ubar0, vbar0, du2, dv2}
                     origin is at normalized (ubar0,vbar0) in iface
                     dirn1  is normal to Face
                     dirn2  is in (du2,dv2) direction
                  elseif csyList contains 5 entries and first is negative
                     {-iedge, tbar, dx2, dy2, dz2}
                     origin is at normalized (tbar) in iedge
                     dirn1  is tangent to Edge
                     dirn2  is part of (dx2,dy2,dz2) that is orthog. to dirn1
                  elseif csysList contains 7 entries
                     {inode, dx1, dy1, dz1, dx2, dy2, dz2}
                     origin is at Node inode
                     dirn1  is in (dx1,dy1,dz1) direction
                     dirn2  is part of (dx2,dy2,dz2) that is orthog. to dirn1
                  else
                     error
                  semicolon-sep lists can instead refer to
                     multi-valued Parameter
                  dirn3 is formed by (dirn1)-cross-(dirn2)
                  does not create a Branch

CYLINDER  xbeg ybeg zbeg xend yend zend radius
          use:    create a cylinder Body
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. xbeg, ybeg, zbeg, xend, yend,
                     zend, radius
                  computes Face, Edge, and Node sensitivities analytically
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is: beg, end, umin, umax
                     if x-aligned: umin=ymin, umax=ymax
                     if y-aligned: umin=xmax, umax=xmin
                     if z-aligned: umin=ymax, umax=ymin
                  signals that may be thrown/caught:
                     $illegal_value

DESPMTR   $pmtrName values
          use:    define a design Parameter
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  statement may not be used in a function-type .udc file
                  pmtrName can be in form 'name' or 'name[irow,icol]'
                  pmtrName must not start with '@'
                  pmtrName must not refer to an LOCALVAR/OUTPMTR/CONPMTR
                  pmtrName will be marked as DESPMTR
                  pmtrName is used directly (without evaluation)
                  irow and icol cannot contain a comma or open bracket
                  if irow is a colon (:), then all rows    are input
                  if icol is a colon (:), then all columns are input
                  pmtrName[:,:] is equivalent to pmtrName
                  values cannot refer to any other Parameter
                  if value already exists, it is not overwritten
                  values are defined across rows, then across columns
                  if values has more entries than needed, extra values
                     are lost
                  if values has fewer entries than needed, last value
                     is repeated
                  does not create a Branch
                  cannot be followed by ATTRIBUTE or CSYSTEM

DIMENSION $pmtrName nrow ncol
          use:    set up or redimensions an array Parameter
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  nrow >= 1
                  ncol >= 1
                  pmtrName must not start with '@'
                  if applied to a DESPMTR or CFGPMTR, must be in either
                      .csm file or top-level include-style .udc file
                  a legacy fourth argument (despmtr) is no longer used
                  old values are not overwritten
                  cannot be followed by ATTRIBUTE or CSYSTEM

DUMP      $filename remove=0 toMark=0 withTess=0
          use:    write a file that contains the Body
          pops:   Body1 (if remove=1)
          pushes: -
          notes:  Solver may not be open
                  if file exists, it is overwritten
                  filename is used directly (without evaluation)
                  if filename starts with '$/', it is prepended with path of
                     the .csm file
                  if remove=1, then Body1 is removed after dumping
                  if toMark=1, all Bodys back to the Mark (or all if no Mark)
                     are combined into a single model
                  if toMark=1, the remove flag is ignored
                  if withTess!=0, add tessellations to .egads file
                  for .ugrid files, toMark must be 0
                  valid filetypes are:
                     .brep   .BREP   --> OpenCASCADE output
                     .bstl   .BSTL   --> binary stl  output
                     .egads  .EGADS  --> EGADS       output
                     .egg    .EGG    --> EGG restart output
                     .iges   .IGES   --> IGES        output
                     .igs    .IGS    --> IGES        output
                     .plot   .PLOT   --> ASCII plot  output
                     .sens   .SENS   --> ASCII sens  output
                     .step   .STEP   --> STEP        output
                     .stl    .STL    --> ASCII stl   output
                     .stp    .STP    --> STEP        output
                     .tess   .TESS   --> ASCII tess  output
                     .ugrid  .UGRID  --> ASCII AFRL3 output
                  if .bstl, use _stlColor from Face, Body, or 0 for color
                  if .egads, set _despmtr_* and _outpmtr_ Attributes on Model
                  signals that may be thrown/caught:
                     $file_not_found
                     $insufficient_bodys_on_stack

ELSE
          use:    execute or skip a Block of Branches
          pops:   -
          pushes: -
          notes:  inner-most Block must be an Ifthen Block
                  must follow an IFTHEN or ELSEIF statment
                  if preceeding (matching) IFTHEN or ELSEIF evaluated true,
                     then skip Branches up to the matching ENDIF
                  cannot be followed by ATTRIBUTE or CSYSTEM

ELSEIF    val1 $op1 val2 $op2=and val3=0 $op3=eq val4=0
          use:    execute or skip a sequence of Branches
          pops:   -
          pushes: -
          notes:  inner-most Block must be an Ifthen Block
                  must follow an IFTHEN or ELSEIF statement
                  if preceeding (matching) IFTHEN or ELSEIF evaluated true,
                     then skip Branches up to matching ENDIF
                  op1 must be one of: lt LT le LE eq EQ ge GE gt GT ne NE
                  op2 must be one of: or OR and AND xor XOR
                  op3 must be one of: lt LT le LE eq EQ ge GE gt GT ne NE
                  if expression evaluates false, skip Branches up to next
                     ELSEIF, ELSE, or ENDIF
                  cannot be followed by ATTRIBUTE or CSYSTEM

END
          use:    signifies end of .csm or .udc file
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  Bodys on Stack are returned last-in-first-out
                  cannot be followed by ATTRIBUTE or CSYSTEM

ENDIF
          use:    terminates an Ifthen Block of Branches
          pops:   -
          pushes: -
          notes:  inner-most Block must be an Ifthen Block
                  must follow an IFTHEN, ELSEIF, or ELSE
                  closes Ifthen Block
                  cannot be followed by ATTRIBUTE or CSYSTEM

EVALUATE  $type arg1 ...
          use:    evaluate coordinates of NODE, EDGE, or FACE
          pops:   -
          pushes: -
          notes:  if     arguments are: "node ibody inode"
                     ibody is Body number (1:nbody)
                     inode is Node number (1:nnode)
                     return in @edata:
                        x, y, z
                  elseif arguments are: "edge ibody iedge t"
                     ibody is Body number (1:nbody)
                     iedge is Edge number (1:nedge)
                     evaluate Edge at given t
                     return in @edata:
                        t (clipped),
                        x,      y,      z,
                        dxdt,   dydt,   dzdt,
                        d2xdt2, d2ydt2, d2zdt2
                  elseif arguments are: "edge ibody iedge $beg"
                     ibody is Body number (1:nbody)
                     iedge is Edge number (1:nedge)
                     evaluate Edge at given t
                     return in @edata:
                        t (clipped),
                        x,      y,      z,
                        dxdt,   dydt,   dzdt,
                        d2xdt2, d2ydt2, d2zdt2
                  elseif arguments are: "edge ibody iedge $end"
                     ibody is Body number (1:nbody)
                     iedge is Edge number (1:nedge)
                     evaluate Edge at given t
                     return in @edata:
                        t (clipped),
                        x,      y,      z,
                        dxdt,   dydt,   dzdt,
                        d2xdt2, d2ydt2, d2zdt2
                  elseif arguments are: "edgerng ibody iedge"
                     ibody is Body number (1:nbody)
                     iedge is Edge number (1:nedge)
                     return in @edata:
                        tmin, tmax
                  elseif arguments are: "edgeinv ibody iedge x y z"
                     ibody is Body number (1:nbody)
                     iedge is Edge number (1:nedge)
                     inverse evaluate Edge at given (x,y,z)
                     return in @edata:
                        t,
                        xclose,  yclose,  zclose
                  elseif arguments are: "face ibody iface u v"
                     ibody is Body number (1:nbody)
                     iface is Face number (1:nface)
                     evaluate Face at given (u,v)
                     return in @edata:
                        u (clipped), v (clipped),
                        x,       y,       z,
                        dxdu,    dydu,    dzdu,
                        dxdv,    dydv,    dzdv,
                        d2xdu2,  d2ydu2,  d2zdu2,
                        d2xdudv, d2ydudv, d2zdudv,
                        d2xdv2,  d2ydv2,  d2zdv2,
                        normx,   normy,   normz
                  elseif arguments are: "facerng ibody iface"
                     ibody is Body number (1:nbody)
                     iface is Face number (1:nface)
                     return in @edata:
                        umin, umax, vmin, vmax
                  elseif arguments are: "faceinv ibody iface x y z"
                     ibody is Body number (1:nbody)
                     iface is Face number (1:nface)
                     inverse evaluate Face at given (x,y,z)
                     return in @edata:
                        u,       v,
                        xclose,  yclose,  zclose
                  cannot be followed by ATTRIBUTE or CSYSTEM
                  signals that may be thrown/caught:
                     $body_not_found
                     $edge_not_found
                     $face_not_found
                     $node_not_found

EXTRACT   entList
          use:    extract Face(s) or Edge(s) from a Body
          pops:   Body1
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  all members of entList must have the same sign
                  Body1 must be a SolidBody or a SheetBody
                  if     entList entries are all positive
                     create SheetBody from entList Face(s) of Body1
                  elseif entList entries are all negative
                     create WireBody from -entList Edge(0) of Body1
                  elseif Body1=SolidBody and entList=0
                     create SheetBody from outer Shell of Body1
                  elseif Body1=SheetBody and entList=0
                     create WireBody from outer Loop of Body1
                  endif
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack
                     $did_not_create_body
                     $illegal_value
                     $edge_not_found
                     $face_not_found

EXTRUDE   dx dy dz
          use:    create a Body by extruding an Xsect
          pops:   Xsect
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if Xsect is a SheetBody, then a SolidBody is created
                  if Xsect is a WireBody, then a SheetBody is created
                  if Xsect is a NodeBody, then a WireBody is created
                  sensitivity computed w.r.t. dx, dy, dz
                  computes Face sensitivities analytically
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  Attributes on Xsect are maintained
                  face-order is: (base), (end), feat1, ...
                  signals that may be thrown/caught:
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

FILLET    radius edgeList=0 listStyle=0
          use:    apply a fillet to a Body
          pops:   Body
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if listStyle==0
                     if previous operation is boolean, apply to all new Edges
                     edgeList=0 is the same as edgeList=[0;0]
                     edgeList is a multi-value Parameter or a semicolon-separated
                        list
                     pairs of edgeList entries are processed in order
                     pairs of edgeList entries are interpreted as follows:
                        col1  col2   meaning
                         =0    =0    add all Edges
                         >0    >0    add    Edges between iford=+icol1
                                                      and iford=+icol2
                         <0    <0    remove Edges between iford=-icol1
                                                      and iford=-icol2
                         >0    =0    add    Edges adjacent to iford=+icol1
                         <0    =0    remove Edges adjacent to iford=-icol1
                  else
                     edgeList contains Edge number(s)
                  sensitivity computed w.r.t. radius
                  sets up @-parameters
                  new Faces all receive the Branch's Attributes
                  face-order is based upon order that is returned from EGADS
                  signals that may be thrown/caught:
                     $illegal_argument
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

GETATTR   $pmtrName attrID global=0
          use:    store an Attribute value(s) in a LOCALVAR
          pops:   -
          pushes: -
          Notes:  pmtrName must be in form 'name', without subscripts
                  pmtrName must not start with '@'
                  pmtrName must not refer to an DESPMTR/CONPMTR Parameter
                  pmtrName will be marked as LOCALVAR (or OUTPMTR)
                  pmtrName is used directly (without evaluation)
                  the type of pmtrName is changed to match the result
                  if global==0, then
                     applies to Attributes on the selected Body
                  else
                     applies to global Attributes
                  if attrID is $_nattr_ then number of Attributes
                     will be retrieved into a scalar or indexed entry
                  if attrID is an integer (i), then the name of the
                     i'th (bias-1) Attribute will be retreived into a
                     string Parameter
                  Attributes are retrieved from last Body or from a Body,
                     Face, or Edge if it follows a SELECT statement
                  signals that may be thrown/caught:
                     $illegal_pmtr_index, $illegal_attribute

GROUP     nbody=0
          use:    create a Group of Bodys since Mark for subsequent
                     transformations
          pops:   Body1 ... Mark  -or-  Body1 ...
          pushes: Body1 ...
          notes:  Sketch may not be open
                  Solver may not be open
                  if nbody>0,   then nbody Bodys on stack are in Group
                  if nbody<0,   then Bodys are ungrouped
                  if no Mark on stack, all Bodys on stack are in Group
                  the Mark is removed from the stack
                  Attributes are set on all Bodys in Group
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

HOLLOW    thick=0 entList=0 listStyle=0
          use:    hollow out a SolidBody or SheetBody
          pops:   Body
          pushes: Body
          notes:  results can be unpredictable due to OpenCASCADE issues
                  Sketch may not be open
                  Solver may not be open
                  if SolidBody (radius is ignored)
                     if thick=0 and entList==0         (case A)
                         convert to SheetBody
                     if thick=0 and entList!=0         (case B)
                        convert to SheetBody without Faces in entList (if connected)
                     if thick>0 and entList==0         (case C)
                        larger offset Body is created
                     if thick<0 and entList==0         (case D)
                        smaller offset Body is created
                     if thick>0 and entList!=0         (case E)
                        hollow (removing entList) with new Faces inside  original Body
                     if thick<0 and entList!=0         (case F)
                        hollow (removing entList) with new Faces outside original Body
                  if a SheetBody with only one Face
                     if thick=0 and entList==0         (case G)
                        convert to WireBody (if connected)
                     if thick=0 and entList!=0         (case H)
                        convert to WireBody without Edges in entList (if connected)
                     if thick>0 and entList==0         (case I)
                        hollow with new Edges inside  original Body
                     if thick<0 and entList==0         (case J)
                        hollow with new Edges outside original Body
                     if thick>0 and entList!=0         (case K)
                        hollow (removing entList) with new Edges inside  original Body
                     if thick<0 and entList!=0         (case L)
                        hollow (removing entList) with new Edges outside original Body
                  if a SheetBody with multiple Faces
                     if thick=0 and entList!=0         (case M)
                        remove Faces in entList (if connected)
                     if thick>0 and entList==0         (case N)
                        hollow all Faces with new Edges inside original Faces
                     if thick>0 and entList!=0         (case P)
                        hollow Faces in entList with new Edges inside original Faces
                  entList is multi-valued Parameter, or a semicolon-separated list
                  if listStyle==0 and a SolidBody
                     pairs of entList entries are processed in order
                        the first  entry in a pair indicates the Body when
                           Face was generated (see first number in _body Attribute)
                        the second entry in a pair indicates the face-order (see
                           second number in _body Attribute)
                  otherwise
                     entries in entList are Edge or Face numbers
                  sensitivity computed w.r.t. thick
                  sets up @-parameters
                  new Faces all receive the Branch's Attributes
                  face-order is based upon order that is returned from EGADS
                  signals that may be thrown/caught:
                     $illegal_argument
                     $did_not_create_body
                     $insufficient_bodys_on_stack

IFTHEN    val1 $op1 val2 $op2=and val3=0 $op3=eq val4=0
          use:    execute or skip a Block of Branches
          pops:   -
          pushes: -
          notes:  works in combination with ELSEIF, ELSE, and ENDIF statements
                  op1 must be one of: lt LT le LE eq EQ ge GE gt GT ne NE
                  op2 must be one of: or OR and AND xor XOR
                  op3 must be one of: lt LT le LE eq EQ ge GE gt GT ne NE
                  if expression evaluates false, skip Block of Branches up
                     to next (matching) ELSEIF, ELSE, or ENDIF are skipped
                  cannot be followed by ATTRIBUTE or CSYSTEM

IMPORT    $filename bodynumber=1
          use:    import from filename
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  filename is used directly (without evaluation)
                  if filename starts with '$$/', use path relative to .csm file
                  if bodynumber=-1, then all Bodys are returned in one Group
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is based upon order in file
                  signals that may be thrown/caught:
                     $did_not_create_body
                     udp-specific code

INTERFACE $argName $argType default=0
          use:    defines an argument for a .udc file
          pops:   -
          pushes: -
          notes:  only allowed in a .udc file
                  must be placed before any executable statement
                  argType must be "in", "out", "dim", or "all"
                  argType="dim" is obsolete, use DIMENSION instead
                  if argType=="all", a new scope is not created (and
                                     $argName is ignored)
                  a string variable can be passed into UDC if default
                     is a string
                  a string varaible can be passed out of UDC
                  cannot be followed by ATTRIBUTE or CSYSTEM
                  signals that may be thrown/caught:
                     $pmtr_is_conpmtr

INTERSECT $order=none index=1 maxtol=0
          use:    perform Boolean intersection (Body2 & Body1)
          pops:   Body1 Body2
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if     Body1=SolidBody and Body2=SolidBody
                     create SolidBody that is common part of Body1 and Body2
                     if index=-1, then all Bodys are returned
                  elseif Body1=SolidBody and Body2=SheetBody
                     create SheetBody that is the part of Body2 that is
                        inside Body1
                     if index=-1, then all Bodys are returned
                  elseif Body1=SolidBody and Body2=WireBody
                     create WireBody that is the part of Body2 that is
                        inside Body1
                     if index=-1, then all Bodys are returned
                  elseif Body1=SheetBody and Body2=SolidBody
                     create SheetBody that is the part of Body1 that is
                        inside Body2
                     if index=-1, then all Bodys are returned
                  elseif Body1=SheetBody and Body2=SheetBody and Bodys are
                        co-planar
                     create SheetBody that is common part of Body1 and Body2
                     CURRENTLY NOT IMPLEMENTED
                  elseif Body1=SheetBody and Body2=SheetBody and Bodys are not
                        co-planar
                     create WireBody at the intersection of Body1 and Body2
                     CURRENTLY NOT IMPLEMENTED
                  elseif Body1=SheetBody and Body2=WireBody
                     create WireBody that is the part of Body2 that is
                        inside Body1
                     CURRENTLY NOT IMPLEMENTED
                  elseif Body1=WireBody and Body2=SolidBody
                     create WireBody that is the part of Body1 that is
                        inside Body2
                     if index=-1, then all Bodys are returned
                  elseif Body1=WireBody and Body2=SheetBody
                     create WireBody that is the part of Body1 that is
                        inside Body2
                     CURRENTLY NOT IMPLEMENTED
                  endif
                  if intersection does not produce at least index Bodys, an
                     error is returned
                  order may be one of:
                     none    same order as returned from geometry engine
                     xmin    minimum xmin   is first
                     xmax    maximum xmax   is first
                     ymin    minimum ymin   is first
                     ymax    maximum ymax   is first
                     zmin    minimum zmin   is first
                     zmax    maximum zmax   is first
                     amin    minimum area   is first
                     amax    maximum area   is first
                     vmin    minimum volume is first
                     vmax    maximum volume is first
                  order is used directly (without evaluation)
                  if maxtol>0, then tolerance can be relaxed until successful
                  if maxtol<0, then use -maxtol as only tolerance to use
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $did_not_create_body
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

JOIN      toler=0 toMark=0
          use:    join two Bodys at a common Edge or Face
          pops:   Body1 Body2
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if toMark=1 and all Bodys to Mark are SheetBodys
                     create SheetBody
                  elseif toMark=1 and all Bodys to Mark are WireBodys
                     create WireBody
                  elseif Body1=SolidBody and Body2=SolidBody
                     create SolidBody formed by joining Body1 and Body2 at
                        common Faces
                  elseif Body1=SheetBody and Body2=SheetBody
                     create SheetBody formed by joining Body1 and Body2 at
                        common Edges
                  elseif Body1=WireBody and Body2=WireBody
                     create WireBody formed by joining Body1 and Body2 at
                        common end Node
                  endif
                  change in v1.19: if common Edges are not found, return
                                   $edge_not_found
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $created_too_many_bodys
                     $did_not_create_body
                     $edge_not_found
                     $face_not_found
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

LBOUND    $pmtrName bounds
          use:    defines a lower bound for a DESPMTR or CFGPMTR
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  statement may not be used in a function-type .udc file
                  if value of Parameter is smaller than bounds, a warning is
                     generated
                  pmtrName must have been defined previously by DESPMTR
                     statement
                  pmtrName can be in form 'name' or 'name[irow,icol]'
                  pmtrName must not start with '@'
                  pmtrName is used directly (without evaluation)
                  irow and icol cannot contain a comma or open bracket
                  if irow is a colon (:), then all rows    are input
                  if icol is a colon (:), then all columns are input
                  pmtrName[:,:] is equivalent to pmtrName
                  bounds cannot refer to any other Parameter
                  bounds are defined across rows, then across columns
                  if bounds has more entries than needed, extra bounds
                     are lost
                  if bounds has fewer entries than needed, last bound
                     is repeated
                  any previous bounds are overwritten
                  does not create a Branch
                  cannot be followed by ATTRIBUTE or CSYSTEM

LINSEG    x y z
          use:    create a new line segment, connecting the previous
                     and specified points
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  sensitivity computed w.r.t. x, y, z
                  signals that may be thrown/caught:

LOFT      smooth
          use:    create a Body by lofting through Xsects since Mark
          pops:   Xsect1 ... Mark
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  all Xsects must have the same number of Segments
                  if Xsect is a SheetBody, then a SolidBody is created
                  if Xsect is a WireBody, then a SheetBody is created
                  the Faces all receive the Branch's Attributes
                  Attributes on Xsects are not maintained
                  face-order is: (base), (end), feat1, ...
                  if NINT(smooth)=1, then sections are smoothed
                  the first and/or last Xsect can be a point

                  LOFT (through OpenCASCADE) is not very robust
                  use BLEND or RULE if possible
                  sets up @-parameters
                  MAY BE DEPRECATED (use RULE or BLEND)
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

MACBEG    imacro
          use:    marks the start of a macro
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  imacro must be between 1 and 100
                  cannot overwrite a previous macro
                  cannot be followed by ATTRIBUTE or CSYSTEM
                  MAY BE DEPRECATED (use UDPRIM)

MACEND
          use:    ends a macro
          pops:   -
          pushes: -
          notes:  cannot be followed by ATTRIBUTE or CSYSTEM
                  MAY BE DEPRECATED (use UDPRIM)

MARK
          use:    used to identify groups such as in RULE, BLEND, or GROUP
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  cannot be followed by ATTRIBUTE or CSYSTEM

MESSAGE   $text $schar=_
          use:    generate a message to be displayed
          pops:   -
          pushes: -
          notes:  schar must contain a single character
                  any character in text that matches schar will
                     be converted to a space
                  cannot be followed by ATTRIBUTE or CSYSTEM

MIRROR    nx ny nz dist=0
          use:    mirrors Group on top of Stack
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  normal of the mirror plane is given by nx,ny,nz
                  mirror plane is dist from origin
                  sensitivity computed w.r.t. nx, ny, nz, dist
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

NAME      $branchName
          use:    names the entry on top of Stack
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  does not create a Branch

OUTPMTR   $pmtrName
          use:    define an OUTPMTR
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  statement may not be used in a .udc file
                  pmtrName must be in form 'name'
                  pmtrName must not start with '@'
                  pmtrName will be marked as OUTPMTR
                  pmtrName is used directly (without evaluation)
                  does not create a Branch
                  cannot be followed by ATTRIBUTE or CSYSTEM

PATBEG    $pmtrName ncopy
          use:    execute a Block of Branches ncopy times
          pops:   -
          pushes: -
          notes:  Solver may not be open
                  Block contains all Branches up to matching PATEND
                  pmtrName must not start with '@'
                  pmtrName takes values from 1 to ncopy (see below)
                  pmtrName is used directly (without evaluation)
                  cannot be followed by ATTRIBUTE or CSYSTEM

PATBREAK  expr
          use:    break out of inner-most Patbeg Block if expr>0
          pops:   -
          pushes: -
          notes:  Solver may not be open
                  must be in a Patbeg Block
                  skip to Branch after matching PATEND if expr>0
                  cannot be followed by ATTRIBUTE or CSYSTEM

PATEND
          use:    designates the end of a Patbeg Block
          pops:   -
          pushes: -
          notes:  Solver may not be open
                  inner-most Block must be a Patbeg Block
                  closes Patbeg Block
                  cannot be followed by ATTRIBUTE or CSYSTEM

POINT     xloc yloc zloc
          use:    create a single point Body
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. xloc, yloc, zloc
                  computes Node sensitivity analytically
                  sets up @-parameters

PROJECT   x y z dx dy dz useEdges=0
          use:    find the first projection from given point in given
                     direction
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  if useEdges!=1
                      look for intersections with Faces and overwrite @iface
                  else
                      look for intersections with Edges and overwrite @iedge
                  endif
                  over-writes the following @-parameters: @xcg, @ycg, and @zcg
                  cannot be followed by ATTRIBUTE or CSYSTEM
                  signals that may be thrown/caught:
                     $face_not_found
                     $insufficient_bodys_on_stack

RECALL    imacro
          use:    recalls copy of macro from a storage location imacro
          pops:   -
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  storage location imacro must have been previously filled by
                     a MACBEG statement
                  MAY BE DEPRECATED (use UDPRIM)

REORDER   ishift iflip=0
          use:    change the order of Edges in a Body
          pops:   Body1
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  it is generally better to use reorder argument in
                     RULE and BLEND than this command
                  Body1 must be either WireBody or SheetBody Body
                  Body1 must contain 1 Loop
                  if the Loop is open, ishift must be 0
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

RESTORE   $name index=0
          use:    restores Body(s) that was/were previously stored
          pops:   -
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  $name is used directly (without evaluation)
                  if $name is . (dot), then duplicate Body on stack
                  if index<0, get all Bodys that match $name
                  sets up @-parameters
                  error results if nothing has been stored in name
                  the Faces all receive the Branch's Attributes
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack
                     $name_not_found

REVOLVE   xorig yorig zorig dxaxis dyaxis dzaxis angDeg
          use:    create a Body by revolving an Xsect around an axis
          pops:   Xsect
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if Xsect is a SheetBody, then a SolidBody is created
                  if Xsect is a WireBody, then a SheetBody is created
                  sensitivity computed w.r.t. xorig, yorig, zorig, dxaxis,
                     dyaxis, dzaxis, andDeg
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  Attributes on Xsect are maintained
                  face-order is: (base), (end), feat1, ...
                  signals that may be thrown/caught:
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

ROTATEX   angDeg yaxis=0 zaxis=0
          use:    rotates Group on top of Stack around an axis that
                  passes through (0, yaxis, zaxis) and is parallel
                  to the x-axis
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. angDeg, yaxis, zaxis
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

ROTATEY   angDeg zaxis=0 xaxis=0
          use:    rotates Group on top of Stack around an axis that
                  passes through (xaxis, 0, zaxis) and is parallel
                  to the y-axis
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. angDeg, zaxis, xaxis
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

ROTATEZ   angDeg xaxis=0 yaxis=0
          use:    rotates Group on top of Stack around an axis that
                  passes tgrough (xaxis, yaxis, 0) and is parallel
                  to the z-axis
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. angDeg, xaxis, yaxis
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

RULE      reorder=0 periodic=0
          use:    create a Body by creating ruled surfaces thru Xsects
                     since Mark
          pops:   Xsect1 ... Mark
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if reorder!=0 then Xsects are reordered to minimize Edge
                     lengths
                  first Xsect is unaltered if reorder>0
                  last  Xsect is unaltered if reorder<0
                  all Xsects must have the same number of Edges
                  if all Xsects are NodeBodys
                     a WireBody is created
                  elseif all Xsects are WireBodys (or a NodeBody at one end)
                     a SheetBody is created
                  else
                     a SolidBody is created
                  the first and/or last Xsect can be a NodeBody
                  if periodic=1 then connect the first and last Xsects
                  computes Face sensitivities analytically
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  Attributes on Xsects are maintained
                  face-order is: (base), (end), feat1:part1,
                     feat1:part2, ... feat2:part1, ...
                  signals that may be thrown/caught:
                     $did_not_create_body
                     $error_in_bodys_on_stack
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

SCALE     fact xcent=0 ycent=0 zcent=0
          use:    scales Group on top of Stack around given point
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  (xcent,ycent,zcent are not yet implemented)
                  sensitivity computed w.r.t. fact
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

SELECT    $type arg1 ...
          use:    selects entity for which @-parameters are evaluated
          pops:   -
          pushes: -
          notes:  if     arguments are: "body"
                     sets @seltype to -1
                     sets @selbody to @nbody
                     sets @sellist to -1
                  elseif arguments are: "body ibody"
                     sets @seltype to -1
                     sets @selbody to ibody
                     sets @sellist to -1
                  elseif arguments are: "body -n"
                     sets @seltype to -1
                     sets @selbody to the nth from the top of the stack
                     sets @sellist to -1
                  elseif arguments are: "body attrName1    attrValue1
                                              attrName2=$* attrValue2=$*
                                              attrName3=$* attrValue3=$*"
                     sets @seltype to -1
                     uses @selbody to Body that match all Attributes
                     sets @sellist to -1
                  elseif arguments are: "face"
                     sets @seltype to 2
                     uses @selbody
                     sets @sellist to all Faces
                  elseif arguments are: "face iface"
                     sets @seltype to 2
                     uses @selbody
                     sets @sellist to iface
//$$$                  elseif arguments are: "face 0 iford1" or
//$$$                                        "face ibody1 0"
//$$$                     sets @seltype to 2
//$$$                     uses @selbody
//$$$                     sets @sellist with Faces in @selbody that matches ibody1/iford1
//$$$                                   (with 0 being treated as a wildcard)
                  elseif arguments are: "face ibody1 iford1 iseq=1"
                     sets @seltype to 2
                     uses @selbody
                     sets @sellist with Face in @selbody that matches ibody1/iford1
                     (note that 0 can be used as a wildcard in any numeric field)
                     (note that if any fields are 0, the default iseq=0)
                  elseif arguments are: "face xmin xmax ymin ymax zmin zmax"
                     sets @seltype to 2
                     uses @selbody
                     sets @sellist to Faces whose bboxs are completely in given range
                  elseif arguments are: "face attrName1    attrValue1
                                              attrName2=$* attrValue2=$*
                                              attrName3=$* attrValue3=$*"
                     sets @seltype to 2
                     uses @selbody
                     sets @sellist to Faces in @selbody that match all Attributes
                  elseif arguments are: "edge"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to all Edges
                  elseif arguments are: "edge iedge"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to iedge
//$$$                  elseif arguments are: "edge 0 iford1 ibody2 iford2" or
//$$$                                        "edge ibody1 0 ibody2 iford2" or
//$$$                                        "edge ibody1 iford1 0 iford2" or
//$$$                                        "edge ibody1 iford1 ibody2 0" or
//$$$                                        "edge ibody1 0 ibody2 0"
//$$$                     sets @seltype to 1
//$$$                     uses @selbody
//$$$                     sets @sellist to Edge in @selbody that adjoins Faces
//$$$                        ibody1/iford1 and ibody2/iford2 (with 0 being treated as wildcard)
                  elseif arguments are: "edge ibody1 iford1 ibody2 iford2 iseq=1"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to Edge in @selbody that adjoins Faces
                        ibody1/iford1 and ibody2/iford2
                     (note that 0 can be used as a wildcard in any numeric field)
                     (note that if any fields are 0, the default iseq=0)
                  elseif arguments are: "edge xmin xmax ymin ymax zmin zmax"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to Edges whose bboxs are completely in given range
                  elseif arguments are: "edge attrName1    attrValue1
                                              attrName2=$* attrValue2=$*
                                              attrName3=$* attrValue3=$*"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to Edges in @selbody that match all Attributes
                  elseif arguments are: "edge x y z"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to Edge whose center is closest to (x,y,z)
                  elseif arguments are: "loop iface iloop"
                     sets @seltype to 1
                     uses @selbody
                     sets @sellist to Edges in order in the Loop
                  elseif arguments are: "node"
                     sets @seltype to 0
                     uses @selbody
                     sets @sellist to all Nodes
                  elseif arguments are: "node inode"
                     sets @seltype to 0
                     uses @selbody
                     sets @sellist to inode
                  elseif arguments are: "node x y z"
                     sets @seltype to 0
                     uses @selbodt
                     sets @sellist to Node closest to (x,y,z)
                  elseif arguments are: "node xmin xmax ymin ymax zmin zmax"
                     sets @seltype to 0
                     uses @selbody
                     sets @sellist to Nodes that are in given range
                  elseif arguments are: "node attrName1    attrValue1
                                              attrName2=$* attrValue2=$*
                                              attrName3=$* attrValue3=$*"
                     sets @seltype to 0
                     uses @selbody
                     sets sellist to Nodes in @selbody that match all Attributes
                  elseif arguments are: "add attrName1    attrValue1
                                             attrName2=$* attrValue2=$*
                                             attrName3=$* attrValue3=$*"
                     uses @seltype
                     uses @selbody
                     appends to @selList the Nodes/Edges/Faces that match all Attributes
                  elseif arguments are: "add ibody1 iford1 iseq=1" and @seltype is 2
                     uses @selbody
                     appends to @sellist the Face in @selbody that matches ibody1/iford1
                     (a 0 matches ibody1=0 amd/or iford1=0)
                  elseif arguments are: "add ibody1 iford1 ibody2 iford2 iseq=1" and @seltype is 1
                     uses @selbody
                     appends to @sellist the Edge in @selbody that adjoins Faces
                     (a 0 matches ibody1=0, iford1=0, ibody2=0, and/or iford2=0)
                  elseif arguments are: "add iface" and @seltype is 2
                     uses @selbody
                     appends to @sellist Face iface in @selbody
                     (wildcarding is not allowed)
                  elseif arguments are: "add iedge" and @seltype is 1
                     uses @selbody
                     appends to @sellist Edge iedge in @selbody
                     (wildcarding is not allowed)
                  elseif arguments are: "add inode" and @seltype is 0
                     uses @selbody
                     appends to @sellist Node inode in @selbody
                     (wildcarding is not allowed)
                  elseif arguments are: "sub attrName1    attrValue1
                                             attrName2=$* attrValue2=$*
                                             attrName3=$* attrValue3=$*"
                     uses @seltype
                     uses @selbody
                     removes from @sellist the Nodes/Edges/Faces that match all Attributes
                  elseif arguments are: "sub ibody1 iford1 iseq=1" and @seltype is 2
                     uses @selbody
                     removes from @sellist the Face in @selbody that matches ibody1/iford1
                     (a 0 matches ibody1=0 amd/or iford1=0)
                  elseif arguments are: "sub ibody1 iford1 ibody2 iford2 iseq=1" and @seltype is 1
                     uses @selbody
                     removes from @sellist the Edge in @selbody that adjoins Faces
                     (a 0 matches ibody1=0, iford1=0, ibody2=0, and/or iford2=0)
                  elseif arguments are: "sub ient" and ient is in @sellist
                     removes from @sellist ient
                     (wildcarding is not allowed)
                  elseif arguments are: "sort $key"
                     sorts @sellist based upon $key which can be: $xmin, $ymin, $zmin,
                        $xmax, $ymax, $zmax, $xcg, $ycg, $zcg, $area, or $length

                  Face specifications are stored in _faceID Attribute
                  Edge specifications are stored in _edgeID Attribute
//                Node specifications are stored in _nodeID Attribute
                  iseq selects from amongst multiple Faces/Edges/Nodes that
                     match the ibody/iford specifications
                  attrNames and attrValues can be wild-carded with $*
                  avoid using forms "SELECT face iface" and "SELECT edge iedge"
                     since iface and iedge are not guaranteed to be the same during
                     rebuilds or on different OpenCASCADE versions or computers
                  sets up @-parameters
                  cannot be followed by CSYSTEM
                  signals that may be thrown/caught:
                     $body_not_found
                     $edge_not_found
                     $face_not_found
                     $node_not_found

SET       $pmtrName exprs
          use:    define or update a (redefinable) local or output variable
          pops:   -
          pushes: -
          notes:  Solver may not be open
                  pmtrName can be in form 'name', 'name[irow]', or 'name[irow,icol]'
                  pmtrName must not start with '@'
                  pmtrName must not refer to an DESPMTR/CONPMTR Parameter
                  pmtrName will be marked as LOCALVAR (or OUTPMTR)
                  pmtrName is used directly (without evaluation)
                  irow and icol cannot contain a comma or open bracket
                  if exprs has multiple values (separated by ;), then
                     any subscripts in pmtrName are ignored
                  if exprs starts with $ (or evaluates to a string), then any
                     subscripts in pmtrName are ignored and a string value is defined
                  if exprs is the name of a multi-valued parameter, it is
                      treated as if its values were listed as a semi-colon-
                      separated list
                  if pmtrName is in the form 'name' or 'name[0]' or 'name[0,0]'
                     if exprs is longer than Parameter size, extra exprs are lost
                     if exprs is shorter than Parameter size, last expr is repeated
                  if pmtrName is in the form 'name[irow]' or 'name[irow,0]', then
                     the irow'th element is defined (where elements are counted
                     across rows)
                  if pmtrName is in the form 'name[irow,icol]' and irow is between
                     1 and nrow and icol is between 1 and ncol, then the
                     [irow,icol]th element is set
                  if no Bodys have been created yet
                     associated ATTRIBUTEs are global Attributes
                  otherwise
                     cannot be followed by ATTRIBUTE
                  cannot be folowed by CSYSTEM

SKBEG     x y z relative=0
          use:    start a new Sketch with the given point
          pops:   -
          pushes: -
          notes:  opens Sketch
                  Solver may not be open
                  if relative=1, then all values in sketch are relative to x,y,z
                  sensitivity computed w.r.t. x, y, z
                  cannot be followed by ATTRIBUTE or CSYSTEM

SKCON     $type index1 index2=-1 $value=0
          use:    creates a Sketch constraint
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  may only follow SKVAR or another SKCON statement
                  $type
                     X  ::x[index1]=value
                     Y  ::y[index1]=value
                     Z  ::z[index1]=value
                     P  segments adjacent to point index1 are perpendicular
                     T  segments adjacent to point index1 are tangent
                     A  segments adjacent to point index1 have
                                                           angle=$value (deg)
                     W  width:  ::x[index2]-::x[index1]=value  if plane==xy
                                ::y[index2]-::y[index1]=value  if plane==yz
                                ::z[index2]-::z[index1]=value  if plane==zx
                     D  depth:  ::y[index2]-::y[index1]=value  if plane==xy
                                ::z[index2]-::z[index1]=value  if plane==zx
                                ::x[index2]-::x[index1]=value  if plane==zx
                     H  segment from index1 and index2 is horizontal
                     V  segment from index1 and index2 is vertical
                     I  segment from index1 and index2 has
                                                     inclination=$value (deg)
                     L  segment from index1 and index2 has length=$value
                     R  cirarc  from index1 and index2 has radius=$value
                     S  cirarc  from index1 and index2 has sweep=$value (deg)
                  index=1 refers to point in SKBEG statement
                  $value can include the following variables
                     ::x[i]  X-coordinate of point i
                     ::y[i]  Y-coordinate of point i
                     ::z[i]  Z-coordinate of point i
                     ::d[i]  dip associated with segment starting at point i
                  $value can include the following shorthands
                     ::L[i]  length      of segment starting at point i
                     ::I[i]  inclination of segment starting at point i  (degrees)
                     ::R[i]  radius of arc          starting at point i
                     ::S[i]  sweep  of rc           starting at point i  (degrees)
                  cannot be followed by ATTRIBUTE or CSYSTEM

SKEND     wireonly=0
          use:    completes a Sketch
          pops:   -
          pushes: Sketch
          notes:  Sketch must be open
                  Solver may not be open
                  if Sketch contains SKVAR/SKCON, then Sketch variables are
                     updated first
                  if wireonly=0, all LINSEGs and CIRARCs must be x-, y-, or
                     z-co-planar
                  if Sketch is     closed and wireonly=0,
                     then a SheetBody is created
                  if Sketch is     closed and wireonly=1,
                     then a WireBody  is created
                  if Sketch is not closed,
                     then a WireBody  is created
                  if SKEND immediately follows SKBEG, then a NODE is created
                     (which can be used at either end of a LOFT or BLEND)
                  closes Sketch
                  new Face receives the Branch's Attributes
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $underconstrained
                     $overconstrained
                     $not_converged
                     $colinear_sketch_points
                     $non_coplnar_sketch_points
                     $self_intersecting

SKVAR     $type valList
          use:    create multi-valued Sketch variables and their initial
                     values
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  may only follow SKBEG statement
                  $type
                     xy valList contains ::x[1]; ::y[1]; ::d[1]; ::x[2]; ...
                     yz valList contains ::y[1]; ::z[1]; ::d[1]; ::y[2]; ...
                     zx valList contains ::z[1]; ::x[1]; ::d[1]; ::z[2]; ...
                  valList is a semicolon-separated list
                  valList must end with a semicolon
                  the number of entries in valList is taken from number of
                     semicolons
                  the number of entries in valList must be evenly divisible by 3
                  enter :d[i] as zero for LINSEGs
                  values of ::x[1], ::y[1], and ::z[1] are overwritten by
                     values in SKBEG
                  cannot be followed by ATTRIBUTE or CSYSTEM

SOLBEG    $varList
          use:    starts a Solver Block
          pops:   -
          pushes: -
          notes:  Solver must not be open
                  opens the Solver
                  varList is a list of semicolon-separated LOCALVARs
                  varList must end with a semicolon
                  cannot be followed by ATTRIBUTE or CSYSTEM

SOLCON    $expr
          use:    constraint used to set Solver parameters
          pops:   -
          pushes: -
          notes:  Sketch must not be open
                  Solver must be open
                  SOLEND will drive expr to zero
                  cannot be followed by ATTRIBUTE or CSYSTEM

SOLEND
          use:    designates the end of a Solver Block
          pops:   -
          pushes: -
          notes:  Sketch must not be open
                  inner-most Block must be a Solver Block
                  adjust parameters to drive constraints to zero
                  closes Solver Block
                  cannot be followed by ATTRIBUTE or CSYSTEM

SPHERE    xcent ycent zcent radius
          use:    create a sphere Body
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. xcent, ycent, zcent, radius
                  computes Face, Edge, and Node sensitivities analytically
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is: ymin, ymax
                  signals that may be thrown/caught:
                     $illegal_value

SPLINE    x y z
          use:    add a point to a spline
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  sensitivity computed w.r.t. x, y, z
                  signals that may be thrown/caught:

SSLOPE    dx dy dz
          use:    define the slope at the beginning or end of a SPLINE
          pops:   -
          pushes: -
          notes:  Sketch must be open
                  Solver may not be open
                  for defining slope at beginning:
                      must not follow a SPLINE statement
                      must    precede a SPLINE statement
                  for definiing slope at end:
                      must      follow a SPLINE statement
                      must not precede a SPLINE statement
                  dx, dy, and dz must not all be zero
                  sensitivity computed w.r.t. x, y, z
                  signals that may be thrown/caught:
                     $illegal_value

STORE     $name index=0 keep=0
          use:    stores Group on top of Stack
          pops:   any
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  $name is used directly (without evaluation)
                  if index<0, use first available index
                  previous Group in name/index is overwritten
                  if $name=.   then Body is popped off stack
                                    but not actually stored
                  if $name=..  then pop Bodys off stack back
                                    to the Mark
                  if $name=... then the stack is cleared
                  if keep==1, the Group is not popped off stack
                  cannot be followed by ATTRIBUTE or CSYSTEM
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

SUBTRACT  $order=none index=1 maxtol=0
          use:    perform Boolean subtraction (Body2 - Body1)
          pops:   Body1 Body2
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if     Body1=SolidBody and Body2=SolidBody
                     create SolidBody that is the part of Body1 that is
                        outside Body2
                     if index=-1, then all Bodys are returned
                  elseif Body1=SolidBody and Body2=SheetBody
                     create SolidBody that is Body1 scribed with Edges at
                        intersection with Body2
                  elseif Body1=SheetBody and Body2=SolidBody
                     create SheetBody that is part of Body1 that is
                        outside Body2
                     if index=-1, then all Bodys are returned
                  elseif Body1=SheetBody and Body2=SheetBody
                     create SheetBody that is Body1 scribed with Edges at
                        intersection with Body2
                  elseif Body1=WireBody and Body2=SolidBody
                     create WireBody that is part of Body1 that is outside Body2
                     CURRENTLY NOT IMPLEMENTED
                  elseif Body1=WireBody and Body2=SheetBody
                     create WireBody that is Body1 scribed with Nodes at
                        intersection with Body2
                     CURRENTLY NOT IMPLEMENTED
                  elseif Body2=NodeBody
                     split Edges in Body1 at Body2
                  endif
                  if subtraction does not produce at least index Bodys,
                     an error is returned
                  order may be one of:
                     none    same order as returned from geometry engine
                     xmin    minimum xmin   is first
                     xmax    maximum xmax   is first
                     ymin    minimum ymin   is first
                     ymax    maximum ymax   is first
                     zmin    minimum zmin   is first
                     zmax    maximum zmax   is first
                     amin    minimum area   is first
                     amax    maximum area   is first
                     vmin    minimum volume is first
                     vmax    maximum volume is first
                  if maxtol>0, then tolerance can be relaxed until successful
                  if maxtol<0, then use -maxtol as only tolerance to use
                  sets up @-parameters
                  order is used directly (without evaluation)
                  signals that may be thrown/caught:
                     $did_not_create_body
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

SWEEP
          use:    create a Body by sweeping an Xsect along an Xsect
          pops:   Xsect1 Xsect2
          pushes: Body
          notes:  results can be unpredictable due to OpenCASCADE issues
                  Sketch may not be open
                  Solver may not be open
                  Xsect1 must be either a SheetBody or WireBody
                  Xsect2 must be a WireBody
                  if Xsect2 is not slope-continuous, result may not be
                     as expected
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  Attributes on Xsect are maintained
                  face-order is: (base), (end), feat1a, feat1b, ...
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack

THROW     sigCode
          use:    set current signal to sigCode
          pops:   -
          pushes: -
          notes:  skip statements until a matching CATBEG Branch is found
                  sigCode>0 are usually user-generated signals
                  sigCode<0 are usually system-generated signals
                  cannot be followed by ATTRIBUTE or CSYSTEM

TORUS     xcent ycent zcent dxaxis dyaxis dzaxis majorRad minorRad
          use:    create a torus Body
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. xcent, ycent, zcent, dxaxis,
                     dyaxis, dzaxis, majorRad, minorRad
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is: xmin/ymin, xmin/ymax, xmax/ymin, xmax/ymax
                  signals that may be thrown/caught:
                     $illegal_value

TRANSLATE dx dy dz
          use:    translates Group on top of Stack
          pops:   any
          pushes: any
          notes:  Sketch may not be open
                  Solver may not be open
                  sensitivity computed w.r.t. dx, dy, dz
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $insufficient_bodys_on_stack

UBOUND    $pmtrName bounds
          use:    defines an upper bound for a DESPMTR or CFGPMTR
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  statement may not be used in a function-type .udc file
                  if value of Parameter is larger than bounds, a warning is
                     generated
                  pmtrName must have been defined previously by DESPMTR
                     statement
                  pmtrName can be in form 'name' or 'name[irow,icol]'
                  pmtrName must not start with '@'
                  pmtrName is used directly (without evaluation)
                  irow and icol cannot contain a comma or open bracket
                  if irow is a colon (:), then all rows    are input
                  if icol is a colon (:), then all columns are input
                  pmtrName[:,:] is equivalent to pmtrName
                  bounds cannot refer to any other Parameter
                  bounds are defined across rows, then across columns
                  if bounds has more entries than needed, extra bounds
                     are lost
                  if bounds has fewer entries than needed, last bound
                     is repeated
                  any previous bounds are overwritten
                  does not create a Branch
                  cannot be followed by ATTRIBUTE or CSYSTEM

UDPARG    $primtype $argName1 argValue1 $argName2 argValue2 ...
                    $argName3 argValue3 $argName4 argValue4
          use:    pre-set arguments for next UDPRIM statement
          pops:   -
          pushes: -
          notes:  Sketch may not be open
                  Solver may not be open
                  there can be no statements except other UDPARGs before the
                     next matching UDPRIM
                  primtype determines the type of primitive
                  primtype must match primtype of next UDPRIM statement
                  primtype is used directly (without evaluation)
                  arguments are specified in name/value pairs and are
                      not positional
                  argName#  is used directly (without evaluation)
                  argValue# is used directly if it starts with '$', otherwise it
                     is evaluated
                  if argValue starts with '$$/', use path relative to .csm file
                  arguments for following UDPRIM statement are evaluated
                     in the order they are encountered (UDPARG first)
                  sensitivity computed w.r.t. argValue1, argValue2, argValue3,
                     argValue4
                  cannot be followed by ATTRIBUTE or CSYSTEM

UDPRIM    $primtype $argName1 argValue1 $argName2 argValue2 ...
                    $argName3 argValue3 $argName4 argValue4
          use:    create a Body by executing a UDP, UDC, or UDF
          pops:   -
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  primtype  determines the type of primitive and the number of
                     argName/argValue pairs
                  if primtype begins with a letter
                     then a compiled udp whose name is primtype.so is used
                  if primtype starts with a /
                     then a .udc file in the current directory will be used
                  if primtype starts with $/
                     then a .udc file in the parent (.csm or .udc)
                     directory will be used
                  if primtype starts with $$/
                     then a .udc file in ESP_ROOT/udc will be used
                  primtype  is used directly (without evaluation)
                  arguments are specified in name/value pairs and are
                      not positional
                  argName#  is used directly (without evaluation)
                  argValue# is used directly if it starts with '$', otherwise it
                     is evaluated
                  if argValue# is <<, use data to matching >> as inline file
                  if argValue# starts with '$$/', use path relative to .csm file
                  extra arguments can be set with UDPARG statement
                  when called to execute a .udc file:
                     the level is incremented
                     LOCALVARs are created for all INTERFACE stmts
                        for "in"  the value is set to its default
                        for "out" the value is set to its default
                        for "dim" an array is created (of size=value) with
                           value=dot=0
                     the associated UDPARG and UDPRIM statements are processed
                        in order
                        if argName matches a Parameter created by an INTERFACE
                           statement
                           if argValueX matches the name of a Parameter at
                              level-1
                              the values are copied into the new Parameter
                           else
                              argValueX is evalued and stored in the new
                                 Parameter
                        else
                           an error is returned
                     the statements in the .udc are executed until an END
                        statement
                        a SET statement either creates a new Parameter or
                           overwrites a value
                     during the execution of the END statement
                        for values associated with an INTERFACE "out" statement
                           the value is copied to the appropriate @@-parameter
                              (at level-1)
                        all Parameters at the current level are destroyed
                        the level is decremented
                  sensitivity computed w.r.t. argValue1, argValue2, argValue3,
                     argValue4
                  computes Face and Edge sensitivities analytically (if supplied
                     by the udp)
                  sets up @-parameters
                  the Faces all receive the Branch's Attributes
                  face-order is based upon order returned from UDPRIM
                  signals that may be thrown/caught:
                     $did_not_create_body
                     $insufficient_bodys_on_stack
                     udp-specific code
                  see udp documentation for full information

UNION     toMark=0 trimList=0 maxtol=0
          use:    perform Boolean union
          pops:   Body1 Body2  -or-  Body1 ... Mark
          pushes: Body
          notes:  Sketch may not be open
                  Solver may not be open
                  if     toMark=1
                     create SolidBody that is combination of SolidBodys
                        since Mark
                  elseif Body1=SolidBody and Body2=SolidBody
                     if trimList=0
                        create SolidBody that us combination of Body1 and Body2
                     else
                        create SolidBody that is trimmed combination of Body1
                           and Body2
                        trimList contains x;y;z;dx;dy;dz where
                           (x,y,z) is inside the Body to be trimmed
                           (dx,dy,dz) is step toward the trimming Body
                     endif
                  elseif Body1=SheetBody and Body2=SheetBody
                     create SheetBody that is the combination of Bodys with
                        possible new Edges
                  endif
                  if maxtol>0, then tolerance can be relaxed until successful
                  if maxtol<0, then use -maxtol as only tolerance to use
                  sets up @-parameters
                  signals that may be thrown/caught:
                     $did_not_create_body
                     $illegal_value
                     $insufficient_bodys_on_stack
                     $wrong_types_on_stack
*/

/*
 ************************************************************************
 *                                                                      *
 * Number and string rules                                              *
 *                                                                      *
 ************************************************************************
 */

/*
Numbers:
    start with a digit or decimal (.)
    followed by zero or more digits and/or decimals (.)
    there can be at most one decimal in a number
    optionally followed by an e, e+, e-, E, E+, or E-
    if there is an e or E, it must be followed by one or more digits
Strings:
    introduced with a dollar sign ($) that is not part of the value
    followed by one to 128 characters from the set
       letter                     a-z or A-Z
       digit                      0-9
       at-sign                    @
       underscore                 _
       colon                      :
       semicolon                  ;
       dollar-sign                $
       period                     .
       escaped comma              ',
       escaped plus               '+
       minus                      -
       star                       *
       slash                      /
       caret                      ^
       question                   ?
       percent                    %
       open-parenthesis           (
       escaped close-parenthesis  ')
       open-bracket               [
       close-bracket              ]
       open-brace                 {
       close-brace                }
       less-than                  <
       greater-than               >
       equal                      =
    the following characters are not allowed in strings
       apostrophe                 '  (except to escape ', '+ or ') )
       quotation                  "
       hashtag                    #
       backslash                  \
       vertical bar               |
       tilde                      ~
       ampersand                  &
       exclamation                !
*/

/*
 ************************************************************************
 *                                                                      *
 * Parameter rules                                                      *
 *                                                                      *
 ************************************************************************
 */

/*
Valid names:
    start with a letter, colon(:), or at-sign(@)
    contains letters, digits, at-signs(@), underscores(_), and colons(:)
    contains fewer than 64 characters
    names that start with an at-sign cannot be set by a CONPMTR, DESPMTR,
       SET, PATBEG, or GETATTR statement
    if a name has a dot-suffix, a property of the name (and not its
        value) is returned
       x.nrow   number of rows     in x or 0 if a string
       x.ncol   number of columns  in x or 0 if a string
       x.size   number of elements in x (=x.nrow*x.ncol) or
                     length of string x
       x.sum    sum of elements    in x
       x.norm   norm of elements   in x (=sqrt(x[1]^2+x[2]^2+...))
       x.min    minimum value      in x
       x.max    maximum value      in x
       x.dot    velocity           of x

Array names:
    basic format is: name[irow,icol] or name[ielem]
    name must follow rules above
    irow, icol, and ielem must be valid expressions
    irow, icol, and ielem start counting at 1
    values are stored across rows ([1,1], [1,2], ..., [2,1], ...)

Types:
    DESPMTR
        declared by a DESPMTR statement only at the top level
        if a scalar, declared and defined by a DESPMTR statement
        if an array, declared by a DIMENSION statement
            values defined by one or more DESPMTR statements
        each value can only be defined in one DESPMTR statement
        values are not over-written by subsequent DESPMTR statements
        can have an optional lower bound
        can have an optional upper bound
        is only available at the top level
        can be set  outside ocsmBuild by a call to ocsmSetValu
        can be read outside ocsmBuild by a call to ocsmGetValu
        can be used to find sensitivities
    CFGPMTR
        declared by a CFGPMTR statement only at the top level
        if a scalar, declared and defined by a CFGPMTR statement
        if an array, declared by a DIMENSION statement
            values defined by one or more CFGPMTR statements
        each value can only be defined in one CFGPMTR statement
        values are not over-written by subsequent CFGPMTR statements
        can have an optional lower bound
        can have an optional upper bound
        is only available at the top level
        can be set  outside ocsmBuild by a call to ocsmSetValu
        can be read outside ocsmBuild by a call to ocsmGetValu
    CONPMTR
        declared by a CONPMTR statement at the top level or
            via ocsmLoadDict
        cannot be declared in a DIMENSION statment
        is available everywhere
        can be set  outside ocsmBuild by a call to ocsmSetValu
        can be read outside ocsmBuild by a call to ocsmGetValu
    OUTPMTR
        declared by a OUTPMTR statement only at the top level
        can only be set at the top level
        can be read outside ocsmBuild by a call to ocsmGetValu
    LOCALVAR
        if a scalar, declared and defined by a SET, PATBEG, or
            GETATTR statement or via an INTERFACE IN statement
        if an array, declared by a DIMENSION statement
            values defined by one or more SET statements
        can have a string value
        values can be overwritten by subsequent statements
        are only available in the .csm or .udc file in which
            it is defined or in an include-type UDC that is
            called by the program unit in which iit was defined
    SOLVER
        is a scalar defined by a SOLBEG statement
        only available between SOLBEG and SOLEND

                                                            L
                                                D  C  C  O  O
                                                E  F  O  U  C
                                                S  G  N  T  A
                                                P  P  P  P  L
                                                M  M  M  M  V
                                                T  T  T  T  A
                                                R  R  R  R  R
                                                -  -  -  -  -
      Can be vector or array of numbers         Y  Y  Y  Y  Y
      Can have a string value                   N  N  N  Y  Y
      Can be restricted by LBOUND or UBOUND     Y  Y  N  N  N
      Scope (T=top-level, G=global, L=local)    T  T  G  L  L
      Defined during ocsmLoad or ocsmLoadDict   Y  Y  Y  N  N
      Can be set via ocsmSetValu(D)             Y  Y  N  N  N
      Defined and set during ocsmBuild          N  N  N  Y  Y
      Can be read via ocsmGetValu(S)            Y  Y  Y  Y  N
      Can find associated sensitivity           Y  N  N  N  N

    @-parameters depend on the last SELECT statement(s).
        each time a new Body is added to the Stack, 'SELECT body' is
            implicitly called
        depending on last SELECT statement, the values of the
             @-parameters are given by:

               body face edge node  <- last SELECT

        @seltype -1    2    1    0   selection type (0=node,1=edge,2=face)
        @selbody  x    -    -    -   current Body
        @sellist -1    x    x    x   list of Nodes/Edges/Faces

        @nbody    x    x    x    x   number of Bodys
        @ibody    x    x    x    x   current   Body
        @nface    x    x    x    x   number of Faces in @ibody
        @iface   -1    x   -1   -1   current   Face  in @ibody (or -2)
        @nedge    x    x    x    x   number of Edges in @ibody
        @iedge   -1   -1    x   -1   current   Edge  in @ibody (or -2)
        @nnode    x    x    x    x   number of Nodes in @ibody
        @inode   -1   -1   -1    x   current   Node  in @ibody (or -2)
        @igroup   x    x    x    x   group of @ibody
        @itype    x    x    x    x   0=NodeBody, 1=WireBody,
                                                 2=SheetBody, 3=SolidBody
        @nbors   -1    x    -    x   number of incident Edges
        @nbors   -1    -    x    -   number of incident Faces

        @ibody1  -1    x    x   -1   1st element of 'Body' Attr in @ibody
        @ibody2  -1    x    x   -1   2nd element of 'Body' Attr in @ibody

        @xmin     x    x    *    x   x-min of bboxes or x at beg of Edge
        @ymin     x    x    *    x   y-min of bboxes or y at beg of Edge
        @zmin     x    x    *    x   z-min of bboxes or z at beg of Edge
        @xmax     x    x    *    x   x-max of bboxes or x at end of Edge
        @ymax     x    x    *    x   y-max of bboxes or y at end of Edge
        @zmax     x    x    *    x   z-max of bboxes or z at end of Edge

        @length   0    0    x    0   length of Edges
        @area     x    x    0    0   area of Faces or surface area of body
        @volume   x    0    0    0   volume of body (if a solid)

        @xcg      x    x    x    x   location of center of gravity
        @ycg      x    x    x    x
        @zcg      x    x    x    x

        @Ixx      x    x    x    0   centroidal moment of inertia
        @Ixy      x    x    x    0
        @Ixz      x    x    x    0
        @Iyx      x    x    x    0
        @Iyy      x    x    x    0
        @Iyz      x    x    x    0
        @Izx      x    x    x    0
        @Izy      x    x    x    0
        @Izz      x    x    x    0

        @toler    x    x    x    0   maximum tolerance (at last SELECT)
        @signal   x    x    x    x   current signal code
        @nwarn    x    x    x    x   number of warnings (at last SELECT)

        @edata                       only set up by EVALUATE statement
        @stack                       Bodys on stack; 0=Mark; -1=none
        @scope                       scoping level (at last SELECT)
        @version                     version number

        in above table:
           x -> value is set
           - -> value is unchanged
           * -> special value is set (if single Edge)
           0 -> value is set to  0
          -1 -> value is set to -1

*/

/*
 ************************************************************************
 *                                                                      *
 * Expression rules                                                     *
 *                                                                      *
 ************************************************************************
 */

/*
Valid operators (in order of precedence):
    ( )            parentheses, inner-most evaluated first
    func(a,b)      function arguments, then function itself
    ^              exponentiation             (evaluated left to right)
    * /            multiply and divide        (evaluated left to right)
    + -            add/concat and subtract    (evaluated left to right)

    An expression that consists of only the name of a Parameter my be
    prepended by a unary + or -

Valid function calls:
    pi(x)                        3.14159...*x
    min(x,y)                     minimum of x and y
    max(x,y)                     maximum of x and y
    sqrt(x)                      square root of x
    abs(x)                       absolute value of x
    int(x)                       integer part of x  (3.5 -> 3, -3.5 -> -3)
                                     produces derivative=0
    nint(x)                      nearest integer to x
                                     produces derivative=0
    ceil(x)                      smallest integer not less than x
                                     produces derivative=0
    floor(x)                     largest integer not greater than x
                                     produces derivative=0
    mod(a,b)                     mod(a/b), with same sign as a and b>=0
    sign(test)                   returns -1, 0, or +1
                                     produces derivative=0
    exp(x)                       exponential of x
    log(x)                       natural logarithm of x
    log10(x)                     common logarithm of x
    sin(x)                       sine of x          (in radians)
    sind(x)                      sine of x          (in degrees)
    asin(x)                      arc-sine of x      (in radians)
    asind(x)                     arc-sine of x      (in degrees)
    cos(x)                       cosine of x        (in radians)
    cosd(x)                      cosine of x        (in degrees)
    acos(x)                      arc-cosine of x    (in radians)
    acosd(x)                     arc-cosine of x    (in degrees)
    tan(x)                       tangent of x       (in radians)
    tand(x)                      tangent of x       (in degrees)
    atan(x)                      arc-tangent of x   (in radians)
    atand(x)                     arc-tangent of x   (in degrees)
    atan2(y,x)                   arc-tangent of y/x (in radians)
    atan2d(y,x)                  arc-tangent of y/x (in degrees)
    hypot(x,y)                   hypotenuse: sqrt(x^2+y^2)
    hypot3(x,y,z)                hypotenuse: sqrt(x^2+y^2+z^2)
    incline(xa,ya,dab,xb,yb)     inclination of chord (in degrees)
                                     produces derivative=0
    Xcent(xa,ya,dab,xb,yb)       X-center of circular arc
                                     produces derivative=0
    Ycent(xa,ya,dab,xb,yb)       Y-center of circular arc
                                     produces derivative=0
    Xmidl(xa,ya,dab,xb,yb)       X-point at midpoint of circular arc
                                     produces derivative=0
    Ymidl(xa,ya,dab,xb,yb)       Y-point at midpoint of circular arc
                                     produces derivative=0
    seglen(xa,ya,dab,xb,yb)      length of segment
                                     produces derivative=0
    radius(xa,ya,dab,xb,yb)      radius of curvature (or 0 for LINSEG)
                                     produces derivative=0
    sweep(xa,ya,dab,xb,yb)       sweep angle of circular arc (in degrees)
                                     produces derivative=0
    turnang(xa,ya,dab,xb,yb,...
                     dbc,xc,yc)  turnnig angle at b (in degrees)
                                     produces derivative=0
    dip(xa,ya,xb,yb,rad)         acute dip between arc and chord
                                     produces derivative=0
    smallang(x)                  ensures -180<=x<=180
    val2str(num,digits)          convert num to string ("%d" if digits=0,
                                     "%f" if digits>0, "%e" if digits<0)
    str2val(string)              convert string to value
    findstr(str1,str2)           find locn of str2 in str1 (bias-1 or 0)
    slice(str,ibeg,iend)         substring of str from ibeg to iend
                                     (bias-1)
    path($pwd)                   returns present working directory
    path($csm)                   returns directory of current .csm,
                                     .cpc, or .udc file
    path($root)                  returns $ESP_ROOT
    path($file)                  returns name of .csm, .cpc, or .udc file
    ifzero(test,ifTrue,ifFalse)  if test=0, return ifTrue, else ifFalse
    ifpos(test,ifTrue,ifFalse)   if test>0, return ifTrue, else ifFalse
    ifneg(test,ifTrue,ifFalse)   if test<0, return ifTrue, else ifFalse
    ifmatch(str,pat,ifTrue,...
                      ifFalse)   if str match pat, return ifTrue,
                                     else ifFalse
                                        ? matches any one character
                                       '+ matches one  or more characters
                                        * matches zero or more characters
    ifnan(test,ifTrue,ifFalse)   if test is NaN, return ifTrue,
                                     else ifFalse
*/

/*
 ************************************************************************
 *                                                                      *
 * Attribute rules (accessible through EGADS)                           *
 *                                                                      *
 ************************************************************************
 */

/*
Attributes assigned to Bodys:

    _body       Body index (1:nbody)

    _brch       Branch index (1:nbrch)

    _tParams    tessellation parameters that were used

    _csys_*     arguments when CSYSTEM was defined

    _hist       history for WireBodys only

    <any>       all global Attributes

    <any>       all Attributes associated with Branch that created Body

    <any>       all Attributes associated with "SELECT $body" statement

                Note: if the Attribute name is ".tParams", then its
                      corresponding values are:
                       .tParams[1] = maximum triangle side length
                       .tParams[2] = maximum sag (distance between
                                                  chord and arc)
                       .tParams[3] = maximum angle between edge
                                                  segments (deg)

                Note: if the Attribute name is ".qParams" and it
                      value is any string, then the tessellation
                      templates are not used

                Note: if the Attribute name is ".qParams", then its
                      corresponding values are:
                      .qParams[1] = Edge matching expressed as the
                                    deviation from alignment
                      .qParams[2] = maximum quad side ratio point
                                    count to allow
                      .qParams[3] = number of smoothing iterations

Special User-defined Attributes for Bodys:

    _makeQuads  to make quads on all Faces in Body

    _name       string used in ESP interface for a Body

    _stlColor   color to use for all Faces in an .stl file

Attributes assigned to Faces:

    _body       non-unique 2-tuple associated with first Face creation
        [0]     Body index in which Face first existed (1:nbody)
        [1]     face-order associated with creation (see above)

    _brch       non-unique even-numbered list associated with Branches
                   that are active when the Face is created (most
                   recent Branch is listed first)
        [2*i  ] Branch index (1:nbrch)
        [2*i+1] (see below)

                Branches that contribute to brch Attribute are
                   primitive  (for which _brch[2*i+1] is face-order)
                   UDPRIM.udc (for which _brch[2*i+1] is 1)
                   grown      (for which _brch[2*i+1] is face-order)
                   applied    (for which _brch[2*i+1] is face-order)
                   sketch     (for which _brch[2*i+1] is Sketch primitive
                               if making WireBody)
                   PATBEG     (for which _brch[2*i+1] is pattern index)
                   IFTHEN     (for which _brch[2*i+1] is -1)
                   RECALL     (for which _brch[2*i+1] is +1)
                   RESTORE    (for which _brch[2*i+1] is Body numr stored)

    _faceID     unique 3-tuple that is assigned automatically
          [0]   _body[0]
          [1]   _body[1]
          [2]   sequence number

                if multiple Faces have same _faceID[0] and _faceID[1],
                   then the sequence number is defined based upon the
                   first rule that applies:
                   * Face with smaller xcg  has lower sequence number
                   * Face with smaller ycg  has lower sequence number
                   * Face with smaller zcg  has lower sequence number
                   * Face with smaller area has lower sequence number

    _hist       list of Bodys that contained this Face (oldest to newest)

    <any>       all Attributes associated with Branch that first
                    created Face
                    (BOX, CONE, CYLINDER, IMPORT, SPHERE, TORUS, UDPRIM)
                    (BLEND, EXTRUDE, LOFT, REVOLVE, RULE, SWEEP)
                    (SKEND)
                    (CHAMFER, CONNECT, FILLET, HOLLOW)

    <any>       all Attributes associated with Branch if a RESTORE
                    statement

    <any>       all Attributes associated with "SELECT FACE" statement

Special User-defined Attributes for Faces:

    _color      color of front of Face in ESP
                either R,G,B in three 0-1 reals
                or $red, $green, $blue, $yellow, $magenta,
                $cyan, $white, or $black

    _bcolor     color of back of Face in ESP (see _color)

    _gcolor     color of grid of Face in ESP (see _color)

    _makeQuads  to make quads for this Face

    _stlColor   color to use for this Face in an .stl file

Attributes assigned to Edges:

    _body       non-unique 2-tuple associated with first Edge creation
        [0]     Body index in which Edge first existed (1:nbody)
        [1]     10000 * min(_body[1][ileft],_body[1][irite])
                      + max(_body[1][ileft],_body[1][irite])
                (or -3 if non-manifold)

    _edgeID     unique 5-tuple that is assigned automatically
          [0]   _faceID[0] of Face 1 (or 0 if non-manifold)
          [1]   _faceID[1] of Face 1 (or 0 if non-manifold)
          [2]   _faceID[0] of Face 2 (or 0 if non-manifold)
          [3]   _faceID[1] of Face 2 (or 0 if non-manifold)
          [4]   sequence number

                _edgeID[0]/[1] swapped with edge[2]/[3]
                   10000*_edgeID[0]+_edgeID[1] > 10000*_edgeID[2]+_edgeID[3]
                if multiple Edges have same _edgeID[0], _edgeID[1],
                   _edgeID[2], and _edgeID[3], then the sequence number
                   is defined based upon the first rule that applies:
                   * Edge with smaller xcg    has lower sequence number
                   * Edge with smaller ycg    has lower sequence number
                   * Edge with smaller zcg    has lower sequence number
                   * Edge with smaller length has lower sequence number

    _nface      number of incident Faces

    <any>       all Attributes associated with "SELECT EDGE" statement

Special User-defined Attributes for Edges:

    _color      color of front of Edge in ESP
                either R,G,B in three 0-1 reals
                or $red, $green, $blue, $yellow, $magenta,
                $cyan, $white, or $black

    _gcolor     color of grid of Edge in ESP (see _color)

Attributes assigned to Nodes:

    _nodeID     unique integer that is assigned automatically

    _nedge      number of incident Edges

    <any>       all Attributes associated with "SELECT FACE" statement

Special User-defined Attributes for Nodes:

    _color      color of Node in ESP
                either R,G,B in three 0-1 reals
                or $red, $green, $blue, $yellow, $magenta,
                $cyan, $white, or $black
*/

/*
 ************************************************************************
 *                                                                      *
 * Structures                                                           *
 *                                                                      *
 ************************************************************************
 */

/* "Attr" is a Branch Attribute */
typedef struct {
    char          *name;                /* Attribute name */
    char          *defn;                /* Attribute definition */
    int           type;                 /* ATTRREAL or ATTRCSYS */
} attr_T;

/* "GRatt" is a graphic Attribute */
typedef struct {
    void          *object;              /* pointer to GvGraphic (or NULL) */
    int           active;               /* =1 if entity should be rendered */
    int           color;                /* entity color in form 0x00rrggbb */
    int           bcolor;               /* back   color in form 0x00rrggbb */
    int           mcolor;               /* mesh   color in form 0x00rrggbb */
    int           lwidth;               /* line width in pixels */
    int           ptsize;               /* point size in pixels */
    int           render;               /* render flags: */
                                        /*    2 GV_FOREGROUND */
                                        /*    4 GV_ORIENTATION */
                                        /*    8 GV_TRANSPARENT */
                                        /*   16 GV_FACETLIGHT */
                                        /*   32 GV_MESH */
                                        /*   64 GV_FORWARD */
    int           dirty;                /* =1 if Attributes have been changed */
} grat_T;

/* "Varg" is the (multi-) value of an argument associated with a Body */
typedef struct {
    int           nval;                 /* number of values (or 0 if string) */
    int           nrow;                 /* number of rows */
    int           ncol;                 /* number of columns */
    union {
        char      *str;                 /* character array if nval==0 */
        double    *val;                 /* array of values if nval>0  */
    };
    double        *dot;                 /* array of velocities if nval>0 */
} varg_T;

/* "Mprp" are the mass properties used by a Branch */
typedef struct {
    char          name[10];             /* name     of the mass property */
    double        val;                  /* value    of the mass property */
} mprp_T;

/* "Node" is a 0-D topological entity in a Body */
typedef struct {
    int           nedge;                /* number of indicent Edges */
    double        x;                    /* x-coordinate */
    double        y;                    /* y-coordinate */
    double        z;                    /* z-coordinate */
    int           ibody;                /* Body index (1:nbody) */
    grat_T        gratt;                /* GRatt of the Node */
    double        *dxyz;                /* tessellation velocity (or NULL) */
    ego           enode;                /* EGADS node object */
} node_T;

/* "Edge" is a 1-D topological entity in a Body */
typedef struct {
    int           itype;                /* Edge type */
    int           ibeg;                 /* Node at beginning */
    int           iend;                 /* Node at end */
    int           ileft;                /* Face on the left */
    int           irite;                /* Face on the rite */
    int           nface;                /* number of incident Faces */
    int           ibody;                /* Body index (1:nbody) */
    int           iford;                /* face-order */
    int           imark;                /* value of "mark" Attribute (or -1) */
    grat_T        gratt;                /* GRatt of the Edge */
    double        *dxyz;                /* tessellation velocity (or NULL) */
    double        *dt;                  /* parametric   velocity (or NULL) */
    int           globid;               /* global ID (bias-1) */
    ego           eedge;                /* EGADS edge object */
} edge_T;

/* "Face" is a 2-D topological entity in a Body */
typedef struct {
    int           ibody;                /* Body index (1:nbody) */
    int           iford;                /* face-order */
    int           imark;                /* value of "mark" Attribute (or -1) */
    grat_T        gratt;                /* GRatt of the Face */
    void          *eggdata;             /* pointer to external grid generator data */
    double        *dxyz;                /* tessellation velocity (or NULL) */
    double        *duv;                 /* parametric   velocity (or NULL) */
    int           globid;               /* global ID (bias-1) */
    ego           eface;                /* EGADS face object */
} face_T;

/* "Body" is a boundary representation */
typedef struct {
    int           ibrch;                /* Branch associated with Body */
    int           brtype;               /* Branch type (see below) */
    int           ileft;                /* left parent Body (or 0) */
    int           irite;                /* rite parent Body (or 0) */
    int           ichld;                /* child Body (or 0 for root) */
    int           igroup;               /* Group number */
    varg_T        arg[10];              /* array of evaluated arguments (actually use 1-9) */

    ego           ebody;                /* EGADS Body         object(s) */
    ego           eebody;               /* EGADS EBody        object(s) */
    ego           etess;                /* EGADS Tessellation object(s) */
    ego           eetess;               /* EGADS Tessellation object(s) associated with eebody */
    int           npnts;                /* total number of unique points */
    int           ntris;                /* total number of triangles */

    int           onstack;              /* =1 if on stack (and returned); =0 otherwise */
    int           hasdots;              /* =1 if an argument has a dot; =2 if UDPARG is changed; =0 otherwise */
    int           hasdxyz;              /* =1 if Body has associated velocities */
    int           botype;               /* Body type (see below) */
    double        CPU;                  /* CPU time (sec) */
    int           nnode;                /* number of Nodes */
    node_T        *node;                /* array  of Nodes */
    int           nedge;                /* number of Edges */
    edge_T        *edge;                /* array  of Edges */
    int           nface;                /* number of Faces */
    face_T        *face;                /* array  of Faces */
    int           hassens;              /* flag for caching sensitivity info */
    grat_T        gratt;                /* GRatt of the Nodes */
} body_T;

/* "Brch" is a Branch in a feature tree */
typedef struct {
    char          *name;                /* name  of Branch */
    int           type;                 /* type  of Branch (see OpenCSM.h) */
    int           bclass;               /* class of Branch (see OpenCSM.h) */
    int           level;                /* =0 if from .csm, >0 if from .udc */
    int           indent;               /* indentation */
    char          *filename;            /* filename where Branch is defined */
    int           linenum;              /* line number in file where Branch is defined */
    int           actv;                 /* activity of Branch (see OpenCSM.h) */
    int           dirty;                /* =1 if dirty */
    int           nattr;                /* number of Attributes and Csystem */
    attr_T        *attr;                /* array  of Attributes and Csystem */
    int           ileft;                /* left parent Branch (or 0)*/
    int           irite;                /* rite parent Branch (or 0)*/
    int           ichld;                /* child Branch (or 0 for root) */
    int           nmprp;                /* number of mass properties */
    mprp_T        *mprp;                /* array  of mass properties */
    int           narg;                 /* number of arguments */
    char          *arg1;                /* definition for args[1] */
    char          *arg2;                /* definition for args[2] */
    char          *arg3;                /* definition for args[3] */
    char          *arg4;                /* definition for args[4] */
    char          *arg5;                /* definition for args[5] */
    char          *arg6;                /* definition for args[6] */
    char          *arg7;                /* definition for args[7] */
    char          *arg8;                /* definition for args[8] */
    char          *arg9;                /* definition for args[9] */
} brch_T;

/* "Pmtr" is a CONPMTR, driving (DESPMTR/CFGPMTR), or driven (LOCALVAR/OUTPMTR) Parameter */
typedef struct {
    char          *name;                /* name of Parameter */
    int           type;                 /* Parameter type (see below) */
    int           scope;                /* associated scope (nominally 0) */
    int           mprop;                /* =1 if associated with mass properties */
    int           nrow;                 /* number of rows    (=0 for string) */
    int           ncol;                 /* number of columns (=0 for string) */
    double        *value;               /* current value(s) */
    double        *dot;                 /* current velocity(s) */
    double        *lbnd;                /* lower Bound(s) */
    double        *ubnd;                /* upper Bound(s) */
    char          *str;                 /* string value */
} pmtr_T;

/* "Stor" is the storage locations used by STORE/RESTORE */
typedef struct {
    char          name[MAX_NAME_LEN+1]; /* name  of Storage */
    int           index;                /* index of Storare */
    int           nbody;                /* number of Bodys stored */
    int           *ibody;               /* array  of Body numbers stored */
    ego           *ebody;               /* array  of EGADS Bodys stored */
} stor_T;

/* "Prof" is profile data */
typedef struct {
    int           ncall;                /* number of calls */
    clock_t       time;                 /* total time */
} prof_T;

/* handle to callback functions */
typedef void (*mesgCB_H)   (char message[]);
typedef void (*sizeCB_H)   (void *modl, int ipmtr, int nrow, int ncol);

/* handles to functions in the external grid generator */
typedef int (*eggGenerate_H) (double[], int[], void**);
typedef int (*eggMorph_H)    (void*, double*, void**);
typedef int (*eggInfo_H)     (void*, int*, int*, const double**, const int**, int*, const int**);
typedef int (*eggDump_H)     (void*, FILE*);
typedef int (*eggLoad_H)     (FILE*, void**);
typedef int (*eggFree_H)     (void*);

/* "Modl" is a constructive solid model consisting of a tree of Branches
         and (possibly) a set of Parameters as well as the associated Bodys */
typedef struct modl_T {
    int           magic;                /* magic number to check for valid *modl */
    int           checked;              /* =1 if successfully passed checks */
    int           ibrch;                /* Branch number being executed */
    int           nextseq;              /* number of next automatcally-numbered item */
    int           ngroup;               /* number of Groups */
    int           recycle;              /* last Body to recycle */
    int           verify;               /* =1 if verification ASSERTs are checked */
    int           cleanup;              /* number of Branches before context is cleaned up */
    int           dumpEgads;            /* =1 if Bodys are dumped during build */
    int           loadEgads;            /* =1 if Bodys are loaded during build */
    int           hasMPs;               /* =1 if mass properties have been calculated */
    int           printStack;           /* =1 to print stack after every command */
    int           tessAtEnd;            /* =1 to tessellate Bodys on stack at end of ocsmBuild */
    int           erepAtEnd;            /* =1 to generate Erep based upon _erepAttr and _erepAngle */
    int           bodyLoaded;           /* Body index of last Body loaded */

    int           seltype;              /* selection type: 0=Node, 1=Edge, 2=Face, or -1 */
    int           selbody;              /* Body selected (or -1)  (1:nbody) */
    int           selsize;              /* number of selected entities */
    int           *sellist;             /* array  of selected entities */

    int           level;                /* level of file (=0 for .csm, >0 for .udc) */
    int           scope[11];            /* variable scope at this level */
    char          filename[MAX_FILENAME_LEN]; /* name of .csm file */
    int           tmpDirNum;            /* the xx of tmp_OpenCSM_xx.yy (or -1 if none) */

    int           ninline;              /* number of bytes used in inline file stream */
    int           minline;              /* size of inline file stream */
    char          *sinline;             /* inline file stream */

    int           nattr;                /* number of global Attributes */
    attr_T        *attr;                /* array  of global Attributes */

    int           nstor;                /* number of storages */
    stor_T        *stor;                /* array  of storages */

    int           nbrch;                /* number of Branches */
    int           mbrch;                /* maximum   Branches */
    brch_T        *brch;                /* array  of Branches */

    int           npmtr;                /* number of Parameters */
    int           mpmtr;                /* maximum   Parameters */
    pmtr_T        *pmtr;                /* array  of Parameters */

    int           nbody;                /* number of Bodys */
    int           mbody;                /* maximum   Bodys */
    body_T        *body;                /* array  of Bodys */

    int           needFDs;              /* =1 if finite differences are needed */
    int           numdots;              /* number of non-zero dots associated with DESPMTRS */
    struct modl_T *perturb;             /* model of perturbed body for sensitivty */
    struct modl_T *basemodl;            /* base MODL while creating perturbation */
    double        dtime;                /* time step in sensitivity */
                                        /*  0.001 = initial value */
                                        /* -2     = problem with previous attempt to create perturb */

    ego           context;              /* EGADS context */
    void          *userdata;            /* pointer for user data */
    mesgCB_H      mesgCB;               /* handle of message function (or NULL) */
    sizeCB_H      sizeCB;               /* handle of size change function (or NULL) */
    char          eggname[256];         /* name of external grid generator (or NULL) */
    eggGenerate_H eggGenerate;          /* handle of egg_generate function */
    eggMorph_H    eggMorph;             /* handle of egg_morph function */
    eggInfo_H     eggInfo;              /* handle of egg_info function */
    eggDump_H     eggDump;              /* handle of egg_dump function */
    eggLoad_H     eggLoad;              /* handle of egg_load function */
    eggFree_H     eggFree;              /* handle of egg_free function */

    int           nwarn;                /* number of warnings */
    int           sigCode;              /* current signal code */
    char          *sigMesg;             /* current signal message */

    prof_T        profile[100];         /* profile data */
} modl_T;

/*
 ************************************************************************
 *                                                                      *
 * Callable routines                                                    *
 *                                                                      *
 ************************************************************************
 */

/* return current version */
int ocsmVersion(int   *imajor,          /* (out) major version number */
                int   *iminor);         /* (out) minor version number */

/* set output level */
int ocsmSetOutLevel(int    ilevel);     /* (in)  output level: */
                                        /*       <0 do not change */
                                        /*       =0 errors only */
                                        /*       =1 nominal (default) */
                                        /*       =2 debug */
                                        /* (out) previous outLevel */

/* create a MODL by reading a .csm file */
int ocsmLoad(char   filename[],         /* (in)  file to be read (with .csm) */
             void   **modl);            /* (out) pointer to MODL */

/* create a MODL from an ego */
int ocsmLoadFromModel(ego    emodel,    /* (in)  egads MODEL */
                      void   **modl);   /* (out) pointer to MODL */

/* load dictionary from dictname */
int ocsmLoadDict(void   *modl,          /* (in)  pointer to MODL */
                 char   dictname[]);    /* (in)  file that contains dictionary */

/* update CFGPMTRs and DESPMTRs from filename */
int ocsmUpdateDespmtrs(void   *modl,    /* (in)  pointer to MODL */
                       char   filename[]); /* (in) file that contains DESPMTRs */

/* get a list of all .csm, .cpc. and .udc files */
int ocsmGetFilelist(void   *modl,       /* (in)  pointer to MODL */
                    char   *filelist[]);/* (out) bar-sepatared list of files */
                                        /*       must be freed by user */

/* save a MODL to a file */
int ocsmSave(void   *modl,              /* (in)  pointer to MODL */
             char   filename[]);        /* (in)  file to be written (with extension) */
                                        /*       .csm -> write outer .csm file */
                                        /*       .cpc -> write checkpointed .csm file */
                                        /*       .udc -> write a .udc file */

/* save CFGPMTRs and DESPMTRs to file */
int ocsmSaveDespmtrs(void   *modl,      /* (in)  pointer to MODL */
                     char   filename[]);/* (in)  file to which DESPMTRS should be written */

/* copy a MODL (but not the Bodys) */
int ocsmCopy(void   *srcModl,           /* (in)  pointer to source MODL */
             void   **newModl);         /* (out) pointer to new    MODL */

/* free up all storage associated with a MODL */
int ocsmFree(
   /*@null@*/void   *modl);             /* (in)  pointer to MODL (or NULL) */

/* get info about a MODL */
int ocsmInfo(void   *modl,              /* (in)  pointer to MODL */
             int    *nbrch,             /* (out) number of Branches */
             int    *npmtr,             /* (out) number of Parameters */
             int    *nbody);            /* (out) number of Bodys */


/* check that Branches are properly ordered */
int ocsmCheck(void   *modl);            /* (in)  pointer to MODL */

/* register a callback function for exporting messages */
int ocsmRegMesgCB(void    *modl,        /* (in)  pointer to MODL */
                  void    (*callback)(char[]));   /* (in)  handle of callback function */

/* register a callback function for DESPMTR size changes */
int ocsmRegSizeCB(void    *modl,        /* (in)  pointer to MODL */
                  void    (*callback)(void*, int, int, int));   /* (in)  handle of callback function */

/* build Bodys by executing the MODL up to a given Branch */
int ocsmBuild(void   *modl,             /* (in)  pointer to MODL */
              int    buildTo,           /* (in)  last Branch to execute (or 0 for all, or -1 for no recycling) */
              int    *builtTo,          /* (out) last Branch executed successfully */
              int    *nbody,            /* (in)  number of entries allocated in body[] */
                                        /* (out) number of Bodys on the stack */
    /*@null@*/int    body[]);           /* (out) array  of Bodys on the stack (LIFO)
                                                 (at least nbody long) */

/* get information about one Body */
int ocsmBodyDetails(void   *modl,       /* (in)  pointer to MODL */
                    char   fiename[],   /* (in)  name of file from which Bdy was created */
                    int    linenum,     /* (in)  line in filename from which Body was created */
                    char   *info[]);    /* (out) info about the Body (freeable) */

/* create a perturbed MODL */
int ocsmPerturb(void   *modl,           /* (in)  pointer to MODL */
                int    npmtrs,          /* (in)  number of perturbed Parameters (or 0 to remove) */
      /*@null@*/int    ipmtrs[],        /* (in)  array of Parameter indices (1:npmtr) */
      /*@null@*/int    irows[],         /* (in)  array of row       indices (1:nrow) */
      /*@null@*/int    icols[],         /* (in)  array of column    indices (1:ncol) */
      /*@null@*/double values[]);       /* (in)  array of perturbed values */

/* update a tessellation from a file */
int ocsmUpdateTess(void   *modl,        /* (in)  pointer to MODL */
                   int    ibody,        /* (in)  Body index (1:nbody) */
                   char   filename[]);  /* (in)  name of .tess file */

/* create a new Branch */
int ocsmNewBrch(void   *modl,           /* (in)  pointer to MODL */
                int    iafter,          /* (in)  Branch index (0-nbrch) after which to add */
                int    type,            /* (in)  Branch type (see below) */
                char   filename[],      /* (in)  filename where Branch is defined */
                int    linenum,         /* (in)  linenumber where Branch is defined (bias-1) */
      /*@null@*/char   arg1[],          /* (in)  Argument 1 (or NULL) */
      /*@null@*/char   arg2[],          /* (in)  Argument 2 (or NULL) */
      /*@null@*/char   arg3[],          /* (in)  Argument 3 (or NULL) */
      /*@null@*/char   arg4[],          /* (in)  Argument 4 (or NULL) */
      /*@null@*/char   arg5[],          /* (in)  Argument 5 (or NULL) */
      /*@null@*/char   arg6[],          /* (in)  Argument 6 (or NULL) */
      /*@null@*/char   arg7[],          /* (in)  Argument 7 (or NULL) */
      /*@null@*/char   arg8[],          /* (in)  Argument 8 (or NULL) */
      /*@null@*/char   arg9[]);         /* (in)  Argument 9 (or NULL) */

/* get info about a Branch */
int ocsmGetBrch(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                int    *type,           /* (out) Branch type (see below) */
                int    *bclass,         /* (out) Branch class (see below) */
                int    *actv,           /* (out) Branch Activity (see below) */
                int    *ichld,          /* (out) ibrch of child (or 0 if root) */
                int    *ileft,          /* (out) ibrch of left parent (or 0) */
                int    *irite,          /* (out) ibrch of rite parent (or 0) */
                int    *narg,           /* (out) number of Arguments */
                int    *nattr);         /* (out) number of Attributes */

/* set activity for a Branch */
int ocsmSetBrch(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                int    actv);           /* (in)  Branch activity (see below) */

/* delete a Branch (or whole Sketch if SKBEG) */
int ocsmDelBrch(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch);          /* (in)  Branch index (1:nbrch) */

/* print Branches to file */
int ocsmPrintBrchs(void   *modl,        /* (in)  pointer to MODL */
                   char   filename[]);  /* (in)  file to which output is appended (or "" for stdout) */

/* get an Argument for a Branch */
int ocsmGetArg(void   *modl,            /* (in)  pointer to MODL */
               int    ibrch,            /* (in)  Branch index (1:nbrch) */
               int    iarg,             /* (in)  Argument index (1:narg) */
               char   defn[],           /* (out) Argument definition (at least MAX_STRVAL_LEN long) */
               double *value,           /* (out) Argument value */
               double *dot);            /* (out) Argument velocity */

/* set an Argument for a Branch */
int ocsmSetArg(void   *modl,            /* (in)  pointer to MODL */
               int    ibrch,            /* (in)  Branch index (1:nbrch) */
               int    iarg,             /* (in)  Argument index (1:narg) */
               char   defn[]);          /* (in)  Argument definition */

/* return an Attribute for a Branch by index */
int ocsmRetAttr(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                int    iattr,           /* (in)  Attribute index (1:nattr) */
                char   aname[],         /* (out) Attribute name  (at least MAX_STRVAL_LEN long) */
                char   avalue[]);       /* (out) Attribute value (at least MAX_STRVAL_LEN long) */

/* get an Attribute for a Branch by name */
int ocsmGetAttr(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) or 0 for global */
                char   aname[],         /* (in)  Attribute name */
                char   avalue[]);       /* (out) Attribute value (at least MAX_STRVAL_LEN long) */

/* set an Attribute for a Branch */
int ocsmSetAttr(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) or 0 for global */
                char   aname[],         /* (in)  Attribute name */
                char   avalue[]);       /* (in)  Attribute value (or blank to delete) */

/* return a Csystem for a Branch by index */
int ocsmRetCsys(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                int    icsys,           /* (in)  Csystem index (1:nattr) */
                char   cname[],         /* (out) Csystem name  (at least MAX_STRVAL_LEN long) */
                char   cvalue[]);       /* (out) Csystem value (at least MAX_STRVAL_LEN long) */

/* get a Csystem for a Branch by name */
int ocsmGetCsys(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                char   cname[],         /* (in)  Csystem name */
                char   cvalue[]);       /* (out) Csystem value (at least MAX_STRVAL_LEN long) */

/* set a Csystem for a Branch */
int ocsmSetCsys(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch)  */
                char   cname[],         /* (in)  Csystem name */
                char   cvalue[]);       /* (in)  Csystem value (or blank to delete) */

/* print global Attributes to file */
int ocsmPrintAttrs(void   *modl,        /* (in)  pointer to MODL */
                   char   filename[]);  /* (in)  file to which output is appended (or "" for stdout) */

/* get the name of a Branch */
int ocsmGetName(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                char   name[]);         /* (out) Branch name (at least MAX_STRVAL_LEN long) */

/* set the name for a Branch */
int ocsmSetName(void   *modl,           /* (in)  pointer to MODL */
                int    ibrch,           /* (in)  Branch index (1:nbrch) */
                char   name[]);         /* (in)  Branch name */

/* get string data associated with a Sketch */
int ocsmGetSketch(void   *modl,         /* (in)  pointer to MODL */
                  int    ibrch,         /* (in)  Branch index (1:nbrch) within Sketch */
                  int    maxlen,        /* (in)  length of begs, vars, cons, and segs */
                  char   begs[],        /* (out) string with SKBEG info */
                                        /*       "xarg;xval;yarg;yval;zarg;zval;" */
                  char   vars[],        /* (out) string with Sketch variables */
                                        /*       "x1;y1;d1;x2; ... dn;" */
                  char   cons[],        /* (out) string with Sketch constraints */
                                        /*       "type1;index1_1;index2_1;value1; ... valuen;" */
                                        /*       index1 and index2 are bias-1 */
                  char   segs[]);       /* (out) string with Sketch segments */
                                        /*       "type1;ibeg1;iend1; ... iendn;" */
                                        /*       ibeg and iend are bias-1 */

/* solve for new Sketch variables */
int ocsmSolveSketch(void   *modl,       /* (in)  pointer to MODL */
                    char   vars_in[],   /* (in)  string with Sketch variables */
                    char   cons[],      /* (in)  string with Sketch constraints */
                    char   vars_out[]); /* (out) string (1024 long) with new Sketch variables */

/* overwrite Branches associated with a Sketch */
int ocsmSaveSketch(void   *modl,        /* (in)  pointer to MODL */
                   int    ibrch,        /* (in)  Branch index (1:nbrch) within Sketch */
                   char   vars[],       /* (in)  string with Sketch variables */
                   char   cons[],       /* (in)  string with Sketch constraints */
                   char   segs[]);      /* (in)  string with Sketch segments */

/* make an EBody from a given Body */
int ocsmMakeEBody(void   *modl,         /* (in)  pointer to MODL */
                  int    ibody,         /* (in)  Body index (1:nbody) */
                  double dihedral,      /* (in)  dihedral angle (deg) for implicit .Keeps */
        /*@null@*/char   entList[]);    /* (in)  entities to use in make EBody (or NULL to delete) */

/* create a new Parameter */
int ocsmNewPmtr(void   *modl,           /* (in)  pointer to MODL */
                char   name[],          /* (in)  Parameter name */
                int    type,            /* (in)  Parameter type */
                int    nrow,            /* (in)  number of rows */
                int    ncol);           /* (in)  number of columns */

/* delete a Parameter */
int ocsmDelPmtr(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr);          /* (in)  Parameter index (1:npmtr) */

/* find (or create) a Parameter */
int ocsmFindPmtr(void   *modl,          /* (in)  pointer to MODL */
                 char   name[],         /* (in)  Parameter name */
                 int    type,           /* (in)  Parameter type    (or 0 to not create) */
                 int    nrow,           /* (in)  number of rows */
                 int    ncol,           /* (in)  number of columns */
                 int    *ipmtr);        /* (out) Parameter index (1:npmtr) (or 0 if not found) */

/* get info about a Parameter */
int ocsmGetPmtr(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr,           /* (in)  Parameter index (1:npmtr) */
                int    *type,           /* (out) Parameter type */
                int    *nrow,           /* (out) number of rows */
                int    *ncol,           /* (out) number of columns */
                char   name[]);         /* (out) Parameter name (at least MAX_NAME_LEN long) */

/* print DESPMTRs and OUTPMTRs to file */
int ocsmPrintPmtrs(void   *modl,        /* (in)  pointer to MODL */
                   char   filename[]);  /* (in)  file to which output is appended (or "" for stdout) */

/* get the Value of a Parameter */
int ocsmGetValu(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr,           /* (in)  Parameter index (1:npmtr) */
                int    irow,            /* (in)  row    index (1:nrow) */
                int    icol,            /* (in)  column index (1:ncol) */
                double *value,          /* (out) Parameter value */
                double *dot);           /* (out) Parameter velocity */

/* get the Value of a string Parameter */
int ocsmGetValuS(void   *modl,          /* (in)  pointer to MODL */
                 int    ipmtr,          /* (in)  Parameter index (1:npmtr) */
                 char   str[]);         /* (out) Parameter value (at least MAX_STRVAL_LEN long) */

/* set a Value for a Parameter */
int ocsmSetValu(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr,           /* (in)  Parameter index (1:npmtr) */
                int    irow,            /* (in)  row    index (1:nrow) */
                int    icol,            /* (in)  column index (1:ncol) or 0 for index*/
                char   defn[]);         /* (in)  definition of Value */

/* set a (double) Value for a Parameter */
int ocsmSetValuD(void   *modl,          /* (in)  pointer to MODL */
                 int    ipmtr,          /* (in)  Parameter index (1:npmtr) */
                 int    irow,           /* (in)  row    index (1:nrow) */
                 int    icol,           /* (in)  column index (1:ncol) or 0 for index*/
                 double value);         /* (in)  value to set */

/* get the Bounds of a Parameter */
int ocsmGetBnds(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr,           /* (in)  Parameter index (1:npmtr) */
                int    irow,            /* (in)  row    index (1:nrow) */
                int    icol,            /* (in)  column index (1:ncol) */
                double *lbound,         /* (out) lower Bound */
                double *ubound);        /* (out) upper Bound */

/* set the Bounds of a Parameter */
int ocsmSetBnds(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr,           /* (in)  Parameter index (1:npmtr) */
                int    irow,            /* (in)  row    index (1:nrow) */
                int    icol,            /* (in)  column index (1:ncol) */
                double lbound,          /* (in)  lower Bound to set */
                double ubound);         /* (in)  upper Bound to set */

/* set sensitivity FD time step (or select analytic) */
int ocsmSetDtime(void   *modl,          /* (in)  pointer to MODL */
                 double dtime);         /* (in)  time step (or 0 to choose analytic) */

/* set the velocity for a Parameter */
int ocsmSetVel(void   *modl,            /* (in)  pointer to MODL */
               int    ipmtr,            /* (in)  Parameter index (1:npmtr) or 0 for all */
               int    irow,             /* (in)  row    index (1:nrow)     or 0 for all */
               int    icol,             /* (in)  column index (1:ncol)     or 0 for index */
               char   defn[]);          /* (in)  definition of Velocity */

/* set the (double) velocity for a Parameter */
int ocsmSetVelD(void   *modl,           /* (in)  pointer to MODL */
                int    ipmtr,           /* (in)  Parameter index (1:npmtr) or 0 for all */
                int    irow,            /* (in)  row    index (1:nrow)     or 0 for all */
                int    icol,            /* (in)  column index (1:ncol)     or 0 for index */
                double dot);            /* (in)  velocity to set */

/* get the EGO associated with a Body or its parts */
int ocsmGetEgo(void   *modl,            /* (in)  pointer to MODL */
               int    ibody,            /* (in)  Body index (1:nbody) */
               int    seltype,          /* (in)  OCSM_BODY, OCSM_NODE, OCSM_EDGE, or OCSM_FACE */
               int    iselect,          /* (in)  if OCSM_BODY, 0 for Body, 1 for Tessellation, 2 for context
                                                 3 for EBody, 4 for Tessellation on EBody */
                                        /*       if OCSM_NODE, Node index (1:nnode) */
                                        /*       if OCSM_EDGE, Edge index (1:nedge) */
                                        /*       if OCSM_FACE, Face index (1:nface) */
               ego    *theEgo);         /* (out) context or associated EGO */

/* set a tessellation or EBody for a Body */
int ocsmSetEgo(void   *modl,            /* (in)  pointer to MODL */
               int    ibody,            /* (in)  Body index (1:nbody) */
               int    iselect,          /* (in)  1 for Tessellation,
                                                 3 for EBody, 4 for Tessellation on EBody */
               ego    theEgo);          /* (in)  associated EGO */

/* find the entity nnumber given FaceID or EdgeID */
int ocsmFindEnt(void   *modl,           /* (in)  pointer to MODL */
                int    ibody,           /* (in)  Body index (1:nbody) */
                int    seltype,         /* (in)  OCSM_EDGE, or OCSM_FACE */
                int    entID[],         /* (in)  entity ID (either faceID or edgeID) */
                int    *ient);          /* (out) entity number */

/* get the parametric coordinates on an Edge or Face */
int ocsmGetUV(void   *modl,             /* (in)  pointer to MODL */
              int    ibody,             /* (in)  Body index (1:nbody) */
              int    seltype,           /* (in)  OCSM_EDGE, or OCSM_FACE */
              int    iselect,           /* (in)  Edge or Face index (1:nent) */
              int    npnt,              /* (in)  number of points */
    /*@null@*/double xyz[],             /* (in)  coordinates (NULL or 3*npnt in length) */
              double uv[]);             /* (out) para coords (1*npnt or 2*npnt in length) */

/* get the coordinates on a Node, Edge, or Face */
int ocsmGetXYZ(void   *modl,            /* (in)  pointer to MODL */
               int    ibody,            /* (in)  Body index (1:nbody) */
               int    seltype,          /* (in)  OCSM_NODE, OCSM_EDGE, or OCSM_FACE */
               int    iselect,          /* (in)  Node, Edge, or Face index (1:nent) */
               int    npnt,             /* (in)  number of points */
     /*@null@*/double uv[],             /* (in)  para coords (NULL, 1*npnt, or 2*npnt) */
               double xyz[]);           /* (out) coordinates (3*npnt in length) */

/* get the unit normals for a Face */
int ocsmGetNorm(void   *modl,           /* (in)  pointer to MODL */
                int    ibody,           /* (in)  Body index (1:nbody) */
                int    iface,           /* (in)  Face index (1:nface) */
                int    npnt,            /* (in)  number of points */
      /*@null@*/double uv[],            /* (in)  para coords (NULL or 2*npnt in length) */
                double norm[]);         /* (out) normals (3*npnt in length) */

/* get the velocities of coordinates on a Node, Edge, or Face */
int ocsmGetVel(void   *modl,            /* (in)  pointer to MODL */
               int    ibody,            /* (in)  Body index (1:nbody) */
               int    seltype,          /* (in)  OCSM_NODE, OCSM_EDGE, or OCSM_FACE */
               int    iselect,          /* (in)  Node, Edge, or Face index (1:nent) */
               int    npnt,             /* (in)  number of points */
     /*@null@*/double uv[],             /* (in)  para coords
                                                    NULL           for OCSM_NODE
                                                    NULL or 1*npnt for OCSM_EDGE
                                                    NULL or 2*npnt for OCSM_FACE */
               double vel[]);           /* (out) velocities (in pre-allocated array)
                                                    3      for OCSM_NODE
                                                    3*npnt for OCEM_EDGE
                                                    3*npnt for OCSM_FACE */

/* set up alternative tessellation by an external grid generator */
int ocsmSetEgg(void   *modl,            /* (in)  pointer to MODL */
               char   *eggname);        /* (in)  name of dynamically-loadable file */

/* get the number of tess points on an Edge or Face */
int ocsmGetTessNpnt(void   *modl,       /* (in)  pointer to MODL */
                    int    ibody,       /* (in)  Body index (1:nbody)) */
                    int    seltype,     /* (in)  OCSM_EDGE , or OCSM_FACE */
                    int    iselect,     /* (in)  Edge  or Face index (1:nent) */
                    int    *npnt);      /* (out) number of tessellation points */

/* get the tessellation velocities on a Node, Edge, or Face */
int ocsmGetTessVel(void   *modl,        /* (in)  pointer to MODL */
                   int    ibody,        /* (in)  Body index (1:nbody)) */
                   int    seltype,      /* (in)  OCSM_NODE, OCSM_EDGE, or OCSM_FACE */
                   int    iselect,      /* (in)  Node, Edge, or Face index (1:nent) */
             const double *dxyz[]);     /* (out) pointer to storage containing velocities */

/* get info about a Body */
int ocsmGetBody(void   *modl,           /* (in)  pointer to MODL */
                int    ibody,           /* (in)  Body index (1:nbody) */
                int    *type,           /* (out) Branch type (see below) */
                int    *ichld,          /* (out) ibody of child (or 0 if root) */
                int    *ileft,          /* (out) ibody of left parent (or 0) */
                int    *irite,          /* (out) ibody of rite parent (or 0) */
                double vals[],          /* (out) array  of Argument values (at least 10 long) */
                int    *nnode,          /* (out) number of Nodes */
                int    *nedge,          /* (out) number of Edges */
                int    *nface);         /* (out) number of Faces */

/* print all Bodys to file */
int ocsmPrintBodys(void   *modl,        /* (in)  pointer to MODL */
                   char   filename[]);  /* (in)  file to which output is appended (or "" for stdout) */

/* print the BRep associated with a specific Body */
int ocsmPrintBrep(void   *modl,         /* (in)  pointer to MODL */
                  int    ibody,         /* (in)  Body index (1:nbody) */
                  char   filename[]);   /* (in)  file to which output is appended (or "" for stdout) */

/* evaluate an expression */
int ocsmEvalExpr(void   *modl,          /* (in)  pointer to MODL */
                 char   expr[],         /* (in)  expression */
                 double *value,         /* (out) value */
                 double *dot,           /* (out) velocity */
                 char   str[]);         /* (out) value if string-valued (w/o leading $) */

/* print the contents of an EGADS ego */
void ocsmPrintEgo(
        /*@null@*/ego    obj);          /* (in)  EGADS ego */

/* convert an OCSM code to text */
/*@observer@*/
char *ocsmGetText(int    icode);        /* (in)  code to look up */

/* convert text to an OCSM code */
int ocsmGetCode(char   *text);          /* (in)  text to look up */

/*
 ************************************************************************
 *                                                                      *
 * Defined constants                                                    *
 *                                                                      *
 ************************************************************************
 */

#define           OCSM_DIMENSION  100   /* not Branch or OCSM_UTILITY */
#define           OCSM_LBOUND     101
#define           OCSM_UBOUND     102
#define           OCSM_NAME       103
#define           OCSM_ATTRIBUTE  104
#define           OCSM_CSYSTEM    105

#define           OCSM_POINT      111   /* OCSM_PRIMITIVE */
#define           OCSM_BOX        112
#define           OCSM_SPHERE     113
#define           OCSM_CONE       114
#define           OCSM_CYLINDER   115
#define           OCSM_TORUS      116
#define           OCSM_IMPORT     117
#define           OCSM_UDPRIM     118
#define           OCSM_RESTORE    119

#define           OCSM_EXTRUDE    121   /* OCSM_GROWN */
#define           OCSM_RULE       122
#define           OCSM_LOFT       123
#define           OCSM_BLEND      124
#define           OCSM_REVOLVE    125
#define           OCSM_SWEEP      126

#define           OCSM_FILLET     131   /* OCSM_APPLIED */
#define           OCSM_CHAMFER    132
#define           OCSM_HOLLOW     133
#define           OCSM_CONNECT    134

#define           OCSM_INTERSECT  141   /* OCSM_BOOLEAN */
#define           OCSM_SUBTRACT   142
#define           OCSM_UNION      143
#define           OCSM_JOIN       144
#define           OCSM_EXTRACT    145
#define           OCSM_COMBINE    146

#define           OCSM_TRANSLATE  151   /* OCSM_TRANSFORM */
#define           OCSM_ROTATEX    152
#define           OCSM_ROTATEY    153
#define           OCSM_ROTATEZ    154
#define           OCSM_SCALE      155
#define           OCSM_MIRROR     156
#define           OCSM_APPLYCSYS  157
#define           OCSM_REORDER    158

#define           OCSM_SKBEG      160   /* OCSM_SKETCH */
#define           OCSM_SKVAR      161
#define           OCSM_SKCON      162
#define           OCSM_LINSEG     163
#define           OCSM_CIRARC     164
#define           OCSM_ARC        165
#define           OCSM_ELLARC     166
#define           OCSM_SPLINE     167
#define           OCSM_SSLOPE     168
#define           OCSM_BEZIER     169
#define           OCSM_SKEND      170

#define           OCSM_SOLBEG     171   /* OCSM_SOLVER */
#define           OCSM_SOLCON     172
#define           OCSM_SOLEND     173

#define           OCSM_INTERFACE  174   /* OCSM_UTILITY */
#define           OCSM_END        175
#define           OCSM_SET        176
#define           OCSM_EVALUATE   177
#define           OCSM_GETATTR    178
#define           OCSM_UDPARG     179
#define           OCSM_SELECT     180
#define           OCSM_PROJECT    181
#define           OCSM_MACBEG     182
#define           OCSM_MACEND     183
#define           OCSM_RECALL     184
#define           OCSM_STORE      185
#define           OCSM_PATBEG     186
#define           OCSM_PATBREAK   187
#define           OCSM_PATEND     188
#define           OCSM_IFTHEN     189
#define           OCSM_ELSEIF     190
#define           OCSM_ELSE       191
#define           OCSM_ENDIF      192
#define           OCSM_THROW      193
#define           OCSM_CATBEG     194
#define           OCSM_CATEND     195
#define           OCSM_MARK       196
#define           OCSM_GROUP      197
#define           OCSM_DUMP       198
#define           OCSM_ASSERT     199
#define           OCSM_MESSAGE    200
#define           OCSM_SPECIAL    299

#define           OCSM_PRIMITIVE  201   /* Branch classes */
#define           OCSM_GROWN      202
#define           OCSM_APPLIED    203
#define           OCSM_BOOLEAN    204
#define           OCSM_TRANSFORM  205
#define           OCSM_SKETCH     206
#define           OCSM_SOLVER     207
#define           OCSM_UTILITY    208

#define           OCSM_ACTIVE     300   /* Branch activities (also in ESP.html) */
#define           OCSM_SUPPRESSED 301
#define           OCSM_INACTIVE   302
#define           OCSM_DEFERRED   303

#define           OCSM_SOLID_BODY 400   /* Body types */
#define           OCSM_SHEET_BODY 401
#define           OCSM_WIRE_BODY  402
#define           OCSM_NODE_BODY  403
#define           OCSM_NULL_BODY  404

#define           OCSM_DESPMTR    500   /* Parameter types (also in ESP.html) */
#define           OCSM_CFGPMTR    501
#define           OCSM_CONPMTR    502
#define           OCSM_LOCALVAR   503
#define           OCSM_OUTPMTR    504
#define           OCSM_UNKNOWN    505

#define           OCSM_NODE       600   /* Selector types */
#define           OCSM_EDGE       601
#define           OCSM_FACE       602
#define           OCSM_BODY       603
#define           OCSM_EEDGE      604
#define           OCSM_EFACE      605
#define           OCSM_EBODY      606

#define           OCSM_UNDEFINED  -123.456

/*
 ************************************************************************
 *                                                                      *
 * Return codes (errors are -201 to -299)                               *
 *                                                                      *
 ************************************************************************
 */

#define SUCCESS                                 0

#define OCSM_FILE_NOT_FOUND                  -201
#define OCSM_ILLEGAL_STATEMENT               -202
#define OCSM_NOT_ENOUGH_ARGS                 -203
#define OCSM_NAME_ALREADY_DEFINED            -204
#define OCSM_NESTED_TOO_DEEPLY               -205
#define OCSM_IMPROPER_NESTING                -206
#define OCSM_NESTING_NOT_CLOSED              -207
#define OCSM_NOT_MODL_STRUCTURE              -208
#define OCSM_PROBLEM_CREATING_PERTURB        -209
#define OCSM_EBODY_NOT_FOUND                 -210

#define OCSM_MISSING_MARK                    -211
#define OCSM_INSUFFICIENT_BODYS_ON_STACK     -212
#define OCSM_WRONG_TYPES_ON_STACK            -213
#define OCSM_DID_NOT_CREATE_BODY             -214
#define OCSM_CREATED_TOO_MANY_BODYS          -215
#define OCSM_TOO_MANY_BODYS_ON_STACK         -216
#define OCSM_ERROR_IN_BODYS_ON_STACK         -217
#define OCSM_MODL_NOT_CHECKED                -218
#define OCSM_NEED_TESSELLATION               -219

#define OCSM_BODY_NOT_FOUND                  -221
#define OCSM_FACE_NOT_FOUND                  -222
#define OCSM_EDGE_NOT_FOUND                  -223
#define OCSM_NODE_NOT_FOUND                  -224
#define OCSM_ILLEGAL_VALUE                   -225
#define OCSM_ILLEGAL_ATTRIBUTE               -226
#define OCSM_ILLEGAL_CSYSTEM                 -227
#define OCSM_NO_SELECTION                    -228

#define OCSM_SKETCH_IS_OPEN                  -231
#define OCSM_SKETCH_IS_NOT_OPEN              -232
#define OCSM_COLINEAR_SKETCH_POINTS          -233
#define OCSM_NON_COPLANAR_SKETCH_POINTS      -234
#define OCSM_TOO_MANY_SKETCH_POINTS          -235
#define OCSM_TOO_FEW_SPLINE_POINTS           -236
#define OCSM_SKETCH_DOES_NOT_CLOSE           -237
#define OCSM_SELF_INTERSECTING               -238
#define OCSM_ASSERT_FAILED                   -239

#define OCSM_ILLEGAL_CHAR_IN_EXPR            -241
#define OCSM_CLOSE_BEFORE_OPEN               -242
#define OCSM_MISSING_CLOSE                   -243
#define OCSM_ILLEGAL_TOKEN_SEQUENCE          -244
#define OCSM_ILLEGAL_NUMBER                  -245
#define OCSM_ILLEGAL_PMTR_NAME               -246
#define OCSM_ILLEGAL_FUNC_NAME               -247
#define OCSM_ILLEGAL_TYPE                    -248
#define OCSM_ILLEGAL_NARG                    -249

#define OCSM_NAME_NOT_FOUND                  -251
#define OCSM_NAME_NOT_UNIQUE                 -252
#define OCSM_PMTR_IS_DESPMTR                 -253
#define OCSM_PMTR_IS_LOCALVAR                -254
#define OCSM_PMTR_IS_OUTPMTR                 -255
#define OCSM_PMTR_IS_CONPMTR                 -256
#define OCSM_WRONG_PMTR_TYPE                 -257
#define OCSM_FUNC_ARG_OUT_OF_BOUNDS          -258
#define OCSM_VAL_STACK_UNDERFLOW             -259  /* probably not enough args to func */
#define OCSM_VAL_STACK_OVERFLOW              -260  /* probably too many   args to func */

#define OCSM_ILLEGAL_BRCH_INDEX              -261  /* should be from 1 to nbrch */
#define OCSM_ILLEGAL_PMTR_INDEX              -262  /* should be from 1 to npmtr */
#define OCSM_ILLEGAL_BODY_INDEX              -263  /* should be from 1 to nbody */
#define OCSM_ILLEGAL_ARG_INDEX               -264  /* should be from 1 to narg  */
#define OCSM_ILLEGAL_ACTIVITY                -265  /* should OCSM_ACTIVE or OCSM_SUPPRESSED */
#define OCSM_ILLEGAL_MACRO_INDEX             -266  /* should be between 1 and 100 */
#define OCSM_ILLEGAL_ARGUMENT                -267
#define OCSM_CANNOT_BE_SUPPRESSED            -268
#define OCSM_STORAGE_ALREADY_USED            -269
#define OCSM_NOTHING_PREVIOUSLY_STORED       -270

#define OCSM_SOLVER_IS_OPEN                  -271
#define OCSM_SOLVER_IS_NOT_OPEN              -272
#define OCSM_TOO_MANY_SOLVER_VARS            -273
#define OCSM_UNDERCONSTRAINED                -274
#define OCSM_OVERCONSTRAINED                 -275
#define OCSM_SINGULAR_MATRIX                 -276
#define OCSM_NOT_CONVERGED                   -277

#define OCSM_UDP_ERROR1                      -281
#define OCSM_UDP_ERROR2                      -282
#define OCSM_UDP_ERROR3                      -283
#define OCSM_UDP_ERROR4                      -284
#define OCSM_UDP_ERROR5                      -285
#define OCSM_UDP_ERROR6                      -286
#define OCSM_UDP_ERROR7                      -287
#define OCSM_UDP_ERROR8                      -288
#define OCSM_UDP_ERROR9                      -289

#define OCSM_OP_STACK_UNDERFLOW              -291
#define OCSM_OP_STACK_OVERFLOW               -292
#define OCSM_RPN_STACK_UNDERFLOW             -293
#define OCSM_RPN_STACK_OVERFLOW              -294
#define OCSM_TOKEN_STACK_UNDERFLOW           -295
#define OCSM_TOKEN_STACK_OVERFLOW            -296
#define OCSM_UNSUPPORTED                     -298
#define OCSM_INTERNAL_ERROR                  -299

/* constants for backward compatability */
#define OCSM_EXTERNAL                        OCSM_DESPMTR
#define OCSM_CONFIG                          OCSM_CFGPMTR
#define OCSM_CONSTANT                        OCSM_CONPMTR
#define OCSM_OUTPUT                          OCSM_OUTPMTR
#define OCSM_INTERNAL                        OCSM_LOCALVAR
#define OCSM_PMTR_IS_EXTERNAL                OCSM_PMTR_IS_DESPMTR
#define OCSM_PMTR_IS_CONSTANT                OCSM_PMTR_IS_CONPMTR
#define OCSM_PMTR_IS_OUTPUT                  OCSM_PMTR_IS_OUTPMTR
#define OCSM_PMTR_IS_INTERNAL                OCSM_PMTR_IS_LOCALVAR

#endif  /* _OPENCSM_H_ */
