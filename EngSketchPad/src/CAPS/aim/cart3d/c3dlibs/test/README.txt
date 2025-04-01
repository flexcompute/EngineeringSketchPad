$Id: README-TEST.txt,v 1.1.1.1 2022/11/08 18:23:13 mnemec Exp $


C3DIO TESTER
------------

To test reading and writing surface triangulations, make sure that the
test_c3dio executable is compiled and installed in '../bin'. See '../
README.txt' for instructions.

To read and write triangulation files in the traditional Cart3D
format:

  ../bin/test_c3dio cube.tri test.tri

To read and write extended triangulation files (VTK/VTU):

  ../bin/test_c3dio cube_vtk.tri test_vtk.tri

In both cases, the generated output should be the same as the input, i.e.,

  diff cube.tri test.tri
  diff cube_vtk.tri test_vtk.tri

should not find any differences.  The sample test code
../c3dio/test_c3dio.c demonstrates basic usage of the libc3dio.a library. 


XDDM TESTER
-----------

To test reading and writing XDDM files, make sure that the test_xddm
executable is compiled and installed in '../bin/test_xddm'. See '../
README.txt' for instructions.

To parse the content of an XDDM file for given XPath expression, use

  ../bin/test_xddm basic.xml /Model

The test generates xtest_echo.xml that should be identical to
basic.xml:

  diff basic.xml xtest_echo.xml

The test also writes xtest_opt.xml that echoes a node created in
../xddm/main.c.  The XDDM specification is provided in '../doc'.


COPYRIGHT

Copyright © 2022 United States Government as represented by the
Administrator of the National Aeronautics and Space Administration.
All Rights Reserved.


DISCLAIMERS

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE
WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT
OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY
OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.
FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES
REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE,
AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF
THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES,
EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM
PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT
SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED
STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE
REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL
TERMINATION OF THIS AGREEMENT.
