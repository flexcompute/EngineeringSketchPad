$Id: README.txt,v 1.2 2022/11/10 19:22:55 mnemec Exp $

Support Libraries for Cart3D I/O Functions and Extensible Design
Description Markup
================================================================

Read COPYRIGHT.txt (NOSA.pdf, INDIVIDUAL-CLA.pdf and
CORPORATE-CLA.pdf) files for legal and license terms.


DESCRIPTION

This is a collection of functions used for various I/O functions of
the Cart3D aerodynamic analysis and optimization package. This
includes reading and writing surface triangulation files, volume mesh
files and files for aerodynamic shape optimization
problems. Description of the classic TRI file format is available
here:
https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/cart3dTriangulations.html#1.%20Component%20file%20forma
Description of the newer VTK/VTU file format is available here:
https://kitware.github.io/vtk-examples/site/VTKFileFormats/#xml-file-formats
All triangulation and grid IO functions are in the 'c3dio'
directory. The XDDM library is in the 'xddm' directory.  It implements
the Extensible Design Description Markup specification, which is
provided in the 'doc' directory.


REQUIREMENTS

  1) Linux, macOS or Windows
  2) C compiler, e.g., gcc
  3) libxml2 (https://gitlab.gnome.org/GNOME/libxml2/-/wikis/home)
  4) For Windows: tested with Cygwin. Requires make, libiconv, zlib, liblzma

  To check availability of libxml2, use:
  
    xml2-config --help


INSTALL

  To build IO libraries for reading Cart3D triangulations and grids:

    cd c3dio
    make

  This should build the libc3dio.a library. To test the library:

    make test

  This should build the test_c3dio executable. To install the library,
  tester and headers:

   make install

  The default install location is the current directory: ./lib, ./bin
  and ./include.

  Additional make commands:

    make clean
    make distclean

  To build XDDM:

    cd xddm
    make
    make test
    make install


TEST and EXAMPLES

  Read test/README.txt file.


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
