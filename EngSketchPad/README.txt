                        ESP: The Engineering Sketch Pad
                             Rev 1.25 -- June 2024

                          https://acdl.mit.edu/ESP/

  ******************************************************************
  * THESE DIRECTIONS ARE ONLY REQUIRED IF YOU WILL BE BUILDING ESP *
  *    THEY ARE NOT NEEDED IF YOU USE A PRE-BUILT DISTRIBUTION     *
  ******************************************************************

0. Preamble

    Windows 7 & 8 are no longer supported, only Windows 10 is tested (we have
    not begun testing against Windows 11). This also means that older versions 
    of MS Visual Studio are no longer being tested either. Only MSVS versions 
    2017 and up are fully supported.

    This ESP release no longer works with Python 2.7. The minimum supported
    version is now Python 3.8. Also, we now only support OpenCASCADE at Rev 
    7.6 or higher. And these must be the versions taken from the ESP website
    (and not from elsewhere). At this point we recommend 7.8.0.

    It is advisable to unblock browser tabs on the web browser in use.

    Training material can be found at the ESP website at 
    http://acdl.mit.edu/ESP/Training, which is in 2 parts. The first is on
    ESP geometry construction, labeled with ESP, and the second on analysis, 
    labeled with CAPS. 

    Apple notes: 
    (1) You CANNOT download the distributions using a browser. For instructions 
        on how to get ESP see MACdownloads.txt on the web site.
    (2) You must have XQuartz at a minimum release of 2.8.1 for some supplied
        executables to function.
    (3) Apple arm64 (M#) computers are natively supported but require Rosetta2 
        for the running of some legacy CAPS apps. Rosetta2 can be installed by 
        executing the following command: "softwareupdate --install-rosetta".
    (4) Apple arm64 (M#) builds must be done in a "native" shell. 
        That is, typing "arch" must return "arm64".
    (5) If Safari blocks a pop-up (for example, the flowchart in ESP),
        you can press the rectangular button in the Smart Search field
        and allow the file to be seen.

    Windows notes: 
    (1) You CANNOT download the distributions using a browser. For instructions 
        on how to get ESP see WINdownloads.txt on the web site.

1. Prerequisites

    The most significant prerequisite for this software is OpenCASCADE.
    This ESP release only supports the prebuilt versions marked 7.7.0 
    and 7.8.0, which are available at http://acdl.mit.edu/ESP. Please DO 
    NOT report any problems with any other versions of OpenCASCADE, much 
    effort has been spent in "hardening" the OpenCASCADE code. It is advised 
    that all ESP users update to 7.7.0/7.8.0 because of better robustness and
    performance. If you are still on a LINUX box with a version of gcc less 
    than 4.8, you will have to upgrade to a newer OS or version of gcc.

    Another prerequisite is a WebGL/Websocket capable Browser. In general 
    these include Mozilla's FireFox, Google Chrome and Apple's Safari. 
    Internet Explorer and legacy versions of Edge are NOT supported because 
    of a problem in their Websockets implementation. The "Chromium" version
    of Microsoft Edge is fully supported. Also, note that there have been some
    reports of problems with Intel Graphics and some WebGL Browsers. For 
    LINUX, "zlib" development is required.

    CAPS has a dependency on UDUNITS2, and potentially on Python and other 
    applications. See Section 2.3.

1.1 Source Distribution Layout

    In the discussions that follow, $ESP_ROOT will be used as the name of the 
    directory that contains:

    README.txt        - this file
    bin               - a directory that will contain executables
    CAPSexamples      - a directory with CAPS examples
    config            - files that allow for automatic configuration
    contributions     - user contributions
    data              - test and example scripts
    doc               - documentation
    ESP               - web client code for the Engineering Sketch Pad
    externApps        - the ESP connections to 3rd party Apps (outside of CAPS)
    include           - location for all ESP header files
    lib               - a directory that will contain libraries, shared objects
                        and DLLs
    LICENSE.txt       - the GNU Lesser General Public license (LGPL 2.1) text
    pyESP             - Python bindings
    SLUGS             - the browser code for Slugs (web Slugs client)
    src               - source files (contains EGADS, CAPS, wvServer & OpenCSM)
    udc               - a collection of User Defined Components
    wvClient          - simple examples of Web viewing used by EGADS

1.2 Release Notes

1.2.1 EGADS

    The significant updates made to EGADS from Rev 1.24 are:

    * Upgrade to OpenCASCADE 7.8
    * Add Face Color attribute for STEP/IGES
    * Improved Name attribute for STEP/IGES
    * Allow copyObject on PCurve without transformation
    * Use consistent CURVE/SURFACE approximation algorithms and preserve periodicity in EG_convertToBspline
    * Add EG_addKnots
    * Add EG\_isIsoPCurve (only works for OCC 7.6 or above)
    * Enable FullAttr for EG_imprintBody
    * EG_matchBodyFaces and EG_matchBodyEdges API change, please see egads.pdf for details
    * Check that all edges are imprinted by OCC in EG\_imprintBody
    * General bug fixes

1.2.2 OpenCSM

    In addition to many big fixes (see $ESP_ROOT/src/OpenCSMnotes.txt
    for a full list), the significant upgrades are documented in section 
    8.1 of ESP-help.html; bug fixes are documented in section 8.2 of the 
    same document.

1.2.3 CAPS

    * Refactor mesh writing fremwork
      - Mesh files may be requested with Mesh_Format, and file extensions may have changed
      - Removed Mesh_ASCII_Flag AIM inputs
    * Upgrade to AFLR 11.5.14
    * Updated config file generator for SU2-8.0.1
    * Add cbaero AIM
    * Add zaero AIM
    * Add ability to link mass properties
    * Add CLAF avl input
    * Add MassProp outputs for nastran.
    * OUTPMTR derivatives computed by AIM's are properly cached
    * Fix Normaliztion misspelling
    * Updates to refine AIM to support 2D meshes
    * General bug fixes

1.2.4 ESP

    * Transition to a new server: serveESP
    * A complete "Integrated Design Environment" now exists in ESP. In the
      help look at Tutorial #6 to get a flavor of what you can now do.

1.2.5 Known issues in v1.25:

    * data/fighter4 does not functoin with Intel macOS 13.3 and gerater
    * AFLR2 does not always generate valid quads


2. Building the Software

    The config subdirectory contains scripts that need to be used to generate 
    the environment both to build and run the software here. There are two 
    different build procedures based on the OS:

    If using Windows, skip to section 2.2.

2.1 Linux and MAC OS

    The configuration is built using the path where the OpenCASCADE runtime 
    distribution can be found.  This path can be located in an OpenCASCADE 
    distribution by looking for a subdirectory that includes an "inc" or 
    "include" directory and either a "lib" or "$ARCH/lib" (where $ARCH is 
    the name of your architecture) directory.  Once that is found, execute 
    the commands:

        % cd $ESP_ROOT/config
        % ./makeEnv **absolute_path_of_OpenCASCADE_directory_containing_inc_and_lib**

    An optional second argument to makeEnv is required if the distribution 
    of OpenCASCADE has multiple architectures. In this case it is the 
    subdirectory name that contains the libraries for the build of interest 
    (CASARCH). Apple arm64 (M1/M2) CPUs should indicate "DARWIN_ARM64" as 
    the architecture.

    This procedure produces 2 files at the top level: ESPenv.sh and
    ESPenv.csh.  These are the environments for both sh (bash/zsh) and csh 
    (tcsh) respectively.  The appropriate file can be "source"d or included 
    in the user's startup scripts. This must be done before either building 
    and/or running the software. 
    
    For example, if using the csh or tcsh:

        % cd $ESP_ROOT
        % source ESPenv.csh

    or if using bash/zsh:

        $ cd $ESP_ROOT
        $ source ESPenv.sh

    Skip to section 2.3.

2.2 Windows Configuration

    IMPORTANT: The ESP distribution and OpenCASCADE MUST be unpackaged 
               into a location ($ESP_ROOT) that has NO spaces in the path!

    The configuration is built from the path where the OpenCASCADE runtime 
    distribution can be found. MS Visual Studio is required and a command 
    shell where the 64bit C/C++ compiler (use x64, not x86) should be opened 
    and the following executed in that window (note that MS VS 2017, 2019, 
    and 2022 are fully supported). The Windows environment is built simply 
    by going to the config subdirectory and executing the script "winEnv" 
    in a bash shell (run from the command window):

        C:\> cd %ESP_ROOT%\config
        C:\> bash winEnv D:\OpenCASCADE7.8.0

    winEnv (like makeEnv) has an optional second argument that is only 
    required if the distribution of OpenCASCADE has multiple architectures. 
    In this case it is the subdirectory name that contains the libraries 
    for the build of interest (CASARCH).

    This procedure produces a single file at the top level: ESPenv.bat.  
    This file needs to be executed before either building and/or running 
    the software. This is done with:

        C:\> cd %ESP_ROOT%
        C:\> ESPenv

    Check that the method that you used to unzip the distribution created 
    directories named %ESP_ROOT%\bin and %ESP_ROOT%\lib. If it did not, 
    create them using the commands:

        C:\> cd %ESP_ROOT%
        C:\> mkdir bin
        C:\> mkdir lib

    Also note that the winEnv script will find the version of Python if in
    the PATH and on the same drive.

2.3 CAPS Options

    CAPS requires the use of the Open Source Project UDUNITS2 for all unit 
    conversions. Since there are no prebuilt package distributions for the 
    MAC and Windows, the CAPS build procedure copies prebuilt DLL/DyLibs 
    to the lib directory of ESP. Because most Linux distributions contain 
    a UDUNITS2 package, another procedure is used. If the UDUNITS2 
    development package is loaded in the OS, then nothing is done. If not 
    loaded, then a pre-built shared object is moved to the ESP lib directory.
    Be careful if you install UDUNITS2 from your OS package manager at some
    later time.

2.3.1 Python with CAPS (pyCAPS)

    Python may be used with CAPS to provide testing, general scripting and
    demonstration capabilities. The execution of pyCAPS requires a single 
    environment variable:

    PYTHONPATH is a Python environment variable that needs to have the path
               $ESP_ROOT/pyESP included.

2.3.2 3rd Party Environment Variables

    CAPS is driven by a plugin technology using AIMs (Analysis Interface  
    Modules). These AIMs allow direct coupling between CAPS and the external 
    meshers and solvers. Many are built by default (where there are no 
    licensing problems or other dependencies). The CAPS build subsystem will 
    respond to the following (if these are not set, then the AIMs for these 
    systems will not be built):

    AFLR (See Section 4.1):
      AFLR      is the path where the AFLR distribution has been deposited
      AFLR_ARCH is the architecture flag to use (MacOSX-arm64, MacOSX-x86-64, 
                Linux-x86-64 or WIN64) -- note that this is set by the config 
                procedure

    AWAVE
      AWAVE is the location to find the FORTRAN source for AWAVE

    TETGEN
      TETGEN is the path where the TetGen distribution has been unpacked

    Some of the AIMs have Python embedded. Building these AIMs with
    Python embedding is enabled by 2 environment variables (the Python 
    development package is required under Linux):
    
    PYTHONINC  is the include path to find the Python includes for building
    PYTHONLIB  is a description of the location of the Python libraries and
               the library to use

    The exact same version of Python that was used to compile the embedding
    must be used when executing Python scripts.

    For MACs and LINUX the configuration procedure inserts these environment 
    variables with the locations it finds by executing the version of Python 
    available in the shell performing the build. If makeEnv emits any errors
    related to Python, the resultant environment file(s) will need to be 
    updated in order to use Python in the AIMs (the automatic procedure has 
    failed).

    For Windows ESPenv.bat may need be edited (unless configured from a 
    command prompt that has both the MSVS and Python environments), the "rem" 
    removed and the appropriate information set (if Python exists on the 
    machine). Also note that the bit size (32 or 64) of Python that gets 
    used on Windows must be consistent with the build of ESP, which is 
    64bit.

    For Example on Windows (after downloading and installing Python on C:):
      set PYTHONINC=C:\Python311\include
      set PYTHONLIB=C:\Python311\Libs\python311.lib
      set PYTHONPATH=%ESP_ROOT%\lib:%ESP_ROOT%\pyESP

2.3.3 The Cart3D Design Framework

    The application ESPxddm is the ESP connection to the Cart3D Design
    Framework. On LINUX this requires that the libxml2 development package
    be installed. If it is, then ESPxddm is built, otherwise it is not.

2.3.4 CAPS AIM Documentation

    The CAPS documentation can be seen in PDF form from within the directory
    $ESP_ROOT/doc/CAPSdoc. Or in html by $ESP_ROOT/doc/CAPSdoc/html/index.html.

2.4 The Build

    For any of the operating systems, after properly setting the environment 
    in the command window (or shell), follow this simple procedure:

        % cd $ESP_ROOT/src
        % make

    or

        C:\> cd $ESP_ROOT\src
        C:\> make

    You can use "make clean" which will clean up all object modules or 
    "make cleanall" to remove all objects, executables, libraries, shared 
    objects and dynamic libraries.


3.0 Running

3.1 serveESP

    To start ESP there are two steps: (1) start the "server" and (2) start 
    the "browser". This can be done in a variety of ways, but the two most 
    common follow.

3.1.1 Procedure 1: have ESP automatically started

    If it exists, the ESP_START environment variable contains the command 
    that should be executed to start the browser once the server has 
    created its scene graph.  On a MAC, you can set this variable with 
    commands such as

        % setenv ESP_START "open -a /Applications/Firefox.app $ESP_ROOT/ESP/ESP.html"

    or

        % export ESP_START="open -a /Applications/Firefox.app $ESP_ROOT/ESP/ESP.html"

    depending on the shell in use.  The commands in other operating systems 
    will differ slightly, depending on how the browser can be started from 
    the command line, for example for Windows it may be:

        % set ESP_START="start /b "C:\Program Files (x86)\Mozilla Firefox\firefox.exe" %ESP_ROOT%\ESP\ESP.html"

    To run the program, use:

         % cd $ESP_ROOT/bin
         % ./serveESP ../data/tutorial1

3.1.2 Procedure 2: start the browser manually

    If the ESP_START environment variable does not exist, issuing the
    commands:

        % cd $ESP_ROOT/bin
        % ./serveESP ../data/tutorial1

    will start the server.  The last lines of output from serveESP tells 
    the user that the server is waiting for a browser to attach to it. 
    This can be done by starting a browser (FireFox and GoogleChrome have 
    been tested) and loading the file:

        $ESP_ROOT/ESP/ESP.html

    Whether you used procedure 1 or 2, as long as the browser stays connected 
    to serveESP, serveESP will stay alive and handle requests sent to it from 
    the browser. Once the last browser that is connected to serveESP exits, 
    serveESP will shut down.

    Note that the default "port" used by serveESP is 7681. One can change 
    the port in the call to serveESP with a command such as:

        % cd $ESP_ROOT/bin
        % ./serveESP ../data/tutorial1 -port 7788

    Once the browser starts, you will be prompted for a "hostname:port".  
    Make the appropriate response depending on the network situation. Once 
    the ESP GUI is functional, press the "help" button in the upper left 
    if you want to execute the tutorial.

3.2 egads2cart

    This example takes an input geometry file and generates a Cart3D "tri" 
    file. The acceptable input is STEP, EGADS or OpenCASCADE BRep files 
    (which can all be generated from an OpenCSM "dump" command).

        % cd $ESP_ROOT/bin
        % ./egads2cart geomFilePath [angle relSide relSag]

3.3 vTess and wvClient

    vTess allows for the examination of geometry through its discrete
    representation. Like egads2cart, the acceptable geometric input is STEP, 
    EGADS or OpenCASCADE BRep files. vTess acts like serveESP and wvClient 
    should be used like ESP in the discussion in Section 3.1 above.

        % cd $ESP_ROOT/bin
        % ./vTess geomFilePath [angle maxlen sag]

3.4 Executing CAPS through Python

    % python pyCAPSscript.py  (Note: many example Python scripts can be 
                                     found in $ESP_ROOT/CAPSexamples/pyCAPS)


4.0 Notes on 3rd Party Analyses

4.1 The AFLR suite

    Building the AFLR AIMs (AFLR2, AFLR3 and AFLR4) requires AFLR_LIB at
    11.5.14. Note that built versions of the so/DLLs are now provided
    with the ESP source. There is no longer need to copy them from the PreBuilt.

4.2 Athena Vortex Lattice

    The interface to AVL is designed for V3.40, and the avl executable
    is provided in $ESP_ROOT/bin.
    
4.3 Astros and mAstros

    Both Astros and a limited version (microAstros or more simply) mAstros
    can be run with the Astros AIM. An mAstros executable is part of this
    distribution and is used with the CAPS training material. The pyCAPS 
    Astros examples use the environment variable ASTROS_ROOT (set to the 
    path where Astros can be found) to locate Astros and its runtime files.
    If not set it defaults to mAstros execution. 

4.4 Cart3D

    The interfaces to Cart3D will only work with V1.5.9.

4.5 Fun3D

    The Fun3D AIM supports Fun3D V12.4 or higher.

4.6 Mystran

    On Windows, follow the install instructions MYSTRAN-Install-Manual.pdf
    carefully. The CAPS examples function only if MYSTRAN.ini in the
    Mystran bin directory is an empty file and the MYSTRAN_directory 
    environment variable points to the Mystran bin directory.

4.7 Pointwise

    The CAPS connection to Pointwise is now handled internally but requires,
    at a minimum Pointwise V18.2 R2, but V18.4 or higher is recommended. This 
    setup performs automatic unstructured meshing. Note that the environment 
    variable CAPS_GLYPH is set by the ESP configure script and points to the 
    Glyph scripts that should be used with CAPS and the current release of 
    Pointwise.

4.8 SU2

    Supported versions are: 
        4.1.1 (Cardinal)
        5.0.0 (Raven)
        6.2.0 (Falcon)
        7.5.1 (Blackbird)
        8.0.1 (Harrier)

    SU2 version 6.0 will work except for the use of 
    displacements in a Fluid/Structure Interaction setting.

4.9 xfoil

    The interface to xfoil is designed for V6.99, and the xfoil executable
    is provided in $ESP_ROOT/bin. Note that multiple 'versions' of xfoil
    6.99 have been released with differing output file formats.
