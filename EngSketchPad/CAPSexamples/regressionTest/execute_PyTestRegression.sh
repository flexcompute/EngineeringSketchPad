#!/bin/bash -e

stat=0
notExpectedTests=""
notRun=""
LOGFILENAME=pyCAPSlog.txt

echo
echo "================================================="
echo "Using python : `which python`"
echo "================================================="
echo


expectPythonSuccess ()
{
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    local status
    set -x
    echo "$1 test;" | tee -a $LOGFILE
    time $VALGRIND_COMMAND python -u $1 -outLevel=0 $2  | tee -a $LOGFILE; status=${PIPESTATUS[0]}
    set +x
    if [[ $status == 0 ]]; then
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo "CAPS verification (using pyCAPS) case $1 passed (as expected)" | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    else
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo "CAPS verification (using pyCAPS) case $1 failed (as NOT expected)" | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    stat=1
    notExpectedTests="$notExpectedTests\n$1"
    fi
}

# function to run a case that is expected to fail
expectPythonFailure ()
{
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    local status
    set -x
    echo "$1 test;" | tee -a $LOGFILE
    time $VALGRIND_COMMAND python -u $1 -outLevel=0 $2 | tee -a $LOGFILE; status=${PIPESTATUS[0]}
    set +x
    if [[ $status == 0 ]]; then
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo "CAPS verification (using pyCAPS) case $1 passed (as NOT expected)" | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    stat=1
    notExpectedTests="$notExpectedTests\n$1"
    else
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo "CAPS verification (using pyCAPS) case $1 failed (as expected)" | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    fi
}

# Get the current directory
curDir=${PWD##*/} 
regTestDir=${PWD}

if [[ "$curDir" != "regressionTest" || ! -d "../csmData" || ! -d "../pyCAPS" ]]; then
    echo "execute_PyTestRegression.sh must be executed in the regressionTest directory!"
    echo "Current directory: $regTestDir"
    exit 1
fi

# make the log file abosolute path to the current directory
LOGFILE=$regTestDir/$LOGFILENAME

# Remove old log file
rm -f $LOGFILE
touch $LOGFILE

# Move into the pyCAPS directory
cd ../pyCAPS

######### Execute pyTest examples ########

TYPE="ALL"
testsRan=0

if [ "$1" != "" ]; then
    TYPE="$1"
else
    echo "No options provided...defaulting to ALL"
fi

if [ "$OS" == "Windows_NT" ]; then
    EXT=dll
else
    EXT=so
fi

#### Some minimal fast running tests ####
if [ "$TYPE" == "MINIMAL" ]; then
    echo "Running.... MINIMAL PyTests"

    if [[ `command -v avl` ]]; then
         expectPythonSuccess "avl_PyTest.py"
    else
        notRun="$notRun\navl"
    fi

    if [[ `command -v xfoil` ]]; then
        expectPythonSuccess "xfoil_PyTest.py" -noPlotData
    else
        notRun="$notRun\nxfoil"
    fi

    if [[ `command -v awave` ]]; then
        expectPythonSuccess "awave_PyTest.py"
    else
        notRun="$notRun\nawave"
    fi

    if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT ]]; then
        expectPythonSuccess "aflr2_PyTest.py"
    else
        notRun="$notRun\nAFLR"
    fi

    if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
        expectPythonSuccess "delaundo_PyTest.py"
    else
        notRun="$notRun\ndelaundo"
    fi

    if [[ "$ASTROS_ROOT" != "" ]]; then
        expectPythonSuccess "astros_ThreeBar_PyTest.py"
        expectPythonSuccess "astros_Flutter_15degree.py"
    else
        notRun="$notRun\nAstros"
    fi

    expectPythonSuccess "masstran_PyTest.py"

    ######  HSM ###### 
    if [[ -f $ESP_ROOT/lib/hsmAIM.$EXT ]]; then
        expectPythonSuccess "hsm_SingleLoadCase_PyTest.py"
    else
        notRun="$notRun\nHSM"
    fi 

    testsRan=1
fi

###### Linear aero. examples ######
if [[ "$TYPE" == "LINEARAERO" || "$TYPE" == "ALL" ]]; then
    echo "Running.... LINEARAERO PyTests"
    
    ###### AVL ######
    if [[ `command -v avl` ]]; then
        expectPythonSuccess "avl_PyTest.py"
        expectPythonSuccess "avl_AutoSpan_PyTest.py"
        expectPythonSuccess "avl_DesignSweep_PyTest.py" -noPlotData
        expectPythonSuccess "avl_ControlSurf_PyTest.py"
        expectPythonSuccess "avl_EigenValue_PyTest.py"
        if ( python -c 'import openmdao' ); then
            expectPythonSuccess "avl_OpenMDAO_3_PyTest.py"
        else
            notRun="$notRun\avl_OpenMDAO_3_PyTest.py"
        fi
    else
        notRun="$notRun\navl"
    fi
    
    ###### awave ######
    if [[ `command -v awave` ]]; then
        expectPythonSuccess "awave_PyTest.py"
    else
        notRun="$notRun\nawave"
    fi
    
    ###### friction ######
    if [[ `command -v friction` ]]; then
        expectPythonSuccess "friction_PyTest.py"
    else
        notRun="$notRun\nfriction"
    fi
    
    ###### tsFoil ######
    if [[ `command -v tsfoil2` && "$OS" != "Windows_NT" ]]; then
        expectPythonSuccess "tsfoil_PyTest.py"
    else
        notRun="$notRun\ntsfoil"
    fi
    
    ###### XFoil ######
    if [[ `command -v xfoil` ]]; then
        expectPythonSuccess "xfoil_PyTest.py" -noPlotData
    else
        notRun="$notRun\nxfoil"
    fi

    ###### AutoLink with awave, friction, and avl ######
    if [[ `command -v awave` && `command -v friction` && `command -v avl` ]]; then
        expectPythonSuccess "autoLink_Pytest.py" -noPlotData
    else
        notRun="$notRun\nautoLink"
    fi

    testsRan=1
fi

###### Meshing examples ######
if [[ "$TYPE" == "MESH" || "$TYPE" == "ALL" ]]; then
    echo "Running.... MESH PyTests"

    # EGADS Tess
    expectPythonSuccess "egadsTess_PyTest.py"
    expectPythonSuccess "egadsTess_Spheres_Quad_PyTest.py"
    expectPythonSuccess "egadsTess_Nose_Quad_PyTest.py"
 
    # Example of controlling quading
    expectPythonSuccess "quading_Pytest.py"
    
    # AFLR
    if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT && \
          -f $ESP_ROOT/lib/aflr3AIM.$EXT && \
          -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
        expectPythonSuccess "aflr2_PyTest.py"
        expectPythonSuccess "aflr4_PyTest.py"
        expectPythonSuccess "aflr4_Symmetry_PyTest.py"
        expectPythonSuccess "aflr4_Generic_Missile.py"
        expectPythonSuccess "aflr4_TipTreat_PyTest.py" -noPlotData
        expectPythonSuccess "aflr4_and_aflr3_PyTest.py"
    
        if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "aflr4_and_Tetgen_PyTest.py"
            expectPythonSuccess "aflr4_tetgen_Regions_PyTest.py"
        fi
    else
        notRun="$notRun\nAFLR"
    fi
    
    # Tetgen
    if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
        expectPythonSuccess "tetgen_PyTest.py"

        # TODO: Sort out why this hangs on windoze...
        if [ "$OS" != "Windows_NT" ]; then
            expectPythonSuccess "tetgen_Holes_PyTest.py"
            expectPythonSuccess "tetgen_Regions_PyTest.py"
        fi
    else
        notRun="$notRun\nTetGen"
    fi
   
    # Delaundo 
    if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
        expectPythonSuccess "delaundo_PyTest.py"
    else
        notRun="$notRun\ndelaundo"
    fi

    # pointwise
    if [[ `command -v pointwise` && -f $ESP_ROOT/lib/pointwiseAIM.$EXT ]]; then
        expectPythonSuccess "pointwise_PyTest.py"
        expectPythonSuccess "pointwise_Symmetry_PyTest.py"
    else
        notRun="$notRun\nPointwise"
    fi

    testsRan=1
fi

###### CFD examples ######
if [[ "$TYPE" == "CFD" || "$TYPE" == "ALL" ]]; then
    echo "Running.... CFD PyTests"
    
    if [ "$SU2_RUN" != "" ] && [ "$OS" != "Windows_NT" ]; then
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT ]]; then
            expectPythonSuccess "su2_and_AFLR2_NodeDist_PyTest.py"
            expectPythonSuccess "su2_and_AFLR2_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/aflr3AIM.$EXT && 
              -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
            expectPythonSuccess "su2_and_AFLR4_AFLR3_PyTest.py"
        fi
        if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
            expectPythonSuccess "su2_and_Delaundo_PyTest.py"
        fi
        if [[ `command -v pointwise` ]]; then
            expectPythonSuccess "su2_and_pointwise_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "su2_and_Tetgen_PyTest.py"
            expectPythonSuccess "su2_X43a_PyTest.py"
        fi
    else
        notRun="$notRun\nSU2"
    fi
   
    if [[ `command -v flowCart` ]]; then
        ulimit -s unlimited || true # Cart3D requires unlimited stack size
        expectPythonSuccess "cart3d_PyTest.py"
        if ( python -c 'import openmdao' ); then
            expectPythonSuccess "cart3d_OpenMDAO_3_alpha_PyTest.py"
            expectPythonSuccess "cart3d_OpenMDAO_3_twist_PyTest.py"
        else
            notRun="$notRun\ncart3d_OpenMDAO_3_alpha_PyTest.py"
            notRun="$notRun\ncart3d_OpenMDAO_3_twist_PyTest.py"
        fi
    else
        notRun="$notRun\nCart3D"
    fi

    if [[ `command -v nodet_mpi` ]]; then
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT &&
              -f $ESP_ROOT/lib/aflr3AIM.$EXT &&
              -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
            expectPythonSuccess "fun3d_and_AFLR2_NodeDist_PyTest.py"
            expectPythonSuccess "fun3d_and_AFLR2_PyTest.py"
            expectPythonSuccess "fun3d_and_AFLR4_AFLR3_PyTest.py"
        fi
        if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
            expectPythonSuccess "fun3d_and_Delaundo_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "fun3d_and_Tetgen_AlphaSweep_PyTest.py"
            expectPythonSuccess "fun3d_and_Tetgen_PyTest.py"
        fi
        expectPythonSuccess "fun3d_PyTest.py"
    else
        notRun="$notRun\nFun3D"
    fi

    testsRan=1
fi

###### Structural examples ######
if [[ "$TYPE" == "STRUCTURE" || "$TYPE" == "ALL" ]]; then
    echo "Running.... STRUCTURE PyTests"

    ###### masstran ######
    expectPythonSuccess "masstran_PyTest.py"
    expectPythonSuccess "masstran_AGARD445_PyTest.py"

    ######  HSM ###### 
    if [[ -f $ESP_ROOT/lib/hsmAIM.$EXT ]]; then
        expectPythonSuccess "hsm_SingleLoadCase_PyTest.py"
        expectPythonSuccess "hsm_CantileverPlate_PyTest.py"
    else
        notRun="$notRun\nHSM"
    fi 

    ###### Astros ######
    if [[ "$ASTROS_ROOT" != "" ]]; then
        expectPythonSuccess "astros_Aeroelastic_PyTest.py"
        expectPythonSuccess "astros_AGARD445_PyTest.py"
        expectPythonSuccess "astros_Composite_PyTest.py"
        expectPythonSuccess "astros_CompositeWingFreq_PyTest.py"
        expectPythonSuccess "astros_Flutter_15degree.py"
        expectPythonSuccess "astros_MultipleLoadCase_PyTest.py"
        expectPythonSuccess "astros_PyTest.py"
        expectPythonSuccess "astros_SingleLoadCase_PyTest.py"
        expectPythonSuccess "astros_ThreeBarFreq_PyTest.py"
        expectPythonSuccess "astros_ThreeBarLinkDesign_PyTest.py"
        expectPythonSuccess "astros_ThreeBarMultiLoad_PyTest.py"
        expectPythonSuccess "astros_ThreeBar_PyTest.py"
        expectPythonSuccess "astros_Trim_15degree.py"
 
        # TODO: Create these examples
        #expectPythonSuccess "astros_CompositeWingDesign_PyTest.py"
        #expectPythonSuccess "astros_OptimizationMultiLoad_PyTest.py"
        #expectPythonSuccess "astros_ThreeBarDesign_PyTest.py"
    else
        notRun="$notRun\nAstros"
    fi

    ###### Nastran ######
    if [[ `command -v nastran` ]]; then
        expectPythonSuccess "nastran_Aeroelastic_PyTest.py"
        expectPythonSuccess "nastran_AGARD445_PyTest.py"
        expectPythonSuccess "nastran_Composite_PyTest.py"
        expectPythonSuccess "nastran_CompositeWingDesign_PyTest.py"
        expectPythonSuccess "nastran_CompositeWingFreq_PyTest.py"
        expectPythonSuccess "nastran_Flutter_15degree.py"
        expectPythonSuccess "nastran_MultiLoadCase_PyTest.py"
        expectPythonSuccess "nastran_OptimizationMultiLoad_PyTest.py"
        expectPythonSuccess "nastran_PyTest.py"
        expectPythonSuccess "nastran_SingleLoadCase_PyTest.py"
        expectPythonSuccess "nastran_ThreeBarDesign_PyTest.py"
        expectPythonSuccess "nastran_ThreeBarFreq_PyTest.py"
        expectPythonSuccess "nastran_ThreeBarLinkDesign_PyTest.py"
        expectPythonSuccess "nastran_ThreeBarMultiLoad_PyTest.py"
        expectPythonSuccess "nastran_ThreeBar_PyTest.py"
        expectPythonSuccess "nastran_Trim_15degree.py"
    else
        notRun="$notRun\nNastran"
    fi

    ###### Mystran ######
    if [[ `command -v mystran.exe` ]]; then
        if [[ "$OS" != "Windows_NT" ]]; then
            # Mystran runs out of heap memory on windows running this example
            expectPythonSuccess "mystran_PyTest.py"
        fi
        expectPythonSuccess "mystran_AGARD445_PyTest.py"
        expectPythonSuccess "mystran_SingleLoadCase_PyTest.py"
        expectPythonSuccess "mystran_MultiLoadCase_PyTest.py"
    else
        notRun="$notRun\nMystran"
    fi

    testsRan=1
fi

###### Aeroelastic examples ######
if [[ "$TYPE" == "AEROELASTIC" || "$TYPE" == "ALL" ]]; then
    echo "Running.... AEROELASTIC PyTests"

    ###### Astros ######
    if [[ "$ASTROS_ROOT" != "" && "$SU2_RUN" != "" && "$OS" != "Windows_NT" ]]; then
        # SU2 6.0.0 on Windows is no longer supported
        expectPythonSuccess "aeroelasticSimple_Pressure_SU2_and_Astros.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Displacement_SU2_and_Astros.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Iterative_SU2_and_Astros.py"

        expectPythonSuccess "aeroelastic_Iterative_SU2_and_Astros.py"
    else
        notRun="$notRun\nAstros"
    fi

    ###### Mystran ######
    if [[ `command -v mystran.exe` && "$SU2_RUN" != "" && "$OS" != "Windows_NT" ]]; then
        # SU2 6.0.0 on Windows does not work with displacements
        expectPythonSuccess "aeroelasticSimple_Pressure_SU2_and_Mystran.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Displacement_SU2_and_Mystran.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Iterative_SU2_and_Mystran.py"
    else
        notRun="$notRun\nMystran and SU2"
    fi

    testsRan=1
fi

echo ""
echo "*************************************************"
echo "*************************************************"
echo ""

# Make sure something actually was executed
if [[ $testsRan == 0 ]]; then
    echo "error: unknown test selection : $TYPE"
    echo
    exit 1
fi

if [[ "$notRun" != "" ]]; then
    echo "================================================="
    echo "Did not run examples for:"
    printf $notRun
    echo ""
    echo "================================================="
fi

if [[ $stat != 0 ]]; then
    echo "================================================="
    echo "Summary of examples NOT finished as expected"
    printf $notExpectedTests
    echo ""
    echo ""
    echo "See $LOGFILENAME for more details."
    echo "================================================="
else
    echo ""
    echo "All tests pass!"
    echo ""
fi

TESSFAIL=`grep EG_fillTris $LOGFILE | wc -l`

if [[ $TESSFAIL != 0 ]]; then
    echo ""
    echo "================================================="    
    awk '/test;$/ { test = $0 }; /EG_fillTris/ { printf("%-35s EGADS Tess Failure\n", test, $0) }' $LOGFILE
    echo "================================================="
fi

# exit with the status
exit $stat
