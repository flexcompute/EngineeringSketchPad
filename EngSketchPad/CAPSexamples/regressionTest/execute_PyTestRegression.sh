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

    ###### abaqus ######
    if [[ `command -v abaqus` ]]; then
         expectPythonSuccess "abaqus_SingleLoadCase_PyTest.py"
    else
        notRun="$notRun\nabaqus"
    fi

    ###### avl ######
    if [[ `command -v avl` ]]; then
         expectPythonSuccess "avl_PyTest.py"
    else
        notRun="$notRun\navl"
    fi

    ###### xfoil ######
    if [[ `command -v xfoil` ]]; then
        expectPythonSuccess "xfoil_PyTest.py" -noPlotData
    else
        notRun="$notRun\nxfoil"
    fi

    ###### mses ######
    if [[ `command -v mses` ]]; then
        expectPythonSuccess "mses_PyTest.py" -noPlotData
    else
        notRun="$notRun\nmses"
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
    
    ###### delaundo ######
    if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
        expectPythonSuccess "delaundo_PyTest.py"
    else
        notRun="$notRun\ndelaundo"
    fi

    ###### EGADS Tess ######
    expectPythonSuccess "egadsTess_PyTest.py"
    
    ###### AFLR ######
    if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT && \
          -f $ESP_ROOT/lib/aflr3AIM.$EXT && \
          -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
        expectPythonSuccess "aflr2_PyTest.py"
        expectPythonSuccess "aflr4_and_aflr3_PyTest.py"
    else
        notRun="$notRun\nAFLR"
    fi
    
    ###### Tetgen ######
    if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
        expectPythonSuccess "tetgen_PyTest.py"
    else
        notRun="$notRun\nTetGen"
    fi

    ######  SU2 ###### 
    if [ "$SU2_RUN" != "" ] && [ "$OS" != "Windows_NT" ]; then
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT ]]; then
            expectPythonSuccess "su2_and_AFLR2_PyTest.py"
        else
            notRun="$notRun\nSU2 and AFLR2"
        fi
    else
        notRun="$notRun\nSU2"
    fi
   
    ######  Cart3D ###### 
    if [[ `command -v flowCart` ]]; then
        ulimit -s unlimited || true # Cart3D requires unlimited stack size
        expectPythonSuccess "cart3d_PyTest.py"
        #notRun="$notRun\nCart3D"
    else
        notRun="$notRun\nCart3D"
    fi

    ######  ASTROS ###### 
    if [[ "$ASTROS_ROOT" != "" ]]; then
        expectPythonSuccess "astros_ThreeBar_PyTest.py"
        expectPythonSuccess "astros_Flutter_15degree.py"
    else
        notRun="$notRun\nAstros"
    fi

    ######  MASSTRAN ###### 
    expectPythonSuccess "masstran_PyTest.py"

    ######  HSM ###### 
    if [[ -f $ESP_ROOT/lib/hsmAIM.$EXT ]]; then
        expectPythonSuccess "hsm_SingleLoadCase_PyTest.py"
    else
        notRun="$notRun\nHSM"
    fi 

    ###### Mystran ######
    if [[ ( `command -v mystran.exe` || `command -v mystran` ) && "$OS" != "Windows_NT" ]]; then
        expectPythonSuccess "mystran_PyTest.py"
    else
        notRun="$notRun\nMystran"
    fi
    
    ######  CORSAIRLITE ###### 
    expectPythonSuccess "../corsairlite/qp.py"
    
    testsRan=1
fi

###### Linear aero. examples ######
if [[ "$TYPE" == "LINEARAERO" || "$TYPE" == "ALL" ]]; then
    echo "Running.... LINEARAERO PyTests"
    echo ""
    
    ###### AVL ######
    if [[ `command -v avl` ]]; then
        echo "avl: `which avl`"
        expectPythonSuccess "avl_PyTest.py"
        expectPythonSuccess "avl_AutoSpan_PyTest.py"
        expectPythonSuccess "avl_DesignSweep_PyTest.py" -noPlotData
        expectPythonSuccess "avl_ControlSurf_PyTest.py"
        expectPythonSuccess "avl_EigenValue_PyTest.py"
        if ( python -c 'import openmdao' ); then
            expectPythonSuccess "avl_OpenMDAO_3_PyTest.py"
        else
            notRun="$notRun\navl_OpenMDAO_3_PyTest.py"
        fi
    else
        notRun="$notRun\navl"
    fi
    
    ###### awave ######
    if [[ `command -v awave` ]]; then
        echo "awave: `which awave`"
        expectPythonSuccess "awave_PyTest.py"
    else
        notRun="$notRun\nawave"
    fi
    
    ###### friction ######
    if [[ `command -v friction` ]]; then
        echo "friction: `which friction`"
        expectPythonSuccess "friction_PyTest.py"
    else
        notRun="$notRun\nfriction"
    fi
    
    ###### tsFoil ######
    if [[ `command -v tsfoil2` && "$OS" != "Windows_NT" ]]; then
        echo "tsfoil2: `which tsfoil2`"
        expectPythonSuccess "tsfoil_PyTest.py"
    else
        notRun="$notRun\ntsfoil"
    fi
    
    ###### MSES ######
    if [[ `command -v mses` ]]; then
        echo "mses: `which mses`"
        expectPythonSuccess "mses_PyTest.py" -noPlotData
        if ( python -c 'import openmdao' ); then
          expectPythonSuccess "mses_OpenMDAO_3_PyTest.py"
        else
            notRun="$notRun\nmses_OpenMDAO_3_PyTest.py"
        fi
    else
        notRun="$notRun\nmses"
    fi
    
    ###### XFoil ######
    if [[ `command -v xfoil` ]]; then
        echo "xfoil: `which xfoil`"
        expectPythonSuccess "xfoil_PyTest.py" -noPlotData
    else
        notRun="$notRun\nxfoil"
    fi

    ###### AutoLink with awave, friction, and avl ######
    if [[ `command -v awave` && `command -v friction` && `command -v avl` ]]; then
        echo "awave: `which awave`"
        expectPythonSuccess "autoLink_Pytest.py" -noPlotData
    else
        notRun="$notRun\nautoLink"
    fi

    testsRan=1
fi

###### Meshing examples ######
if [[ "$TYPE" == "MESH" || "$TYPE" == "ALL" ]]; then
    echo "Running.... MESH PyTests"
    echo ""

    # EGADS Tess
    expectPythonSuccess "egadsTess_PyTest.py"
    expectPythonSuccess "egadsTess_Box_Quad_Pytest.py"
    expectPythonSuccess "egadsTess_Spheres_Quad_PyTest.py"
    expectPythonSuccess "egadsTess_Nose_Quad_PyTest.py"

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
        expectPythonSuccess "aflr4_and_aflr3_Symmetry_PyTest.py"
        expectPythonSuccess "aflr4_QuadMesh_PyTest.py"

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
        echo "delaundo: `which delaundo`"
        expectPythonSuccess "delaundo_PyTest.py"
    else
        notRun="$notRun\ndelaundo"
    fi

    # pointwise
    if [[ `command -v pointwise` && -f $ESP_ROOT/lib/pointwiseAIM.$EXT ]]; then
        echo "pointwise: `which pointwise`"
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
    echo ""
    
    # SU2
    if [ "$SU2_RUN" != "" ] && [ "$OS" != "Windows_NT" ]; then
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT ]]; then
            expectPythonSuccess "su2_and_AFLR2_NodeDist_PyTest.py"
            expectPythonSuccess "su2_and_AFLR2_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/aflr3AIM.$EXT && 
              -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
            expectPythonSuccess "su2_and_AFLR4_AFLR3_PyTest.py"
            expectPythonSuccess "su2_Morph_PyTest.py"
        fi
        if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
            echo "delaundo: `which delaundo`"
            expectPythonSuccess "su2_and_Delaundo_PyTest.py"
        fi
        if [[ `command -v pointwise` ]]; then
            echo "pointwise: `which pointwise`"
            expectPythonSuccess "su2_and_pointwise_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "su2_and_Tetgen_PyTest.py"
            expectPythonSuccess "su2_X43a_PyTest.py"
        fi
    else
        notRun="$notRun\nSU2"
    fi
   
    # Cart3D
    if [[ `command -v flowCart` ]]; then
        ulimit -s unlimited || true # Cart3D requires unlimited stack size
        echo "flowCart: `which flowCart`"
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

    # Fun3D
    if [[ `command -v nodet_mpi` ]]; then
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT &&
              -f $ESP_ROOT/lib/aflr3AIM.$EXT &&
              -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
            expectPythonSuccess "fun3d_and_AFLR2_NodeDist_PyTest.py"
            expectPythonSuccess "fun3d_and_AFLR2_PyTest.py"
            expectPythonSuccess "fun3d_and_AFLR4_AFLR3_PyTest.py"
            expectPythonSuccess "fun3d_refine_PyTest.py"
            expectPythonSuccess "fun3d_Morph_PyTest.py"
        fi
        if [[ `command -v delaundo` || `command -v delaundo.exe` ]]; then
            echo "delaundo: `which delaundo`"
            expectPythonSuccess "fun3d_and_Delaundo_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "fun3d_and_Tetgen_AlphaSweep_PyTest.py"
            expectPythonSuccess "fun3d_and_Tetgen_PyTest.py"
        fi
        expectPythonSuccess "fun3d_PyTest.py"
    else
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT &&
              -f $ESP_ROOT/lib/aflr3AIM.$EXT &&
              -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
            expectPythonSuccess "fun3d_Morph_PyTest.py" -noAnalysis
        fi
        notRun="$notRun\nFun3D"
    fi

    testsRan=1
fi

###### Structural examples ######
if [[ "$TYPE" == "STRUCTURE" || "$TYPE" == "ALL" ]]; then
    echo "Running.... STRUCTURE PyTests"
    echo ""

    ###### abaqus ######
    if [[ `command -v abaqus` ]]; then
         expectPythonSuccess "abaqus_SingleLoadCase_PyTest.py"
    else
        notRun="$notRun\nabaqus"
    fi

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
        echo "ASTROS_ROOT: $ASTROS_ROOT"
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
        echo "nastran: `which nastran`"
    else
       NASTRAN_NOANALYSIS=-noAnalysis
       notRun="$notRun\nNastran"
    fi
    # Nastran has saved output files to allow for -noAnalysis runs
    expectPythonSuccess "nastran_AGARD445_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_Aeroelastic_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_CompositeWingDesign_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_CompositeWingFreq_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_Composite_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_Flutter_15degree.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_MultiLoadCase_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_OptimizationMultiLoad_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_PlateMisc_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_SingleLoadCase_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_ThreeBarDesign_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_ThreeBarFreq_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_ThreeBarLinkDesign_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_ThreeBarMultiLoad_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_ThreeBar_PyTest.py" $NASTRAN_NOANALYSIS
    expectPythonSuccess "nastran_Trim_15degree.py" $NASTRAN_NOANALYSIS


    ###### Mystran ######
    if [[ `command -v mystran.exe` || `command -v mystran` ]]; then
        echo "mystran: `which mystran`"
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

    ###### Plato ######
    if [[ -f $ESP_ROOT/lib/exodusWriter.$EXT ]]; then
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT && \
              -f $ESP_ROOT/lib/aflr3AIM.$EXT && \
              -f $ESP_ROOT/lib/aflr4AIM.$EXT ]]; then
            expectPythonSuccess "plato_aflr_airfoil_PyTest.py"
            expectPythonSuccess "plato_aflr_table_PyTest.py"
            expectPythonSuccess "plato_aflr_table_morph_PyTest.py"
            expectPythonSuccess "plato_aflr_cyli_box_PyTest.py"
            expectPythonSuccess "plato_aflr_airfoil_morph_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "plato_tetgen_table_PyTest.py"
            expectPythonSuccess "plato_tetgen_cyli_box_PyTest.py"
        fi
        if [[ -f $ESP_ROOT/lib/aflr2AIM.$EXT && \
              -f $ESP_ROOT/lib/aflr3AIM.$EXT && \
              -f $ESP_ROOT/lib/aflr4AIM.$EXT && \
              -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
            expectPythonSuccess "plato_aflr_tetgen_table_PyTest.py"
            expectPythonSuccess "plato_aflr_tetgen_cyli_box_PyTest.py"
        fi
    else
        notRun="$notRun\nPlato"
    fi

    testsRan=1
fi

###### Aeroelastic examples ######
if [[ "$TYPE" == "AEROELASTIC" || "$TYPE" == "ALL" ]]; then
    echo "Running.... AEROELASTIC PyTests"
    echo ""

    ###### SU2 and Astros ######
    if [[ "$ASTROS_ROOT" != "" && "$SU2_RUN" != "" && "$OS" != "Windows_NT" ]]; then
        # SU2 6.0.0 on Windows is no longer supported
        expectPythonSuccess "aeroelasticSimple_Pressure_SU2_and_Astros.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Displacement_SU2_and_Astros.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Iterative_SU2_and_Astros.py"

        expectPythonSuccess "aeroelastic_Iterative_SU2_and_Astros.py"
    else
        notRun="$notRun\nAstros"
    fi

    ###### SU2 and Mystran ######
    if [[ (`command -v mystran.exe` || `command -v mystran`) && "$SU2_RUN" != "" && "$OS" != "Windows_NT" ]]; then
        echo "mystran: `which mystran`"
        # SU2 6.0.0 on Windows does not work with displacements
        expectPythonSuccess "aeroelasticSimple_Pressure_SU2_and_Mystran.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Displacement_SU2_and_Mystran.py" -noPlotData
        expectPythonSuccess "aeroelasticSimple_Iterative_SU2_and_Mystran.py"
    else
        notRun="$notRun\nMystran and SU2"
    fi

    ###### SU2 and TACS ######
    if [[ "$SU2_RUN" != "" && "$OS" != "Windows_NT" ]]; then
        expectPythonSuccess "aeroelasticSimple_Temperature_SU2_and_TACS.py" -noPlotData
    else
        notRun="$notRun\nTACS and SU2"
    fi

    ###### Cart3D and Mystran ######
    if [[ (`command -v mystran.exe` || `command -v mystran`) && `command -v flowCart` ]]; then
        echo "mystran: `which mystran`"
        echo "flowCart: `which flowCart`"
        expectPythonSuccess "aeroelasticSimple_Iterative_Cart3D_and_Mystran.py"
    else
        notRun="$notRun\nMystran and Cart3D"
    fi

    testsRan=1
fi

###### CorsairLite examples ######
if [[ "$TYPE" == "CORSAIRLITE" || "$TYPE" == "ALL" ]]; then
    echo "Running.... CORSAIRLITE PyTests"
    echo ""
    
    cd $regTestDir/../corsairlite
    expectPythonSuccess "gp.py"
    expectPythonSuccess "qp.py"
    expectPythonSuccess "sp.py"
    expectPythonSuccess "boydBox.py"
    expectPythonSuccess "hoburg_gp.py"
    if ( python -c 'import scipy' ); then
      expectPythonSuccess "hoburg_blackbox.py"
    fi

    cd $regTestDir/../corsairlite/capsPhase
    rm -rf ostrich
    expectPythonSuccess "pyCAPS/sizeWing.py"
    if [[ `command -v avl` ]]; then
        echo "avl: `which avl`"
        expectPythonSuccess "pyCAPS/optTaperAvl.py"
        expectPythonSuccess "pyCAPS/optWingAvl.py"
    fi

    cd $regTestDir/../corsairlite/multifidelity
    rm -rf hoburg
    expectPythonSuccess "pyCAPS/GPSize.py"
    if [[ `command -v mses` ]]; then
        echo "mses: `which mses`"
        expectPythonSuccess "pyCAPS/MSES.py"
        expectPythonSuccess "pyCAPS/Camber.py"
        expectPythonSuccess "pyCAPS/CMConstraint.py"
        if ( python -c 'import torch' ); then
            expectPythonSuccess "pyCAPS/Kulfan2.py"
            expectPythonSuccess "pyCAPS/FlowTrip.py"
            expectPythonSuccess "pyCAPS/Kulfan4.py"
        fi
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
