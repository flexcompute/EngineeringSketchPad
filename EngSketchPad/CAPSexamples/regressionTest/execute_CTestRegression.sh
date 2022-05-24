#!/bin/bash -e

stat=0
notExpectedTests=""
notRun=""
LOGFILENAME=cCAPSlog.txt

expectCSuccess ()
{
    # $1 = program
    # $2 = outLevel
     
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    local status
    set -x
    echo "$1 test;" | tee -a $LOGFILE
    time $VALGRIND_COMMAND $1 $2 | tee -a $LOGFILE; status=${PIPESTATUS[0]}
    set +x
    if [[ $status == 0 ]]; then
    echo | tee -a $LOGFILE
    echo "=============================================" | tee -a $LOGFILE
    echo "CAPS verification (using c-executables) case $1 passed (as expected)" | tee -a $LOGFILE
    echo "=============================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    else
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo "CAPS verification (using c-executables) case $1 failed (as NOT expected)" | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    stat=1
    notExpectedTests="$notExpectedTests\n$1"
    fi
}

# function to run a case that is expected to fail
expectCFailure ()
{
    # $1 = program
    # $2 = outLevel
    
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    local status
    echo "$1 test;" | tee -a $LOGFILE
    set -x
    time $VALGRIND_COMMAND $1 $2 | tee -a $LOGFILE; status=${PIPESTATUS[0]}
    set +x
    if [[ $status == 0 ]]; then
    echo | tee -a $LOGFILE
    echo "=============================================" | tee -a $LOGFILE
    echo "CAPS verification (using c-executables) case $1 passed (as NOT expected)" | tee -a $LOGFILE
    echo "=============================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    stat=1
    notExpectedTests="$notExpectedTests\n$1"
    else
    echo | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo "CAPS verification (using c-executables) case $1 failed (as expected)" | tee -a $LOGFILE
    echo "=================================================" | tee -a $LOGFILE
    echo | tee -a $LOGFILE
    fi
}

# Get the current directory
curDir=${PWD##*/} 
regTestDir=${PWD}


if [[ "$curDir" != "regressionTest" || ! -d "../csmData" || ! -d "../cCAPS" ]]; then
  echo "execute_CTestRegression.sh must be executed in the regressionTest directory!"
  echo "Current directory: $regTestDir"
  exit 1
fi

# make the log file abosolute path to the current directory
LOGFILE=$regTestDir/$LOGFILENAME

# Remove old log file
rm -f $LOGFILE
touch $LOGFILE

# Move into the cCAPS directory
cd ../cCAPS

######### Execute c-Test examples ########

TYPE="ALL"
testsRan=0

if [[ "$1" != "" ]]; then
    TYPE="$1"
else
    echo "No options provided...defaulting to ALL"
fi

if [[ "$OS" == "Windows_NT" ]]; then
    EXT=dll
else
    EXT=so
fi

###################################################
if [[ "$TYPE" == "MINIMAL" || "$TYPE" == "ALL" ]]; then
    echo "Running.... MINIMAL c-Tests"
    testsRan=1
fi

###### Linear aero. examples ######
if [[ "$TYPE" == "LINEARAERO" || "$TYPE" == "ALL" ]]; then
    echo "Running.... LINEARAERO c-Tests"
    
    ###### AVL ###### 
    if [[ "`which avl`" != "" ]]; then
        expectCSuccess "./avlTest" 0
    else
        notRun="$notRun\navl"
    fi

    ###### MSES ###### 
    if [[ "`which mses`" != "" ]]; then
        expectCSuccess "./msesTest" 0
    else
        notRun="$notRun\nmses"
    fi

    ###### AWAVE ###### 
    if [[ "`which awave`" != "" ]]; then
        expectCSuccess "./awaveTest" 0
    else
        notRun="$notRun\nawave"
    fi

    ######  FRICTION ###### 
    if [[ "`which friction`" != "" ]]; then
        expectCSuccess "./frictionTest" 0
    else
        notRun="$notRun\nfriction"
    fi

    testsRan=1
fi

###### Meshing examples ######
if [[ "$TYPE" == "MESH" || "$TYPE" == "ALL" ]]; then
    echo "Running.... MESH c-Tests"

   # if [[ "$AFLR" != "" ]]; then
   #     
   #     if [[ "$TETGEN" != "" ]]; then
   #     fi
   # fi
    
   # if [[ "$TETGEN" != "" ]]; then
   # fi

    if [[ "$OS" != "Windows_NT" ]]; then
        POINTWISE_AIM=pointwiseAIM.so
    fi
    
    # pointwise
    if [[ "`which pointwise`" != "" && "$VALGRIND_COMMAND" == "" && -f $ESP_ROOT/lib/$POINTWISE_AIM ]]; then
        expectCSuccess "./pointwiseTest" 0
    else
        notRun="$notRun\npointwise"
    fi

    testsRan=1
fi

###### CFD examples ######
if [[ "$TYPE" == "CFD" || "$TYPE" == "ALL" ]]; then
    echo "Running.... CFD c-Tests"
    
    ###### FUN3D ######

    if [[  -f $ESP_ROOT/lib/aflr2AIM.$EXT ]]; then
        expectCSuccess "./fun3dAFLR2Test" 0
    else
        notRun="$notRun\nfun3d AFLR"
    fi

    if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
        expectCSuccess "./fun3dTetgenTest" 0
    else
        notRun="$notRun\nfun3d TetGen"
    fi

    testsRan=1
fi

###### Structural examples ######
if [[ "$TYPE" == "STRUCTURE" || "$TYPE" == "ALL" ]]; then
    echo "Running.... STRUCTURE c-Tests"
    
    ######  Mystran ###### 
    if [[ "`which mystran.exe`" != "" && "$OS" != "Windows_NT" ]]; then
        expectCSuccess "./mystranTest" 0
    else
        notRun="$notRun\nMystran"
    fi 

    ######  Interference ###### 
    expectCSuccess "./interferenceTest"
    
    ######  HSM ###### 
    if [[ -f $ESP_ROOT/lib/hsmAIM.$EXT  ]]; then
        expectCSuccess "./hsmCantileverPlateTest" 0
        expectCSuccess "./hsmSimplePlateTest" 0
        expectCSuccess "./hsmJoinedPlateTest" 0
    else
        notRun="$notRun\nHSM"
    fi 

    testsRan=1
fi

###### Aeroelastic examples ######
if [[ "$TYPE" == "AEROELASTIC" || "$TYPE" == "ALL" ]]; then
    echo "Running.... AEROELASTIC c-Tests"
    
    if [[ -f $ESP_ROOT/lib/tetgenAIM.$EXT ]]; then
    
        if [[ "`which mystran.exe`" != "" && "`which nodet_mpi`" != "" ]]; then
           expectCSuccess "./aeroelasticTest" 0
        fi

            # SU2 6.0.0 on Windows does not work with displacements
        if [[ "`which mystran.exe`" != "" && "$OS" != "Windows_NT" ]]; then
           expectCSuccess "./aeroelasticSimple_Iterative_SU2_and_MystranTest" 0
        fi
        
    else
        notRun="$notRun\naeroelastic"
    fi

    testsRan=1
fi


# Cleanup if we have success
if [[ $stat == 0 ]]; then
    rm -fr $cRegDir
fi

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
