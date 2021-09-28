                 CAPS: Computational Aircraft Prototype Syntheses
                              Rev 1.15 -- May 2019
                                   Examples
                                

1. Prerequisites

    All the examples assume that the appropriate executables can be invoked 
    via a system call (i.e. they are in the PATH environment variable). The 
    one exception is the Astros examples, which use use the environment 
    variable ASTROS_ROOT (set to the path where the Astros executable can 
    be found). The analysis tools xfoil and avl are provided with the 
    distribution, and all other software must be obtained elsewhere.
    
    The examples are provided in the cCAPS and pyCAPS sub-directories.
    File names indicate the analysis tools executed within the example.
    
    The bash scripts regressionTest/execute_CTestRegression.sh and
    regressionTest/execute_PyTestRegression.sh can be used to test the 
    CAPS installation on OSX and Linux. Only examples where the analysis 
    software is found in the path are executed. These scripts must be 
    executed from within the regressionTest directory.


2.  C-code examples

    2.1 Low Fidelity CFD
        avlTest                          - Simple execution of avl
 
    2.2 Mesh Generation
         
    2.3 CFD and Mesh Generation
        fun3dTest                        - Execute fun3d without mesh generation
        fun3dAFLR2Test                   - 2D mesh generation and execute Fun3d
        fun3dTetegenTest                 - 3D mesh generation and execute Fun3d

    2.4 Structural Analysis
        mystranTest                      - Simple execution of mystran

    2.5 Aeroelastic data transfer
        aeroelasticSimple_Iterative_SU2_and_MystranTest - Iterative aeroelastic


3.  pyCAPS examples

    3.1 Low Fidelity CFD
        avl_*                            - Execution of avl
        xfoil_PyTest                     - Execution of xfoil
        autoLink_PyTest                  - Example of linking outputs and inputs
       
    3.2 Mesh Generation
        aflr2_PyTest                     - 2D mesh generation
        delaundo_PyTest                  - 2D mesh generation
        aflr4_PyTest                     - Surface mesh generation
        quading_PyTest                   - Surface quad mesh generation
        aflr3_PyTest                     - 3D mesh generation
        aflr4_and_*_PyTest               - Surface and 3D mesh generation

    3.3 CFD and Mesh Generation
        su2_and_AFLR2_*_PyTest           - 2D mesh generation and execute SU2
        su2_and_Delaundo_PyTest          - 2D mesh generation and execute SU2
        su2_and_*AFLR3_PyTest            - 3D mesh generation and execute SU2
        su2_and_Tetgen_PyTest            - 3D mesh generation and execute SU2
        su2_X43a_PyTest                  - 3D mesh generation and execute SU2

        fun3d_PyTest                     - Execute Fun3D
        fun3d_and_AFLR2_*_PyTest         - 2D mesh generation and execute Fun3D
        fun3d_and_Delaundo_PyTest        - 2D mesh generation and execute Fun3D
        fun3d_and_egadsTess_*_PyTest     - 2D mesh generation and execute Fun3D
        fun3d_and_AFLR4_AFLR3_PyTest     - Surface, 3D mesh generation and execute Fun3D
        fun3d_and_Tetgen_*_PyTest        - 3D mesh generation and execute Fun3D

    3.4 Structural Analysis
        astros_ThreeBar*_PyTest          - Beam element analysis
        astros_PyTest                    - Eigen analysis
        astros_Composite_PyTest          - Eigen analysis
        astros_CompositeWingFreq_PyTest  - Eigen analysis
        astros_SingleLoadCase_PyTest     - Static analysis
        astros_MultipleLoadCase_PyTest   - Static analysis
        astros_AGARD445_PyTest           - Modal analysis
        astros_Aeroelastic_PyTest        - Astros internal aeroelastic analysis
        astros_Flutter_15degree          - Astros internal aeroelastic analysis

        nastran_ThreeBar*_PyTest         - Beam element analysis
        nastran_PyTest                   - Eigen analysis
        nastran_Composite_PyTest         - Eigen analysis
        nastran_CompositeWingFreq_PyTest - Eigen analysis
        nastran_SingleLoadCase_PyTest    - Static analysis
        nastran_MultipleLoadCase_PyTest  - Static analysis
        nastran_AGARD445_PyTest          - Modal analysis
        nastran_Aeroelastic_PyTest       - Nastran internal aeroelastic analysis
        nastran_Flutter_15degree         - Nastran internal aeroelastic analysis

        mystran_PyTest                   - Eigen analysis
        mystran_AGARD445_PyTest          - Modal analysis
        mystran_SingleLoadCase_PyTest    - Static analysis
        mystran_MultipleLoadCase_PyTest  - Static analysis

    3.5 Aeroelastic data transfer
        aeroelasticSimple_Pressure_*     - Transfer pressure load from CFD to FEA
        aeroelasticSimple_Displacement_* - Transfer displacements from FEA to CFD
        aeroelasticSimple_Iterative_*    - Iterative data transfer between CFD and FEA
                                         
        aeroelastic_Iterative_*          - Iterative data transfer between CFD and FEA

