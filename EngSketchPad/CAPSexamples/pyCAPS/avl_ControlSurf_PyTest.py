# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse
import argparse
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AVL Contral Surface Pytest Example',
                                 prog = 'avl_ControlSurf_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AVLControlAnalsyisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","avlWingTail.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.area = 10.0

# Load desired aim
print ("Loading AIM")
myAnalysis = myProblem.analysis.create(aim = "avlAIM")

print ("Setting design parameters")

# Set new Mach/Alt parameters
myAnalysis.input.Mach      = 0.5
myAnalysis.input.Alpha     = 1.0
myAnalysis.input.Beta      = 0.0
myAnalysis.input.RollRate  = 0.0
myAnalysis.input.PitchRate = 0.0
myAnalysis.input.YawRate   = 0.0

wing = {"groupName"    : "Wing",
        "numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 24,
        "spaceSpan"    : 1.0}

tail = {"numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 12,
        "spaceSpan"    : 1.0}

myAnalysis.input.AVL_Surface = {"Wing": wing, "Vertical_Tail": tail}


# Set control surface parameters
flapLE = {"controlGain"    : 0.5,
          "deflectionAngle" : 25.0}

flapTE = {"controlGain"    : 1.0,
          "deflectionAngle" : 15.0}


myAnalysis.input.AVL_Control = {"WingRightLE": flapLE,
                                "WingRightTE": flapTE,
                               #"WingLeftLE": flapLE, # Don't change anything
                               #"WingLeftTE": flapTE
                                "Tail"      : flapTE}


# Get Output Data from AVL
print ("Alpha  " + str(myAnalysis.output["Alpha" ].value))
print ("Beta   " + str(myAnalysis.output["Beta"  ].value))
print ("Mach   " + str(myAnalysis.output["Mach"  ].value))
print ("pb/2V  " + str(myAnalysis.output["pb/2V" ].value))
print ("qc/2V  " + str(myAnalysis.output["qc/2V" ].value))
print ("rb/2V  " + str(myAnalysis.output["rb/2V" ].value))
print ("p'b/2V " + str(myAnalysis.output["p'b/2V"].value))
print ("r'b/2V " + str(myAnalysis.output["r'b/2V"].value))
print ("CXtot  " + str(myAnalysis.output["CXtot" ].value))
print ("CYtot  " + str(myAnalysis.output["CYtot" ].value))
print ("CZtot  " + str(myAnalysis.output["CZtot" ].value))
print ("Cltot  " + str(myAnalysis.output["Cltot" ].value))
print ("Cmtot  " + str(myAnalysis.output["Cmtot" ].value))
print ("Cntot  " + str(myAnalysis.output["Cntot" ].value))
print ("Cl'tot " + str(myAnalysis.output["Cl'tot"].value))
print ("Cn'tot " + str(myAnalysis.output["Cn'tot"].value))
print ("CLtot  " + str(myAnalysis.output["CLtot" ].value))
print ("CDtot  " + str(myAnalysis.output["CDtot" ].value))
print ("CDvis  " + str(myAnalysis.output["CDvis" ].value))
print ("CLff   " + str(myAnalysis.output["CLff"  ].value))
print ("CYff   " + str(myAnalysis.output["CYff"  ].value))
print ("CDind  " + str(myAnalysis.output["CDind" ].value))
print ("CDff   " + str(myAnalysis.output["CDff"  ].value))
print ("e      " + str(myAnalysis.output["e"     ].value))

Cl = myAnalysis.output["CLtot"].value
Cd = myAnalysis.output["CDtot"].value

# Check assertation
assert abs(Cl-0.61999) <= 1E-3 and abs(Cd-0.02984) <= 1E-3
