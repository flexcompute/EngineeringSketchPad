# Import pyCAPS and os module
import pyCAPS
import os

# Import argparse
import argparse
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'autoLink Pytest Link Example',
                                 prog = 'autoLink_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot data")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Define units
m    = pyCAPS.Unit("meter")
ft   = pyCAPS.Unit("ft")
kg   = pyCAPS.Unit("kg")
s    = pyCAPS.Unit("s")
K    = pyCAPS.Unit("Kelvin")
deg  = pyCAPS.Unit("degree")


def combineDrag(myProblem):

    drags = ["CDtot", "CDwave", "CDtotal"] # CDtotal = Cdform+Cdfric
    CDTotal = 0
    for myAnalysis in myProblem.analysis.values():

        # Skip analysis for AVL for Mach numbers > .75 - We wont have any analysis
        if myProblem.parameter["Mach"].value > 0.75 and myAnalysis.name == "avl":
            continue
        # Skip analysis for AWave for Mach < 1.0 - We wont have any analysis
        if myProblem.parameter["Mach"].value < 1.0 and myAnalysis.name == "awave":
            continue

        for dragType in drags:
            if dragType in myAnalysis.output:
                CDTotal += myAnalysis.output[dragType].value

    return CDTotal

# Create working directory variable
problemName = str(args.workDir[0]) + "LinkTest"

geometryScript = os.path.join("..","csmData","linearAero.csm")

# Initialize Problem object
myProblem = pyCAPS.Problem(problemName, capsFile=geometryScript, outLevel=0)

# Load desired aim
fric = myProblem.analysis.create(aim = "frictionAIM",
                                 name = "friction")

awave = myProblem.analysis.create(aim = "awaveAIM",
                                  name = "awave",
                                  unitSystem = {"angle": deg})

avl = myProblem.analysis.create(aim = "avlAIM",
                                name = "avl",
                                unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

# Create mission parameters
mach     = myProblem.parameter.create("Mach", 0.1) # Mach number
alpha    = myProblem.parameter.create("Alpha", 0.0 * deg) # Angle of attack
altitude = myProblem.parameter.create("Altitude", 5000.0 * ft) # Flying altitude

machRange     = [mach.value + i*0.1 for i in range(7)]
altitudeRange = [altitude.value + i*5000 * ft for i in range(12)]

# Autolink the mission parameters to inputs of the loaded AIMs
myProblem.autoLinkParameter()

# Set geometry inputs
# despmtr   thick     0.12      frac of local chord
# despmtr   camber    0.04      frac of local chord
# despmtr   tlen      5.00      length from wing LE to Tail LE
# despmtr   toff      0.5       tail offset
#
# despmtr   area      10.0
# despmtr   aspect    6.00
# despmtr   taper     0.60
# despmtr   sweep     20.0      deg (of c/4)
#
# despmtr   washout  -5.00      deg (down at tip)
# despmtr   dihedral  4.00      deg
myProblem.geometry.despmtr.area = 7.5

# Set additional information for AVL
avlSurfaces = {}
for i in avl.geometry.attrList("capsGroup"):

     if "tail" in i.lower():
         surf = {"numChord"     : 8,
                 "spaceChord"   : 1.0,
                 "numSpanTotal" : 6,
                 "spaceSpan"    : 1.0}
     elif "wing" in i.lower():
         surf = {"numChord"     : 10,
                 "spaceChord"   : 0.5,
                 "numSpanTotal" : 40,
                 "spaceSpan"    : 1.0}
     else:
         continue # Not interested in other capsGroups

     avlSurfaces[i] = surf

# Add surface information to AVL
avl.input.AVL_Surface = avlSurfaces

#Create a list of lists for storage
Cd = [] # altitude index

# Loop through ranges
for alt in altitudeRange:

    Cd.append([]) # mach index

    for ma in machRange:

        # Set mission parameters - AIMs will be updated automatically due to link
        mach.value = ma
        altitude.value = alt

        # Append the combined drag (analysis is executed automatically just-in-time)
        Cd[-1].append(combineDrag(myProblem))

# Print out drag values
print("Cd = ", Cd)

if (args.noPlotData == False):
    # Create a contour plot in matplotlib and numpy is available
    try:
        import matplotlib.pyplot as plt

        import numpy

        X = numpy.arange(min(machRange), max(machRange)+0.0001, machRange[1] - machRange[0])
        Y = numpy.arange(min(altitudeRange)/ft, max(altitudeRange)/ft+0.0001, (altitudeRange[1] - altitudeRange[0])/ft)
        Z = numpy.array(Cd)

        CS = plt.contourf(X, Y, Z, 30)

        plt.title('Coefficient of Drag')
        plt.xlabel('Mach number')
        plt.ylabel('Altitude (ft)')

        # Make a colorbar for the ContourSet returned by the contourf call.
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('Cd')

        plt.show()

    except:
        print("Unable to plot data")
