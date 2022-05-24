# Import other need modules
from __future__ import print_function

import unittest
import os
import glob
import shutil

# Import pyCAPS class file
import pyCAPS


# Geometry is verified with avl by plotting the camber using the commands:
#
# avl caps
# oper
# g
# ca

class TestAVL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # Create working directory variable
        cls.problemName = "workDir_avlAnalysisTest"

        cls.cleanUp()

        # Initialize Problem object
        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile="../csmData/avlSections.csm", outLevel=0)

    @classmethod
    def tearDownClass(cls):
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove problemName directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_numSpan(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanTotal"      : 24, # Can specify both Total and PerSection
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        wing = {"groupName"  : "Wing",
                "numChord"   : 8,
                "spaceChord" : 1.0,
                "numSpan"    : 24, # numSpan is depricated
                "spaceSpan"  : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_alpha_custom_increment(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        Alpha      = [0.0, 3.0, 9.0]
        CLtotTrue  = [ 0.23858,  0.41257,  0.74994]
        CDtotTrue  = [ 0.01285,  0.02492,  0.06719]
        CmtotTrue  = [-0.10501, -0.09945, -0.08685]

        for i in range(0,len(Alpha)):
            # Set custom AoA
            avl.input.Alpha = Alpha[i]

            # Retrieve results
            CLtot = avl.output.CLtot
            CDtot = avl.output.CDtot
            Cmtot = avl.output.Cmtot
            #print("Alpha = ", Alpha[i])
            #print("CLtot = ", CLtot)
            #print("CDtot = ", CDtot)
            #print("Cmtot = ", Cmtot)

            self.assertAlmostEqual(CLtotTrue[i], CLtot, 4)
            self.assertAlmostEqual(CDtotTrue[i], CDtot, 4)
            self.assertAlmostEqual(CmtotTrue[i], Cmtot, 4)

            # Check derivatives

            # Stability axis angles
            self.assertAlmostEqual(avl.output["CLa"].value, avl.output["CLtot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["CLb"].value, avl.output["CLtot"].deriv("Beta"), 8)

            self.assertAlmostEqual(avl.output["Cl'a"].value, avl.output["Cl'tot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["Cl'b"].value, avl.output["Cl'tot"].deriv("Beta"), 8)

            self.assertAlmostEqual(avl.output["Cma"].value, avl.output["Cmtot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["Cmb"].value, avl.output["Cmtot"].deriv("Beta"), 8)

            self.assertAlmostEqual(avl.output["Cn'a"].value, avl.output["Cn'tot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["Cn'b"].value, avl.output["Cn'tot"].deriv("Beta"), 8)

            # Body axis rotation rates
            self.assertAlmostEqual(avl.output["CXp"].value, avl.output["CXtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["CXq"].value, avl.output["CXtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["CXr"].value, avl.output["CXtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["CYp"].value, avl.output["CYtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["CYq"].value, avl.output["CYtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["CYr"].value, avl.output["CYtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["CZp"].value, avl.output["CZtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["CZq"].value, avl.output["CZtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["CZr"].value, avl.output["CZtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["Clp"].value, avl.output["Cltot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["Clq"].value, avl.output["Cltot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["Clr"].value, avl.output["Cltot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["Cmp"].value, avl.output["Cmtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["Cmq"].value, avl.output["Cmtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["Cmr"].value, avl.output["Cmtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["Cnp"].value, avl.output["Cntot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["Cnq"].value, avl.output["Cntot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["Cnr"].value, avl.output["Cntot"].deriv("YawRate"), 8)

            # Control surfaces in the stability axis
            ControlStability = avl.output["ControlStability"].value

            self.assertAlmostEqual(ControlStability['LeftAileron']["CLtot"] , avl.output["CLtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["CYtot"] , avl.output["CYtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["Cl'tot"], avl.output["Cl'tot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["Cmtot"] , avl.output["Cmtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["Cn'tot"], avl.output["Cn'tot"].deriv("LeftAileron"), 8)

            self.assertAlmostEqual(ControlStability['RightAileron']["CLtot"] , avl.output["CLtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["CYtot"] , avl.output["CYtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["Cl'tot"], avl.output["Cl'tot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["Cmtot"] , avl.output["Cmtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["Cn'tot"], avl.output["Cn'tot"].deriv("RightAileron"), 8)


            # Control surfaces in the body axis
            ControlBody = avl.output["ControlBody"].value

            self.assertAlmostEqual(ControlBody['LeftAileron']["CXtot"], avl.output["CXtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["CYtot"], avl.output["CYtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["CZtot"], avl.output["CZtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["Cltot"], avl.output["Cltot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["Cmtot"], avl.output["Cmtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["Cntot"], avl.output["Cntot"].deriv("LeftAileron"), 8)

            self.assertAlmostEqual(ControlBody['RightAileron']["CXtot"], avl.output["CXtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["CYtot"], avl.output["CYtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["CZtot"], avl.output["CZtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["Cltot"], avl.output["Cltot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["Cmtot"], avl.output["Cmtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["Cntot"], avl.output["Cntot"].deriv("RightAileron"), 8)


#==============================================================================
    def test_wing_tail(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        htail = {"groupName"         : "hTail",
                 "numChord"          : 8,
                 "spaceChord"        : 0.7,
                 "numSpanTotal"      : 10}

        vtail = {"numChord"          : 8,
                 "spaceChord"        : 1.0,
                 "numSpanTotal"      : 8,
                 "spaceSpan"         : 0.5}

        avl.input.AVL_Surface = {"Wing" : wing,
                                 "hTail": htail,
                                 "vTail": vtail}

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.15827
        CDtotTrue  = 0.02262
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)

#==============================================================================
    def test_wing_Vtail(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        vtail = {"numChord"          : 8,
                 "spaceChord"        : 1.0,
                 "numSpanTotal"      : 8,
                 "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing"  : wing,
                                 "VTail1": vtail,
                                 "VTail2": vtail}

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.14681
        CDtotTrue  = 0.02022
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)

#==============================================================================
    def test_geom_change(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.29683
        CDtotTrue  = 0.01617
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)

        # change the geometry and rerun
        camber = self.myProblem.geometry.despmtr.camber
        self.myProblem.geometry.despmtr.camber = 1.1*camber

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.32434
        CDtotTrue  = 0.01801
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)

        # reset the geometry
        self.myProblem.geometry.despmtr.camber = camber

#==============================================================================
    def test_MassProp_noUnits(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        mass = 0.1773
        x =  0.02463
        y = 0.
        z = 0.2239
        Ixx = 1.350
        Iyy = 0.7509
        Izz = 2.095

        avl.input.MassProp = {"Aircraft": {"mass":mass, "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz]}}

        # check there are errors if information is missing
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.Gravity  = 9.8         # m/s
        avl.input.Density  = 1.22557083  # kg/m^3
        avl.input.Velocity = 19.67167008 # "m/s"

        # make sure there are no errsos
        avl.runAnalysis()

        avl.input.MassProp = {"Aircraft":{"mass":mass, "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz, 1.0, 2.0, 3.0]},
                              "Engine"  :{"mass":mass, "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz]}}

        # again should not cause errors
        avl.runAnalysis()

        # test error handling of the mass properties parsing

        avl.input.MassProp = {"Aircraft": "1"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": "foo"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":mass, "CG":[x,y,z]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

        avl.input.MassProp = {"Aircraft": {"mass":"mass", "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")


#==============================================================================
    def test_MassProp_Units(self):

        m    = pyCAPS.Unit("meter")
        kg   = pyCAPS.Unit("kg")
        s    = pyCAPS.Unit("s")
        K    = pyCAPS.Unit("Kelvin")
        deg  = pyCAPS.Unit("degree")
        ft   = pyCAPS.Unit("ft")
        slug = pyCAPS.Unit("slug")


        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM",
                                             unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0 * deg
        avl.input.Beta  = 0.0 * deg

        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        mass = 0.1773
        x =  0.02463
        y = 0.
        z = 0.2239
        Ixx = 1.350
        Iyy = 0.7509
        Izz = 2.095

        avl.input.MassProp = {"Aircraft": {"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz] * kg*m**2}}

        # check there are errors if information is missing
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.Gravity  = 32.18 * ft/s**2
        avl.input.Density  = 0.002378 * slug/ft**3
        avl.input.Velocity = 64.5396 * ft/s

        # make sure there are no errsos
        avl.runAnalysis()

        avl.input.MassProp = {"Aircraft":{"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz, 1.0, 2.0, 3.0] * kg*m**2},
                              "Engine"  :{"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz] * kg*m**2}}

        # again should not cause errors
        avl.runAnalysis()

        # test error handling of the mass properties parsing

        avl.input.MassProp = {"Aircraft": "1"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": "foo"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":[mass,"kg"], "CG":[[x,y,z],"m"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

        avl.input.MassProp = {"Aircraft": {"mass":(mass), "CG":[[x,y,z],"m"], "massInertia":[[Ixx, Iyy, Izz], "kg*m^2"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":[mass,"kg"], "CG":[x,y,z], "massInertia":[[Ixx, Iyy, Izz], "kg*m^2"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":[mass,"kg"], "CG":[[x,y,z],"m"], "massInertia":[Ixx, Iyy, Izz]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":("mass","kg"), "CG":[[x,y,z],"m"], "massInertia":[[Ixx, Iyy, Izz], "kg*m^2"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")
        
#==============================================================================
    def test_phase(self):
        
        # Initialize Problem object
        problemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile="../csmData/avlSections.csm", outLevel=0)

        # Load avl aim
        avl = myProblem.analysis.create(aim = "avlAIM", name = "avl")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}
        
        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", prevPhaseName="Phase0", outLevel=0)

        avl = myProblem.analysis["avl"]

        # Retrieve results
        self.assertAlmostEqual(CLtot, avl.output.CLtot)
        self.assertAlmostEqual(CDtot, avl.output.CDtot)

        avl.input.Alpha  = 3.0
        self.assertAlmostEqual(0.41257, avl.output.CLtot, 4)


if __name__ == '__main__':
    unittest.main()
