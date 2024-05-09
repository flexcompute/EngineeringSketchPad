import unittest

import os
import glob
import shutil
import __main__

from sys import version_info as pyVersion

import pyCAPS

class TestJournal(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "journalTest"
        cls.iProb = 1

        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove project directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def run_xfoil_AutoExec(self, myProblem, line_exit):

        line = 0
        if line < line_exit and line_exit > 0: self.assertTrue(myProblem.journaling())
        if line == line_exit: return line

        # Set the limits on a design parameter
        myProblem.geometry.despmtr["camber"].limits = [0., 0.5]; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Get the limits on a design parameter
        limits = myProblem.geometry.despmtr["camber"].limits; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Remove the limits on a design parameter
        myProblem.geometry.despmtr["camber"].limits = None; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Change a design parameter - camber in the geometry
        myProblem.geometry.despmtr.camber = 0.1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load xfoil aim
        xfoil = myProblem.analysis.create(aim = "xfoilAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set Mach number
        xfoil.input.Mach = 0.5; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set Reynolds number
        xfoil.input.Re   = 1.0e6; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set custom AoA
        xfoil.input.Alpha = 5.0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Append the polar file if it already exists - otherwise the AIM will delete the file
        xfoil.input.Append_PolarFile = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Retrieve Cl and Cd
        Cl = xfoil.output.CL; line += 1
        self.assertAlmostEqual(Cl, 1.5658, 3)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        Cd = xfoil.output.CD; line += 1
        self.assertAlmostEqual(Cd, 0.02005, 4)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Transition location
        TranX = xfoil.output.Transition_Top; line += 1
        self.assertAlmostEqual(TranX, 0.4424, 3)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # make sure the last call journals everything
        return line+2


#==============================================================================
    # Try journaling with xfoil and auto execution
    def test_journal_xfoil_AutoExec(self):

        file = "airfoilSection.csm"
        myProblem = pyCAPS.Problem(self.problemName, capsFile=file, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_xfoil_AutoExec(myProblem, -1)

        myProblem.close()
        shutil.rmtree(self.problemName)

        #print(80*"=")
        #print(80*"=")
        myProblem = pyCAPS.Problem(self.problemName, phaseName="phase0", capsFile=file, outLevel=0)
        myProblem.close()

        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(self.problemName, phaseName="phase0", outLevel=0)
            self.run_xfoil_AutoExec(myProblem, line_exit)
            myProblem.close()

#==============================================================================
    def run_xfoil_ManualExec(self, myProblem, line_exit):

        line = 0
        if line == line_exit: return line

        # Change a design parameter - camber in the geometry
        myProblem.geometry.despmtr.camber = 0.1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load xfoil aim
        xfoil = myProblem.analysis.create(aim = "xfoilAIM", autoExec = False); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set Mach number
        xfoil.input.Mach = 0.5; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set Reynolds number
        xfoil.input.Re   = 1.0e6; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set custom AoA
        xfoil.input.Alpha = 5.0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Append the polar file if it already exists - otherwise the AIM will delete the file
        xfoil.input.Append_PolarFile = False; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run AIM pre-analysis
        xfoil.preAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        ####### Run xfoil ####################
        #print ("\n\nRunning xFoil......")
        xfoil.system("xfoil < xfoilInput.txt > xfoilOutput.txt"); line += 1 # Run xfoil via system call
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run AIM post-analysis
        xfoil.postAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Retrieve Cl and Cd
        Cl = xfoil.output.CL; line += 1
        #print("Cl = ", Cl)
        self.assertAlmostEqual(Cl, 1.5658, 3)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        Cd = xfoil.output.CD; line += 1
        #print("Cd = ", Cd)
        self.assertAlmostEqual(Cd, 0.02005, 4)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Transition location
        TranX = xfoil.output.Transition_Top; line += 1
        #print("TranX = ", TranX)
        self.assertAlmostEqual(TranX, 0.4424, 3)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # make sure the last call journals everything
        return line+2

#==============================================================================
    # Try journaling with xfoil and auto execution
    def test_journal_xfoil_ManualExec(self):

        file = "airfoilSection.csm"
        myProblem = pyCAPS.Problem(self.problemName, capsFile=file, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_xfoil_ManualExec(myProblem, -1)

        myProblem.close()
        shutil.rmtree(self.problemName)

        #print(80*"=")
        #print(80*"=")
        myProblem = pyCAPS.Problem(self.problemName, phaseName="phase0", capsFile=file, outLevel=0)
        myProblem.close()

        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(self.problemName, phaseName="phase0", outLevel=0)
            self.run_xfoil_ManualExec(myProblem, line_exit)
            myProblem.close()


#==============================================================================
    def run_geometry(self, myProblem, line_exit):

        myGeometry = myProblem.geometry

        line = 0
        if line == line_exit: return line

        #---------
        # Change a design parameter
        myGeometry.despmtr.wing.dihedral = 8; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Check the design parameter
        self.assertEqual( 8, myGeometry.despmtr.wing.dihedral); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        #---------

        #---------
        # compute some outpmtr values consistent with equations in unitGeom.csm
        aspect = myGeometry.despmtr.aspect; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        area   = myGeometry.despmtr.area; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        cmean  = (area/aspect)**0.5
        span   = cmean*aspect

        # Access outpmtr class
        self.assertAlmostEqual(span, myGeometry.outpmtr.span, 5); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        #---------

        #---------
        # Check despmtr matrix access
        self.assertEqual(myGeometry.cfgpmtr.nrow, 3); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.cfgpmtr.ncol, 2); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.despmtr.despMat, [[11.0, 12.0], [13.0, 14.0], [15.0, 16.0]]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        myGeometry.despmtr.despMat = [[11.0, 12.0], [13.0, 14.0], [30.0, 32.0]]; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.despmtr.despMat, [[11.0, 12.0], [13.0, 14.0], [30.0, 32.0]]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        #---------

        #---------
        myGeometry.cfgpmtr.nrow = 5; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.cfgpmtr.nrow, 5); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        myGeometry.despmtr.despMat = [[11.,12.],
                                      [13.,14.],
                                      [15.,16.],
                                      [17.,18.],
                                      [19.,20.]]; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.despmtr.despMat, [[11.,12.],
                                                      [13.,14.],
                                                      [15.,16.],
                                                      [17.,18.],
                                                      [19.,20.]]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        myGeometry.despmtr.despCol = [11.,12.,13.,14.,15.]; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.despmtr.despCol, [11.,12.,13.,14.,15.]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        myGeometry.despmtr.despRow = [11.,12.,13.,14.,15.]; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertEqual(myGeometry.despmtr.despRow, [11.,12.,13.,14.,15.]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        #---------

        #---------
        # Test with generating and catching an error
        try:
            value = myProblem.parameter["foo"].value
        except KeyError:
            value = myProblem.parameter.create("foo", 4).value
        self.assertEqual(value, 4); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # make sure the last call journals everything
        return line+2

#==============================================================================
    # Try journaling with geometry modifications
    def test_journal_Geometry(self):

        file = "unitGeom.csm"
        myProblem = pyCAPS.Problem(self.problemName, capsFile=file, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_geometry(myProblem, -1)

        myProblem.close()
        shutil.rmtree(self.problemName)

        #print(80*"=")
        #print(80*"=")
        myProblem = pyCAPS.Problem(self.problemName, phaseName="Phase0", capsFile=file, outLevel=0)
        myProblem.close()

        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(self.problemName, phaseName="Phase0", outLevel=0)
            self.run_geometry(myProblem, line_exit)
            myProblem.close()

if __name__ == '__main__':
    unittest.main()
