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

        cls.file = "airfoilSection.csm"
        cls.probemName = "journalTest"
        cls.iProb = 1

        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove project directories
        dirs = glob.glob( cls.probemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)
                
#==============================================================================
    def run_xfoil_AutoExec(self, myProblem, line_exit):

        line = 0
        if line < line_exit and line_exit > 0: self.assertTrue(myProblem.journaling())
        if line == line_exit: return line

        # Change a design parameter - area in the geometry
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
        self.assertAlmostEqual(Cl, 1.5266, 3)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        
        Cd = xfoil.output.CD; line += 1
        self.assertAlmostEqual(Cd, 0.01986, 4)
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

        myProblem = pyCAPS.Problem(self.probemName, capsFile=self.file, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_xfoil_AutoExec(myProblem, -1)
        
        myProblem.close()
        shutil.rmtree(self.probemName)
        
        #print(80*"=")
        #print(80*"=")
        myProblem = pyCAPS.Problem(self.probemName, capsFile=self.file, outLevel=0)
        myProblem.close()
        
        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(self.probemName, capsFile=self.file, outLevel=0, useJournal=True)
            self.run_xfoil_AutoExec(myProblem, line_exit)
            myProblem.close()

#==============================================================================
    def run_xfoil_ManualExec(self, myProblem, line_exit):

        line = 0
        if line == line_exit: return line

        # Change a design parameter - area in the geometry
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
        self.assertAlmostEqual(Cl, 1.5266, 3)
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        
        Cd = xfoil.output.CD; line += 1
        #print("Cd = ", Cd)
        self.assertAlmostEqual(Cd, 0.01986, 4)
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

        myProblem = pyCAPS.Problem(self.probemName, capsFile=self.file, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_xfoil_ManualExec(myProblem, -1)
        
        myProblem.close()
        shutil.rmtree(self.probemName)
        
        #print(80*"=")
        #print(80*"=")
        myProblem = pyCAPS.Problem(self.probemName, capsFile=self.file, outLevel=0)
        myProblem.close()
        
        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(self.probemName, capsFile=self.file, outLevel=0, useJournal=True)
            self.run_xfoil_ManualExec(myProblem, line_exit)
            myProblem.close()

if __name__ == '__main__':
    unittest.main()
