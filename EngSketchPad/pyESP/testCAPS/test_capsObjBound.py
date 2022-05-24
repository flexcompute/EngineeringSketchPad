import unittest

import os
import glob
import shutil
import __main__

from sys import version_info as pyVersion
from sys import version_info

from pyCAPS import caps

class TestBound(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.analysisDir = "UnitTest"
        cls.projectName = "testProblem"
        cls.iDir = 1
        cls.iProb = 1
        cls.cleanUp()

        cls.boundName = "Upper_Left"

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.projectName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove default projectName
        base = os.path.basename(__main__.__file__)
        projectName = os.path.splitext(base)[0]
        if os.path.isdir(projectName):
            shutil.rmtree(projectName)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

#=============================================================================-
    # Test bound creation
    def test_boundInit(self):

        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        fun3d = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        astros = problem.makeAnalysis("astrosAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        bound = problem.makeBound(2, self.boundName)
        
        vset_fun3d  = bound.makeVertexSet(fun3d)
        vset_astros = bound.makeVertexSet(astros)

        dset_src = vset_fun3d.makeDataSet("Pressure", caps.fType.FieldOut)
        dset_dst = vset_astros.makeDataSet("Pressure", caps.fType.FieldIn)

        dset_dst.linkDataSet(dset_src, caps.dMethod.Interpolate)

        dset_src = vset_astros.makeDataSet("Displacement", caps.fType.FieldOut)
        dset_dst = vset_fun3d.makeDataSet("Displacement", caps.fType.FieldIn)

        dset_dst.linkDataSet(dset_src, caps.dMethod.Conserve)

        bound.closeBound()

        dset_dst.initDataSet([0,0,0])

if __name__ == '__main__':
    unittest.main() 
