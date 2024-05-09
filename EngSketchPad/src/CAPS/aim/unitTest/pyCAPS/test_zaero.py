import unittest

import os
import glob
import shutil

import sys

import pyCAPS

class TestZAERO(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "workDir_zaero"
        cls.iProb = 1
        cls.cleanUp()

        geomFile = os.path.join("..","csmData","feaWingBEMAero.csm")
        cls.capsProblem = pyCAPS.Problem(cls.problemName, capsFile=geomFile, outLevel=0)

        cls.egadsTess = cls.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                        name = "egadsTess",
                                                        capsIntent = "Structure")

        # Modify local mesh sizing parameters
        cls.egadsTess.input.Edge_Point_Min = 4
        cls.egadsTess.input.Edge_Point_Max = 8

        m    = pyCAPS.Unit("meter")
        mm   = pyCAPS.Unit("mm")
        kg   = pyCAPS.Unit("kg")
        gram = pyCAPS.Unit("gram")
        s    = pyCAPS.Unit("s")
        K    = pyCAPS.Unit("Kelvin")

        # Load masstran aim
        cls.masstran = cls.capsProblem.analysis.create(aim = "masstranAIM",
                                                       capsIntent = "Structure",
                                                       unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

        cls.masstran.input["Surface_Mesh"].link(cls.egadsTess.output["Surface_Mesh"])


        # Set materials
        madeupium    = {"materialType" : "isotropic",
                        "density"      : 10 * gram/mm**3}

        cls.masstran.input.Material = {"madeupium" : madeupium}

        # Set properties
        shell = {"propertyType"      : "Shell",
                 "membraneThickness" : 2.0*mm,
                 "material"          : "madeupium"}

        pointMass = {"propertyType" : "ConcentratedMass",
                     "mass"         : 3.0 * kg,
                     "massInertia"  : [10.0, 11.0, 12.0, 13.0, 14.0, 15.0]*kg*m**2}

        cls.masstran.input.Property = {"Ribs"     : shell,
                                       "Rib_Root" : shell,
                                       "Spar1"    : shell,
                                       "Spar2"    : shell,
                                       "Skin"     : shell,
                                       "Rib_Root_Point": pointMass}

    @classmethod
    def tearDownClass(cls):
        del cls.capsProblem
        del cls.masstran
        del cls.egadsTess
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    # Test inputs
    def test_inputs(self):

        m   = pyCAPS.Unit("meter")
        kg  = pyCAPS.Unit("kg")
        s   = pyCAPS.Unit("s")
        K   = pyCAPS.Unit("Kelvin")
        kPa   = pyCAPS.Unit("kPa")

        zaero = self.capsProblem.analysis.create(aim = "zaeroAIM",autoExec=False,
                                                 unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

        zaero.input.Proj_Name = "zaero_CAPS"
        
        
        FlutterCruise = {"discipline":"flutter", "function":"FixMachDensity", "uaic":"cruise"}
        
        # Set trim analysis
        trim = {
            # case control
            # "label": "Trim at Mach 0.10",
        
            # case unsteady aerodynamics
            "uaic": "climb",
        
            # case analysis definition
            "discipline": "Trim",
            # "dynamicPressure": 50000,
            "dynamicPressure": 10 * kPa,
            # "vectorToCG": [0.0, 0.0, 0.0], # refers to RHOX, RHOY, RHOZ fields in TRIM card
            "gravityAcceleration": 9.81 * m/s**2,
            "accUnit": "G",
            "variables": ["ANGLEA", "URDD3"],
            }
        
        zaero.input.Analysis = {"FlutterCruise":FlutterCruise, "Trim":trim}
        
        zaero.input.File_Format = "Large"
        zaero.input.FEM_1 = {"fem": "nastran_1.f06",
                             "form": "msc",
                             "boundary": "sym",
                             "print": 0}

        zaero.input.FEM_2 = {"fem": "nastran_2.f06",
                             "form": "msc",
                             "boundary": "asym",
                             "print": 0}
        
        zaero.input.CPU = 2
        zaero.input.Memory = "1800MB"
        zaero.input.Smart_Restart = True
        zaero.input.Echo = "sort"

        zaero.input.HFG = {"XZSymmetric" : "yes",
                            "flip": "no",
                          }
 
        uaicClimb = {
            "machNumber": 0.20,
            "method": 0,
            "reducedFreq": [0.02, 0.1, 0.12, 0.14, 0.16, 0.18, 0.25],
            "print": 0,
        }
        uaicCruise = {
            "machNumber": 0.70,
            "method": 0,
            "reducedFreq": [0.04, 0.2, 0.22, 0.24, 0.26, 0.28, 0.35],
            "print": 0,
        }
        zaero.input.UAIC = {"climb": uaicClimb, "cruise": uaicCruise}
        
        zaero.input.Spline = {"method": 1, "attachFlex":0.1, "eps": 1.0e-6}

        wing = {"groupName"    : "Wing",
                "numChord"     : 4,
                "numSpanTotal" : 12}
        zaero.input.VLM_Surface = {"WingSurface": wing}

        flapLE = {"surfaceSymmetry": "SYM"}
        flapTE = {"surfaceSymmetry": "SYM"}
        zaero.input.VLM_Control = {"CntrlLE": flapLE, "CntrlTE": flapTE}
        
        zaero.input.Trim_Variable = {"trimVar": {"label": "CntrlTE", "value": 0.4}}

        angleOfAttack = {
            "label": "ALPHA",
            "value": "free"
        }
        verticalAccel = {
            "label": "NZ",
            "value": 1.0
        }
        zaero.input.Trim_Variable = {"ANGLEA": angleOfAttack,
                                     "URDD3": verticalAccel}

        zaero.input.ReferenceArea = 2 * m*m
        zaero.input.ReferenceChord = 0.6 * m
        zaero.input.ReferenceSpan = 6 * m
        zaero.input.Moment_Center = [1,0.1,0.2] * m
        
        zaero.input["MassPropLink"].link(self.masstran.output["MassPropLink"])

        # Set output file generation
        
        plotAero = {
            "filename": "aero.plt",
            "form": "TECPLOT"
        }

        numModes = 4
        modes = list(range(1, numModes + 1))
        plotMode = {
            "mode": modes,
            "boundary": ["sym" for m in modes],
            "filename": ["aeroMode_{}.plt".format(m) for m in modes],
            "form": ["TECPLOT" for m in modes]
        }
        
        plotTrim = {
            "analysis": "Trim",
            "type": "CP",
            "form": "TECPLOT",
            "filename": "trimcp.plt"
        }
        
        zaero.input.Output = {"Aero": plotAero,
                              "Mode": plotMode,
                              "Trim": plotTrim}


        # Write input files
        zaero.preAnalysis()

#==============================================================================
    def test_Output(self):

        m      = pyCAPS.Unit("meter")
        cm     = pyCAPS.Unit("centimeter")
        watt   = pyCAPS.Unit("watt")
        kelvin = pyCAPS.Unit("kelvin")
        deg    = pyCAPS.Unit("degree")
        bar    = pyCAPS.Unit("bar")

if __name__ == '__main__':
    unittest.main()
