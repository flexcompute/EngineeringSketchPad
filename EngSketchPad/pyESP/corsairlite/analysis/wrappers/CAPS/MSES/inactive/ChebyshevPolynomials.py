import os
import math
import shutil
import numpy as np
import pyCAPS
# from kulfan import Kulfan
# from readSensx import readSensx
# from corsairlite.core.dataTypes.kulfan import Kulfan
# from corsairlite.analysis.wrappers.CAPS.MSES.readSensx import readSensx

def msesRunChebyshev_Alpha(  Cheby,
                             Mach, Re,        
                             Alpha,
                             baselineAirfoil = None,
                             Acrit = None, GridAlpha = None, 
                             Coarse_Iteration = None, Fine_Iteration = None,
                             xGridRange = None, yGridRange = None, cleanup = True ):

    return msesRunChebyshev( 
                        Cheby=Cheby,
                         Mach=Mach, Re = Re,        
                         Alpha = Alpha,
                         CL = None,
                         baselineAirfoil = baselineAirfoil,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )

    
def msesRunChebyshev_CL( Cheby,
                         Mach, Re,        
                         CL,
                         baselineAirfoil = None,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = 199, Fine_Iteration = 205,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunChebyshev( 
                         Cheby=Cheby,
                         Mach=Mach, Re = Re,        
                         Alpha = None,
                         CL = CL,
                         baselineAirfoil = baselineAirfoil,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )



def msesRunChebyshev(
                 Cheby,
                 Mach, Re,        
                 Alpha = None,
                 CL = None,
                 baselineAirfoil = None,
                 Acrit = None, GridAlpha = None, 
                 Coarse_Iteration = None, Fine_Iteration = None,
                 xGridRange = None, yGridRange = None, cleanup = True ):
    
    Coarse_Iteration = 30
    Fine_Iteration = 50
    geometryScript = "Cheby_airfoilSection.csm"
    runDirectory = "msesAnalysisTest"

    if baselineAirfoil is None:
        pstr = ""
        pstr += "attribute capsIntent $LINEARAERO \n"
        pstr += "attribute capsAIM $xfoilAIM;tsfoilAIM;msesAIM \n"
        pstr += "despmtr   camber       %f \n"%(0.0)
        pstr += "despmtr   maxloc       %f \n"%(0.4)
        pstr += "despmtr   thickness    %f \n"%(0.12)
        pstr += "udparg naca thickness  thickness \n"
        pstr += "udparg naca camber     camber \n"
        pstr += "udparg naca maxloc     maxloc \n"
        pstr += "udparg naca sharpte    1 \n"
        pstr += "udprim naca \n"
        pstr += "extract 0 \n"
    else:
        Aupper = baselineAirfoil.upperCoefficients
        Alower = baselineAirfoil.lowerCoefficients

        pstr = ""
        pstr += 'attribute capsAIM $xfoilAIM;tsfoilAIM;msesAIM \n'
        pstr += 'dimension class 1 2 \n'
        pstr += 'dimension ztail 1 2 \n'
        pstr += 'dimension aupper 1 %d \n'%(len(Aupper))
        pstr += 'dimension alower 1 %d \n'%(len(Alower))
        pstr += 'despmtr class "0.5; 1.0;" \n'
        pstr += 'despmtr ztail "0.001; -0.001;" \n'
        pstr += 'despmtr aupper " '
        for vl in Aupper:
            pstr += " %f;"%(vl)
        pstr += '" \n'
        pstr += 'despmtr alower " '
        for vl in Alower:
            pstr += " %f;"%(vl)
        pstr += '" \n'
        pstr += 'udparg kulfan class class \n'
        pstr += 'udparg kulfan ztail ztail \n'
        pstr += 'udparg kulfan aupper aupper \n'
        pstr += 'udparg kulfan alower alower \n'
        pstr += 'udprim kulfan \n'
        pstr += 'extract 0 \n'

    workDir = os.path.join(".",runDirectory)
    if os.path.exists(workDir) and os.path.isdir(workDir):
        shutil.rmtree(workDir)
    os.mkdir(workDir)
    geometryScript_path = os.path.join(workDir,geometryScript)
    f = open(geometryScript_path,'w')
    f.write(pstr)
    f.close()

    myProblem = pyCAPS.Problem(problemName=workDir, capsFile=geometryScript_path, outLevel=0)
    mses = myProblem.analysis.create(aim = "msesAIM")
    
    mses.input.Mach = Mach
    mses.input.Re   = Re
    
    if CL is not None and Alpha is not None:
        raise ValueError('Conflict Detected:  Alpha and CL have both been set')
        
    if CL is None and Alpha is None:
        raise ValueError('Must specifiy either Alpha or CL')
        
    if CL is not None:
        if CL >= 3.0:
            raise ValueError('CL exceeds 3.0, did you mean to set Alpha?')
        mses.input.CL = CL
        
    if Alpha is not None:
        mses.input.Alpha = Alpha
        
    if Acrit is not None:
        mses.input.Acrit            = Acrit # Default 9.0
    if GridAlpha is not None:
        mses.input.GridAlpha        = GridAlpha # Default 0.0
    elif Alpha is not None:
        mses.input.GridAlpha        = Alpha # Default 0.0
    else:
        mses.input.GridAlpha        = math.ceil(180.0 * CL / (4*np.pi**2))  # Approximating alpha wiht thin airfoil theory, best thought I have
    if Coarse_Iteration is not None:
        mses.input.Coarse_Iteration = Coarse_Iteration # Default 20
    if Fine_Iteration is not None:
        mses.input.Fine_Iteration   = Fine_Iteration   # Default 200
    if xGridRange is not None:
        mses.input.xGridRange       = xGridRange       # Default [-1.75, 2.75]
    if yGridRange is not None:
        mses.input.yGridRange       = yGridRange       # Default [-2.0, 2.5]

    mses.input.Airfoil_Points = 231 # points
        
    if Mach <= 0.1 :
        mses.input.ISMOM = 2 # see mses documentation for M<=0.1, though this version seems to still need M>=0.001
    
    mses.input.Cheby_Modes = Cheby
    
    mses.input.xTransition_Upper = 0.05
    mses.input.xTransition_Lower = 0.05
    # mses.input.Acrit = 7

    outputValues = {}
    outputDerivatives = {}
    
    # try:
    firstOutput = mses.output.CD # this runs mses
    # except Exception as e:
    msesOutput = os.path.join(".",runDirectory,"Scratch","msesAIM0","msesOutput.txt")
    f = open(msesOutput,'r')
    outputText = f.read()
    f.close()
    outputLines = outputText.split('\n')
    lastSection = '\n'.join(outputLines[-32:])
    if 'Converged on tolerance' not in lastSection:
        print('Case may not be converged')
    # else:
    #     print('Case Converged')
    
    if Alpha is not None:
        outputValues["CL"]                = mses.output.CL
    if CL is not None:
        outputValues["Alpha"]             = mses.output.Alpha

    outputValues["CD"]                = mses.output.CD
    outputValues["CD_p"]              = mses.output.CD_p
    outputValues["CD_v"]              = mses.output.CD_v
    outputValues["CD_w"]              = mses.output.CD_w
    outputValues["CM"]                = mses.output.CM

    outputDerivatives["CD"] = {}
    outputDerivatives["CD"]["Mach"]      = mses.output["CD"].deriv("Mach")
    outputDerivatives["CD"]["Re"]        = mses.output["CD"].deriv("Re")
    outputDerivatives["CD"]["Cheby_Modes"]    = mses.output["CD"].deriv("Cheby_Modes")

    outputDerivatives["CD_p"] = {}
    outputDerivatives["CD_p"]["Mach"]      = mses.output["CD_p"].deriv("Mach")
    outputDerivatives["CD_p"]["Re"]        = mses.output["CD_p"].deriv("Re")
    outputDerivatives["CD_p"]["Cheby_Modes"]    = mses.output["CD_p"].deriv("Cheby_Modes")

    outputDerivatives["CD_v"] = {}
    outputDerivatives["CD_v"]["Mach"]      = mses.output["CD_v"].deriv("Mach")
    outputDerivatives["CD_v"]["Re"]        = mses.output["CD_v"].deriv("Re")
    outputDerivatives["CD_v"]["Cheby_Modes"]    = mses.output["CD_v"].deriv("Cheby_Modes")

    outputDerivatives["CD_w"] = {}
    outputDerivatives["CD_w"]["Mach"]      = mses.output["CD_w"].deriv("Mach")
    outputDerivatives["CD_w"]["Re"]        = mses.output["CD_w"].deriv("Re")
    outputDerivatives["CD_w"]["Cheby_Modes"]    = mses.output["CD_w"].deriv("Cheby_Modes")

    outputDerivatives["CM"] = {}
    outputDerivatives["CM"]["Mach"]      = mses.output["CM"].deriv("Mach")
    outputDerivatives["CM"]["Re"]        = mses.output["CM"].deriv("Re")
    outputDerivatives["CM"]["Cheby_Modes"]    = mses.output["CM"].deriv("Cheby_Modes")

    if Alpha is not None:
        outputDerivatives["CL"] = {}
        outputDerivatives["CL"]["Alpha"]     = mses.output["CL"].deriv("Alpha")
        outputDerivatives["CL"]["Mach"]      = mses.output["CL"].deriv("Mach")
        outputDerivatives["CL"]["Re"]        = mses.output["CL"].deriv("Re")
        outputDerivatives["CL"]["Cheby_Modes"]    = mses.output["CL"].deriv("Cheby_Modes")

        outputDerivatives["CD"]["Alpha"]     = mses.output["CD"].deriv("Alpha")
        outputDerivatives["CD_p"]["Alpha"]   = mses.output["CD_p"].deriv("Alpha")
        outputDerivatives["CD_v"]["Alpha"]   = mses.output["CD_v"].deriv("Alpha")
        outputDerivatives["CD_w"]["Alpha"]   = mses.output["CD_w"].deriv("Alpha")
        outputDerivatives["CM"]["Alpha"]     = mses.output["CM"].deriv("Alpha")

    if CL is not None:
        outputDerivatives["Alpha"] = {}
        outputDerivatives["Alpha"]["CL"]        = mses.output["Alpha"].deriv("CL")
        outputDerivatives["Alpha"]["Mach"]      = mses.output["Alpha"].deriv("Mach")
        outputDerivatives["Alpha"]["Re"]        = mses.output["Alpha"].deriv("Re")
        outputDerivatives["Alpha"]["Cheby_Modes"]    = mses.output["Alpha"].deriv("Cheby_Modes")

        outputDerivatives["CD"]["CL"]     = mses.output["CD"].deriv("CL")
        outputDerivatives["CD_p"]["CL"]   = mses.output["CD_p"].deriv("CL")
        outputDerivatives["CD_v"]["CL"]   = mses.output["CD_v"].deriv("CL")
        outputDerivatives["CD_w"]["CL"]   = mses.output["CD_w"].deriv("CL")
        outputDerivatives["CM"]["CL"]     = mses.output["CM"].deriv("CL")

    myProblem.close()

    # sensxFile = os.path.join(workDir,'Scratch','msesAIM0','sensx.airfoil')
    # sensDict = readSensx(sensxFile)
    
    if cleanup:
        if os.path.exists(workDir) and os.path.isdir(workDir):
            shutil.rmtree(workDir)

    return [outputValues, outputDerivatives]

