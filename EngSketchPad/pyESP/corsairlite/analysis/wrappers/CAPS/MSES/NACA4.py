import pyCAPS
import os
import numpy as np
import shutil
import math
# from corsairlite.analysis.wrappers.CAPS.MSES.readSensx import readSensx

def msesRunNACA4_Alpha(camber,maxloc,thickness,Mach,Re,Alpha,
                         Acrit = None, Airfoil_Points = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xTransition_Upper = None, xTransition_Lower = None,
                         xGridRange = None, yGridRange = None, cleanup = True, problemObj = None ):
    
    
    return msesRunNACA4( camber, maxloc, thickness,
                         Mach, Re,        
                         Alpha = Alpha,
                         CL = None,
                         Acrit = Acrit, Airfoil_Points = Airfoil_Points, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xTransition_Upper = xTransition_Upper, xTransition_Lower = xTransition_Lower,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup, problemObj = problemObj )
    
def msesRunNACA4_CL(camber,maxloc,thickness,Mach,Re,CL,
                         Acrit = None, Airfoil_Points = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xTransition_Upper = None, xTransition_Lower = None,
                         xGridRange = None, yGridRange = None, cleanup = True, problemObj = None ):
    
    
    return msesRunNACA4( camber, maxloc, thickness,
                         Mach, Re,        
                         Alpha = None,
                         CL = CL,
                         Acrit = Acrit, Airfoil_Points = Airfoil_Points, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xTransition_Upper = xTransition_Upper, xTransition_Lower = xTransition_Lower,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup, problemObj = problemObj )



def msesRunNACA4(camber, maxloc, thickness,
                 Mach, Re,        
                 Alpha = None,
                 CL = None,
                 Acrit = None, Airfoil_Points = None, GridAlpha = None, 
                 Coarse_Iteration = None, Fine_Iteration = None,
                 xTransition_Upper = None, xTransition_Lower = None,
                 xGridRange = None, yGridRange = None, cleanup = True, problemObj = None ):
    
    geometryScript = "NACA4_airfoilSection.csm"
    runDirectory = "msesAnalysisTest"
    workDir = os.path.join(".",runDirectory)
    if os.path.exists(workDir) and os.path.isdir(workDir):
        shutil.rmtree(workDir)
    os.mkdir(workDir)

    if problemObj is not None:
        myProblem = problemObj
        if "mses" in problemObj.analysis:
            mses = problemObj.analysis["mses"]
        else:
            mses = problemObj.analysis.create(aim="msesAIM", name="mses")
        problemObj.geometry.despmtr["thickness"].value = thickness
        problemObj.geometry.despmtr["camber"].value = camber
        problemObj.geometry.despmtr["maxloc"].value = maxloc
    else:
        pstr = ""
        pstr += "attribute capsIntent $LINEARAERO \n"
        pstr += "attribute capsAIM $xfoilAIM;tsfoilAIM;msesAIM \n"
        pstr += "despmtr   camber       %.16f \n"%(camber)
        pstr += "despmtr   maxloc       %.16f \n"%(maxloc)
        pstr += "despmtr   thickness    %.16f \n"%(thickness)
        pstr += "udparg naca thickness  thickness \n"
        pstr += "udparg naca camber     camber \n"
        pstr += "udparg naca maxloc     maxloc \n"
        pstr += "udparg naca sharpte    1 \n"
        pstr += "udprim naca \n"

        print(pstr)
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
    if Airfoil_Points is not None:
        mses.input.Airfoil_Points   = Airfoil_Points # default is 181
    if GridAlpha is not None:
        mses.input.GridAlpha        = GridAlpha # Default 0.0
    elif Alpha is not None:
        mses.input.GridAlpha        = Alpha # Default 0.0
    else:
        mses.input.GridAlpha        = float(math.ceil(180.0 * CL / (4*np.pi**2)))  # Approximating alpha wiht thin airfoil theory, best thought I have
    if Coarse_Iteration is not None:
        mses.input.Coarse_Iteration = Coarse_Iteration # Default 20
    if Fine_Iteration is not None:
        mses.input.Fine_Iteration   = Fine_Iteration   # Default 200
    if xGridRange is not None:
        mses.input.xGridRange       = xGridRange       # Default [-1.75, 2.75]
    if yGridRange is not None:
        mses.input.yGridRange       = yGridRange       # Default [-2.5, 2.5]
    if xTransition_Upper is not None:
        mses.input.xTransition_Upper = xTransition_Upper  # Default is 1.0
    if xTransition_Lower is not None:    
        mses.input.xTransition_Lower = xTransition_Lower   #Default is 1.0

    if Mach <= 0.1 :
        mses.input.ISMOM = 2 # see mses documentation for M<=0.1
    
    mses.input.Design_Variable = {"camber": {}, "maxloc":{}, "thickness": {}} # should correspond to CSM despmtr values

    outputValues = {}
    outputDerivatives = {}

    firstOutput = mses.output.CD # this runs mses
    # print("CD:", firstOutput)
    
    if problemObj is None:
        msesOutput = os.path.join(".",runDirectory,"Scratch","msesAIM0","msesOutput.txt")
        f = open(msesOutput,'r')
        outputText = f.read()
        f.close()
        outputLines = outputText.split('\n')
        lastSection = '\n'.join(outputLines[-32:])
        if 'Converged on tolerance' not in lastSection:
            print('Case may not be converged')
            # print("Mach: %f CL: %f Re: %e camber: %f maxloc: %f thickness: %f "%(Mach,CL,Re,camber,maxloc,thickness))
    
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
    outputDerivatives["CD"]["camber"]    = mses.output["CD"].deriv("camber")
    outputDerivatives["CD"]["maxloc"]    = mses.output["CD"].deriv("maxloc")
    outputDerivatives["CD"]["thickness"] = mses.output["CD"].deriv("thickness")

    outputDerivatives["CD_p"] = {}
    outputDerivatives["CD_p"]["Mach"]      = mses.output["CD_p"].deriv("Mach")
    outputDerivatives["CD_p"]["Re"]        = mses.output["CD_p"].deriv("Re")
    outputDerivatives["CD_p"]["camber"]    = mses.output["CD_p"].deriv("camber")
    outputDerivatives["CD_p"]["maxloc"]    = mses.output["CD_p"].deriv("maxloc")
    outputDerivatives["CD_p"]["thickness"] = mses.output["CD_p"].deriv("thickness")

    outputDerivatives["CD_v"] = {}
    outputDerivatives["CD_v"]["Mach"]      = mses.output["CD_v"].deriv("Mach")
    outputDerivatives["CD_v"]["Re"]        = mses.output["CD_v"].deriv("Re")
    outputDerivatives["CD_v"]["camber"]    = mses.output["CD_v"].deriv("camber")
    outputDerivatives["CD_v"]["maxloc"]    = mses.output["CD_v"].deriv("maxloc")
    outputDerivatives["CD_v"]["thickness"] = mses.output["CD_v"].deriv("thickness")

    outputDerivatives["CD_w"] = {}
    outputDerivatives["CD_w"]["Mach"]      = mses.output["CD_w"].deriv("Mach")
    outputDerivatives["CD_w"]["Re"]        = mses.output["CD_w"].deriv("Re")
    outputDerivatives["CD_w"]["camber"]    = mses.output["CD_w"].deriv("camber")
    outputDerivatives["CD_w"]["maxloc"]    = mses.output["CD_w"].deriv("maxloc")
    outputDerivatives["CD_w"]["thickness"] = mses.output["CD_w"].deriv("thickness")

    outputDerivatives["CM"] = {}
    outputDerivatives["CM"]["Mach"]      = mses.output["CM"].deriv("Mach")
    outputDerivatives["CM"]["Re"]        = mses.output["CM"].deriv("Re")
    outputDerivatives["CM"]["camber"]    = mses.output["CM"].deriv("camber")
    outputDerivatives["CM"]["maxloc"]    = mses.output["CM"].deriv("maxloc")
    outputDerivatives["CM"]["thickness"] = mses.output["CM"].deriv("thickness")

    if Alpha is not None:
        outputDerivatives["CL"] = {}
        outputDerivatives["CL"]["Alpha"]     = mses.output["CL"].deriv("Alpha")
        outputDerivatives["CL"]["Mach"]      = mses.output["CL"].deriv("Mach")
        outputDerivatives["CL"]["Re"]        = mses.output["CL"].deriv("Re")
        outputDerivatives["CL"]["camber"]    = mses.output["CL"].deriv("camber")
        outputDerivatives["CL"]["maxloc"]    = mses.output["CL"].deriv("maxloc")
        outputDerivatives["CL"]["thickness"] = mses.output["CL"].deriv("thickness")
        
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
        outputDerivatives["Alpha"]["camber"]    = mses.output["Alpha"].deriv("camber")
        outputDerivatives["Alpha"]["maxloc"]    = mses.output["Alpha"].deriv("maxloc")
        outputDerivatives["Alpha"]["thickness"] = mses.output["Alpha"].deriv("thickness")
        
        outputDerivatives["CD"]["CL"]     = mses.output["CD"].deriv("CL")
        outputDerivatives["CD_p"]["CL"]   = mses.output["CD_p"].deriv("CL")
        outputDerivatives["CD_v"]["CL"]   = mses.output["CD_v"].deriv("CL")
        outputDerivatives["CD_w"]["CL"]   = mses.output["CD_w"].deriv("CL")
        outputDerivatives["CM"]["CL"]     = mses.output["CM"].deriv("CL")
    
    if problemObj is None:
        myProblem.close()
    
    if cleanup:
        if os.path.exists(workDir) and os.path.isdir(workDir):
            shutil.rmtree(workDir)
        
    return [outputValues, outputDerivatives]#, sensDict]

# import matplotlib.pyplot as plt
# import numpy as np
      
# sweepParamList = ["Mach","Re","Alpha","camber","maxloc","thickness"]

# for sweepParam in sweepParamList:
#     Mach = 0.3
#     Re = 1e7
#     Alpha = 1.0
#     camber = 0.02
#     maxloc = 0.40
#     thickness = 0.12

#     pltTitle = "Plots varying %s from a baseline NACA %d%d%d at %.1f deg AOA, Mach %.2f, and Re %.1e"%(sweepParam.lower(),int(camber*100),int(maxloc*10),int(thickness*100),Alpha,Mach,Re)

#     if sweepParam == "Mach":
#         vals = np.linspace(0.2,0.6,9)
#     elif sweepParam == "Re":
#         vals = 10**np.linspace(6,7.5,11)
#     elif sweepParam == "Alpha":
#         vals = np.linspace(0,10,13)
#     elif sweepParam == "camber":
#         vals = np.linspace(0,.06,7)
#     elif sweepParam == "maxloc":
#         vals = np.linspace(.1,.6,6)
#     elif sweepParam == "thickness":
#         vals = np.linspace(0.07,0.25,10)
#     else:
#         cody
        
#     outputData = []
#     for vl in vals:
#         if sweepParam == "Mach":
#             Mach = vl
#         elif sweepParam == "Re":
#             Re = vl
#         elif sweepParam == "Alpha":
#             Alpha = vl
#         elif sweepParam == "camber":
#             camber = vl
#         elif sweepParam == "maxloc":
#             maxloc = vl
#         elif sweepParam == "thickness":
#             thickness = vl
#         else:
#             cody

#         opt = msesRunNACA4_Alpha(camber, maxloc, thickness, Mach, Re, Alpha)
#         outputData.append(opt)

#     fig , ax = plt.subplots(figsize=(24,12), constrained_layout=True)
#     plt.rcParams['font.size'] = 18
#     fig.suptitle(pltTitle)
#     ylabs = ['CL','CD','CM']
#     for i in range(0,len(ylabs)):
#         plt.subplot(2,3,i+1)
#         resVals = np.array([rw[0][ylabs[i]] for rw in outputData])
#         plt.plot(vals,resVals)
#         plt.xlabel("%s"%(sweepParam.lower()))
#         plt.ylabel(ylabs[i])
#         plt.grid(1)

#         plt.subplot(2,3,i+4)
#         derApp = (resVals[1:]-resVals[0:-1])/(vals[1]-vals[0])
#         vals_avg = (vals[1:]+vals[0:-1])/2
#         plt.plot(vals_avg,derApp)
#         plt.plot(vals, np.array([rw[1][ylabs[i]][sweepParam] for rw in outputData]))

#         plt.xlabel("%s"%(sweepParam.lower()))
#         plt.ylabel("d(%s)/d(%s)"%(ylabs[i],sweepParam.lower()))
#         plt.legend(['FD Approx','CAPS Reported'])
#         plt.grid(1)

#     plt.savefig("NACA4_%s_Study.png"%(sweepParam))
