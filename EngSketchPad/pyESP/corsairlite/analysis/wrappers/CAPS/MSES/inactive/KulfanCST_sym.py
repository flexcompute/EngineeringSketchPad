import pyCAPS
import os
import numpy as np
import shutil
import math
from corsairlite.analysis.wrappers.CAPS.MSES.readSensx import readSensx

def msesRunKulfan_Alpha(Avar, Mach,Re,Alpha,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = 30, Fine_Iteration = 50,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunKulfan( Avar, 
                         Mach, Re,        
                         Alpha = Alpha,
                         CL = None,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )
    
def msesRunKulfan_CL(Avar, Mach,Re,CL,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = 50, Fine_Iteration = None,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunKulfan( Avar,
                         Mach, Re,        
                         Alpha = None,
                         CL = CL,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )



def msesRunKulfan(Avar,
                 Mach, Re,        
                 Alpha = None,
                 CL = None,
                 Acrit = None, GridAlpha = None, 
                 Coarse_Iteration = None, Fine_Iteration = None,
                 xGridRange = None, yGridRange = None, cleanup = True ):

    geometryScript = "Kulfan_airfoilSection.csm"
    runDirectory = "msesAnalysisTest"

    pstr = ""
    pstr += 'attribute capsAIM $xfoilAIM;tsfoilAIM;msesAIM \n'
    pstr += 'dimension class 1 2 \n'
    pstr += 'dimension ztail 1 2 \n'
    pstr += 'dimension avar 1 %d \n'%(len(Avar))
    pstr += 'dimension aupper 1 %d \n'%(len(Avar))
    pstr += 'dimension alower 1 %d \n'%(len(Avar))
    pstr += 'despmtr class "0.5; 1.0;" \n'
    pstr += 'despmtr ztail "0.00; 0.00;" \n'
    pstr += 'despmtr avar "'
    for vl in Avar:
        pstr += '%f;'%(vl)
    pstr += '"\n'
    pstr += 'set aupper "'
    for i in range(0,len(Avar)):
        pstr += 'avar[%d];'%(i+1)
    pstr += '"\n'
    pstr += 'set alower "'
    for i in range(0,len(Avar)):
        pstr += '-avar[%d];'%(i+1)
    pstr += '"\n'
    pstr += 'udparg kulfan class class \n'
    pstr += 'udparg kulfan ztail ztail \n'
    pstr += 'udparg kulfan aupper aupper \n'
    pstr += 'udparg kulfan alower alower \n'
    pstr += 'udprim kulfan \n'
    pstr += 'extract 0\n'

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
        mses.input.GridAlpha        = float(math.ceil(180.0 * CL / (4*np.pi**2))) # Approximating alpha wiht thin airfoil theory, best thought I have
    if Coarse_Iteration is not None:
        mses.input.Coarse_Iteration = Coarse_Iteration # Default 20
    if Fine_Iteration is not None:
        mses.input.Fine_Iteration   = Fine_Iteration   # Default 200
    if xGridRange is not None:
        mses.input.xGridRange       = xGridRange       # Default [-1.75, 2.75]
    if yGridRange is not None:
        mses.input.yGridRange       = yGridRange       # Default [-2.0, 2.5]
        
    # mses.input.Airfoil_Points = 231 # points
    mses.input.xTransition_Upper = 0.05
    mses.input.xTransition_Lower = 0.05
    mses.input.Acrit = 1

    if Mach <= 0.1 :
        mses.input.ISMOM = 2 # see mses documentation for M<=0.1

    mses.input.Design_Variable = {"avar": {}}

    outputValues = {}
    outputDerivatives = {}

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
    outputDerivatives["CD"]["avar"]      = mses.output["CD"].deriv("avar")  

    outputDerivatives["CD_p"] = {}
    outputDerivatives["CD_p"]["Mach"]      = mses.output["CD_p"].deriv("Mach")
    outputDerivatives["CD_p"]["Re"]        = mses.output["CD_p"].deriv("Re")
    outputDerivatives["CD_p"]["avar"]      = mses.output["CD_p"].deriv("avar")  

    outputDerivatives["CD_v"] = {}
    outputDerivatives["CD_v"]["Mach"]      = mses.output["CD_v"].deriv("Mach")
    outputDerivatives["CD_v"]["Re"]        = mses.output["CD_v"].deriv("Re")
    outputDerivatives["CD_v"]["avar"]      = mses.output["CD_v"].deriv("avar")  

    outputDerivatives["CD_w"] = {}
    outputDerivatives["CD_w"]["Mach"]      = mses.output["CD_w"].deriv("Mach")
    outputDerivatives["CD_w"]["Re"]        = mses.output["CD_w"].deriv("Re")
    outputDerivatives["CD_w"]["avar"]      = mses.output["CD_w"].deriv("avar")  

    outputDerivatives["CM"] = {}
    outputDerivatives["CM"]["Mach"]      = mses.output["CM"].deriv("Mach")
    outputDerivatives["CM"]["Re"]        = mses.output["CM"].deriv("Re")
    outputDerivatives["CM"]["avar"]      = mses.output["CM"].deriv("avar")  


    if Alpha is not None:
        outputDerivatives["CL"] = {}
        outputDerivatives["CL"]["Alpha"]     = mses.output["CL"].deriv("Alpha")
        outputDerivatives["CL"]["Mach"]      = mses.output["CL"].deriv("Mach")
        outputDerivatives["CL"]["Re"]        = mses.output["CL"].deriv("Re")
        outputDerivatives["CL"]["avar"]      = mses.output["CL"].deriv("avar")  

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
        outputDerivatives["Alpha"]["avar"]      = mses.output["Alpha"].deriv("avar")  

        outputDerivatives["CD"]["CL"]     = mses.output["CD"].deriv("CL")
        outputDerivatives["CD_p"]["CL"]   = mses.output["CD_p"].deriv("CL")
        outputDerivatives["CD_v"]["CL"]   = mses.output["CD_v"].deriv("CL")
        outputDerivatives["CD_w"]["CL"]   = mses.output["CD_w"].deriv("CL")
        outputDerivatives["CM"]["CL"]     = mses.output["CM"].deriv("CL")

    myProblem.close()
    
    # sensxFile = os.path.join(workDir,'Scratch','msesAIM0','sensx.airfoil')
    # sensDict = readSensx(sensxFile)

    # if cleanup:
    #     if os.path.exists(workDir) and os.path.isdir(workDir):
    #         shutil.rmtree(workDir)

    return [outputValues, outputDerivatives]#, sensDict]




# import matplotlib.pyplot as plt
# import numpy as np

# sweepParams = ["Mach","Re","Alpha", "Aupper1", "Aupper2", "Aupper3", "Aupper4", "Alower1", "Alower2", "Alower3", "Alower4"]

# for sweepParam in sweepParams:
#     Mach = 0.3
#     Re = 1e7
#     Alpha = 1.0
#     Aupper = [0.18, 0.23, 0.18, 0.22]
#     Alower = [-0.16, -0.07, -0.10, -0.06]

#     pltTitle = "Plots varying %s from a baseline NACA 2412 at %.1f deg AOA, Mach %.2f, and Re %.1e using the Kulfan CST"%(sweepParam.lower(),Alpha,Mach,Re)

#     if sweepParam == "Mach":
#         vals = np.linspace(0.2,0.6,9)
#     elif sweepParam == "Re":
#         vals = 10**np.linspace(6,7.5,11)
#     elif sweepParam == "Alpha":
#         vals = np.linspace(0,10,13)
#     elif sweepParam == "Aupper1":
#         vals = Aupper[0] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Aupper2":
#         vals = Aupper[1] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Aupper3":
#         vals = Aupper[2] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Aupper4":
#         vals = Aupper[3] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Alower1":
#         vals = Alower[0] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Alower2":
#         vals = Alower[1] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Alower3":
#         vals = Alower[2] + np.linspace(-.02,.02,10)
#     elif sweepParam == "Alower4":
#         vals = Alower[3] + np.linspace(-.02,.02,10)
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
#         elif 'Aupper' in sweepParam:
#             Aupper[int(sweepParam[-1])-1] = vl
#         elif 'Alower' in sweepParam:
#             Alower[int(sweepParam[-1])-1] = vl
#         else:
#             cody
            
#         opt = msesRunKulfan_Alpha(Aupper, Alower, Mach, Re, Alpha)
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
#         if 'Aupper' in sweepParam:
#             plt.plot(vals, np.array([rw[1][ylabs[i]]['aupper'][int(sweepParam[-1])-1] for rw in outputData]))
#         elif 'Alower' in sweepParam:
#             plt.plot(vals, np.array([rw[1][ylabs[i]]['alower'][int(sweepParam[-1])-1] for rw in outputData]))
#         else:
#             plt.plot(vals, np.array([rw[1][ylabs[i]][sweepParam] for rw in outputData]))
        
#         plt.xlabel("%s"%(sweepParam.lower()))
#         plt.ylabel("d(%s)/d(%s)"%(ylabs[i],sweepParam.lower()))
#         plt.legend(['FD Approx','CAPS Reported'])
#         plt.grid(1)

#     plt.savefig("Kulfan_%s_Study.png"%(sweepParam))
