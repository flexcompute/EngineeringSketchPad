#
# Written by Dr. Ryan Durscher AFRL/RQVC
#
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

## \namespace pyCAPS
# Python extension module for CAPS
#
#\mainpage Introduction
# \tableofcontents
# \section overviewpyCAPS Overview
# pyCAPS is a Python extension module to interact with Computational
# Aircraft Prototype Syntheses (CAPS) routines in the
# Python environment. Written in Cython, pyCAPS natively handles all type
# conversions/casting, while logically grouping CAPS function calls together to simplify
# a user's experience. Additional functionality not directly available through the
# CAPS API (such has saving a geometric view) is also provided.
#
# An overview of the basic pyCAPS functionality is provided in \ref gettingStarted.
#
# \section differencespyCAPS Key differences between pyCAPS and CAPS
# - Manipulating the "owner" information for CAPS objects isn't currently supported
#
# - CAPS doesn't natively support an array of string values (an array of strings is viewed by CAPS as a single
# concatenated string), however pyCAPS does. If a list of strings is provided this list is concatenated, separated by a ';' and
# provided to CAPS as a single string. The number of rows and columns are correctly set to match the original list.
# If a string is received from CAPS by pyCAPS and the rows and columns are set correctly it will be unpacked correctly considering
# entries are separated by a ';'. (Important) If the rows and columns aren't set correctly and the string contains a ';',
#  the data will likely be unpacked incorrectly or raise an indexing error. (Note: not available when setting attributes on objects)
#
# \section clearancepyCAPS Clearance Statement
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

from . import caps
from pyEGADS import egads
from .attributeUtils import _createAttributeValList, _createAttributeMap, _reverseAttributeMap
from .treeUtils import writeTreeHTML
from .geometryUtils import _viewGeometryMatplotlib

# Python version
try:
    from sys import version_info
    PY_MAJOR_VERSION = version_info.major
except:
    print("Error: Unable to import Python.version")
    raise ImportError

# JSON string conversion
try:
    from json import dumps as jsonDumps
    from json import loads as jsonLoads
except:
    print("Error:  Unable to import json\n")
    raise ImportError

# OS module
try:
    import os
except:
    print("Error: Unable to import os\n")
    raise ImportError

# Time module
try:
    import time
except:
    print("Error:  Unable to import time\n")
    raise ImportError

# Math module
try:
    from math import sqrt
except:
    print("Error: Unable to import math\n")
    raise ImportError

import shutil

# Set version
__version__ = "2.2.1"

# ------------------------- #

# Handles all uni-code and and byte conversion. Originally intended to just convert
#     JSON unicode returns to str - taken from StackOverFlow question - "How to get string
#     objects instead of Unicode ones from JSON in Python?"
# def _byteify(data,ignore_dicts = False):
#
#         if isinstance(data, bytes):
#             return data
#
#         if PY_MAJOR_VERSION < 3:
#             if isinstance(data,unicode):
#                 return data.encode('utf-8')
#         else:
#             if isinstance(data,str):
#                 return data.encode('utf-8')
#
#         if isinstance(data,list):
#             return [_byteify(item,ignore_dicts=True) for item in data]
#         if isinstance(data, dict) and not ignore_dicts:
#             return {_byteify(key,ignore_dicts=True): _byteify(value, ignore_dicts=True) for key, value in data.items()}
#         return data

# Handles all byte to uni-code conversion.
# def _strify(data,ignore_dicts = False):
#
#     if PY_MAJOR_VERSION >= 3:
#         if isinstance(data, bytes):
#             return data.decode("utf-8")
#
#         if isinstance(data,list):
#             return [_strify(item,ignore_dicts=True) for item in data]
#         if isinstance(data, dict) and not ignore_dicts:
#             return {_strify(key,ignore_dicts=True): _strify(value, ignore_dicts=True) for key, value in data.items()}
#     return data

# Perform unit conversion
def capsConvert(value, fromUnits, toUnits):
    toUnits = caps.Unit(toUnits)
    return caps.Quantity(value,fromUnits).convert(toUnits)/toUnits

# Create the tree for an object
def createTree(dataObj, showAnalysisGeom, showInternalGeomAttr, reverseMap):

    def checkValue(data):

        if not data:
            return [{"name" : "None"}]

        elif isinstance(data,tuple):

            return [{"name" : data[0], "children" : checkValue(data[1])}]

        elif isinstance(data,list):

            specialType = False

            for i in data:
                if isinstance(i, (tuple, dict)):
                    specialType = True
                    break

            if specialType is True:
                temp = []

                for i in data:
                    temp.append(checkValue(i)[0])

                return temp
            else:
                return [{"name" : str(data)}]

        elif isinstance(data, dict):
            temp = []

            for key, value in data.items():
                temp.append({"name" : key, "children" : checkValue(value)})

            return temp

        else:
            return [{"name" : str(data)}]

    tree = {}

    # Analyis object
    if isinstance(dataObj, capsAnalysis):
        tree = {"name" : dataObj.aimName + " (" + dataObj.officialName + ")",
                "children" : [{"name" : "Directory", "children" : checkValue(dataObj.analysisDir)},
                              {"name" : "Parents"  , "children" : checkValue(dataObj.parents)},
                              {"name" : "Intent"   , "children" : checkValue(dataObj.capsIntent)},
                              {"name" : "Inputs"   , "children" : checkValue(None)},
                              {"name" : "Outputs"  , "children" : checkValue(None)},
                             ]
                }
        inputs = dataObj.getAnalysisVal()
        if inputs:

            tree["children"][3]["children"] = [] # Clear Inputs

            # Loop through analysis inputs and add current value
            for key, value in inputs.items():
                tree["children"][3]["children"].append( {"name" : key, "children" : checkValue(value)} )
        try:
            outputs = dataObj.getAnalysisOutVal()
        except caps.CAPSError as e:
            if e.errorName == "CAPS_DIRTY":
                temp = dataObj.getAnalysisOutVal(namesOnly=True)
                outputs = {}
                for i in temp:
                    outputs[i] = None
            else:
                raise(e)
        if outputs:

            tree["children"][4]["children"] = [] # Clear Outputs

            # Loop through analysis outputs and add current value
            for key, value in outputs.items():
                tree["children"][4]["children"].append( {"name" : key, "children" : checkValue(value)} )

        # Are we supposed to show geometry with the analysis obj?
        if showAnalysisGeom == True:

            attributeMap = dataObj.getAttributeMap(getInternal = showInternalGeomAttr, reverseMap = reverseMap)

            if attributeMap: # If attribute is populated

                tree["children"].append({"name" : "Geometry", "children" : []})

                for key, value in attributeMap.items():
                    tree["children"][5]["children"].append( {"name" : key, "children" : checkValue(value)} )

    # Data set object
    if isinstance(dataObj, _capsDataSet):
        tree = {"name" : dataObj.dataSetName,
                "children" : [{"name" : "Vertex Set", "children" : checkValue(dataObj.capsVertexSet.vertexSetName)},
                              {"name" : "Data Rank", "children"  : checkValue(dataObj.dataRank)},
                              {"name" : "Data Method", "children" : checkValue(dataObj.dataSetMethod)}
                              ]
                }

    # Bound object
    if isinstance(dataObj, capsBound):
        tree = {"name" : dataObj.boundName,
                "children" : [{"name" : "Source Data Set"      , "children" : checkValue(None)},
                              {"name" : "Destination Data Set" , "children" : checkValue(None)},
                             ]
                }

        # Source data sets
        if dataObj.dataSetSrc:

            tree["children"][0]["children"] = [] # Clear Source data set

            for key, value in dataObj.dataSetSrc.items():
                tree["children"][0]["children"].append(jsonLoads(createTree(value,
                                                                           showAnalysisGeom,
                                                                           showInternalGeomAttr,
                                                                           reverseMap
                                                                           )
                                                                )
                                                      )

        # Destination data sets - if there are any
        if dataObj.dataSetDest:

            tree["children"][1]["children"] = [] # Clear Destination data set

            for key, value in dataObj.dataSetDest.items():
                tree["children"][1]["children"].append(jsonLoads(createTree(value,
                                                                            showAnalysisGeom,
                                                                            showInternalGeomAttr,
                                                                            reverseMap
                                                                            )
                                                                 )
                                                       )

    # Geometry object
    if isinstance(dataObj, capsGeometry):
        tree = {"name" : dataObj.geomName,
                "children" : [{"name" : "Design Parameters", "children" : checkValue(None)},
                              {"name" : "Local Variables", "children" : checkValue(None)},
                              {"name" : "Attributes", "children" : checkValue(None)},
                             ]
                }

        # Design variables
        designVars = dataObj.getGeometryVal()
        if designVars:

            tree["children"][0]["children"] = [] # Clear design parameters

            # Loop through design variables
            for key, value in designVars.items():
                tree["children"][0]["children"].append( {"name" : key, "children" : checkValue(value)} )

        # Local variables
        localVars = dataObj.getGeometryOutVal()
        if localVars:

            tree["children"][1]["children"] = [] # Clear local variables

            # Loop through local variables
            for key, value in localVars.items():
                tree["children"][1]["children"].append( {"name" : key, "children" : checkValue(value)} )

        # Attributes
        attributeMap = dataObj.getAttributeMap(getInternal = showInternalGeomAttr, reverseMap = reverseMap)
        if attributeMap:
            tree["children"][2]["children"] = [] # Clear attribute map
            for key, value in attributeMap.items():
                tree["children"][2]["children"].append( {"name" : key, "children" : checkValue(value)} )

    # Problem object
    if isinstance(dataObj, capsProblem):
#         latest = {"CAPS Intent" : dataObj.capsIntent,
#                   "AIM Count"   : dataObj.aimGlobalCount,
#                   "Analysis Directory" : dataObj.analysisDir
#                   }

        tree = {"name" : "myProblem",
                "children" : [{"name" : "Analysis", "children" : checkValue(None)},
                              {"name" : "Geometry", "children" : checkValue(None)},
                              {"name" : "Bound",    "children" : checkValue(None)},
#                               {"name" : "Latest Defaults" , "children" : checkValue(latest)},
                             ]
                }

        # Do analysis objects - if there are any
        if dataObj.analysis:

            tree["children"][0]["children"] = [] # Clear Analysis

            for key, value in dataObj.analysis.items():
                tree["children"][0]["children"].append(jsonLoads(createTree(value,
                                                                            showAnalysisGeom,
                                                                            showInternalGeomAttr,
                                                                            reverseMap)
                                                                 )
                                                       )

        # Do geometry
        if dataObj.geometry: # There should be geometry

            tree["children"][1]["children"] = [jsonLoads(createTree(dataObj.geometry,
                                                                    showAnalysisGeom,
                                                                    showInternalGeomAttr,
                                                                    reverseMap
                                                                    )
                                                         )
                                               ]

        # Do data bounds - if there are any
        if dataObj.dataBound: #

            tree["children"][2]["children"] = [] # Clear Data Bounds

            for key, value in dataObj.dataBound.items():
                tree["children"][2]["children"].append(jsonLoads(createTree(value,
                                                                            showAnalysisGeom,
                                                                            showInternalGeomAttr,
                                                                            reverseMap
                                                                            )
                                                                 )
                                                       )


    return jsonDumps(tree, indent= 2)

def _storeIteration(initial, analysisDir, noFileDir=True):

    def initialInstance(analysisDir):

        # Get current directory
        currentDir = os.getcwd()

        # Change to analysis directory
        os.chdir(analysisDir)

        # Determine instance number
        instance = 0
        key = "Instance_"

        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]

        for i in currentFiles:

            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path) and key in i:
                keepFiles.remove(i)

                if int(i.replace(key, "")) >= instance:
                    instance = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return

        # Make new instance directory
        instance = os.path.join(os.getcwd(), key + str(instance))
        os.mkdir(instance)

        for i in keepFiles:
            # print("Change file - ", os.path.join(instance,i))
            os.rename(i, os.path.join(instance,i))

        # Change back to current directory
        os.chdir(currentDir)

    if initial:
        initialInstance(analysisDir)
        return

    # Get current directory
    currentDir = os.getcwd()

    # Change to analysis directory
    os.chdir(analysisDir)

    # Determine instance number
    interation = 0
    key = "Iteration_"
    key2 = "Instance_"

    currentFiles = os.listdir(os.getcwd());
    keepFiles = currentFiles[:]

    for i in currentFiles:
        path = os.path.join(os.getcwd(), i)

        if os.path.isdir(path):

            if key in i or key2 in i:
                keepFiles.remove(i)

            if key in i:
                if int(i.replace(key, "")) >= interation:
                    interation = int(i.replace(key, "")) + 1

    if noFileDir == False:
        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return

    # Make new interation directory
    interation = os.path.join(os.getcwd(), key + str(interation))
    if not os.path.isdir(interation):
        os.mkdir(interation)

    # Move files into new interation folder
    for i in keepFiles:
        # print"Change file - ", os.path.join(interation,i)
        os.rename(i, os.path.join(interation,i))

    # Change back to current directory
    os.chdir(currentDir)

def createOpenMDAOComponent(analysisObj,
                            inputParam,
                            outputParam,
                            executeCommand = None,
                            inputFile = None,
                            outputFile = None,
                            sensitivity = {},
                            changeDir = None,
                            saveIteration = False):

    try:
        from openmdao import __version__ as version
    except:
        print("Error: Unable to determine OpenMDAO version!")
        return caps.CAPS_NOTFOUND
    try:
        from numpy import ndarray
    except:
        print("Error: Unable to import numpy!")

    print("OpenMDAO version - ", version)
    versionMajor = int(version.split(".")[0])
    versionMinor = int(version.split(".")[1])

    if versionMajor >=2:
        try:
            from  openmdao.api import ExternalCodeComp, ExplicitComponent
        except:
            print("Error: Unable to import ExternalCodeComp or ExplicitComponent from openmdao.api!")
            return caps.CAPS_NOTFOUND

    elif versionMajor == 1 and versionMinor == 7:
        try:
            from  openmdao.api import ExternalCode, Component
        except:
            print("Error: Unable to import ExternalCode or Component from openmdao.api!")
            return caps.CAPS_NOTFOUND
    else:
        print("Error: Unsupported OpenMDAO version!")
        return caps.CAPS_BADVALUE

    # Change to lists if not lists already
    if not isinstance(inputParam, list):
        inputParam = [inputParam]

    if not isinstance(outputParam, list):
        outputParam = [outputParam]

    def filterInputs(varName, validInput, validGeometry):

        value = None

        if ":" in varName: # Is the input trying to modify a Tuple?

            analysisTuple = varName.split(":")[0]
            if analysisTuple not in validInput.keys():
                if varName not in validGeometry.keys(): #Check geometry
                    print("It appears the input parameter", varName, "is trying to modify a Tuple input, but", analysisTuple,
                          "this is not a valid ANALYSISIN or GEOMETRYIN variable! It will not be added to the OpenMDAO component")
                    return None
                else:
                    value = validGeometry[varName]
            else:

                if ":" in varName:
                    tuples = validInput[analysisTuple]
                    if not isinstance(tuples,list):
                        tuples = [tuples]

                    tupleValue = varName.split(":")[1]

                    for i in tuples:
                        if tupleValue == i[0]:
                            if varName.count(":") == 2:
                                dictValue = varName.split(":")[2]
                                value = i[1][dictValue]
                            elif varName.count(":") == 1:
                                value = i[1]
        else:
            # Check to make sure variable is in analysis
            if varName not in validInput.keys():
                if varName not in validGeometry.keys():
                    print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!",
                          "It will not be added to the OpenMDAO component.")
                    return None
                else:
                    value = validGeometry[varName]
            else:
                value = validInput[varName]

        # Check for None
        if value is None:
            value = 0.0

        return value

    def filterOutputs(varName, validOutput):

        value = None

        # Check to make sure variable is in analysis
        #if varName not in validOutput.keys():
        if varName not in validOutput:
            print("Output parameter", varName, "is not a valid ANALYSISOUT variable!",
                  "It will not be added to the OpenMDAO component.")
            return None
        #else:
        #    value = validOutput[varName]

        # Check for None
        if value is None:
            value = 0.0

        return value

    def setInputs(varName, varValue, analysisObj, validInput, validGeometry):

        if isinstance(varValue, ndarray):
            varValue = varValue.tolist()

        if isinstance(varValue,list):
            numRow, numCol = getPythonObjShape(varValue)
            if numRow == 1 and numCol == 1:
                    varValue = varValue[0]
            elif numRow == 1 and numCol > 1:
                varValue = varValue[0]

        if ":" in varName: # Is the input trying to modify a Analysis Tuple?

            analysisTuple = varName.split(":")[0]

            if analysisTuple not in validInput.keys():
                if varName in validGeometry.keys(): #Check geometry
                    analysisObj.capsProblem.geometry.setGeometryVal(varName, varValue)
                else:
                    print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!",
                          "It will not be modified to the OpenMDAO component. This should have already been reported!")
            else:
                if varName.count(":") == 2:

                    value = setElementofCAPSTuple(analysisObj.getAnalysisVal(analysisTuple), # Validity has already been checked
                                                  varName.split(":")[1], # Tuple key
                                                  varValue, # Value to set
                                                  varName.split(":")[2]) # Dictionary key
                else:

                    value = setElementofCAPSTuple(analysisObj.getAnalysisVal(analysisTuple), # Validity has already been checked
                                                  varName.split(":")[1],  # Tuple key
                                                  varValue) # Value to set

                analysisObj.setAnalysisVal(analysisTuple, value)

        else:

            if varName in validInput.keys():

                analysisObj.setAnalysisVal(varName, varValue)

            elif varName in validGeometry.keys():

                analysisObj.capsProblem.geometry.setGeometryVal(varName, varValue)
            else:

                print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!",
                      "It will not be modified to the OpenMDAO component. This should have already been reported!")

#     def gatherParents(analysis, analysisParents):
#
#         analysisParents += analysis.parents
#
#         if analysis.parents:
#             for parent in analysis.parents:
#                 analysisParents = gatherParents(analysis.capsProblem.analysis[parent], analysisParents)
#
#         return set(analysisParents)
#
    def checkParentExecution(analysis, checkInfo = False):

        def infoDictCheck(analysis):

            infoDict = analysis.getAnalysisInfo(printInfo = False, infoDict = True)

            if infoDict["executionFlag"] == True:
                return True
            else:
                print(analysis.aimName, "can't self execute!")
                return False

        if analysis.parents: # If we have parents

            # Check analysis parents
            for parent in analysis.parents:

                if not checkParentExecution(analysis.capsProblem.analysis[parent], True):
                    return False

                if checkInfo:
                    return infoDictCheck(analysis)

        if checkInfo:
            return infoDictCheck(analysis)

    def executeParent(analysis, runAnalysis = False):

        def executeAnalysis(analysis):
            infoDict = analysis.getAnalysisInfo(printInfo = False, infoDict = True)

            # Only execute analsysis if it isn't up to date
            if infoDict["status"] != 0: # 0 == Up to date
                analysis.preAnalysis()
                analysis.postAnalysis()

            return

        if analysis.parents:
            for parent in analysis.parents:
                executeParent(analysis.capsProblem.analysis[parent], runAnalysis = True)

            if runAnalysis:
                executeAnalysis(analysis)

            return

        if runAnalysis:
            executeAnalysis(analysis)

        return

    def initialInstance(analysisDir):

         # Get current directory
        currentDir = os.getcwd()

        # Change to analysis directory
        os.chdir(analysisDir)

        # Determine instance number
        instance = 0
        key = "Instance_"

        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]

        for i in currentFiles:

            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path) and key in i:
                keepFiles.remove(i)

                if int(i.replace(key, "")) >= instance:
                    instance = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return

        # Make new instance directory
        instance = os.path.join(os.getcwd(), key + str(instance))
        os.mkdir(instance)

        for i in keepFiles:
            # print"Change file - ", os.path.join(instance,i)
            os.rename(i, os.path.join(instance,i))

        # Change back to current directory
        os.chdir(currentDir)

    def storeIteration(analysisDir):

        # Get current directory
        currentDir = os.getcwd()

        # Change to analysis directory
        os.chdir(analysisDir)

        # Determine instance number
        interation = 0
        key = "Iteration_"
        key2 = "Instance_"

        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]

        for i in currentFiles:
            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path):

                if key in i or key2 in i:
                    keepFiles.remove(i)

                if key in i:
                    if int(i.replace(key, "")) >= interation:
                        interation = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return

        # Make new interation directory
        interation = os.path.join(os.getcwd(), key + str(interation))
        if not os.path.isdir(interation):
            os.mkdir(interation)

        # Move files into new interation folder
        for i in keepFiles:
            # print"Change file - ", os.path.join(interation,i)
            os.rename(i, os.path.join(interation,i))

        # Change back to current directory
        os.chdir(currentDir)

    if versionMajor >=2:
        class openMDAOExternal_V2(ExternalCodeComp):
            def __init__(self,
                         analysisObj,
                         inputParam, outputParam, executeCommand,
                         inputFile = None, outputFile = None,
                         sensitivity = {},
                         changeDir = True,
                         saveIteration = False):

                super(openMDAOExternal_V2, self).__init__()

                self.analysisObj = analysisObj
                self.inputParam = inputParam
                self.outputParam = outputParam
                self.executeCommand = executeCommand
                self.inputFile = inputFile
                self.outputFile = outputFile
                self.sensitivity = sensitivity
                self.changeDir = changeDir
                self.saveIteration = saveIteration

                # Get check to make sure
                self.validInput    = self.analysisObj.getAnalysisVal()
                self.validOutput   = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()

            def setup(self):

                for i in self.inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue

                    self.add_input(i, val=value)

                for i in self.outputParam:
                    value = filterOutputs(i, self.validOutput)
                    if value is None:
                        continue

                    self.add_output(i, val=value)

                if self.inputFile is not None:
                    if isinstance(self.inputFile, list):
                        self.options['external_input_files'] = self.inputFile
                    else:
                        self.options['external_input_files'] = [self.inputFile]

                if self.outputFile is not None:
                    if isinstance(self.outputFile, list):
                        self.options['external_output_files'] = self.outputFile
                    else:
                        self.options['external_output_files'] = [self.outputFile]

                if isinstance(self.executeCommand, list):
                    self.options['command'] = self.executeCommand
                else:
                    self.options['command'] = [self.executeCommand]

                method = "fd"
                step= None
                form = None
                step_calc = None

                for i in self.sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = self.sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = self.sensitivity[i]
                    elif "form" in i:
                        form = self.sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = self.sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)

                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return caps.CAPS_BADVALUE

#                  self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                self.declare_partials('*', '*', method=method, step=step, form=form, step_calc=step_calc)

                #self.changeDir = changeDir

                #self.saveIteration = saveIteration

                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)

            def compute(self, params, unknowns):

                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)

                # Run parents
                executeParent(self.analysisObj, runAnalysis = False)

                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)

                self.analysisObj.preAnalysis()

                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysisObj.analysisDir)

                #Parent compute function actually runs the external code
                super(openMDAOExternal_V2, self).compute(params, unknowns)

                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None:
                    os.chdir(currentDir)

                self.analysisObj.postAnalysis()

                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)

        class openMDAOExternal_V2ORIGINAL(ExternalCodeComp):
            def __init__(self,
                         analysisObj,
                         inputParam, outputParam, executeCommand,
                         inputFile = None, outputFile = None,
                         sensitivity = {},
                         changeDir = True,
                         saveIteration = False):

                super(openMDAOExternal_V2, self).__init__()

                self.analysisObj = analysisObj

                # Get check to make sure
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()

                for i in inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue

                    self.add_input(i, val=value)

                for i in outputParam:
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue

                    self.add_output(i, val=value)

                if inputFile is not None:
                    if isinstance(inputFile, list):
                        self.options['external_input_files'] = inputFile
                    else:
                        self.options['external_input_files'] = [inputFile]

                if outputFile is not None:
                    if isinstance(outputFile, list):
                        self.options['external_output_files'] = outputFile
                    else:
                        self.options['external_output_files'] = [outputFile]

                if isinstance(executeCommand, list):
                    self.options['command'] = executeCommand
                else:
                    self.options['command'] = [executeCommand]

                method = "fd"
                step= None
                form = None
                step_calc = None

                for i in sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = sensitivity[i]
                    elif "form" in i:
                        form = sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)

                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return caps.CAPS_BADVALUE

#                 self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                self.declare_partials('*', '*', method=method, step=step, form=form, step_calc=step_calc)

                self.changeDir = changeDir

                self.saveIteration = saveIteration

                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)

            def compute(self, params, unknowns):

                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)

                # Run parents
                executeParent(self.analysisObj, runAnalysis = False)

                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)

                self.analysisObj.preAnalysis()

                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysisObj.analysisDir)

                #Parent compute function actually runs the external code
                super(openMDAOExternal_V2, self).compute(params, unknowns)

                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None:
                    os.chdir(currentDir)

                self.analysisObj.postAnalysis()

                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)

            # TO add sensitivities
    #         def compute_partials(self, inputs, partials):
    #             outputs = {}
    #      the parent compute function actually runs the external code
    #         super(ParaboloidExternalCodeCompDerivs, self).compute(inputs, outputs)
    #
    #         # parse the derivs file from the external code and set partials
    #         with open(self.derivs_file, 'r') as derivs_file:
    #             partials['f_xy', 'x'] = float(derivs_file.readline())
    #             partials['f_xy', 'y'] = float(derivs_file.readline())

        class openMDAOComponent_V2(ExplicitComponent):

            def __init__(self,
                         analysisObj,
                         inputParam, outputParam,
                         sensitivity = {},
                         saveIteration = False):

                super(openMDAOComponent_V2, self).__init__()

                self.analysisObj = analysisObj

                # Get check to make sure
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()

                for i in inputParam:

                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue

                    self.add_input(i, val=value)

                for i in outputParam:

                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue

                    self.add_output(i, val=value)

                method = "fd"
                step= None
                form = None
                step_calc = None

                for i in sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = sensitivity[i]
                    elif "form" in i:
                        form = sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)

                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return caps.CAPS_BADVALUE

                self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)

                self.saveIteration = saveIteration

                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)

            def compute(self, params, unknowns):

                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)

                # Run parents
                executeParent(self.analysisObj, runAnalysis = False)

                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)

                self.analysisObj.preAnalysis()
                self.analysisObj.postAnalysis()

                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)

            # TO add sensitivities
    #         def compute_partials(self, inputs, partials):
    #             partials['area', 'length'] = inputs['width']
    #             partials['area', 'width'] = inputs['length']

    else:
        class openMDAOExternal(ExternalCode):
            def __init__(self,
                         analysisObj,
                         inputParam, outputParam, executeCommand,
                         inputFile = None, outputFile = None,
                         sensitivity = {},
                         changeDir = True,
                         saveIteration = False):

                super(openMDAOExternal, self).__init__()

                self.analysisObj = analysisObj

                # Get check to make sure
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()

                for i in inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue

                    self.add_param(i, val=value)

                for i in outputParam:
                    value = filterOutputs(i, validOutput)

                    if value is None:
                        continue

                    self.add_output(i, val=value)

                if inputFile is not None:
                    if isinstance(inputFile, list):
                        self.options['external_input_files'] = inputFile
                    else:
                        self.options['external_input_files'] = [inputFile]

                if outputFile is not None:
                    if isinstance(outputFile, list):
                        self.options['external_output_files'] = outputFile
                    else:
                        self.options['external_output_files'] = [outputFile]

                if isinstance(executeCommand, list):
                    self.options['command'] = executeCommand
                else:
                    self.options['command'] = [executeCommand]

                for i in sensitivity.keys():
                    self.deriv_options[i] = sensitivity[i]

                self.changeDir = changeDir

                self.saveIteration = saveIteration

                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)

            def solve_nonlinear(self, params, unknowns, resids):

                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)

                # Run parents
                executeParent(self.analysisObj, runAnalysis = False)

                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)

                self.analysisObj.preAnalysis()

                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysisObj.analysisDir)

                #Parent solve_nonlinear function actually runs the external code
                super(openMDAOExternal, self).solve_nonlinear(params, unknowns, resids)

                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None:
                    os.chdir(currentDir)

                self.analysisObj.postAnalysis()

                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)

            # TO add sensitivities
            #def linearize(self, params, unknowns, resids):
            #    J = {}
            #    return J

        class openMDAOComponent(Component):

            def __init__(self,
                         analysisObj,
                         inputParam, outputParam,
                         sensitivity = {},
                         saveIteration = False):

                super(openMDAOComponent, self).__init__()

                self.analysisObj = analysisObj

                # Get check to make sure
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()

                for i in inputParam:

                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue

                    self.add_param(i, val=value)

                for i in outputParam:

                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue

                    self.add_output(i, val=value)

                for i in sensitivity.keys():
                    self.deriv_options[i] = sensitivity[i]

                self.saveIteration = saveIteration

                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)

            def solve_nonlinear(self, params, unknowns, resids):

                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)

                # Run parents
                executeParent(self.analysisObj, runAnalysis = False)

                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)

                self.analysisObj.preAnalysis()
                self.analysisObj.postAnalysis()

                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)

            # TO add sensitivities
            #def linearize(self, params, unknowns, resids):
            #    J = {}
            #    return J

    if checkParentExecution(analysisObj) == False: #Parents can't all self execute
        return caps.CAPS_BADVALUE

    infoDict = analysisObj.getAnalysisInfo(printInfo = False, infoDict = True)

    if executeCommand is not None and infoDict["executionFlag"] == True:

        if versionMajor >=2:
            print("An execution command was provided, but the AIM says it can run itself! Switching ExternalCodeComp to ExplicitComponent")
        else:
            print("An execution command was provided, but the AIM says it can run itself! Switching ExternalCode to Component")

    if executeCommand is not None and infoDict["executionFlag"] == False:

        if versionMajor >=2:
            return openMDAOExternal_V2(analysisObj, inputParam, outputParam,
                                    executeCommand, inputFile, outputFile, sensitivity, changeDir,
                                    saveIteration)

        elif versionMajor == 1 and versionMinor == 7:
            return openMDAOExternal(analysisObj, inputParam, outputParam,
                                    executeCommand, inputFile, outputFile, sensitivity, changeDir,
                                    saveIteration)
    else:

        if versionMajor >=2:
           return openMDAOComponent_V2(analysisObj, inputParam, outputParam, sensitivity, saveIteration)

        elif versionMajor == 1 and versionMinor == 7:
            return openMDAOComponent(analysisObj, inputParam, outputParam, sensitivity, saveIteration)


## Defines a CAPS problem object.
# A capsProblem is the top-level object for a single mission/problem. It maintains a single set
# of interrelated geometric models (see \ref pyCAPS.capsGeometry),
# analyses to be executed (see \ref pyCAPS.capsAnalysis),
# connectivity and data (see \ref pyCAPS.capsBound)
# associated with the run(s), which can be both multi-fidelity and multi-disciplinary.
class capsProblem:

    ## Initialize the problem.
    # See \ref problem.py for a representative use case.
    # \param raiseException Raise an exception after a CAPS error is encountered (default - True). See \ref raiseException.
    def __init__(self, raiseException=True):
        ## \example problem.py
        # Basic example use case of the initiation (pyCAPS.capsProblem.__init__) and
        # termination (pyCAPS.capsProblem.closeCAPS)

        ## Raise an exception after a CAPS error is found (default - True). Disabling (i.e. setting to False) may have unexpected
        # consequences; in general the value should be set to True.
        self.raiseException = raiseException

        # The Problem Object
        self.problemObj = None

        # Geometry
        ## Geometry object loaded into the problem. Set via \ref loadCAPS.
        self.geometry = None

        # Analysis
        # Number of AIMs loaded into the problem.
        #self.aimGlobalCount = 0

        # Current analysis directory which was used to load the latest AIM.
        #self.analysisDir = None

        # Default intent used to load the AIM if capsIntent is not provided.
        #self.capsIntent = None

        ## Dictionary of analysis objects loaded into the problem. Set via \ref loadAIM.
        self.analysis = {}

        ## Dictionary of data transfer/bound objects loaded into the problem. Set
        # via \ref createDataBound or \ref createDataTransfer.
        self.dataBound = {}

        ## Dictionary of value objects loaded into the problem. Set via \ref createValue.
        self.value = {}

    def __del__(self):
        del self.geometry
        del self.analysis
        del self.dataBound
        del self.value
        del self.problemObj

    ## Loads a *.csm, *.caps, or *.egads file into the problem.
    # See \ref problem1.py, \ref problem2.py, and \ref problem8.py for example use cases.
    #
    # \param capsFile CAPS file to load. Options: *.csm, *.caps, or *.egads. If the filename has a *.caps extension
    # the pyCAPS analysis, bound, and value objects will be re-populated (see remarks).
    #
    # \param projectName CAPS project name. projectName=capsFile if not provided.
    # \param verbosity Level of output verbosity. See \ref setVerbosity .
    #
    # \return Optionally returns the reference to the \ref geometry class object (\ref pyCAPS.capsGeometry).
    #
    # \remark
    # Caveats of loading an existing CAPS file:
    # - Can currently only load *.caps files generated from pyCAPS originally.
    # - OpenMDAO objects won't be re-populated for analysis objects
    def loadCAPS(self, capsFile, projectName=None, verbosity=1):
        ## \example problem1.py
        # Example use case for the pyCAPS.capsProblem.loadCAPS() function in which are multiple
        # problems with geometry are created.

        ## \example problem2.py
        # Example use case for the pyCAPS.capsProblem.loadCAPS() function - compare "geometry" attribute with
        # returned geometry object.

        ## \example problem8.py
        # Example use case for the pyCAPS.capsProblem.loadCAPS() and pyCAPS.capsProblem.saveCAPS() functions - using a CAPS file
        # to initiate a new problem.

        capsGeometry(self, capsFile, projectName, verbosity)

        return self.geometry

    ## Set the verbosity level of the CAPS output.
    # See \ref problem5.py for a representative use case.
    #
    # \param verbosityLevel Level of output verbosity. Options: 0 (or "minimal"),
    # 1 (or "standard") [default], and 2 (or "debug").
    def setVerbosity(self, verbosityLevel):
        ## \example problem5.py
        # Basic example for setting the verbosity of a problem using pyCAPS.capsProblem.setVerbosity() function.

        verbosity = {"minimal"  : 0,
                     "standard" : 1,
                     "debug"    : 2}

        if not isinstance(verbosityLevel, (int,float)):
            if verbosityLevel.lower() in verbosity.keys():
                verbosityLevel = verbosity[verbosityLevel.lower()]

        if int(verbosityLevel) not in verbosity.values():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while setting verbosity, invalid verbosity level!")

        self.problemObj.setOutLevel(verbosityLevel)

    ## Save a CAPS problem (deprecated and does nothing).
    # See \ref problem8.py for example use case.
    # \param filename File name to use when saving CAPS problem.
    def saveCAPS(self, filename="saveCAPS.caps"):
        pass

    ## Close a CAPS problem.
    # See \ref problem1.py for a representative use case.
    def closeCAPS(self):

        # delete all capsObj instances
        for analysis in self.analysis.values():
            del analysis.analysisObj

        for value in self.value.values():
            del value.valueObj

        for bound in self.dataBound.values():
            for vset in bound.vertexSet.values():
                del vset.vertexSetObj
                
            for dset in bound.dataSetSrc.values():
                del dset.dataSetObj

            for dset in bound.dataSetDest.values():
                del dset.dataSetObj

            del bound.boundObj

        del self.problemObj

        # re-initialize everyhing to an empty state
        self.__init__()

    ## Report what analyses loaded into the problem are dirty.
    #
    # \return Optionally returns a list of names of the dirty analyses. An empty list is returned if no
    # analyses are dirty.
    def dirtyAnalysis(self):

        numDirty = 0
        dirtyAnalysis = []

        for i in self.analysis:
            if self.analysis[i].getAnalysisInfo(printInfo=False) > 0:
                numDirty = numDirty + 1
                dirtyAnalysis.append(self.analysis[i].aimName)

        print("Number of dirty analyses = ", numDirty)
        if len(dirtyAnalysis) > 0:
            print("Dirty analyses = ", dirtyAnalysis)

        return dirtyAnalysis

    ## Load an AIM (Analysis Interface Module) into the problem. See
    # examples \ref problem3.py and  \ref problem4.py for typical representative use cases.
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param aim Name of the requested AIM.
    #
    # \param altName Alternative name to use when referencing AIM inside
    # the problem (dictionary key in \ref analysis). The name of the AIM, aim, will be used if no
    # \ref capsProblem.loadAIM$altName is provided (see remarks).
    #
    # \param analysisDir  Directory for AIM analysis. If none is provided the directory of the last loaded AIM
    # will be used; if no AIMs have been load the current working directory is used. CAPS requires that an unique
    # directory be specified for each instance of AIM.
    #
    # \param capsIntent Analysis intention in which to invoke the AIM.
    #
    # \param parents Single or list of parent AIM names to initilize the AIM with.
    #
    # \param copyAIM Name of AIM to copy. Creates the new AIM instance by duplicating an existing AIM. Analysis directory
    # (\ref loadAIM$analysisDir) and \ref loadAIM$altName should be provided and different
    # from the AIM being copied. See example \ref analysis2.py for a representative use case.
    #
    # \param copyParents When copying an AIM, should the same parents also be used (default - True).
    #
    # \return Optionally returns the reference to the analysis dictionary (\ref analysis) entry created
    # for the analysis class object (\ref pyCAPS.capsAnalysis).
    #
    # \remark If no \ref capsProblem.loadAIM$altName is provided and an AIM with the name, \ref capsProblem.loadAIM$aim,
    # has already been loaded, an alternative
    # name will be automatically specified with the syntax \ref capsProblem.loadAIM$aim_[instance number]. If an
    # \ref capsProblem.loadAIM$altName is provided it must be unique compared to other instances of the loaded AIM.
    def loadAIM(self, **kwargs):
        ## \example problem3.py
        # Example use cases for the pyCAPS.capsProblem.loadAIM() function.

        ## \example problem4.py
        # Example use cases for the pyCAPS.capsProblem.loadAIM() function - compare analyis dictionary entry
        # with returned analyis object.

        ## \example analysis2.py
        # Duplicate an AIM using the \ref pyCAPS.capsProblem.loadAIM()$copyAIM keyword argument of the pyCAPS.capsProblem.loadAIM() function.

        # Load the aim
        return capsAnalysis(self, **kwargs)

    ## Alteranative to \ref createDataBound.
    # Enforces that at least 2 AIMs must be already loaded into the problem. See \ref createDataBound
    # for details.
    def createDataTransfer(self, **kwargs):

        if len(self.analysis.keys()) < 2:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating a data transfer"
                                                         + "At least two AIMs need to be loaded before a data transfer can be setup")

        return self.createDataBound(**kwargs)

    ## Create a CAPS data bound/transfer into the problem.
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param capsBound Name of capsBound to use for the data bound.
    # \param variableName Single or list of variables names to add.
    # \param aimSrc Single or list of AIM names that will be the data sources for the bound.
    # \param aimDest Single or list of AIM names that will be the data destinations during the transfer.
    # \param transferMethod Single or list of transfer methods to use during the transfer.
    # \param initValueDest Single or list of initial values for the destication data.
    #
    # \return Optionally returns the reference to the data bound dictionary (\ref dataBound) entry created
    # for the bound class object (\ref pyCAPS.capsBound).
    def createDataBound(self, **kwargs):
        return capsBound(self, **kwargs)


    ## Create a CAPS value object. Only a value of subtype in PARAMETER is currently supported.
    # See \ref value.py for a representative use case. Value objects are stored the \ref value dictionary.
    #
    # \param name Name used to define the value. This will be used as the keyword entry in the \ref value
    # dictionary.
    #
    # \param data Data value(s) for the variable. Note that type casting in done automatically based on the
    # determined type of the Python object. Be careful with the Integers and Floats/Reals, for example 10 would be
    # type casted as an Integer, while 10.0 would be a float - this small discrepancy may lead to type errors when
    # linking values to analysis inputs.
    #
    # \param units Units associated with the value (default - None).
    #
    # \param limits Valid/acceptable range for the value (default - None).
    #
    # \param fixedLength Should the length of the value object be fixed (default - True)? For example if
    # the object is initialized with a value of [1, 2] it can not be changed to [1, 2, 3].
    #
    # \param fixedShape Should the shape of the value object be fixed (default - True)? For example if
    # the object is initialized with a value of 1 it can not be changed to [1, 2] or [[1, 2, 3], [4, 5, 6]] can
    # not be changed to [[1, 2], [4, 5]].
    #
    # \return Optionally returns the reference to the value dictionary (\ref value) entry created
    # for the value (\ref pyCAPS._capsValue).
    def createValue(self, name, data, units=None, limits=None, fixedLength=True, fixedShape=True):
        ## \example value.py
        # Example use cases for pyCAPS.capsProblem.createValue() function.

        if name in self.value.keys():
            raise caps.CAPSError(caps.CAPS_BADNAME, msg="A value with the name," + name + "already exists!")

        self.value[name] = _capsValue(self, name, data,
                                      units=units, limits=limits,
                                      subType="Parameter",
                                      fixedLength=fixedLength, fixedShape=fixedShape)

        return self.value[name]

    ## Create a link between a created CAPS value parameter and analyis inputs of <b>all</b> loaded AIMs, automatically.
    # Valid CAPS value, parameter objects must be created with \ref createValue(). Note, only links to
    # ANALYSISIN inputs are currently made at this time.
    # See \ref value6.py for a representative use case.
    #
    # \param value Value to use when creating the link (default - None). A combination (i.e. a single or list)
    # of \ref value dictionary entries and/or value object instances (returned from a
    # call to \ref createValue()) can be used. If no value is provided, all entries in the
    # value dictionary (\ref value) will be used.
    def autoLinkValue(self, value=None):
        ## \example value6.py
        # Example use cases for pyCAPS.capsProblem.autoLinkValue() function.

        # If none we want to try to link all values
        if value is None:
            value = list(self.value.keys())

        # If not a list convert it to one
        if not isinstance(value, list):
            value = [value]

        # Loop through the values
        for i in value:

            # If it is a capsValue instance grab the name
            if isinstance(i, _capsValue):
               i = i.name

            # Check to make sure name is in value dictionary
            if i not in self.value.keys():
                raise caps.CAPSError(caps.CAPS_NOTFOUND, msg = "while autolinking values in problem!"
                                                             + "Unable to find" + str(i) + "value dictionary!")

            varname = self.value[i].name

            found = False
            # Loop through analysis dictionary
            for analysisName in self.analysis.keys():

                # Look to see if the value name is present as an analysis input variable
                try:
                    tempObj = self.analysis[analysisName]._analysisObj().childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, varname)
                except caps.CAPSError as e:
                    if e.errorName == "CAPS_NOTFOUND":
                        continue
                    else:
                        raise

                vObj = self.value[i]._valueObj()

                # Try to make the link if we made it this far
                tempObj.linkValue(vObj, caps.tMethod.Copy)

                print("Linked", str(self.value[i].name), "to analysis", str(analysisName), "input", varname)
                found = True

            if not found:
                print("No linkable data found for", str(self.value[i].name))

    ## Add an attribute (that is meta-data) to the problem object. See example
    # \ref problem7.py for a representative use case.
    #
    # \param name Name used to define the attribute.
    #
    # \param data Data value(s) for the attribute. Note that type casting in done
    # automatically based on the determined type of the Python object.
    def addAttribute(self, name, data):
        ## \example problem7.py
        # Example use cases for interacting the
        # pyCAPS.capsProblem.addAttribute() and
        # pyCAPS.capsProblem.getAttribute() functions.

        valueObj = _capsValue(self, name, data, subType="User")
        self.problemObj.setAttr(valueObj.valueObj)

    ## Get an attribute (that is meta-data) that exists on the problem object. See example
    # \ref problem7.py for a representative use case.
    # \param name Name of attribute to retrieve.
    # \return Value of attribute "name"
    def getAttribute(self, name):

        valueObj = self.problemObj.attrByName(name)

        return _capsValue(self, None, None)._setupObj(valueObj).value

    def __getAIMParents(self, parentList):
        """
        Check the requested parent list against the analysis dictionary
        """

        for i in parentList:
            if i in self.analysis.keys():
                pass
            else:
                print("Requested parent", i, "has not been initialized/loaded")
                return False

        return True

    def _repopulateProblem(self):

        # Repopulate analysis class
        numObj = self.problemObj.size(caps.oType.ANALYSIS, caps.sType.NONE)

        for i in range(numObj):
            tempObj = self.problemObj.childByIndex(caps.oType.ANALYSIS, caps.sType.NONE, i+1)

            tempAnalysis = capsAnalysis(self, __fromSave=True)._setupObj(tempObj)

            self.analysis[tempAnalysis.aimName] = tempAnalysis

            # Set other analysis related variable in the problem
            #self.aimGlobalCount += 1
            #self.capsIntent  = self.analysis[tempAnalysis.aimName].capsIntent
            #self.analysisDir = self.analysis[tempAnalysis.aimName].analysisDir

        # Repopulate value class
        numObj = self.problemObj.size(caps.oType.VALUE, caps.sType.PARAMETER)

        for i in range(numObj):
            tempObj = self.problemObj.childByIndex(caps.oType.VALUE, caps.sType.PARAMETER, i+1)

            tempValue =  _capsValue(self, None, None)._setupObj(tempObj)

            self.value[tempValue.name] = tempValue


        # Repopulate databound class
        numObj = self.problemObj.size(caps.oType.BOUND, caps.sType.NONE)

        for i in range(numObj):
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "Unable to currently repopulate dataBound from a *.caps file")

            #     BOUND, VERTEXSET, DATASET
            #tempObj = self.problemObj.childByIndex(self.problemObj,
            #                                       caps.BOUND, caps.NONE,
            #                                       i+1)
            #self.dataBound[tempValue.name] = tempObj

    ## Create a HTML dendrogram/tree of the current state of the problem. See
    # example \ref problem6.py for a representative use case. The HTML file relies on
    # the open-source JavaScript library, D3, to visualize the data.
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file.
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet
    # access will be necessary to view the tree.
    #
    # \param filename Filename to use when saving the tree (default - "myProblem"). Note an ".html" is
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default - True)? If set to True internet access
    # will be necessary to view the tree.
    #
    # \param analysisGeom Show the geometry for each analysis entity (default - False).
    #
    # \param internalGeomAttr Show the internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    #
    # \param reverseMap Reverse the geometry attribute map (default - False).
    def createTree(self, filename = "myProblem", **kwargs):
        ## \example problem6.py
        # Example use case for the pyCAPS.capsProblem.createTree() function.

        embedJSON = kwargs.pop("embedJSON", True)

        if embedJSON is False:
            try:
                # Write JSON data file
                fp = open(filename + ".json",'w')
                fp.write( createTree(self,
                                     kwargs.pop("analysisGeom", False),
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                    )
                         )

                fp.close()
            except:
                raise CAPSError(caps.CAPS_IOERR, msg = "while creating HTML tree")

            # Write HTML file
            writeTreeHTML(filename +".html",
                          filename + ".json",
                          None,
                          kwargs.pop("internetAccess", True))
        else:
            writeTreeHTML(filename+".html",
                          None,
                          createTree(self,
                                     kwargs.pop("analysisGeom", False),
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                    ),
                          kwargs.pop("internetAccess", True)
                         )

    # Sanitize a filename
    def _createFilename(self, filename, directory, extension):

        if filename is None:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating a filename. Filename is None!")

        # Check to see if directory exists
        if not os.path.isdir(directory):
            print("Directory does not currently exist while creating a filename - it will" +
                    " be made automatically")
            os.makedirs(directory)

        filename = os.path.join(directory, filename)

        if "." not in extension:
            extension = "." + extension

        # only add an extension if the filename is missing one
        file_path, file_extension = os.path.splitext(filename)
        if file_extension == "":
            filename = filename+extension

        # Check to see if file exists
        try:
            os.remove(filename)
            print("File already exists - file " + filename + " will be deleted")
        except OSError:
           pass

        return filename

    def _getBodyName(self, ibody, nameList, body):

        name = body.attributeRet("_name")
        if name is None:
            name = "Body_" + str(ibody+1)

        if name in nameList:
            name = "Body_" + str(ibody+1)
            print("Body", ibody, ": Name '", name, "' found more than once. Changing name to: '", name, "'")

        nameList.append(name)

        return name, nameList

## Defines a CAPS analysis object.
# May be created directly by providing a problem object or indirectly through \ref capsProblem.loadAIM .
class capsAnalysis:

    ## Initializes an analysis object. See \ref capsProblem.loadAIM for additional information and caveats.
    #
    # \param problemObject capsProblem object in which the analysis object will be added to.
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param aim Name of the requested AIM.
    #
    # \param altName Alternative name to use when referencing AIM inside
    # the problem (dictionary key in \ref capsProblem.analysis). The name of the AIM, aim, will be used if no
    # alternative name is provided (see remarks).
    #
    # \param analysisDir  Directory for AIM analysis. If none is provided the directory of the last loaded AIM
    # will be used; if no AIMs have been load the current working directory is used. CAPS requires that an unique
    # directory be specified for each instance of AIM.
    #
    # \param capsIntent Analysis intention in which to invoke the AIM.
    #
    # \param parents Single or list of parent AIM names to initilize the AIM with.
    #
    # \param copyAIM Name of AIM to copy. Creates the new AIM instance by duplicating an existing AIM. Analysis directory
    # (\ref capsProblem.loadAIM$analysisDir) and \ref capsProblem.loadAIM$altName should be provided and different
    # from the AIM being copied. See example \ref analysis2.py for a representative use case.
    #
    # \param copyParents When copying an AIM, should the same parents also be used (default - True).
    #
    # \remark If no $altName is provided and an AIM with the name, aim,
    # has already been loaded into the problem, an alternative
    # name will be automatically specified with the syntax $aim_[instance number]. If an
    # $altName is provided it must be unique compared to other instances of the loaded AIMs in the problem.
    def __init__(self, problemObject, **kwargs):

        # Check problem object
        if not isinstance(problemObject, capsProblem) and not isinstance(problemObject, CapsProblem):
            raise TypeError("Provided problemObject is not a capsProblem or CapsProblem class")

        ## Reference to the problem object that loaded the AIM during a call to \ref capsProblem.loadAIM .
        self.capsProblem = problemObject

        # Check the requested parent list against the analysis dictionary
        def checkAIMParents(parentList):
            for i in parentList:
                if i in self.capsProblem.analysis.keys():
                    pass
                else:
                    print("Requested parent", i, "has not been initialized/loaded")
                    return False

            return True

        # Check the analysis dictionary for redundant names - need unique names
        def checkAnalysisDict(name, printInfo=True):
            if name in self.capsProblem.analysis.keys():
                if printInfo:
                    print("\nAnalysis object with the name, " + str(name) +  ", has already been initilized")
                return True

            else:
                return False

        # Check the analysis dictionary for redundant directories if official aim names are the same
        def checkAnalysisDir(name, analysidDir):

            for aim in self.capsProblem.analysis.values():

                if name == aim.officialName and analysidDir == aim.analysisDir:
                    print("\nAnalysis object with the name, ", str(name),  ", and directory, ", str(analysidDir),
                          ", has already been initilized - Please use a different directory!!\n")
                    return True

            return False

        # Check to see if an AIM has been fully specified
        def checkAIMStatus(analysisDir, currentAIM, copyAIM):
            if copyAIM is not None:
                if analysisDir is None:
                    print("\nError: Unable to copy AIM, analysisDir (analysis directory) have not been set!\n")
                    return False
                else:
                    return True

            if analysisDir is None or currentAIM is None:
                print("\nError: Unable to load AIM, either the analysisDir (analysis directory) and/or aim variables have not been set!\n")
                return False
            else:
                return True

        if "capsFidelity" in kwargs.keys():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while loading aim"
                                                         + "\nIMPORTANT - The term capsFidelity is no longer supported please change to capsIntent in both your *.csm and *.py files")

        if self.capsProblem.geometry == None:
            raise caps.CAPSError(caps.CAPS_NULLOBJ, msg = "while loading aim"
                                                        + "\nError: A *.csm, *.caps, or *.egads file has not been loaded!")

        currentAIM = kwargs.pop("aim", None)
        capsIntent = kwargs.pop("capsIntent", None)
        analysisDir = kwargs.pop("analysisDir", None)
        parents = kwargs.pop("parents", [])
        copyAIM = kwargs.pop("copyAIM", None)
        copyParents = kwargs.pop("copyParent", True)
        altName = kwargs.pop("altName", currentAIM)

        __fromSave = kwargs.pop("__fromSave", False)

        if capsIntent:
            if not isinstance(capsIntent, list):
                capsIntent = [capsIntent]
            capsIntent = ';'.join(capsIntent)

        if not analysisDir or analysisDir == ".":
            print("No analysis directory provided - defaulting to altName")
            analysisDir = altName

        if parents:
            if not isinstance(parents, list):
                parents = [parents]

            temp = []
            for i in parents:

                if isinstance(i, capsAnalysis) or isinstance(i,CapsAnalysis):
                    temp.append(i.aimName)
                elif isinstance(i, str):
                    temp.append(i)
                else:
                    raise TypeError("Provided parent is not a capsAnalysis or CapsAnalysis class or a string name")

            # Update parents
            parents = temp[:]

            # Get parent index from anaylsis dictionary
            if not checkAIMParents(parents):
                raise caps.CAPSError(caps.CAPS_NULLOBJ, msg = "while loading aim"
                                                            + "\nUnable to find parents")


        if copyAIM:
            if copyAIM not in self.capsProblem.analysis.keys():
                raise caps.CAPSError(caps.CAPS_BADNAME, msg = "while loading aim"
                                                            + "\nError: AIM " + str(copyAIM) + " to duplicate has not been loaded into the problem!")
            else:
                currentAIM = self.capsProblem.analysis[copyAIM].officialName
                capsIntent = self.capsProblem.analysis[copyAIM].capsIntent

            if copyParents == True:
                parents = self.capsProblem.analysis[copyAIM].parents

        if not __fromSave:
            # Check aim status
            if not checkAIMStatus(analysisDir, currentAIM, copyAIM):
                raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while loading aim - " + currentAIM)

            # Check analysis dictionary for redunant names
            if checkAnalysisDict(altName):

                if altName == currentAIM:
                    i = 0
                    while i < 1000:
                        i +=1

                        altName = currentAIM + "_" + str(i)
                        if not checkAnalysisDict(altName, printInfo=False):
                            break
                    else:
                        raise caps.CAPSError(caps.CAPS_BADNAME, msg = "while loading aim - " + currentAIM
                                                                    + "Unable to auto-specify altName, are there more than 1000 AIM instances loaded?")

                    print("An altName of", altName, "will be used!")

                else:
                    raise caps.CAPSError(caps.CAPS_BADNAME, msg = "while loading aim - " + currentAIM
                                                                + "The altName," + altName + ", has already been used!")

            # Check analysis for redunant directories if aim has already been loaded
            if checkAnalysisDir(currentAIM, analysisDir):
                raise caps.CAPSError(caps.CAPS_DIRERR, msg = "while loading aim - " + currentAIM)

        self.__createObj(altName, currentAIM, analysisDir, capsIntent,
                         unitSystem=kwargs.pop("unitSystem", None),
                         parentList=parents,
                         copyAIM=copyAIM)

        # Load the aim
        self.capsProblem.analysis[altName] = self

    def __createObj(self,
                    name, officialName, analysisDir, capsIntent, # Use if creating an analysis from scratch, else set to None and call _setupObj
                    unitSystem=None,
                    parentList=None,
                    copyAIM=None):

        ## Reference name of AIM loaded for the analysis object (see \ref capsProblem.loadAIM$altName).
        self.aimName = name

        ## Name of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$aim).
        self.officialName = officialName

        if analysisDir is not None:
            ## Adjust legacy directory names
            analysisDir = analysisDir.replace(".","_")
            analysisDir = analysisDir.replace("/","_")
            analysisDir = analysisDir.replace("\\","_")
    
            ## Analysis directory of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$analysis).
            # If the directory does not exist it will be made automatically.
            self.analysisDir = analysisDir + "_" + officialName
        else:
            self.analysisDir = None

        ## Unit system the AIM was loaded with (if applicable).
        self.unitSystem = unitSystem

        ## List of parents of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$parents).
        self.parents = parentList

        ## Analysis intent of the AIM loaded for the analysis object (see \ref capsProblem.loadAIM$capsIntent
        self.capsIntent = capsIntent

        ## OpenMDAO "component" object (see \ref createOpenMDAOComponent for instantiation).
        self.openMDAOComponent = None

        # It is assumed we are going to call _setupObj
        if self.aimName is None:
            return

        if self.analysisDir:
            root = self.capsProblem.problemObj.getRootPath()
            fullPath = os.path.join(root,self.analysisDir)
            
            # Check to see if directory exists
            if not os.path.isdir(fullPath):
                #print("Analysis directory does not currently exist - it will be made automatically")
                os.makedirs(fullPath)
    
            # Remove the top directory as cCAPS will create it
            if os.path.isdir(fullPath):
                shutil.rmtree(fullPath)

        # Allocate a temporary list of analysis object for the parents
        if len(parentList) != 0:

            aobjTemp = [None]*len(self.parents)

            for i, parent in enumerate(self.parents):
                aobjTemp[i] = self.capsProblem.analysis[parent]._analysisObj()
        else:
            aobjTemp = None

        if copyAIM is not None:
            self.analysisObj = self.capsProblem.analysis[copyAIM].analysisObj.dupAnalysis(self.analysisDir)
        else:
            self.analysisObj = self.capsProblem.problemObj.makeAnalysis(officialName,
                                                                        self.analysisDir,
                                                                        unitSys=unitSystem,
                                                                        intent=capsIntent)

        if len(parentList) != 0:

            inNames = self.getAnalysisVal(namesOnly=True)

            inValObj = None

            if "Mesh" in inNames:
                inValObj = self.analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, "Mesh")
            elif "Surface_Mesh" in inNames:
                inValObj = self.analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, "Surface_Mesh")
            elif "Volume_Mesh" in inNames:
                inValObj = self.analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, "Volume_Mesh")

            if inValObj is not None:
                for parent in self.parents:
                    outNames = self.capsProblem.analysis[parent].getAnalysisOutVal(namesOnly=True)
                    
                    if "Surface_Mesh" in outNames:
                        outValObj = self.capsProblem.analysis[parent].analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISOUT, "Surface_Mesh")
                    elif "Area_Mesh" in outNames:
                        outValObj = self.capsProblem.analysis[parent].analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISOUT, "Area_Mesh")
                    elif "Volume_Mesh" in outNames:
                        outValObj = self.capsProblem.analysis[parent].analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISOUT, "Volume_Mesh")
                    else:
                        continue

                    inValObj.linkValue(outValObj, caps.tMethod.Copy)


        #if copyAIM is not None:
        #    raise caps.CAPSError(msg = "while loading aim - " + str(self.aimName)
        #                                     + " (during a call to caps_dupAnalysis)")
        #else:
        #    raise caps.CAPSError(msg = "while loading aim - " + str(self.aimName)
        #                                     + " (during a call to caps_makeAnalysis)")

        # Add attribute information for saving the problem
        self.addAttribute("__aimName"     , self.aimName)
        self.addAttribute("__officialName", self.officialName)
        if self.parents:
            self.addAttribute("__parents", self.parents)

        analysisInfo = self.getAnalysisInfo(printInfo=False, infoDict=True)

        self.analysisDir = analysisInfo["analysisDir"]

    def _setupObj(self, analysisObj):

        self.analysisObj = analysisObj

        self.aimName = self.getAttribute("__aimName")
        self.officialName = self.getAttribute("__officialName")
        try:
            self.parents = self.getAttribute("__parents")

        except caps.CAPSError as e:
            if e.errorName != "CAPS_NOTFOUND":
                raise

            self.parents = None

        analysisInfo = self.getAnalysisInfo(printInfo=False, infoDict=True)

        self.analysisDir = analysisInfo["analysisDir"]
        self.capsIntent = analysisInfo["intent"]
        self.unitSystem = analysisInfo["unitSystem"]

        return self

    def _analysisObj(self):
        return self.analysisObj

    ## Sets an ANALYSISIN variable for the analysis object.
    # \param varname Name of CAPS value to set in the AIM.
    # \param value Value to set. Type casting is automatically done based on the CAPS value type
    # \param units Applicable units of the current variable (default - None). Only applies to real values.
    # specified in the AIM.
    def setAnalysisVal(self, varname, value, units=None):

        tempObj = self.analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, varname)
        
        # This is a complete HACK, but consistent with previous poor use of units...
        if units is None:
            _value = tempObj.getValue()
            if isinstance(_value, caps.Quantity):
                units = _value._units

        tempObj.setValue( caps._withUnits(value, units) )

    ## Gets an ANALYSISIN variable for the analysis object.
    # \param varname Name of CAPS value to retrieve from the AIM. If no name is provided a dictionary
    # containing all ANALYSISIN values is returned. See example \ref analysis4.py for a representative
    # use case.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param units When set to True returns the units along with the value as specified by the analysis. When
    # set to a string (e.g. units="ft") the returned value is converted into the specified units. (default - False).
    #
    # \param namesOnly Return only a list of variable names (no values) if creating a dictionary (default - False).
    #
    # \return Value of "varname" or dictionary of all values. Units are also returned if applicable
    # based on the "units" keyword (does not apply to dictionary returns).
    def getAnalysisVal(self, varname = None, **kwargs):
        ## \example analysis4.py
        # Example use cases for interacting the
        # pyCAPS.capsAnalysis.getAnalysisVal() and
        # pyCAPS.capsAnalysis.getAnalysisOutVal() functions.

        analysisIn = {}

        namesOnly = kwargs.pop("namesOnly", False)
        units = kwargs.pop("units", None)

        if varname is not None: # Get just the 1 variable

            tempObj = self.analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, varname)
            value = tempObj.getValue()

            if value is None:
                if units:
                    if isinstance(units, str):
                        raise caps.CAPSError(caps.CAPS_UNITERR, msg="Analysis variable, " + varname + ", is None; unable to convert units to " + units )
                else:
                    return None
            else:
                if units:
                    if not isinstance(value, caps.Quantity):
                        raise caps.CAPSError(caps.CAPS_UNITERR, msg = "No units assigned to analysis variable, " + varname + ", unable to convert units to " + units )

                    return value.convert(units)._value
                else:
                    if isinstance(value, caps.Quantity):
                        return value._value
                    else:
                        return value

        else: # Build a dictionary

            numAnalysisInput = self.analysisObj.size(caps.oType.VALUE, caps.sType.ANALYSISIN)

            for i in range(1, numAnalysisInput+1):
                tempObj = self.analysisObj.childByIndex(caps.oType.VALUE, caps.sType.ANALYSISIN, i)

                name, otype, stype, link, parent, last = tempObj.info()

                if namesOnly is True:
                    analysisIn[name] = None
                else:
                    value = tempObj.getValue()
                    if isinstance(value,caps.Quantity):
                        value = value._value
                    analysisIn[name] = value

            if namesOnly is True:
                return analysisIn.keys()
            else:
                return analysisIn

    ## Gets an ANALYSISOUT variable for the analysis object.
    # \param varname Name of CAPS value to retrieve from the AIM. If no name is provided a dictionary
    # containing all ANALYSISOUT values is returned. See example \ref analysis4.py for a
    # representative use case.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param units When set to True returns the units along with the value as specified by the analysis. When
    # set to a string (e.g. units="ft") the returned value is converted into the specified units. (default - False).
    #
    # \param namesOnly Return only a list of variable names (no values) if creating a dictionary (default - False).
    #
    # \return Value of "varname" or dictionary of all values. Units are also returned if applicable
    # based on the "units" keyword (does not apply to dictionary returns).
    def getAnalysisOutVal(self, varname = None, **kwargs):
        ## \example analysis4.py
        # Example use cases for interacting the
        # pyCAPS.capsAnalysis.getAnalysisVal() and
        # pyCAPS.capsAnalysis.getAnalysisOutVal() functions.

        analysisOut = {}

        namesOnly = kwargs.pop("namesOnly", False)
        units = kwargs.pop("units", None)

        stateOfAnalysis = self.getAnalysisInfo(printInfo=False)

        if stateOfAnalysis != 0 and namesOnly == False:
            raise caps.CAPSError(caps.CAPS_DIRTY, msg = "while getting analysis out variables"
                                                      + "Analysis state isn't clean (state = " + str(stateOfAnalysis) + ") can only return names of ouput variables; set namesOnly keyword to True")

        if varname is not None:

            tempObj = self.analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISOUT, varname)
            value = tempObj.getValue()

            if value is None:
                if units:
                    if isinstance(units, str):
                        raise caps.CAPSError(caps.CAPS_UNITERR, msg="Analysis variable, " + varname + ", is None; unable to convert units to " + units )
                else:
                    return None
            else:
                if units:
                    if not isinstance(value, caps.Quantity):
                        raise caps.CAPSError(caps.CAPS_UNITERR, msg = "No units assigned to analysis variable, " + varname + ", unable to convert units to " + units )

                    return value.convert(units)._value
                else:
                    if isinstance(value, caps.Quantity):
                        return value._value
                    else:
                        return value

        else: # Build a dictionary

            numAnalysisOutput = self.analysisObj.size(caps.oType.VALUE, caps.sType.ANALYSISOUT)

            for i in range(1, numAnalysisOutput+1):

                tempObj = self.analysisObj.childByIndex(caps.oType.VALUE, caps.sType.ANALYSISOUT, i)

                name, otype, stype, link, parent, last = tempObj.info()
                if namesOnly is True:
                    analysisOut[name] = None
                else:
                    value = tempObj.getValue()
                    if isinstance(value,caps.Quantity):
                        value = value._value
                    analysisOut[name] = value

            if namesOnly is True:
                return analysisOut.keys()
            else:
                return analysisOut

    ## Gets analysis information for the analysis object. See example \ref analysis6.py for
    # a representative use case.
    #
    # \param printInfo Print information to sceen if True.
    # \param **kwargs See below.
    #
    # \return Cleanliness state of analysis object or a dictionary containing analysis
    # information (infoDict must be set to True)
    #
    # Valid keywords:
    # \param infoDict Return a dictionary containing analysis information instead
    # of just the cleanliness state (default - False)
    def getAnalysisInfo(self, printInfo=True, **kwargs):
        ## \example analysis6.py
        # Example use cases for interacting the
        # pyCAPS.capsAnalysis.getAnalysisInfo() function.

        dirtyInfo = {0 : "Up to date",
                     1 : "Dirty analysis inputs",
                     2 : "Dirty geometry inputs",
                     3 : "Both analysis and geometry inputs are dirty",
                     4 : "New geometry",
                     5 : "Post analysis required",
                     6 : "Execution and Post analysis required"}

        analysisPath, unitSystem, major, minor, capsIntent, fnames, franks, fInOut, execution, cleanliness = self.analysisObj.analysisInfo()

        if printInfo:
            numField = len(fnames)
            print("Analysis Info:")
            print("\tName           = ", self.aimName)
            print("\tIntent         = ", capsIntent)
            print("\tAnalysisDir    = ", analysisPath)
            print("\tVersion        = ", str(major) + "." + str(minor))
            print("\tExecution Flag = ", execution)
            print("\tNumber of Fields = ", numField)

            if numField > 0:
                print("\tField name (Rank):")

            for i in range(numField):
                print("\t ", fnames[i], "(", franks[i], ")")

            if cleanliness in dirtyInfo:
                print("\tDirty state    = ", cleanliness, " - ", dirtyInfo[cleanliness] )
            else:
                print("\tDirty state    = ", cleanliness, " ",)

        if execution == 1:
            executionFlag = True
        else:
            executionFlag = False

        if kwargs.pop("infoDict", False):
            return {"name"          : self.aimName,
                    "intent"        : capsIntent,
                    "analysisDir"   : analysisPath,
                    "version"       : str(major) + "." + str(minor),
                    "executionFlag" : executionFlag,
                    "fieldName"     : fnames,
                    "fieldRank"     : franks,
                    "unitSystem"    : unitSystem,
                    "status"        : cleanliness
                    }

        return cleanliness

    ## Retrieve a list of geometric attribute values of a given name
    # ("attributeName") for the bodies loaded
    # into the analysis. Level in which to search the bodies is
    # determined by the attrLevel keyword argument.
    # See \ref analysis3.py for a representative use case.
    #
    # \param attributeName Name of attribute to retrieve values for.
    # \param **kwargs See below.
    #
    # \return A list of attribute values.
    #
    # Valid keywords:
    # \param bodyIndex Specific body in which to retrieve attribute information from.
    # \param attrLevel Level to which to search the body(ies). Options: \n
    #                  0 (or "Body") - search just body attributes\n
    #                  1 (or "Face") - search the body and all the faces [default]\n
    #                  2 (or "Edge") - search the body, faces, and all the edges\n
    #                  3 (or "Node") - search the body, faces, edges, and all the nodes
    #
    def getAttributeVal(self, attributeName, **kwargs):
        ## \example analysis3.py
        # Example use case of the pyCAPS.capsAnalysis.getAttributeVal() function

        # Initiate attribute value list
        attributeList = []

        bodyIndex = kwargs.pop("bodyIndex", None)
        attrLevel = kwargs.pop("attrLevel", "face")

        attribute = {"body" : 0,
                     "face" : 1,
                     "edge" : 2,
                     "node" : 3}

        if not isinstance(attrLevel, (int,float)):
            if attrLevel.lower() in attribute.keys():
                attrLevel = attribute[attrLevel.lower()]

        if int(attrLevel) not in attribute.values():
            print("Error: Invalid attribute level! Defaulting to 'Face'")
            attrLevel = attribute["face"]

        # How many bodies do we have
        numBody = self.analysisObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while getting attribute " + str(attributeName) + " values for aim "
                                                         + str(self.aimName)
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i

            if bodyIndex is not None:
                if int(bodyIndex) != i+1:
                    continue

            body, units = self.analysisObj.bodyByIndex(i+1)

            attributeList += _createAttributeValList(body,
                                                     attrLevel,
                                                     attributeName)

        try:
            return list(set(attributeList))
        except TypeError:
            print("Unhashable type in attributelist")
            return attributeList

    ## Create geometric attribution map (embeded dictionaries) for the bodies loaded
    # into the analysis.
    # Dictionary layout:
    # - Body 1
    #   - Body : Body level attributes
    #   - Faces
    #     - 1 : Attributes on the first face of the body
    #     - 2 : Attributes on the second face of the body
    #     - " : ...
    #   - Edges
    #     - 1 : Attributes on the first edge of the body
    #     - 2 : Attributes on the second edge of the body
    #     - " : ...
    #   - Nodes :
    #     - 1 : Attributes on the first node of the body
    #     - 2 : Attributes on the second node of the body
    #     - " : ...
    # - Body 2
    #   - Body : Body level attributes
    #   - Faces
    #     - 1 : Attributes on the first face of the body
    #     - " : ...
    #   - ...
    # - ...
    #
    # Dictionary layout (reverseMap = True):
    # - Body 1
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not
    #       - Faces : Face numbers at which the attribute exist
    #       - Edges : Edge numbers at which the attribute exist
    #       - Nodes : Node numbers at which the attribute exist
    #     -  Value : Next value of attribute with the same name
    #       - Body : True if value exist at body level, None if not
    #       -  "   : ...
    #     - ...
    #   - Atribute : Attribute name
    #     - Value : Value of attribute
    #       -   "  : ...
    #     - ...
    # - Body 2
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not
    #       -   "   : ...
    #     - ...
    #   - ...
    # - ...
    #
    # \param getInternal Get internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param reverseMap Reverse the attribute map (default - False). See above table for details.
    #
    # \return Dictionary containing attribution map
    def getAttributeMap(self, getInternal = False, **kwargs):

        attributeMap = {}

        # How many bodies do we have
        numBody = self.analysisObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating attribute map"
                                                         + " for aim " + str(self.aimName)
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i

            body, units = self.analysisObj.bodyByIndex(i+1)

            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            attributeMap["Body " + str(i + 1)] = _createAttributeMap(body, getInternal)

        if kwargs.pop("reverseMap", False):
            return _reverseAttributeMap(attributeMap)
        else:
            return attributeMap


    ## Alternative to \ref preAnalysis(). Warning: May be deprecated in future versions.
    def aimPreAnalysis(self):
        print("\n\nWarning: aimPreAnalysis is a deprecated. Please use preAnalysis.\n",
                   "Further use will result in an error in future releases!!!!!!\n")
        self.preAnalysis()

    ## Run the pre-analysis function for the AIM. If the specified analysis directory doesn't exist it
    # will be made automatically.
    def preAnalysis(self):
        self.analysisObj.preAnalysis()

    ## Alternative to \ref postAnalysis(). Warning: May be deprecated in future versions.
    def aimPostAnalysis(self):
        print("\n\nWarning: aimPostAnalysis is a deprecated. Please use postAnalysis.\n",
                   "Further use will result in an error in future releases!!!!!!\n")
        self.postAnalysis()

    ## Run the post-analysis function for the AIM.
    def postAnalysis(self):
        self.analysisObj.postAnalysis()

    ## Create an OpenMDAO Component[1.7.3]/ExplicitComponent[2.8+] object; an external code component
    # (ExternalCode[1.7.2]/ExternalCodeComp[2.8+]) is created if the
    # executeCommand keyword arguement is provided. This functionality should work with either
    # verison 1.7.3 or >=2.8 of OpenMDAO.
    #
    # \param inputVariable Input variable(s)/parameter(s) to add to the OpenMDAO component.
    # Variables may be either
    # analysis input variables or geometry design parameters. Note, that the setting of analysis inputs
    # supersedes the setting of geometry design parameters; issues may arise if analysis input and
    # geometry design variables have the same name.
    # If the analysis parameter wanting to be added to the OpenMDAO component is part of a
    # capsTuple the following notation should be used:
    #    "AnalysisInput:TupleKey:DictionaryKey", for example "AVL_Control:ControlSurfaceA:deflectionAngle"
    # would correspond to the AVL_Control input variable, the ControlSurfaceA element of the input values
    #  (that is the name of the control surface being created) and finally deflectionAngle corresponds
    # to the name of the dictionary entry that is to be used as the component parameter.
    # If the tuple's value isn't a dictionary just "AnalysisInput:TupleKey" is needed.
    #
    # \param outputVariable Output variable(s)/parameter(s) to add to the OpenMDAO component. Only scalar output variables
    # are currently supported
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param changeDir Automatically switch into the analysis directory set for the AIM when executing an
    # external code (default - True).
    #
    # \param saveIteration If the generated OpenMDAO component is going to be called multiple times,
    # the inputs and outputs from the analysis and the AIM will be
    # automatically bookkept ( = True) by moving the files to a folder within the AIM's analysis directory
    # ( \ref analysisDir ) named "Iteration_#" were # represents the iteration number (default - False). By default ( = False)
    # input and output files will be continously overwritten.
    # Notes:
    # - If the AIM has 'parents' their genertated files will not be bookkept.
    # - If previous iteration folders already exist, the iteration folders and any other files in the directory will be moved to a folder named
    # "Instance_#".
    # - This bookkeeping method will likely fail if the iterations are run concurrently!
    #
    # \param executeCommand Command to be executed when running an external code. Command must be a list of command
    # line arguements (see OpenMDAO documentation). If provided an ExternalCode[1.7.2]/ExternalCodeComp[2.8+]
    # object is created; if not provided or set to None a Component[1.7.3]/ExplicitComponent[2.8+] object is created (default - None).
    #
    # \param inputFile Optional list of input file names for OpenMDAO to check the existence of before
    # OpenMDAO excutes the "solve_nonlinear"[1.7.3]/"compute"[2.8+] (default - None). This is redundant as the AIM automatically
    # does this already.
    #
    # \param outputFile Optional list of output names for OpenMDAO to check the existence of before
    # OpenMDAO excutes the "solve_nonlinear"[1.7.3]/"compute"[2.8+]  (default - None). This is redundant as the AIM automatically
    # does this already.
    #
    # \param stdin Set I/O connection for the standard input of an ExternalCode[1.7.2]/ExternalCodeComp[2.8+] component.
    # The use of this depends on the expected AIM execution.
    #
    # \param stdout Set I/O connection for the standard ouput of an ExternalCode[1.7.2]/ExternalCodeComp[2.8+] component.
    # The use of this depends on the expected AIM execution.
    #
    # \param setSensitivity Optional dictionary containing sensitivity/derivative settings/parameters. Currently only
    # Finite difference is is supported!. See OpenMDAO
    # documentation for addtional details of "deriv_options"(version 1.7) or "declare_partials"(version 2.8).
    # Common values for a finite difference calculation would
    # be setSensitivity['type'] = "fd" (Note in the version 2.8 documentation this varibale has been changed to "method"
    # both variations will work when using version 2.8+), setSensitivity['form'] =  "forward" or "backward" or "central", and
    # setSensitivity['step_size'] = 1.0E-6 (Note in the version 2.8 documentation this varibale has been changed to "step"
    # both variations will work when using version 2.8+).
    #
    # \return Optionally returns the reference to the OpenMDAO component object created; a reference to the component is stored
    # in \ref openMDAOComponent . "None" is returned if a failure occurred during object creation.
    def createOpenMDAOComponent(self, inputVariable, outputVariable, **kwargs):

        print("Creating an OpenMDAO component for analysis, " + str(self.aimName))

        self.openMDAOComponent = createOpenMDAOComponent(self,
                                                         inputVariable,
                                                         outputVariable,
                                                         kwargs.pop("executeCommand", None),
                                                         kwargs.pop("inputFile", None),
                                                         kwargs.pop("outputFile", None),
                                                         kwargs.pop("setSensitivity", {}),
                                                         kwargs.pop("changeDir", True),
                                                         kwargs.pop("saveIteration", False))
        if isinstance(self.openMDAOComponent, int):
            raise caps.CAPSError(self.openMDAOComponent, msg = "while getting creating OpenMDAO component for analysis - "
                                                             + str(self.aimName))
        else:

            if kwargs.get("stdin", None) is not None:
                self.openMDAOComponent.stdin = kwargs.get("stdin", None)

            if kwargs.get("stdout", None) is not None:
                self.openMDAOComponent.stdout = kwargs.get("stdout", None)

            return self.openMDAOComponent

    ## Create a HTML dendrogram/tree of the current state of the analysis. See
    # example \ref analysis5.py for a representative use case. The HTML file relies on
    # the open-source JavaScript library, D3, to visualize the data.
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file.
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet
    # access will be necessary to view the tree.
    #
    # \param filename Filename to use when saving the tree (default - "aimName"). Note an ".html" is
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default - True)? If set to True internet access
    # will be necessary to view the tree.
    #
    # \param analysisGeom Show the geometry currently load into the analysis in the tree (default - False).
    #
    # \param internalGeomAttr Show the internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False). Note: "analysisGeom" must also be set to True.
    #
    # \param reverseMap Reverse the attribute map (default - False). See \ref getAttributeMap for details.
    def createTree(self, filename = "aimName", **kwargs):
        ## \example analysis5.py
        # Example use cases for the pyCAPS.capsAnalysis.createTree() function.

        if filename == "aimName":
            filename = self.aimName

        embedJSON = kwargs.pop("embedJSON", True)

        if embedJSON is False:
            try:
                # Write the JSON data file
                fp = open(filename + ".json",'w')
                fp.write( createTree(self,
                                     kwargs.pop("analysisGeom", False),
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                    )
                         )
                fp.close()
            except IOError:
                raise caps.CAPSError(caps.CAPS_IOERR, msg = "while creating HTML tree")

            # Write the HTML file
            writeTreeHTML(filename +".html",
                          filename + ".json",
                          None,
                          kwargs.pop("internetAccess", True)
                         )
        else:
            writeTreeHTML(filename+".html",
                          None,
                          createTree(self,
                                     kwargs.pop("analysisGeom", False),
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                     ),
                          kwargs.pop("internetAccess", True)
                         )

    ## Add an attribute (that is meta-data) to the analysis object. See example
    # \ref analysis7.py for a representative use case.
    #
    # \param name Name used to define the attribute.
    #
    # \param data Data value(s) for the attribute. Note that type casting in done
    # automatically based on the determined type of the Python object.
    def addAttribute(self, name, data):
        ## \example analysis7.py
        # Example use cases for interacting the
        # pyCAPS.capsAnalysis.addAttribute() and
        # pyCAPS.capsAnalysis.getAttribute() functions.

        valueObj = self.capsProblem.problemObj.makeValue(name, caps.sType.USER, data)

        self.analysisObj.setAttr(valueObj)

    ## Get an attribute (that is meta-data) that exists on the analysis object. See example
    # \ref analysis7.py for a representative use case.
    # \param name Name of attribute to retrieve.
    # \return Value of attribute "name".
    def getAttribute(self, name):

        valueObj = self.analysisObj.attrByName(name)
        value = valueObj.getValue()

        if isinstance(value,caps.Quantity):
            value = value._value

        return value

    ## Get sensitivity values. Note the AIM must have aimBackdoor function to interact with.
    # \param inputVar Input variable to retrieve the sensitivity with respect to.
    # \param outputVar Output variable to retrieve the sensitivity value for.
    # \return Sensitivity value.
    #
    # Note for AIM developers:
    #
    # This function makes use of the caps_AIMbackdoor function. The JSONin variable provided is a JSON dictionary with
    # the following form - {"mode": "Sensitivity", "outputVar": outputVar, "inputVar": inputVar}. Similarly, a JSON dictionary
    # is expected to return (JSONout) with the following form - {"sensitivity": sensitivityVal}; furthermore it is
    # assumed that JSON string returned (JSONout) is freeable.
    def getSensitivity(self, inputVar, outputVar):

        JSONin = jsonDumps({"mode": "Sensitivity",
                            "outputVar": str(outputVar),
                            "inputVar" : str(inputVar)})

        return self.backDoor(JSONin)

    ## Alternative to \ref backDoor(). Warning: May be deprecated in future versions.
    def aimBackDoor(self, JSONin):
        print("\n\nWarning: aimBackDoor is a deprecated. Please use backDoor.\n",
                   "Further use will result in an error in future releases!!!!!!\n")
        self.backDoor(JSONin)

    ## Make a call to the AIM backdoor function. Important: it is assumed that the JSON string returned (JSONout) by the AIM
    # is freeable.
    #
    # \param JSONin JSON string input to the backdoor function. If the value isn't a string (e.g. a Python dictionary)
    # it will automatically be converted to JSON string.
    #
    # \return Returns a Python object of the converted JSON out string from the back door function if any.
    def backDoor(self, JSONin):

        if not isinstance(JSONin, str):
            JSONin = jsonDumps(JSONin)

        JSONout = self.analysisObj.AIMbackdoor(JSONin)

        if JSONout:
            temp =  jsonLoads( JSONout )
            return float(temp["sensitivity"])

        return None


    ## View the geometry associated with the analysis.
    #  If the analysis produces a surface tessellation, then that is shown.
    #  Otherwise the bodies are shown with default tessellation parameters.
    #  Note that the geometry must be built and will not autoamtically be built by this function.
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param portNumber Port number to start the server listening on (default - 7681).
    def viewGeometry(self, **kwargs):

        print("viewGeometry is temporarily dissabled...")
        return

        # Initialize the viewer class
        viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))

        # Get the bodies from the analysis
        bodies = self.analysisObj.getBodies()

        if bodies is None:
            print("No bodies available for display.")
            return

        #viewTess = False
        #for i in range(numBody):
        #    if bodies[i+numBody] != <egads.ego> NULL:
        #        viewTess = True
        #        break

        #if viewTess:
        #    print("Viewing analysis tessellation...")
        #    # Show the tessellation used by the aim
        #    viewer.addEgadsTess(numBody, bodies+numBody)
        #else:
        # Show the bodies
        print("Viewing analysis bodies...")
        viewer.addEgadsBody(numBody, bodies)

        # Start up server
        viewer.startServer()

        del viewer

    ## Save the current geometry used by the AIM to a file.
    # \param filename File name to use when saving geometry file.
    # \param directory Directory where to save file. Default current working directory.
    # \param extension Extension type for file if filename does not contain an extension.
    def saveGeometry(self, filename="myGeometry", directory=os.getcwd(), extension=".egads"):

        filename = self.capsProblem._createFilename(filename, directory, extension)

        self.analysisObj.writeGeometry(filename)

    ## Get bounding boxes for the geometry used by the AIM.
    #
    # \return Returns a dictionary of the bounding boxes, {"BodyName", [ x_min, y_min, z_min, x_max, y_max, z_max]}.
    def getBoundingBox(self):

        boundingBox = {}
        nameList = []

        bodies = self.analysisObj.getBodies()

        if bodies is None:
            raise caps.CAPSError(caps.CAPS_SOURCEERR, msg = "No bodies available to retrieve the bounding box for!")

        for i in range(len(bodies)):
            body = bodies[i]
            box  = body.getBoundingBox()

            name, nameList = self.capsProblem._getBodyName(i, nameList, body)

            boundingBox[name] = [box[0][0], box[0][1], box[0][2],
                                 box[1][0], box[1][1], box[1][2]]

        return boundingBox

## Defines a CAPS geometry object.
# May be created directly by providing a problem object or indirectly through \ref capsProblem.loadCAPS .
class capsGeometry:

    ## Initializes a geometry object. See \ref capsProblem.loadCAPS for additional information and caveats.
    #
    # \param problemObject capsProblem object in which the geometry object will be added to.
    #
    # \param capsFile CAPS file to load. Options: *.csm, *.caps, or *.egads. If the filename has a *.caps extension
    # the pyCAPS analysis, bound, and value objects will be re-populated (see remarks in \ref capsProblem.loadCAPS ).
    #
    # \param projectName CAPS project name. projectName=basename(__main__.__file__) without extension if not provided.
    def __init__(self, problemObject, capsFile, projectName=None, verbosityLevel=1):

        # Check problem object
        if not isinstance(problemObject, capsProblem) and not isinstance(problemObject, CapsProblem):
            raise TypeError

        self.capsProblem = problemObject

        if self.capsProblem.geometry is not None:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while loading file - " + str(capsFile) + ". Can NOT load multiple files into a problem."
                                                         + "\nWarning: Can NOT load multiple files into a problem. A CAPS file has already been loaded!\n")

        ## Geometry file loaded into problem. Note that the directory path has been removed.
        self.geomName = os.path.basename(capsFile)

        if projectName is None:
            import __main__
            try:
                base = os.path.basename(__main__.__file__)
                base = os.path.splitext(base)[0]
            except AttributeError:
                base = "__main__"
            iProj = 0
            projectName = base
            while os.path.exists( os.path.join(projectName, "Scratch", "capsLock") ):
                iProj += 1
                projectName = base + str(iProj)

        verbosity = {"minimal"  : 0,
                     "standard" : 1,
                     "debug"    : 2}

        if not isinstance(verbosityLevel, (int,float)):
            if verbosityLevel.lower() in verbosity.keys():
                verbosityLevel = verbosity[verbosityLevel.lower()]

        if int(verbosityLevel) not in verbosity.values():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while setting verbosity, invalid verbosity level!")

        self.capsProblem.problemObj = caps.open(projectName, None, capsFile, verbosityLevel)

        self.capsProblem.geometry = self

        if ".caps" in capsFile:
            self.capsProblem._repopulateProblem()

    ## Manually force a build of the geometry. The geometry will only be rebuilt if a
    # design parameter (see \ref setGeometryVal) has been changed.
    def buildGeometry(self):
        self.capsProblem.problemObj.preAnalysis()

    ## Save the current geometry to a file.
    # \param filename File name to use when saving geometry file.
    # \param directory Directory where to save file. Default current working directory.
    # \param extension Extension type for file if filename does not contain an extension.
    def saveGeometry(self, filename="myGeometry", directory=os.getcwd(), extension=".egads"):

        filename = self.capsProblem._createFilename(filename, directory, extension)

        self.capsProblem.problemObj.writeGeometry(filename)

    ## Sets a GEOMETRYIN variable for the geometry object.
    # \param varname Name of geometry design parameter to set in the geometry.
    # \param value Value of geometry design parameter. Type casting is automatically done based
    # on value indicated by CAPS.
    def setGeometryVal(self, varname, value=None):

        if varname is None:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "No variable name provided")

        valueObj = self.capsProblem.problemObj.childByName(caps.oType.VALUE, caps.sType.GEOMETRYIN, varname)
        valueObj.setValue(value)

    ## Gets a GEOMETRYIN variable for the geometry object.
    # \param varname Name of geometry design parameter to retrieve from the geometry.
    # If no name is provided a dictionary containing all design parameters is returned.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param namesOnly Return only a list of parameter names (no values) if creating a
    # dictionary (default - False).
    #
    # \return Value of "varname" or dictionary of all values.
    def getGeometryVal(self, varname=None, **kwargs):

        geometryIn = {}

        namesOnly = kwargs.pop("namesOnly", False)

        if varname is not None: # Get just the 1 variable

            valueObj = self.capsProblem.problemObj.childByName(caps.oType.VALUE, caps.sType.GEOMETRYIN, varname)

            value = valueObj.getValue()
            if isinstance(value,caps.Quantity):
                value = value._value
            return value

        else: # Build a dictionary

            numGeometryInput = self.capsProblem.problemObj.size(caps.oType.VALUE, caps.sType.GEOMETRYIN)

            for i in range(1, numGeometryInput+1):

                tempObj = self.capsProblem.problemObj.childByIndex(caps.oType.VALUE, caps.sType.GEOMETRYIN, i)

                name, otype, stype, link, parent, last = tempObj.info()

                if namesOnly:
                    geometryIn[name] = None
                else:
                    value = tempObj.getValue()
                    if isinstance(value,caps.Quantity):
                        value = value._value
                    geometryIn[name] = value

            if namesOnly is True:
                return geometryIn.keys()
            else:
                return geometryIn

    ## Gets a GEOMETRYOUT variable for the geometry object.
    # \param varname Name of geometry (local) parameter to retrieve from the geometry.
    # If no name is provided a dictionary containing all local variables is returned.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param namesOnly Return only a list of parameter names (no values) if creating a dictionary (default - False).
    # \param ignoreAt Ignore @ geometry variables when creating a dictionary (default - True).
    #
    # \return Value of "varname" or dictionary of all values.
    def getGeometryOutVal(self, varname=None, **kwargs):

        geometryOut = {}

        namesOnly = kwargs.pop("namesOnly", False)
        ignoreAt  = kwargs.pop("ignoreAt", True)

        if varname is not None: # Get just the 1 variable

            valueObj = self.capsProblem.problemObj.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, varname)

            value = valueObj.getValue()
            if isinstance(value,caps.Quantity):
                value = value._value
            return value

        else: # Build a dictionary

            numGeometryOutput = self.capsProblem.problemObj.size(caps.oType.VALUE, caps.sType.GEOMETRYOUT)

            for i in range(1, numGeometryOutput+1):

                tempObj = self.capsProblem.problemObj.childByIndex(caps.oType.VALUE, caps.sType.GEOMETRYOUT, i)

                name, otype, stype, link, parent, last = tempObj.info()

                varname = name
                if ignoreAt:
                    if "@" in varname:
                        continue

                if namesOnly:
                    geometryOut[name] = None
                else:
                    value = tempObj.getValue()
                    if isinstance(value,caps.Quantity):
                        value = value._value
                    geometryOut[name] = value

            if namesOnly is True:
                return geometryOut.keys()
            else:
                return geometryOut

    ## View or take a screen shot of the geometry configuration. The use of this function to save geometry requires the \b matplotlib module.
    # \b <em>Important</em>: If both showImage = True and filename is not None, any manual view changes made by the user
    # in the displayed image will be reflected in the saved image.
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param viewerType What viewer should be used (default - "capsViewer").
    #                   Options: "capsViewer" or "matplotlib" (options are case insensitive).
    #                   Important: if $filename is not None, the viewer is changed to matplotlib.
    # \param portNumber Port number to start the server listening on (default - 7681).
    # \param title Title to add to each figure (default - None).
    # \param filename Save image(s) to file specified (default - None). Note filename should
    # not contain '.' other than to indicate file type extension (default type = *.png).
    # 'file' - OK, 'file2.0Test' - BAD, 'file2_0Test.png' - OK, 'file2.0Test.jpg' - BAD.
    # \param directory Directory path were to save file. If the directory doesn't exist
    # it will be made. (default - current directory).
    # \param viewType  Type of view for the image(s). Options: "isometric" (default),
    # "fourview", "top" (or "-zaxis"), "bottom" (or "+zaxis"), "right" (or "+yaxis"),
    # "left" (or "-yaxis"), "front" (or "+xaxis"), "back" (or "-xaxis").
    # \param combineBodies Combine all bodies into a single image (default - False).
    # \param ignoreBndBox Ignore the largest body (default - False).
    # \param showImage Show image(s) (default - False).
    # \param showAxes Show the xyz axes in the image(s) (default - False).
    # \param showTess Show the edges of the tessellation (default - False).
    # \param dpi Resolution in dots-per-inch for the figure (default - None).
    # \param tessParam Custom tessellation paremeters, see EGADS documentation for makeTessBody function.
    # values will be scaled by the norm of the bounding box for the body (default - [0.0250, 0.0010, 15.0]).
    def viewGeometry(self, **kwargs):

        viewDefault = "capsviewer"

        viewer = kwargs.pop("viewerType", viewDefault).lower()
        filename = kwargs.get("filename", None)

        # How many bodies do we have
        numBody = self.capsProblem.problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while saving viewing geometry"
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        bodies = [None]*numBody
        for i in range(numBody):
            bodies[i], units = self.capsProblem.problemObj.bodyByIndex(i+1)

        # Are we going to use the viewer?
        if filename is None and viewer.lower() == "capsviewer":
            self._viewGeometryCAPSViewer(**kwargs)
        else:
            _viewGeometryMatplotlib(bodies, **kwargs)

    def _viewGeometryCAPSViewer(self, **kwargs):

        print("viewGeometry with capsViewer is temporarily dissabled...")
        return

        # Initialize the viewer class
        viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))

        # Get the bodies
        bodies = self.capsProblem.problemObj.getBodies()

        # Load the bodies into the viewer
        viewer.addEgadsBody(bodies)

        # Start up server
        viewer.startServer()

    ## Retrieve a list of attribute values of a given name ("attributeName") for the bodies in the current geometry.
    # Level in which to search the bodies is
    # determined by the attrLevel keyword argument.
    #
    # \param attributeName Name of attribute to retrieve values for.
    # \param **kwargs See below.
    #
    # \return A list of attribute values.
    #
    # Valid keywords:
    # \param bodyIndex Specific body in which to retrieve attribute information from.
    # \param attrLevel Level to which to search the body(ies). Options: \n
    #                  0 (or "Body") - search just body attributes\n
    #                  1 (or "Face") - search the body and all the faces [default]\n
    #                  2 (or "Edge") - search the body, faces, and all the edges\n
    #                  3 (or "Node") - search the body, faces, edges, and all the nodes
    #
    def getAttributeVal(self, attributeName, **kwargs):

        # Initiate attribute value list
        attributeList = []

        bodyIndex = kwargs.pop("bodyIndex", None)
        attrLevel = kwargs.pop("attrLevel", "face")

        attribute = {"body" : 0,
                     "face" : 1,
                     "edge" : 2,
                     "node" : 3}

        if not isinstance(attrLevel, (int,float)):
            if attrLevel.lower() in attribute.keys():
                attrLevel = attribute[attrLevel.lower()]

        if int(attrLevel) not in attribute.values():
            print("Error: Invalid attribute level! Defaulting to 'Face'")
            attrLevel = attribute["face"]

        # How many bodies do we have
        numBody = self.capsProblem.problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while getting attribute " + str(attributeName) + " values on geometry"
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i

            if bodyIndex is not None:
                if int(bodyIndex) != i+1:
                    continue

            body, units = self.capsProblem.problemObj.bodyByIndex(i+1)

            attributeList += _createAttributeValList(body,
                                                     attrLevel,
                                                     attributeName)

        try:
            return list(set(attributeList))
        except TypeError:
            print("Unhashable type in attributelist")
            return attributeList

    ## Create attribution map (embeded dictionaries) of each body in the current geometry.
    # Dictionary layout:
    # - Body 1
    #   - Body : Body level attributes
    #   - Faces
    #     - 1 : Attributes on the first face of the body
    #     - 2 : Attributes on the second face of the body
    #     - " : ...
    #   - Edges
    #     - 1 : Attributes on the first edge of the body
    #     - 2 : Attributes on the second edge of the body
    #     - " : ...
    #   - Nodes :
    #     - 1 : Attributes on the first node of the body
    #     - 2 : Attributes on the second node of the body
    #     - " : ...
    # - Body 2
    #   - Body : Body level attributes
    #   - Faces
    #     - 1 : Attributes on the first face of the body
    #     - " : ...
    #   - ...
    # - ...
    #
    # Dictionary layout (reverseMap = True):
    # - Body 1
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not
    #       - Faces : Face numbers at which the attribute exist
    #       - Edges : Edge numbers at which the attribute exist
    #       - Nodes : Node numbers at which the attribute exist
    #     -  Value : Next value of attribute with the same name
    #       - Body : True if value exist at body level, None if not
    #       -  "   : ...
    #     - ...
    #   - Atribute : Attribute name
    #     - Value : Value of attribute
    #       -   "  : ...
    #     - ...
    # - Body 2
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not
    #       -   "   : ...
    #     - ...
    #   - ...
    # - ...
    #
    # \param getInternal Get internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param reverseMap Reverse the attribute map (default - False). See above table for details.
    #
    # \return Dictionary containing attribution map
    def getAttributeMap(self, getInternal = False, **kwargs):

        attributeMap = {}

        # How many bodies do we have
        numBody = self.capsProblem.problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating attribute map"
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        for bodyIndex in range(numBody, 0, -1):

            # No plus +1 on bodyIndex since we are going backwards
            body, units = self.capsProblem.problemObj.bodyByIndex(bodyIndex)

            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            attributeMap["Body " + str(bodyIndex)] = _createAttributeMap(body, getInternal)

        if kwargs.pop("reverseMap", False):
            return _reverseAttributeMap(attributeMap)
        else:
            return attributeMap

    ## Create a HTML dendrogram/tree of the current state of the geometry. The HTML file relies on
    # the open-source JavaScript library, D3, to visualize the data.
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file.
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet
    # access will be necessary to view the tree.
    #
    # \param filename Filename to use when saving the tree (default - "myGeometry"). Note an ".html" is
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default - True)? If set to True internet access
    # will be necessary to view the tree.
    #
    # \param internalGeomAttr Show the internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    #
    # \param reverseMap Reverse the attribute map (default - False). See \ref getAttributeMap for details.
    def createTree(self, filename = "myGeometry", **kwargs):

        embedJSON = kwargs.pop("embedJSON", True)

        if embedJSON is False:
            try:
                # Write the JSON data file
                fp = open(filename + ".json",'w')
                fp.write(createTree(self,
                                    False,
                                    kwargs.pop("internalGeomAttr", False),
                                    kwargs.pop("reverseMap", False)
                                    )
                         )
                fp.close()
            except IOError:
                raise caps.CAPSError(caps.CAPS_IOERR, msg = "while creating HTML tree")

            # Write the HTML file
            writeTreeHTML(filename +".html",
                          filename + ".json",
                          None,
                          kwargs.pop("internetAccess", True))
        else:
            writeTreeHTML(filename+".html",
                          None,
                          createTree(self,
                                     False,
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                     ),
                          kwargs.pop("internetAccess", True)
                         )

    ## Get bounding boxes for all geometric bodies.
    #
    # \return Returns a dictionary of the bounding boxes, {"BodyName", [ x_min, y_min, z_min, x_max, y_max, z_max]}.
    def getBoundingBox(self):

        boundingBox = {}
        nameList = []

        numBody = self.capsProblem.problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "The number of bodies in the problem is less than 1!!!" )

        # Make copies of the bodies to hand to the model
        for i in range(numBody):
            body, units = self.capsProblem.problemObj.bodyByIndex(i+1)
            box = body.getBoundingBox()

            name, nameList = self.capsProblem._getBodyName(i, nameList, body)

            boundingBox[name] = [box[0][0], box[0][1], box[0][2],
                                 box[1][0], box[1][1], box[1][2]]

        return boundingBox

## Defines a CAPS bound object.
# May be created directly by providing a problem object or indirectly through \ref capsProblem.createDataBound or
# \ref capsProblem.createDataTransfer .
class capsBound:

#     @property
#     def dateSet(self):
#         return self.dataSetSrc

    ## Initializes a bound object. See \ref capsProblem.createDataBound for additional information and caveats.
    #
    # \param problemObject capsProblem object in which the analysis object will be added to.
    #
    # \param **kwargs See below.
   #
    # Valid keywords:
    # \param capsBound Name of capsBound to use for the data bound.
    # \param variableName Single or list of variables names to add.
    # \param aimSrc Single or list of AIM names that will be the data sources for the bound.
    # \param aimDest Single or list of AIM names that will be the data destinations during the transfer.
    # \param transferMethod Single or list of transfer methods to use during the transfer.
    # \param initValueDest Single or list of initial values for the destication data.
    def __init__(self, problemObject, **kwargs):

        # Check problem object
        if not isinstance(problemObject, capsProblem) and not isinstance(problemObject, CapsProblem):
            raise TypeError("Provided problemObject is not a capsProblem or CapsProblem class")

        ## Reference to the problem object that loaded the AIM during a call to \ref capsProblem.loadAIM .
        self.capsProblem = problemObject

        # Look for wild card fieldname and compare variable
        def checkFieldWildCard(variable, fieldName):
            #variable = _byteify(variable)
            for i in fieldName:
                try:
                    i_ind = i.index("#")
                except ValueError:
                    i_ind = -1

                try:
                    if variable[:i_ind] == i[:i_ind]:
                        return True
                except:
                    pass

            return False

        def checkAnalysisName(aimList):
            for i in aimList:
                if i not in self.capsProblem.analysis.keys():
                    raise caps.CAPSError(caps.CAPS_BADVALUE, msg="AIM name = " + i + " not found in analysis dictionary")

        if "capsTransfer" in kwargs.keys():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="The keyword, capsTransfer, should be changed to capsBound")

        if "variableRank" in kwargs.keys():
            print("WARNING: Keyword variableRank is not necessary. Please remove it from your script. May raise an error in future releases.")


        capsBound = kwargs.pop("capsBound", None)
        varName = kwargs.pop("variableName", None)
        aimSrc = kwargs.pop("aimSrc", None)
        aimDest = kwargs.pop("aimDest", [])
        transferMethod = kwargs.pop("transferMethod", [])

        initValueDest = kwargs.pop("initValueDest", None)

        if not capsBound:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="No transfer name (capsBound) provided")

        if not varName:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="No variable name (variableName) provided")

        if not aimSrc:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="No source AIM(s) (aimSrc) provided")

        if not isinstance(varName, list):
            varName = [varName] # Convert to list

        numDataSet = len(varName)

        if not isinstance(aimSrc, list):
            aimSrc = [aimSrc]

        temp = []
        for i in aimSrc:
            if isinstance(i, capsAnalysis) or isinstance(i,CapsAnalysis):
                temp.append(i.aimName)
            elif isinstance(i, str):
                temp.append(i)
            else:
                raise TypeError("Provided aimSrc is not a capsAnalysis or CapsAnalysis class or a string name")
        aimSrc = temp[:]

        checkAnalysisName(aimSrc)

        if len(aimSrc) != numDataSet:
            raise caps.CAPSError(caps.CAPS_MISMATCH, msg="aimSrc dimensional mismatch between inputs!")

        # Check the field names in source analysis
        for var, src in zip(varName, aimSrc):

            fieldName = self.capsProblem.analysis[src].getAnalysisInfo(printInfo = False, infoDict = True)["fieldName"]

            if var not in fieldName:
                if not checkFieldWildCard(var, fieldName):
                    raise caps.CAPSError(caps.CAPS_BADVALUE, msg="Variable (field name) " + var + " not found in aimSrc "+ src)

        if aimDest:
            if not isinstance(aimDest, list):
                aimDest = [aimDest]

            temp = []
            for i in aimDest:
                if isinstance(i, capsAnalysis) or isinstance(i,CapsAnalysis):
                    temp.append(i.aimName)
                elif isinstance(i, str):
                    temp.append(i)
                else:
                    raise TypeError("Provided aimDest is not a capsAnalysis or CapsAnalysis class or a string name")

            aimDest = temp[:]

            checkAnalysisName(aimDest)

            if len(aimDest) != numDataSet:
                raise caps.CAPSError(caps.CAPS_MISMATCH, msg="aimDest dimensional mismatch between inputs!")

            if transferMethod:

                if not isinstance(transferMethod, list):
                    transferMethod = [transferMethod]

                temp = []
                for i in transferMethod:
                    if "conserve" == i.lower():
                        temp.append(caps.dMethod.Conserve)

                    elif "interpolate" == i.lower():
                        temp.append(caps.dMethod.Interpolate)
                    else:
                        print("Unreconized data transfer method - defaulting to \"Interpolate\"")
                        temp.append(caps.dMethod.Interpolate)

                transferMethod = temp[:]

            else:
                print("No data transfer method (transferMethod) provided - defaulting to \"Interpolate\"")
                transferMethod =  [caps.dMethod.Interpolate]*numDataSet

            if len(transferMethod) != numDataSet:
                    raise caps.CAPSError(caps.CAPS_MISMATCH, msg="transferMethod dimensional mismatch between inputs!")

        if initValueDest:
            if not isinstance(initValueDest, list):
                initValueDest = [initValueDest]

            if len(initValueDest) != numDataSet:
                raise caps.CAPSError(caps.CAPS_MISMATCH, msg="initValueDest dimensional mismatch between inputs!")
        else:
            if aimDest:
                initValueDest = [None]*numDataSet
            else:
                initValueDest = []

        # Create the object
        self.__createObj(capsBound, varName, aimSrc, aimDest, transferMethod, initValueDest)

        # Set the bound in the problem
        self.capsProblem.dataBound[capsBound] = self

    def __createObj(self, boundName, varName, aimSrc, aimDest, transferMethod, initValueDest):

        ## Bound/transfer name used to set up the data bound.
        self.boundName = boundName

        ## List of variables in the bound.
        self.variables = varName

        ## Dictionary of vertex set object (\ref pyCAPS._capsVertexSet) in the bound.
        self.vertexSet = {}

        ## Dictionary of "source" data set objects (\ref pyCAPS._capsDataSet) in the bound.
        self.dataSetSrc = {}

        ## Dictionary of "destination" data set objects (\ref pyCAPS._capsDataSet) in the bound.
        self.dataSetDest = {}

        if aimDest == None:
            aimDest = []
        if transferMethod == None:
            transferMethod = []
        if initValueDest == None:
            initValueDest = [None]*len(aimDest)

        # Open bound
        self.boundObj = self.capsProblem.problemObj.makeBound(2, boundName)

        # Make all of our vertex sets
        for i in set(aimSrc+aimDest):
            #print"Vertex Set - ", i
            self.vertexSet[i] = _capsVertexSet(i,
                                               self.capsProblem,
                                               self,
                                               self.capsProblem.analysis[i])


        if len(varName) != len(aimSrc):
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="ERROR: varName and aimSrc must all be the same length")

        # Make all of our source data sets
        for var, src in zip(varName, aimSrc):
            #print"DataSetSrc - ", var, src

            rank = self.__getFieldRank(src, var)

            self.dataSetSrc[var] = _capsDataSet(var,
                                                self.capsProblem,
                                                self,
                                                self.vertexSet[src],
                                                caps.fType.FieldOut,
                                                rank)
        # Make destination data set
        if len(aimDest) > 0:

            if len(varName) != len(aimDest) or \
               len(varName) != len(transferMethod) or \
               len(varName) != len(aimSrc) or \
               len(varName) != len(initValueDest):
                raise caps.CAPSError(caps.CAPS_BADVALUE, msg="ERROR: varName, aimDest, transferMethod, aimSrc, and initValueDest must all be the same length")

            for var, dest, method, src, initVal in zip(varName, aimDest, transferMethod, aimSrc, initValueDest):
                #print"DataSetDest - ", var, dest, method

                rank = self.__getFieldRank(src, var)

                self.dataSetDest[var] = _capsDataSet(var,
                                                     self.capsProblem,
                                                     self,
                                                     self.vertexSet[dest],
                                                     method,
                                                     rank,
                                                     initVal)

        # Close bound
        self.boundObj.closeBound()

    def __boundObj(self):
        return self.boundObj

    # Get the field rank based on source and field name
    def __getFieldRank(self, src, var):

        # Look for wild card fieldname and compare variable
        def __checkFieldWildCard(variable, fieldName):
            #variable = _byteify(variable)
            for i in fieldName:
                try:
                    i_ind = i.index("#")
                except ValueError:
                    i_ind = -1

                try:
                    if variable[:i_ind] == i[:i_ind]:
                        return fieldName.index(i)
                except:
                    pass

            return -1

        fieldName = self.capsProblem.analysis[src].getAnalysisInfo(printInfo=False, infoDict=True)["fieldName"]

        try:

           index = fieldName.index(var)
        except:
            index  = __checkFieldWildCard(var, fieldName)
            if  index < 0:
                raise caps.CAPSError(caps.CAPS_BADVALUE, msg="ERROR: Variable (field name)" + str(var) + "not found in aimSrc" + str(src) )

        return self.capsProblem.analysis[src].getAnalysisInfo(printInfo=False, infoDict=True)["fieldRank"][index]

    ## Gets bound information for the bound object.
    #
    # \param printInfo Print information to sceen if True.
    # \param **kwargs See below.
    #
    # \return State of bound object.
    #
    # Valid keywords:
    # \param infoDict Return a dictionary containing bound information instead
    # of just the state (default - False)
    def getBoundInfo(self, printInfo=True, **kwargs):

        boundInfo = {-2 : "multiple BRep entities - Error in reparameterization!y",
                     -1 : "Open",
                      0 : "Empty and Closed",
                      1 : "Single BRep entity",
                      2 : "Multiple BRep entities"}

        state, dimension, limits = self.boundObj.boundInfo()

        if printInfo:
            print("Bound Info:")
            print("\tName           = ", self.boundName)
            print("\tDimension      = ", dimension)
            # Add something for limits

            if state in boundInfo:
                print("\tState          = ", state, " -> ", boundInfo[state] )
            else:
                print("\tState          = ", state, " ",)

        if kwargs.pop("infoDict", False):
            return {"name"     : self.boundName,
                    "dimension": dimension,
                    "status"   : state
                    }

        return state

    ## Populates VertexSets for the bound. Must be called to finalize
    # the bound after all mesh generation aim's have been executed
    def fillVertexSets(self):
        pass

    ## Execute data transfer for the bound.
    # \param variableName Name of variable to implement the data transfer for. If no name is provided
    # the first variable in bound is used.
    def executeTransfer(self, variableName = None):

        if not self.dataSetDest:
            print("No destination data exists in bound!")
            return

        self.getDestData(variableName)


    def getDestData(self, variableName = None):

        if not self.dataSetDest:
            print("No destination data exists in bound!")
            return

        if variableName is None:
            print("No variable name provided - defaulting to variable - " + self.variables[0])
            variableName = self.variables[0]

        elif variableName not in self.variables:
            print("Unreconized varible name " + str(variableName))
            return

        self.dataSetDest[variableName].getData()

    ## Visualize data in the bound.
    # The function currently relies on matplotlib or the capViewer class (webviewer) to plot the data.
    # \param variableName Name of variable to visualize. If no name is provided
    # the first variable in the bound is used.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param viewerType What viewer should be used (default - capsViewer). Options: capsViewer or matplotlib (options are
    # case insensitive). Important: if $filename isn't set to None, the default to changed to matplotlib.
    # \param portNumber Port number to start the server listening on (default - 7681).
    # \param filename Save image(s) to file specified (default - None). Not available when using the webviewer
    # \param colorMap  Valid string for a, matplotlib::cm,  colormap (default - 'Blues'). Not as options are
    # available when using the webviewer (see \ref capsViewer for additional details).
    # \param showImage Show image(s) (default - True).
    def viewData(self, variableName=None, **kwargs):

        filename = kwargs.pop("filename", None)
        if (filename is not None and "." not in filename):
            filename += ".png"

        viewer = kwargs.pop("viewerType", "capsViewer").lower()

        if variableName is None:
            print("No variable name provided - defaulting to variable - " + self.variables[0])
            variableName = self.variables[0]

        elif variableName not in self.variables:
            raise caps.CAPSError(caps.CAPS_NOTFOUND, msg = "Unrecognized varible name " + str(variableName))

        # Are we going to use the viewer?
        if filename is None and viewer == "capsviewer":

            # Initialize the viewer class
            viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))

            # Add src data set to viewer
            viewer.addDataSet(self.dataSetSrc[variableName])

            # Add dest data set to viewer
            if self.dataSetDest:
                viewer.addDataSet(self.dataSetDest[variableName])

            # Start up server
            viewer.startServer()

            del viewer

        else:

            try:
                import matplotlib.pyplot as plt
            except:
                print("Error: Unable to import matplotlib - viewing the data is not possible")
                raise ImportError



            showImage = kwargs.pop("showImage", True)
            kwargs["showImage"] =  False

            kwargs["filename"] = None

            fig = plt.figure() #figsize=(12,15)


            if self.dataSetDest:

                self.dataSetSrc[variableName].viewData(fig,  2, 0, **kwargs)
                self.dataSetDest[variableName].viewData(fig, 2, 1, **kwargs)

            else:
                self.dataSetSrc[variableName].viewData(fig,  1, 0, **kwargs)

            if filename is not None:
                print("Saving figure - ", filename)
                fig.savefig(filename)

            if showImage:
                plt.show()

            # Cleanup
            fig.clf()
            del fig

    ## Write a Tecplot compatiable file for the data in the bound.
    # \param filename Name of file to save data to.
    # \param variableName Single or list of variables to write data for. If no name is provided
    # all variables in the bound are used.
    def writeTecplot(self, filename, variableName=None):

        if filename is None:
            print("No filename provided")
            return

        if ("." not in filename):
            filename += ".dat"

        try:
            file = open(filename, "w")
        except:
            raise caps.CAPSError(caps.CAPS_IOERR, msg = "while open file for writing Tecplot file for data set - "
                                                      + str(self.dataSetName))

        if variableName is None:
            print("No variable name provided - all variables will be writen to file")
            variableName = self.variables

        if not isinstance(variableName,list):
            variableName = [variableName]

        for var in variableName:
            if var not in self.variables:
                print("Unreconized varible name " + str(var))
                continue

            self.dataSetSrc[var].writeTecplot(file)

            if self.dataSetDest:

                self.dataSetDest[var].writeTecplot(file)

        # Close the file
        file.close()

    ## Create a HTML dendrogram/tree of the current state of the bound. The HTML file relies on
    # the open-source JavaScript library, D3, to visualize the data.
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file.
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet
    # access will be necessary to view the tree.
    #
    # \param filename Filename to use when saving the tree (default - "boundName"). Note an ".html" is
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default True)? If set to True internet access
    # will be necessary to view the tree.
    def createTree(self, filename = "boundName", **kwargs):

        if filename == "boundName":
            filename = self.boundName

        embedJSON = kwargs.pop("embedJSON", True)

        if embedJSON is False:
            try:
                # Write JSON data file
                fp = open(filename + ".json",'w')
                fp.write(createTree(self, False, False, False))
                fp.close()
            except IOError:
                raise caps.CAPSError(caps.CAPS_IOERR, msg = "while creating HTML tree")

            # Write the HTML file
            writeTreeHTML(filename +".html",
                          filename + ".json",
                          None,
                          kwargs.pop("internetAccess", True))
        else:
            writeTreeHTML(filename+".html",
                          None,
                          createTree(self, False, False, False),
                          kwargs.pop("internetAccess", True))

## Functions to interact with a CAPS vertexSet object.
# Should be initiated within \ref pyCAPS.capsBound (not a standalone class)
class _capsVertexSet:

    def __init__(self, name, problemParent, boundParent, analysisParent):

        ## Vertex set name (analysis name the vertex belongs to).
        self.vertexSetName = name

        # Problem object the vertex set belongs to.
        self.capsProblem  = problemParent

        ## Bound object (\ref pyCAPS.capsBound)the vertex set belongs to.
        self.capsBound    = boundParent

        ## Analysis object (\ref pyCAPS.capsAnalysis) the vertex set belongs to.
        self.capsAnalysis = analysisParent

        self.vertexSetObj = self.capsBound.boundObj.makeVertexSet(self.capsAnalysis.analysisObj)

    def __vertexSetObj(self):
        return self.vertexSetObj

## Functions to interact with a CAPS dataSet object.
# Should be initiated within \ref pyCAPS.capsBound (not a standalone class)
class _capsDataSet:

    def __init__(self, name, problemParent, boundParent, vertexSetParent, dataMethod, dataRank, initValue = None):

        ## Data set name (variable name).
        self.dataSetName = name

        # Reference to the problem object that data set pertains to.
        self.capsProblem = problemParent

        ## Reference to the bound object (\ref pyCAPS.capsBound) that data set pertains to.
        self.capsBound = boundParent

        ## Reference to the vertex set object (\ref pyCAPS._capsVertexSet) that data set pertains to.
        self.capsVertexSet = vertexSetParent

        ## Data method: Analysis, Interpolate, Conserve.
        self.dataSetMethod = None

        self.dataMethod = dataMethod

        if self.dataMethod == caps.dMethod.Conserve:
            self.dataSetMethod = "Conserve"
        elif self.dataMethod == caps.dMethod.Interpolate:
            self.dataSetMethod = "Interpolate"
        elif self.dataMethod == caps.dMethod.Analysis:
            self.dataSetMethod = "Analysis"

        ## Rank of data set.
        self.dataRank = dataRank

        self.dataSetObj = self.capsVertexSet.vertexSetObj.makeDataSet(self.dataSetName, self.dataMethod, self.dataRank)

        # Set an initial constant value for the data set
        if initValue is not None:
            self.initDataSet(initValue)

    def __dataSetObj(self):
        return self.dataSetObj

    ## Executes caps_initDataSet on data set object to set an inital constant value needed for cyclic data transfers.
    def initDataSet(self, initValue):
        self.dataSetObj.initDataSet(initValue)

    ## Executes caps_getData on data set object to retrieve data set variable, \ref dataSetName.
    # \return Optionally returns a list of data values. Data with a rank greater than 1 returns a list of lists (e.g.
    # data representing a displacement would return [ [Node1_xDisplacement, Node1_yDisplacement, Node1_zDisplacement],
    # [Node2_xDisplacement, Node2_yDisplacement, Node2_zDisplacement], etc. ]
    def getData(self):
        data = self.dataSetObj.getData()
        if isinstance(data,caps.Quantity):
            data = data._value
        return data

    ## Executes caps_getData on data set object to retrieve XYZ coordinates of the data set.
    # \return Optionally returns a list of lists of x,y, z values (e.g. [ [x2, y2, z2], [x2, y2, z2],
    # [x3, y3, z3], etc. ] )
    def getDataXYZ(self):
        tempDataObj = self.capsVertexSet.vertexSetObj.childByName(caps.oType.DATASET, caps.sType.NONE, "xyz")
        xyz = tempDataObj.getData()
        if isinstance(xyz,caps.Quantity):
             xyz = xyz._value
        return xyz

    ## Executes caps_getTriangles on data set's vertex set to retrieve the connectivity (triangles only) information
    # for the data set.
    # \return Optionally returns a list of lists of connectivity values
    # (e.g. [ [node1, node2, node3], [node2, node3, node7], etc. ] ) and a list of lists of data connectivity (not this is
    # an empty list if the data is node-based) (eg. [ [node1, node2, node3], [node2, node3, node7], etc. ]
    def getDataConnect(self):
        triConn, dataConn = self.capsVertexSet.vertexSetObj.triangulate()
        return (triConn, dataConn)

    ## Visualize data set.
    # The function currently relies on matplotlib to plot the data.
    # \param fig Figure object (matplotlib::figure) to append image to.
    # \param numDataSet Number of data sets in \ref _capsDataSet.viewData$fig.
    # \param dataSetIndex Index of data set being added to \ref _capsDataSet.viewData$fig.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param filename Save image(s) to file specified (default - None).
    # \param colorMap  Valid string for a, matplotlib::cm,  colormap (default - 'Blues').
    # \param showImage Show image(s) (default - True).
    # \param title Set a custom title on the plot (default - VertexSet= 'name', DataSet = 'name', (Var. '#') ).
    def viewData(self, fig=None, numDataSet=1, dataSetIndex=0, **kwargs):

        #viewer = kwargs.pop("viewerType", "capsViewer").lower()

        # Are we going to use the viewer?
        #if viewer == "capsviewer":

            # Initialize the viewer class
            #viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))

            # Add src data set to viewer
            #viewer.addDataSet(self)

            # Start up server
            #viewer.startServer()

            #del viewer

        #else:
            #self._viewDataMatplotLib(self,fig,numDataSet,dataSetIndex,**kwargs)

    #def _viewDataMatplotLib(self, fig=None, numDataSet=1, dataSetIndex=0, **kwargs):

        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection

            from matplotlib import colorbar
            from matplotlib import colors
            from matplotlib import cm as colorMap
        except:
            print("Error: Unable to import matplotlib - viewing the data is not possible")
            raise ImportError

        filename = kwargs.pop("filename", None)
        if (filename is not None and "." not in filename):
            filename += ".png"

        showImage =  kwargs.pop("showImage", True)
        title = kwargs.pop("title", None)

        if (filename is None and fig is None and showImage is False):
            print("No 'filename' or figure instance provided and 'showImage' is false. Nothing to do here!")
            return

        cMap = kwargs.pop("colorMap", "Blues")

        try:
            colorMap.get_cmap(cMap)
        except ValueError:
            print("Colormap ",  cMap, "is not recognized. Defaulting to 'Blues'")
            cMap = "Blues"

        # Initialize values
        colorArray = []
        triConn = NULL
        dataConn = NULL

        self.getData()
        self.getDataXYZ()

        triConn, dataConn = self.capsVertexSet.vertexSetObj.triangulate()

        if not ((self.dataNumPoint == self.numNode) or (self.dataNumPoint == numTri)):
            raise caps.CAPSError(caps.CAPS_MISMATCH, msg = "viewData only supports node-centered or cell-centered data!")

        if numData != 0:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "viewData does NOT currently support cell-centered based data sets!")

        if fig is None:
            fig = plt.figure() #figsize=(9,9)

        cellCenter = (self.dataNumPoint == numTri)

        ax = []
        for j in range(self.dataRank):

            if self.dataRank > 1 or numDataSet > 1:
                ax.append(fig.add_subplot(self.dataRank, numDataSet, numDataSet*j + 1 + dataSetIndex, projection='3d'))
            else:
                ax.append(fig.gca(projection='3d'))

            colorArray = []
            for i in range(self.dataNumPoint):
                colorArray.append(self.dataValue[self.dataRank*i + j])

            xLimits = [9E6, -9E6] # [minX, maxX] - initial values are arbitary
            yLimits = [9E6, -9E6]
            zLimits = [9E6, -9E6]

            minColor = min(colorArray)
            maxColor = max(colorArray)
            dColor   = maxColor - minColor

            p = 0.0
            if (maxColor == 0.0 and minColor == 0.0):
                zeroP = True
            else:
                zeroP = False

            for i in range (numTri):
                tri = []

                index = triConn[i][0]-1
                tri.append([self.xyz[index][0], self.xyz[index][1], self.xyz[index][2]])

                if zeroP == False and not cellCenter:
                    p = (colorArray[index] - minColor)/dColor

                index = triConn[3*i+1]-1
                tri.append([self.xyz[index][0], self.xyz[index][1], self.xyz[index][2]])

                if zeroP == False and not cellCenter:
                    p += (colorArray[index] - minColor)/dColor

                index = triConn[3*i+2]-1
                tri.append([self.xyz[index][0], self.xyz[index][1], self.xyz[index][2]])

                if zeroP == False and not cellCenter:
                    p += (colorArray[index] - minColor)/dColor

                if cellCenter:
                    p = 255*(colorArray[i] - minColor)/dColor
                else:
                    p = p/3.0*255

                xLimits[0] = min(xLimits[0], tri[-1][0], tri[-2][0], tri[-3][0])
                xLimits[1] = max(xLimits[1], tri[-1][0], tri[-2][0], tri[-3][0])

                yLimits[0] = min(yLimits[0], tri[-1][1], tri[-2][1], tri[-3][1])
                yLimits[1] = max(yLimits[1], tri[-1][1], tri[-2][1], tri[-3][1])

                zLimits[0] = min(zLimits[0], tri[-1][2], tri[-2][2], tri[-3][2])
                zLimits[1] = max(zLimits[1], tri[-1][2], tri[-2][2], tri[-3][2])

                tri = Poly3DCollection([tri])

                facecolor = colorMap.get_cmap(cMap)(int(p))

                tri.set_facecolor(facecolor)
                tri.set_edgecolor(facecolor)
                ax[-1].add_collection3d(tri)

            #ax[-1].autoscale_view()
            if not title:
                title = "VertexSet = " + str(self.capsVertexSet.vertexSetName) + \
                        "\nDataSet = " + str(self.dataSetName)
                if self.dataRank > 1:
                    title += " (Var. " + str(j+1) + ")"

            ax[-1].set_title(title, size="small")
            ax[-1].set_axis_off() # Turn the axis off
            #ax[-1].invert_xaxis() # Invert x- axis
            ax[-1].view_init(azim=120, elev=30) # Isometric

            center = []
            center.append((xLimits[1] + xLimits[0])/2.0)
            center.append((yLimits[1] + yLimits[0])/2.0)
            center.append((zLimits[1] + zLimits[0])/2.0)

            scale = max(abs(xLimits[1] - xLimits[0]),
                        abs(yLimits[1] - yLimits[0]),
                        abs(zLimits[1] - zLimits[0]))

            scale = 0.6*(scale/2.0)

            ax[-1].set_xlim(-scale + center[0], scale + center[0])
            ax[-1].set_ylim(-scale + center[1], scale + center[1])
            ax[-1].set_zlim(-scale + center[2], scale + center[2])

            cax, kw = colorbar.make_axes(ax[-1], location="right")
            norm = colors.Normalize(vmin=min(colorArray), vmax=max(colorArray))
            colorbar.ColorbarBase(cax, cmap=cMap, norm=norm)

            del colorArray[:]

        if triConn:
            egads.EG_free(triConn)

        if dataConn:
            egads.EG_free(dataConn)

        #fig.subplots_adjust(hspace=0.1) #left=0.0, bottom=0.0, right=1.0, top=.93, wspace=0.0 )

        if showImage:
            plt.show()

        if filename is not None:
            print("Saving figure - ", filename)
            fig.savefig(filename)


    ## Write data set to a Tecplot compatible data file. A triagulation of the data set will be used
    # for the connectivity.
    # \param file Optional open file object to append data to. If not provided a
    # filename must be given via the keyword arguement \ref _capsDataSet.writeTecplot$filename.
    # \param filename Write Tecplot file with the specified name.
    def writeTecplot(self, file=None,  filename=None):

        # If no file
        if file == None:

            if filename == None:
                raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while writing Tecplot file for data set - " + str(self.dataSetName)
                                                             + " (during a call to caps_triangulate)"
                                                             + "\nNo file name or open file object provided" )

            if "." not in filename:
                filename += ".dat"

            try:
                fp = open(filename, "w")
            except:
                self.capsProblem.status = caps.CAPS_IOERR
                raise caps.CAPSError(msg = "while open file for writing Tecplot file for data set - "
                                             + str(self.dataSetName))
        else:
            fp = file

        data = self.getData()
        if not data: # Do we have empty data
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "the data set is empty for - "
                                                         + str(self.dataSetName))
        else:
            if not isinstance(data[0], list): # Only returns a list if rank 1
                temp = []
                for d in data:
                    temp.append([d])
                data = temp

        xyz = self.getDataXYZ()

        connectivity, dconnectivity = self.getDataConnect()

        #TODO: dconnectivity does not appear to indicate a cell centered data set?
        if dconnectivity:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "writeTecplot does NOT currently support dconnectivity based data sets!")

        if not ((len(data) == len(xyz)) or (len(data) == len(connectivity))):
            raise caps.CAPSError(caps.CAPS_MISMATCH, msg = "writeTecplot only supports node-centered or cell-centered data!")

        title = 'TITLE = "VertexSet: ' + str(self.capsVertexSet.vertexSetName) + \
                ', DataSet = ' + str(self.dataSetName) + '"\n'
        fp.write(title)

        title = 'VARIABLES = "x", "y", "z"'
        for i in range(self.dataRank):
            title += ', "Variable_' + str(i) + '"'
        title += "\n"
        fp.write(title)

        #numNode = self.numNode if len(data) == len(xyz) else 3*len(connectivity)
        numNode = len(xyz)

        title = "ZONE T = " + \
                '"VertexSet: ' + str(self.capsVertexSet.vertexSetName) + \
                ', DataSet = ' + str(self.dataSetName) + '"' +\
                 ", N = " + str(numNode) + \
                 ", E = " + str(len(connectivity)) + \
                 ", DATAPACKING = POINT, ZONETYPE= FETRIANGLE\n"
        fp.write(title)

        if len(data) == len(xyz): # Node centered data set
            # Write coordinates + data values
            for i, j in zip(xyz, data):

                dataString = str(i[0]) +" "+ str(i[1]) +" "+ str(i[2])
                for k in j:
                    dataString += " " + str(k)

                fp.write(dataString + "\n")

            # Write connectivity
            for i in connectivity:

                dataString = str(i[0]) +" "+ str(i[1]) +" "+ str(i[2])

                fp.write(dataString + "\n")

        elif len(data) == len(connectivity): # Cell centered data set
            # Write coordinates + data values
            for i, j in zip(connectivity, data):

                for t in i:
                    dataString = str(xyz[t-1][0]) +" "+ str(xyz[t-1][1]) +" "+ str(xyz[t-1][2])
                    for k in j:
                        dataString += " " + str(k)

                    fp.write(dataString + "\n")

            # Write connectivity
            for i in xrange(len(connectivity)):

                dataString = str(3*i+1) +" "+ str(3*i+2) +" "+ str(3*i+3)

                fp.write(dataString + "\n")

        if file == None:
            fp.close()

## Functions to interact with a CAPS value object.
# Should be created with \ref capsProblem.createValue (not a standalone class)
class _capsValue(object):

    def __init__(self, problemParent,
                  name, value, # Use if creating a value from scratch, else set to None and call _setupObj
                  subType="parameter",
                  units=None, limits=None,
                  fixedLength=True, fixedShape=True):

        ## Reference to the problem object that loaded the value.
        self.capsProblem = problemParent

        # If the CAPS value object is provided we don't need to make it!
        if name is not None or value is not None :

            ## Variable name
            self.name = name

            ## Units of the variable
            self.units = units

            # Save subtype for quick access
            if subType.lower() == "parameter":
                self.subType = caps.sType.PARAMETER
            elif subType.lower() == "user":
                self.subType = caps.sType.USER
            else:
                raise caps.CAPSError(caps.CAPS_BADTYPE, msg = "while creating value object - " + str(self.name)
                                                            + "\nUnreconized subType for value object!")

            self.valueObj = self.capsProblem.problemObj.makeValue(self.name, self.subType, caps._withUnits(value, units) )

            ## Acceptable limits for the value. Limits may be set directly.
            #See \ref value3.py for a representative use case.
            self.limits = limits

            ## Value of the variable. Value may be set directly.
            #See \ref value2.py for a representative use case.
            self.value = value # redundant - here just for the documentation purposes

            # Get shape parameters
            dim, pmtr, lfixed, sfixed, ntype = self.valueObj.getValueProps()

            # Possibly modify parameters
            if fixedLength:
                lfixed = caps.Fixed.Fixed
            else:
                lfixed = caps.Fixed.Change

            if fixedShape:
                sfixed = caps.Fixed.Fixed
            else:
                sfixed = caps.Fixed.Change

            # Set the new parameters
            self.valueObj.setValueProps(dim, lfixed, sfixed, ntype)

    def _setupObj(self, valueObj):

        self.valueObj = valueObj

        self.name, objType, self.subType, link, parent, owner = self.valueObj.info()

        # Get rest of value
        self.value

        return self

    def _valueObj(self):
        return self.valueObj

    @property
    def value(self):

        self._value = self.valueObj.getValue()
        if isinstance(self._value, caps.Quantity):
            self._value = self._value._value

        if isinstance(self._value, str) and ";" in self._value:
            self._value = self._value.split(";")

        return self._value

    @value.setter
    def value(self, dataValue):

        self._value = dataValue
        if not isinstance(dataValue, caps.Quantity):
            dataValue = caps._withUnits(dataValue, self.units)

        self.valueObj.setValue(dataValue)

    @property
    def limits(self):
        self._limits = self.valueObj.getLimits()

        if isinstance(self._limits, caps.Quantity):
            self._limits = self._limits._value

        return self._limits

    @limits.setter
    def limits(self, limitsValue):

        self._limits = limitsValue
        if not isinstance(limitsValue, caps.Quantity):
            limitsValue = caps._withUnits(limitsValue, self.units)

        self.valueObj.setLimits(limitsValue)

    ## \example value2.py
    # Example use cases for interacting the a pyCAPS._capsValue object for setting the value.

    ## \example value3.py
    # Example use cases for interacting the a pyCAPS._capsValue object for setting the limits.


    ## Change the value of the object. See \ref value2.py for a representative use case.
    # \param data Data value(s) for the variable. Note that data will be type casted to match the type
    # used to original create the capsValue object.
    def setVal(self, data):
        self.value = data

    ## Get the current value of object. See \ref value2.py for a representative use case.
    # \return Current value of set for object.
    def getVal(self):
        return self.value

    ## Set new limits. See \ref value3.py for a representative use case.
    # \param newLimits New values to set for the limits. Should be 2 element list - [min value, max value].
    def setLimits(self, newLimits):
        self.limits = newLimits

    ## Get the current value for the limits. See \ref value3.py for a representative use case.
    # \return Current value for the limits.
    def getLimits(self):
        return self.limits

    ## Return the current value of the object in the desired, specified units. Note that this neither
    # changes the value or units of the object, only returns a converted value.
    # See \ref value4.py for a representative use case.
    # \return Current value of the object in the specified units.
    def convertUnits(self, toUnits):
        ## \example value4.py
        # Example use cases for interacting the with pyCAPS._capsValue.convertUnits() function.

        if self.units == None:
            raise caps.CAPSError(caps.CAPS_UNITERR, msg = "while converting units for value object - " + str(self.name)
                                                        + "\nValue doesn't have defined units")

        return capsConvert(self.value, self.units, toUnits)

## Pythonic alias to capsProblem
class CapsProblem(capsProblem):
    pass

## Pythonic alias to capsAnalysis
class CapsAnalysis(capsAnalysis):
    pass

## Pythonic alias to capsGeometry
class CapsGeometry(capsGeometry):
    pass

## Pythonic alias to capsBound
class CapsBound(capsBound):
    pass
