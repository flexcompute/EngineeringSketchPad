
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
# \section clearancepyCAPS Clearance Statement

from . import caps
from .attributeUtils import _createAttributeValList, _createAttributeMap, _reverseAttributeMap
from .treeUtils import writeTreeHTML
from .geometryUtils import _viewGeometryMatplotlib
from .openmdao import createOpenMDAOComponent

from pyEGADS import egads
import os
import functools
import re
import shutil
from json import dumps as jsonDumps
from json import loads as jsonLoads

from sys import version_info as pyVersion

# Swith for dissabling deprecation warnings
LEGACY = False

#==============================================================================
_deprecated = {}
class PrintDeprecated:
    def __init__(self, msg):
        self.msg = msg
    def __del__(self):
        if not LEGACY:
            print(self.msg)

def deprecated(newname):
    """Add a deprecation warning"""
    def decorator_deprecated(func):
        @functools.wraps(func)
        def wrapper_deprecated(*args, **kwargs):
            if pyVersion.major > 2:
                funcname = func.__qualname__
            else:
                funcname = func.__name__
            if funcname not in _deprecated:
                _deprecated[funcname] = PrintDeprecated(("WARNING: {!r} is deprecated. Please use {}!").format(funcname, newname))
            return func(*args, **kwargs)
        return wrapper_deprecated
    return decorator_deprecated

def deprecate(oldname, newname):
    if oldname not in _deprecated:
        _deprecated[oldname] = PrintDeprecated(("WARNING: {!r} is deprecated. Please use {}!").format(oldname, newname))

#==============================================================================
# Sanitize a filename
def _createFilename(filename, directory, extension):

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

#==============================================================================
def _getBodyName(ibody, nameList, body):

    name = body.attributeRet("_name")
    if name is None:
        name = "Body_" + str(ibody+1)

    if name in nameList:
        name = "Body_" + str(ibody+1)
        print("Body", ibody, ": Name '", name, "' found more than once. Changing name to: '", name, "'")

    nameList.append(name)

    return name

#==============================================================================
def _ocsmPath(filename):

    relname = os.path.relpath(filename)
    if len(relname) < len(filename): filename = relname
    return filename

#==============================================================================
# Perform unit conversion
@deprecated("'pyCAPS.Unit'")
def capsConvert(value, fromUnits, toUnits):
    toUnits = caps.Unit(toUnits)
    return caps.Quantity(value,fromUnits).convert(toUnits)/toUnits

#==============================================================================
# Create the tree for an object
def createTree(dataItem, showAnalysisGeom, showInternalGeomAttr, reverseMap):

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
    if isinstance(dataItem, Analysis):
        tree = {"name" : dataItem.name, # + " (" + dataItem.officialName + ")",
                "children" : [{"name" : "Directory", "children" : checkValue(dataItem.analysisDir)},
                              {"name" : "Intent"   , "children" : checkValue(None)}, #checkValue(dataItem.capsIntent)},
                              {"name" : "Inputs"   , "children" : checkValue(None)},
                              {"name" : "Outputs"  , "children" : checkValue(None)},
                             ]
                }

        if dataItem.input:

            tree["children"][2]["children"] = [] # Clear Inputs

            # Loop through analysis inputs and add current value
            for key, value in dataItem.input.items():
                tree["children"][2]["children"].append( {"name" : key, "children" : checkValue(value)} )

        if dataItem.output:

            tree["children"][3]["children"] = [] # Clear Outputs

            # Loop through analysis outputs and add current value
            for key, value in dataItem.output.items():
                tree["children"][3]["children"].append( {"name" : key, "children" : checkValue(value)} )

        # Are we supposed to show geometry with the analysis obj?
        if showAnalysisGeom == True:

            attributeMap = dataItem.geometry.attrMap(getInternal = showInternalGeomAttr, reverseMap = reverseMap)

            if attributeMap: # If attribute is populated

                tree["children"].append({"name" : "Geometry", "children" : []})

                for key, value in attributeMap.items():
                    tree["children"][4]["children"].append( {"name" : key, "children" : checkValue(value)} )

    # Data set object
    if isinstance(dataItem, DataSet):
        tree = {"name" : dataItem.name,
                "children" : [{"name" : "Vertex Set", "children" : checkValue(dataItem.capsVertexSet.vertexSetName)},
                              {"name" : "Data Rank", "children"  : checkValue(dataItem.dataRank)},
                              {"name" : "Data Method", "children" : checkValue(dataItem.dataSetMethod)}
                              ]
                }

    # Bound object
    if isinstance(dataItem, Bound):
        tree = {"name" : dataItem.boundName,
                "children" : [{"name" : "Source Data Set"      , "children" : checkValue(None)},
                              {"name" : "Destination Data Set" , "children" : checkValue(None)},
                             ]
                }

        # Source data sets
        if dataItem.dataSetSrc:

            tree["children"][0]["children"] = [] # Clear Source data set

            for key, value in dataItem.dataSetSrc.items():
                tree["children"][0]["children"].append(jsonLoads(createTree(value,
                                                                           showAnalysisGeom,
                                                                           showInternalGeomAttr,
                                                                           reverseMap
                                                                           )
                                                                )
                                                      )

        # Destination data sets - if there are any
        if dataItem.dataSetDest:

            tree["children"][1]["children"] = [] # Clear Destination data set

            for key, value in dataItem.dataSetDest.items():
                tree["children"][1]["children"].append(jsonLoads(createTree(value,
                                                                            showAnalysisGeom,
                                                                            showInternalGeomAttr,
                                                                            reverseMap
                                                                            )
                                                                 )
                                                       )

    # Geometry object
    if isinstance(dataItem, ProblemGeometry):
        tree = { #"name" : dataItem.name,
                "children" : [{"name" : "Design Parameters", "children" : checkValue(None)},
                              {"name" : "Local Variables", "children" : checkValue(None)},
                              {"name" : "Attributes", "children" : checkValue(None)},
                             ]
                }

        # Design variables
        if dataItem.despmtr:

            tree["children"][0]["children"] = [] # Clear design parameters

            # Loop through design variables
            for key, value in dataItem.despmtr.items():
                tree["children"][0]["children"].append( {"name" : key, "children" : checkValue(value)} )

        # Local variables
        if dataItem.outpmtr:

            tree["children"][1]["children"] = [] # Clear local variables

            # Loop through local variables
            for key, value in dataItem.outpmtr.items():
                tree["children"][1]["children"].append( {"name" : key, "children" : checkValue(value)} )

        # Attributes
        attributeMap = dataItem.attrMap(getInternal = showInternalGeomAttr, reverseMap = reverseMap)
        if attributeMap:
            tree["children"][2]["children"] = [] # Clear attribute map
            for key, value in attributeMap.items():
                tree["children"][2]["children"].append( {"name" : key, "children" : checkValue(value)} )

    # Problem object
    if isinstance(dataItem, Problem):
#         latest = {"CAPS Intent" : dataItem.capsIntent,
#                   "AIM Count"   : dataItem.aimGlobalCount,
#                   "Analysis Directory" : dataItem.analysisDir
#                   }

        tree = {"name" : "myProblem",
                "children" : [{"name" : "Analysis", "children" : checkValue(None)},
                              {"name" : "Geometry", "children" : checkValue(None)},
                              {"name" : "Bound",    "children" : checkValue(None)},
#                               {"name" : "Latest Defaults" , "children" : checkValue(latest)},
                             ]
                }

        # Do analysis objects - if there are any
        if dataItem.analysis:

            tree["children"][0]["children"] = [] # Clear Analysis

            for key, value in dataItem.analysis.items():
                tree["children"][0]["children"].append(jsonLoads(createTree(value,
                                                                            showAnalysisGeom,
                                                                            showInternalGeomAttr,
                                                                            reverseMap)
                                                                 )
                                                       )

        # Do geometry
        tree["children"][1]["children"] = [jsonLoads(createTree(dataItem.geometry,
                                                                showAnalysisGeom,
                                                                showInternalGeomAttr,
                                                                reverseMap
                                                                )
                                                     )
                                           ]

        # Do data bounds - if there are any
        if dataItem.bound: #

            tree["children"][2]["children"] = [] # Clear Data Bounds

            for key, value in dataItem.bound.items():
                tree["children"][2]["children"].append(jsonLoads(createTree(value,
                                                                            showAnalysisGeom,
                                                                            showInternalGeomAttr,
                                                                            reverseMap
                                                                            )
                                                                 )
                                                       )


    return jsonDumps(tree, indent= 2)


#==============================================================================
## Defines a CAPS Problem Object.
#
# The Problem Object is the top-level object for a single mission/problem. It maintains a single set
# of interrelated geometric models (see \ref ProblemGeometry),
# analyses to be executed (see \ref Analysis),
# connectivity and data (see \ref Bound)
# associated with the run(s), which can be both multi-fidelity and multi-disciplinary.
#
# \param Problem.geometry  \ref ProblemGeometry instances representing the CSM geometry
# \param Problem.analysis  \ref AnalysisSequence of \ref Analysis instances
# \param Problem.parameter \ref ParamSequence of \ref ValueIn parameters
# \param Problem.bound     \ref BoundSequence of \ref Bound instances
# \param Problem.attr      \ref AttrSequence of \ref ValueIn attributes
#
class Problem(object):

    __slots__ = ['_problemObj', 'geometry', 'parameter', 'analysis', 'bound', 'attr', '_name']

    ## Initialize the problem.
    # \param problemName CAPS problem name that serves as the root directory for all file I/O.
    # \param phaseName the current phase name (None is equivalent to 'Scratch')
    # \param prevPhaseName the previous phase name
    # \param capsFile CAPS file to load. Options: *.csm or *.egads.
    # \param outLevel Level of output verbosity. See \ref setOutLevel .
    # \param useJournal Use Journaling to continue execution of an interrupted script.
    #
    def __init__(self, problemName, phaseName=None, prevPhaseName=None, capsFile=None, outLevel=1, useJournal=False):

        verbosity = {"minimal"  : 0,
                     "standard" : 1,
                     "debug"    : 2}

        if not isinstance(outLevel, (int,float)):
            if outLevel.lower() in verbosity.keys():
                outLevel = verbosity[outLevel.lower()]

        if int(outLevel) not in verbosity.values():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "invalid verbosity level! outLevel={!r}".format(outLevel))

        if isinstance(problemName,caps.capsObj):
            problemObj = problemName
            problemName, otype, stype, link, parent, last = problemObj.info()
        else:
            if capsFile is not None:
                ptr = capsFile
                flag = caps.oFlag.oFileName
                if useJournal:
                    flag = caps.oFlag.oContinue
            elif prevPhaseName is not None:
                flag = caps.oFlag.oPhaseName
                ptr = prevPhaseName


            problemObj = caps.open(problemName, phaseName, flag, ptr, outLevel)

        super(Problem, self).__setattr__("_problemObj", problemObj)

        super(Problem, self).__setattr__("geometry",  ProblemGeometry(problemObj))
        super(Problem, self).__setattr__("parameter", ParamSequence(problemObj))
        super(Problem, self).__setattr__("analysis",  AnalysisSequence(self))
        super(Problem, self).__setattr__("bound",     BoundSequence(problemObj))
        super(Problem, self).__setattr__("attr",      AttrSequence(problemObj))
        super(Problem, self).__setattr__("_name",     problemName)

    def __setattr__(self, name, data):
        raise AttributeError("Cannot set attribute: {!r}".format(name))

    def __delattr__(self, name):
        raise AttributeError("Cannot del attribute: {!r}".format(name))


    # delete all references to _problemObj after closing the problem
    def _clean(self):
        super(Problem, self).__delattr__("geometry")
        super(Problem, self).__delattr__("parameter")
        super(Problem, self).__delattr__("analysis")
        super(Problem, self).__delattr__("bound")
        super(Problem, self).__delattr__("attr")

    ## Exlicitly closes CAPS Problem Object
    #
    # This method is mainly useful for testing purposes
    #
    def close(self):
        self._problemObj.close()
        self._clean()

    ## Completes the Phase and closes the CAPS Problem Object
    # \param phaseName Phase Name of the Scratch phase is closed as complete
    #
    def closePhase(self, phaseName=None):
        self._problemObj.close(1, phaseName)
        self._clean()

    ## Property returns the name of the CAPS Problem Object
    @property
    def name(self):
        return self._name

    ## Boolean indicator if the CAPS Problem Object is currently journaling
    def journaling(self):
        return self._problemObj.journalState()

    def outputObjects(self):
        self._problemObj.outputObjects()

    ## Set the verbosity level of the CAPS output.
    # See \ref problem5.py for a representative use case.
    #
    # \param outLevel Level of output verbosity. Options: 0 (or "minimal"),
    # 1 (or "standard") [default], and 2 (or "debug").
    def setOutLevel(self, outLevel):
        ## \example problem5.py
        # Basic example for setting the verbosity of a problem using pyCAPS.Problem.setOutLevel() function.

        verbosity = {"minimal"  : 0,
                     "standard" : 1,
                     "debug"    : 2}

        if not isinstance(outLevel, (int,float)):
            if outLevel.lower() in verbosity.keys():
                outLevel = verbosity[outLevel.lower()]

        if int(outLevel) not in verbosity.values():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "invalid verbosity level! outLevel={!r}".format(outLevel))

        self._problemObj.setOutLevel(outLevel)


    ## Create a link between a created CAPS parameter and analyis inputs of <b>all</b> loaded AIMs, automatically.
    # Valid CAPS value, parameter objects must be created with Problam.parameter.create(). Note, only links to
    # ANALYSISIN inputs are currently made at this time.
    #
    # \param param Parameter to use when creating the link (default - None). A combination (i.e. a single or list)
    # of \ref ValueIn dictionary entries and/or value object instances (returned from a
    # call to Problam.parameter.create()) can be used. If no value is provided, all entries in the
    # ValueIn dictionary (\ref ValueIn) will be used.
    def autoLinkParameter(self, param=None):

        # If none we want to try to link all values
        if param is None:
            param = list(self.parameter.keys())

        # If not a list convert it to one
        if not isinstance(param, list):
            param = [param]

        # Loop through the values
        for i in param:

            # If it is a ValueIn instance grab the name
            if isinstance(i, (ValueIn, _capsValue)):
               i = i.name

            # Check to make sure name is in parameter Sequence
            if i not in self.parameter.keys():
                raise caps.CAPSError(caps.CAPS_NOTFOUND, msg = "while autolinking values in problem! "
                                                             + "Unable to find parameter " + str(i))

            varname = self.parameter[i].name

            found = False
            # Loop through analysis dictionary
            for analysis in self.analysis.values():

                # Look to see if the value name is present as an analysis input variable
                if varname in analysis.input:

                    analysis.input[varname].link(self.parameter[i])

                    print("Linked Parameter", str(self.parameter[i].name), "to analysis", str(analysis.name), "input", varname)
                    found = True

            if not found:
                print("No linkable data found for", str(self.parameter[i].name))


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

#==============================================================================
# deprecated
class capsProblem(Problem):

    __slots__ = ['value','dataBound','capsIntent']

    @deprecated("'Problem.__init__'")
    def __init__(self, raiseException=True):

        super(Problem, self).__setattr__("_problemObj", None)
        super(Problem, self).__setattr__("geometry",  None)
        super(Problem, self).__setattr__("parameter", None)
        super(Problem, self).__setattr__("analysis",  None)
        super(Problem, self).__setattr__("bound",     None)
        super(Problem, self).__setattr__("attr",      None)

        super(Problem, self).__setattr__("value", dict())
        super(Problem, self).__setattr__("dataBound", dict())
        super(Problem, self).__setattr__("capsIntent", None)
        pass

    def __setattr__(self, name, data):
        super(Problem, self).__setattr__(name, data)
        if name == "capsIntent":
            deprecate("capsProblem.capsIntent", "'None'")

    @deprecated("'Problem.__init__'")
    def loadCAPS(self, capsFile, projectName=None, verbosity=1):

        capsGeometry(self, capsFile, projectName, verbosity)

        for name, val in self.parameter.items():
            self.value[name] = self.parameter[name]

        return self.geometry

    @deprecated("'Problem.analysis.create'")
    def loadAIM(self, aim=None, altName=None, unitSystem=None, capsIntent=None, analysisDir=None, parents=None, copyAIM=None):
        if self.capsIntent and capsIntent is None:
            capsIntent = self.capsIntent
        return capsAnalysis(self, aim, altName, unitSystem, capsIntent, analysisDir, parents, copyAIM)

    @deprecated("'Problem.setOutLevel'")
    def setVerbosity(self, outlevel):
        return self.setOutLevel(outlevel)

    @deprecated("'Problem.parameter.create'")
    def createValue(self, name, data, units=None, limits=None, fixedLength=True, fixedShape=True):
        if name in self.value.keys():
            raise caps.CAPSError(caps.CAPS_BADNAME, msg="A value with the name," + name + "already exists!")

        if isinstance(units, str):
            unit = caps.Unit(units)
            data = data*unit
            if limits: limits = limits*unit

        self.parameter.create(name, data, limits=limits, fixedLength=fixedLength, fixedShape=fixedShape)
        self.value[name] = _capsValue(self.parameter[name]._valObj, units)
        return self.value[name]

    @deprecated("'Problem.analysis.dirty'")
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

    @deprecated("'close'")
    def closeCAPS(self):

        # delete all references to _problemObj
        super(Problem, self).__delattr__("geometry")
        super(Problem, self).__delattr__("parameter")
        super(Problem, self).__delattr__("analysis")
        super(Problem, self).__delattr__("bound")
        super(Problem, self).__delattr__("attr")

        super(Problem, self).__delattr__("value")
        super(Problem, self).__delattr__("dataBound")
        super(Problem, self).__delattr__("capsIntent")

        super(Problem, self).__delattr__('_problemObj')

        # re-initialize everyhing to an empty state
        self.__init__()

    @deprecated("'None'")
    def saveCAPS(self, filename="saveCAPS.caps"):
        pass

    @deprecated("'Problem.attr.create'")
    def addAttribute(self, name, data):
        self.attr.create(name,data,True)

    @deprecated("'Problem.attr[\"name\"].value'")
    def getAttribute(self, name):
        value = self.attr[name].value

        if isinstance(value,caps.Quantity):
            value = value._value

        return value

    @deprecated("'Problem.autoLinkParameter'")
    def autoLinkValue(self, value=None):
        self.autoLinkParameter(value)

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

## Pythonic alias to capsProblem (deprecated)
CapsProblem = capsProblem

#=============================================================================
## Defines Problem Geometry Object
#
# \param ProblemGeometry.despmtr \ref ValueInSequence of \ref ValueIn CSM design parameters
# \param ProblemGeometry.cfgpmtr \ref ValueInSequence of \ref ValueIn CSM configuration parameters
# \param ProblemGeometry.conpmtr \ref ValueInSequence of \ref ValueIn CSM constant parameters
# \param ProblemGeometry.outpmtr \ref ValueOutSequence of \ref ValueOut CSM outputs

class ProblemGeometry(object):

    __slots__ = ['_problemObj', 'despmtr', 'cfgpmtr', 'conpmtr', 'outpmtr']

    def __init__(self, problemObj):
        super(ProblemGeometry, self).__setattr__('_problemObj', problemObj)

        super(ProblemGeometry, self).__setattr__("despmtr", ValueInSequence (_values(self._problemObj, caps.sType.GEOMETRYIN)))
        super(ProblemGeometry, self).__setattr__("cfgpmtr", ValueInSequence (_values(self._problemObj, caps.sType.GEOMETRYIN, 1)))
        super(ProblemGeometry, self).__setattr__("conpmtr", ValueInSequence (_values(self._problemObj, caps.sType.GEOMETRYIN, 2)))
        super(ProblemGeometry, self).__setattr__('outpmtr', ValueOutSequence(_values(self._problemObj, caps.sType.GEOMETRYOUT)))

    def __setattr__(self, name, data):
        raise AttributeError("Cannot add attribute: {!r}".format(name))

    def __delattr__(self, name):
        raise AttributeError("Cannot del attribute: {!r}".format(name)) #from None

    ## Exlicitly build geometry.
    def build(self):
        self._problemObj.execute()

    ## Save the current geometry to a file.
    # \param filename File name to use when saving geometry file.
    # \param directory Directory where to save file. Default current working directory.
    # \param extension Extension type for file if filename does not contain an extension.
    def save(self, filename="myGeometry", directory=os.getcwd(), extension=".egads"):
        filename = _createFilename(filename, directory, extension)
        flag = 1

        self._problemObj.writeGeometry(filename, flag)

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
    def view(self, **kwargs):

        viewDefault = "capsviewer"

        viewer = kwargs.pop("viewerType", viewDefault).lower()
        filename = kwargs.get("filename", None)

        # Get the bodies
        bodies, units = self.bodies()

        # Are we going to use the viewer?
        if filename is None and viewer.lower() == "capsviewer":
            self._viewGeometryCAPSViewer(**kwargs)
        else:
            _viewGeometryMatplotlib(bodies, **kwargs)

    def _viewGeometryCAPSViewer(self, **kwargs):

        egadsFile = os.path.join(self._problemObj.getRootPath(), "viewProblemGeometry.egads")
        relFile = os.path.relpath(egadsFile, os.getcwd())

        # use the shortest filename
        if len(relFile) < len(egadsFile):
            egadsFile = relFile

        self.save( egadsFile )

        # dissable the viewer (used for testing)
        if os.getenv("CAPS_BATCH"):
            os.remove(egadsFile)
            return

        egadsFile = _ocsmPath(egadsFile)
        os.system("serveCSM -port " + str(kwargs.pop("portNumber", 7681)) + " -outLevel 0 " + egadsFile)
        os.remove(egadsFile)
        os.remove("autoEgads.csm")

        return

        # Initialize the viewer class
        viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))

        # Get the bodies
        bodies = self._problemObj.getBodies()

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
    def attrList(self, attributeName, **kwargs):

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
            raise caps.CAPSError(CAPS_BADVALUE, "Error: Invalid attribute level!")

        # How many bodies do we have
        numBody = self._problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while getting attribute " + str(attributeName) + " values on geometry"
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i

            if bodyIndex is not None:
                if int(bodyIndex) != i+1:
                    continue

            body, units = self._problemObj.bodyByIndex(i+1)

            attributeList += _createAttributeValList(body,
                                                     attrLevel,
                                                     attributeName)

        try:
            return list(set(attributeList))
        except TypeError:
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
    def attrMap(self, getInternal = False, **kwargs):

        attributeMap = {}
        nameList = []

        # Get the bodies
        bodies, units = self.bodies()

        if len(bodies) < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating attribute map"
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        for name, body in bodies.items():

            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            attributeMap[name] = _createAttributeMap(body, getInternal)

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
    # \param reverseMap Reverse the attribute map (default - False). See \ref attrMap for details.
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

    ## Get dict of geometric bodies.
    #
    # \return Returns a dictionary of the bodies and the capsLength unit.
    # Keys use the body "_name" attribute or "Body_#".
    def bodies(self):

        bodies = {}
        nameList = []

        numBody = self._problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "The number of bodies in the problem is less than 1!!!" )

        for i in range(numBody):
            body, units = self._problemObj.bodyByIndex(i+1)

            name = _getBodyName(i, nameList, body)

            bodies[name] = body

        return bodies, units

    ## Get the lenght Unit of geometric bodies.
    #
    # \return Returns the length unit defined by capsLength attribute.
    def lengthUnit(self):

        numBody = self._problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "The number of bodies in the problem is less than 1!!!" )

        body, units = self._problemObj.bodyByIndex(1)

        return units

    ## Write an OpenCSM Design Parameter file to disk
    #
    # \param filename Filename of the OpenCSM Design Parameter file
    def writeParameters(self, filename):
        self._problemObj.writeParameters(filename)

    ## Read an OpenCSM Design Parameter file from disk and and overwrites (makes dirty)
    #   the current state of the geometry
    #
    # \param filename Filename of the OpenCSM Design Parameter file
    def readParameters(self, filename):
        self._problemObj.readParameters(filename)

#==============================================================================
# deprecated
class capsGeometry(ProblemGeometry):

    __slots__ = []

    @deprecated("'Problem.__init__'")
    def __init__(self, problem, capsFile, projectName=None, verbosity=1):

        # Check problem object
        if not isinstance(problem, capsProblem) and not isinstance(problem, CapsProblem):
            raise TypeError

        if problem.geometry is not None:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while loading file - " + str(capsFile) + ". Can NOT load multiple files into a problem."
                                                         + "\nWarning: Can NOT load multiple files into a problem. A CAPS file has already been loaded!\n")

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

        super(capsProblem, problem).__init__(projectName, None, None, capsFile, verbosity)
        super(capsGeometry, self).__init__(problem._problemObj)

        super(Problem, problem).__setattr__("geometry", self)

    @deprecated("'Problem.geometry.build'")
    def buildGeometry(self):
        self.build()

    @deprecated("'Problem.geometry.despmtr[\"varname\"].value'")
    def setGeometryVal(self, varname, value=None):
        if varname in self.despmtr:
            self.despmtr[varname].value = value
        elif varname in self.cfgpmtr:
            self.cfgpmtr[varname].value = value
        else:
            raise caps.CAPSError(caps.CAPS_BADNAME, "No variable name {!r}".format(varname))

    @deprecated("'Problem.geometry.despmtr[\"varname\"].value'")
    def getGeometryVal(self, varname=None, **kwargs):

        geometryIn = {}

        namesOnly = kwargs.pop("namesOnly", False)

        if varname is not None: # Get just the 1 variable

            valueObj = self._problemObj.childByName(caps.oType.VALUE, caps.sType.GEOMETRYIN, varname)

            value = valueObj.getValue()
            if isinstance(value,caps.Quantity):
                value = value._value
            return value

        else: # Build a dictionary

            numGeometryInput = self._problemObj.size(caps.oType.VALUE, caps.sType.GEOMETRYIN)

            for i in range(1, numGeometryInput+1):

                tempObj = self._problemObj.childByIndex(caps.oType.VALUE, caps.sType.GEOMETRYIN, i)

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

    @deprecated("'Problem.geometry.outpmtr[\"varname\"].value'")
    def getGeometryOutVal(self, varname=None, **kwargs):

        geometryOut = {}

        namesOnly = kwargs.pop("namesOnly", False)
        ignoreAt  = kwargs.pop("ignoreAt", True)

        if varname is not None: # Get just the 1 variable

            valueObj = self._problemObj.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, varname)

            value = valueObj.getValue()
            if isinstance(value,caps.Quantity):
                value = value._value
            return value

        else: # Build a dictionary

            numGeometryOutput = self._problemObj.size(caps.oType.VALUE, caps.sType.GEOMETRYOUT)

            for i in range(1, numGeometryOutput+1):

                tempObj = self._problemObj.childByIndex(caps.oType.VALUE, caps.sType.GEOMETRYOUT, i)

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

    @deprecated("'Problem.geometry.attrList'")
    def getAttributeVal(self, attributeName, **kwargs):

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
            kwargs["attrLevel"] = attrLevel

        return self.attrList(attributeName, **kwargs)

    @deprecated("'Problem.geometry.attrMap'")
    def getAttributeMap(self, getInternal = False, **kwargs):

        attributeMap = {}

        # How many bodies do we have
        numBody = self._problemObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating attribute map"
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        for bodyIndex in range(numBody, 0, -1):

            # No plus +1 on bodyIndex since we are going backwards
            body, units = self._problemObj.bodyByIndex(bodyIndex)

            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            attributeMap["Body " + str(bodyIndex)] = _createAttributeMap(body, getInternal)

        if kwargs.pop("reverseMap", False):
            return _reverseAttributeMap(attributeMap)
        else:
            return attributeMap

    @deprecated("'Problem.geometry.save'")
    def saveGeometry(self, filename="myGeometry", directory=os.getcwd(), extension=".egads"):
        self.save(filename, directory, extension)

    @deprecated("'Problem.geometry.view'")
    def viewGeometry(self, **kwargs):
        return self.view(**kwargs)

    @deprecated("'Problem.geometry.bodies'")
    def getBoundingBox(self):

        boundingBox = {}

        bodies, units = self.bodies()

        if bodies is None:
            raise caps.CAPSError(caps.CAPS_SOURCEERR, msg = "No bodies available to retrieve the bounding box!")

        for name, body in bodies.items():
            box  = body.getBoundingBox()

            boundingBox[name] = [box[0][0], box[0][1], box[0][2],
                                 box[1][0], box[1][1], box[1][2]]

        return boundingBox

## Pythonic alias to capsProblem (deprecated)
CapsGeometry = capsGeometry

#=============================================================================
# deprecated
class _capsValue(object):

    @deprecated("'Problem.parameter.create'")
    def __init__(self, valueObj, units):
        self.valueObj = valueObj
        self.units = units

        self.name, objType, self.subType, link, parent, owner = self.valueObj.info()

    @property
    def value(self):

        self._value = self.valueObj.getValue()
        if isinstance(self._value, caps.Quantity):
            self._value = self._value._value

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

    def setVal(self, data):
        self.value = data

    def getVal(self):
        return self.value

    def setLimits(self, newLimits):
        self.limits = newLimits

    def getLimits(self):
        return self.limits

    def convertUnits(self, toUnits):
        if self.units == None:
            raise caps.CAPSError(caps.CAPS_UNITERR, msg = "while converting units for value object - " + str(self.name)
                                                        + "\nValue doesn't have defined units")

        return caps.Quantity(self.value, self.units)/caps.Unit(toUnits)

#==============================================================================
## Base class for all CAPS Sequence classes
#
# A CAPS Sequence only contains instances of a single type
# Items are added to the Seqence via the 'create' method in derived classes
# Items cannot be removed from the sequence (except for CAPS Attributes)
class Sequence(object):

    __slots__ = ["_capsItems", "_typeName"]

    def __init__(self, typeName):
        super(Sequence, self).__setattr__('_capsItems', dict())
        super(Sequence, self).__setattr__('_typeName', typeName)

    def __getitem__(self, key):
        if key not in self._capsItems:
            raise KeyError("No {}: {!r}".format(self._typeName, key))
        return self._capsItems[key]

    def __setitem__(self, key, value):
        if key in self._capsItems:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot set {}: {!r}".format(self._typeName, key))
        else:
            raise KeyError("No {}: {!r}".format(self._typeName, key))

    def __delitem__(self, key):
        if key in self._capsItems:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot del {}: {!r}".format(self._typeName,key))
        else:
            raise KeyError("No {}: {!r}".format(self._typeName,key))

    def __contains__(self, item):
        return item in self._capsItems

    def __len__(self):
        return len(self._capsItems)

    def __iter__(self):
        return self._capsItems.__iter__()

    def __next__(self):
        return self._capsItems.__next__()

    def __reversed__(self):
        return self._capsItems.__reversed__()

    def __bool__(self):
        return bool(self._capsItems)

    ## Returns the keys of the Sequence
    def keys(self):
        return self._capsItems.keys()

    ## Returns the values of the Sequence
    def values(self):
        return self._capsItems.values()

    ## Returns the items of the Sequence
    def items(self):
        return self._capsItems.items()


#==============================================================================
## Defines a CAPS input Value Object
class ValueIn(object):

    __slots__ = ["_valObj", "_name"]

    def __init__(self, valObj):
        self._valObj = valObj
        self._name, otype, stype, link, parent, last = valObj.info()

    ## Property getter returns a copy the values stored in the CAPS Value Object
    @property
    def value(self):
        return self._valObj.getValue()

    ## Property setter sets the value in the CAPS Value Object
    @value.setter
    def value(self, val):
        self._valObj.setValue(val)

    ## Property getter returns a copy the limits of the CAPS Value Object
    @property
    def limits(self):
        return self._valObj.getLimits()

    ## Property setter sets the limits in the CAPS Value Object (if changable)
    @limits.setter
    def limits(self, limit):
        self._valObj.setLimits(limit)

    ## Property returns the name of the CAPS Value Object
    @property
    def name(self):
        return self._name

    def props(self):
        dim, pmtr, lfix, sfix, ntype = self._valObj.getValueProps()
        return {"dim":dim, "pmtr":pmtr, "lfix":lfix, "sfix":sfix, "ntype":ntype}

    ## Link this input value to an other CAPS Value Object
    # \param source The source Value Object
    # \param tmethod Transfter method: tMethod.Copy or "Copy", tMethod.Integrate or "Integrate", tMethod.Average or "Average"
    def link(self, source, tmethod = caps.tMethod.Copy):

        tMethods = {"copy"      : caps.tMethod.Copy,
                    "integrate" : caps.tMethod.Integrate,
                    "average"   : caps.tMethod.Average}

        if not isinstance(tmethod, (int,float)):
            if tmethod.lower() in tMethods.keys():
                tmethod = tMethods[tmethod.lower()]

        if not (isinstance(source, ValueIn) or isinstance(source, ValueOut)):
            raise caps.CAPSError(caps.CAPS_BADVALUE, "source must be a Value Object")

        self._valObj.linkValue(source._valObj, tmethod)

    ## Remove an existing link
    def unlink(self):
        self._valObj.removeLink()

    ## Transfer values from src to self
    # \param tmethod 0 - copy, 1 - integrate, 2 - weighted average  -- (1 & 2 only for DataSet src)
    # \param source the source value object
    def transferValue(self, tmethod, source):
        if not (isinstance(source, ValueIn) or isinstance(source, ValueOut)):
            raise caps.CAPSError(caps.CAPS_BADVALUE, "source must be a Value Object")

        self._valObj.transferValues(tmethod, source._valObj)

#==============================================================================

def _values(obj, stype, usepmtr=0):
    nVal = obj.size(caps.oType.VALUE, stype)
    valObjs = []

    for i in range(nVal):
        valObj = obj.childByIndex(caps.oType.VALUE, stype, i+1)
        name, otype, stype, link, parent, last = valObj.info()
        dim, pmtr, lfix, sfix, ntype = valObj.getValueProps()
        if usepmtr != pmtr: continue
        valObjs.append(valObj)

    return valObjs

#==============================================================================
## Defines a Sequence of CAPS input Value Objects
class ValueInSequence(Sequence):

    __slots__ = ["_subVals"]

    def __init__(self, valObjs, level=0):
        super(self.__class__, self).__init__('In Value')

        subVals = {}

        for valObj in valObjs:
            name, otype, stype, link, parent, last = valObj.info()
            if ":" in name:
                names = name.split(":")
                self._capsItems[":".join(names[level:])] = ValueIn(valObj)
                if level+1 == len(names):
                    continue

                if names[level] not in subVals:
                    subVals[names[level]] = []
                subVals[names[level]].append(valObj)
            else:
                self._capsItems[name] = ValueIn(valObj)

        super(self.__class__, self).__setattr__('_subVals', dict())
        for name, vals in subVals.items():
            self._subVals[name] = ValueInSequence(vals, level+1)

    def __getattr__(self, name):
        if name in self._capsItems:
            return self._capsItems[name].value
        elif name in self._subVals:
            return self._subVals[name]

        raise AttributeError("No In Value: {!r}".format(name))

    def __setattr__(self, name, value):
        if name in self._capsItems:
            self._capsItems[name].value = value
            return
        elif name in self._subVals:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot set: {!r}".format(name))

        raise AttributeError("No In Value: {!r}".format(name))

    def __delattr__(self, name):
        if name in self._capsItems or name in self._subVals:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot del In Value: {!r}".format(name))
        else:
           raise AttributeError("No In Value: {!r}".format(name))

    def __getitem__(self, name):
        if name in self._capsItems:
            return self._capsItems[name]
        elif name in self._subVals:
            return self._subVals[name]

        raise KeyError("No In Value: {!r}".format(name))

#==============================================================================
## Defines a Sequence of CAPS Attribute Value Objects
class AttrSequence(Sequence):

    __slots__ = ["_capsObj"]

    def __init__(self, obj):
        super(self.__class__, self).__init__('CAPS Attribute')

        super(self.__class__, self).__setattr__('_capsObj', obj)

        nAttr = obj.size(caps.oType.ATTRIBUTES, caps.sType.NONE)

        for i in range(nAttr):
            valObj = obj.attrByIndex(i+1)
            name, otype, stype, link, parent, last = valObj.info()
            self._capsItems[name] = ValueIn(valObj)

    ## Create an attribute (that is meta-data) to the CAPS Object. See example
    #
    # \param name Name used to define the attribute.
    #
    # \param data Initial data value(s) for the attribute. Note that type casting in done
    # automatically based on the determined type of the Python object.
    #
    # \param overwrite Flag to overwrite any existing attribute with the same 'name'
    #
    # \return The new Value Object is added to the sequence and returned
    def create(self, name, data, overwrite=False):

        if not overwrite and name in self._capsItems:
            raise KeyError("Parameter {!r} already exists! Use overwrite=True to overwrite it.".format(name))

        valObj = self._capsObj.problemObj().makeValue(name, caps.sType.USER, data)

        self._capsObj.setAttr(valObj)
        self._capsItems[name] = ValueIn(valObj)

        return self._capsItems[name]

    def __getattr__(self, name):
        if name not in self._capsItems:
            raise AttributeError("No CAPS Attribute: {!r}".format(name))
        return self._capsItems[name].value

    def __setattr__(self, name, value):
        if name not in self._capsItems:
            raise AttributeError("No CAPS Attribute: {!r}".format(name))
        self._capsItems[name].value = value

    def __delattr__(self, name):
        if name not in self._capsItems:
            raise AttributeError("No CAPS Attribute: {!r}".format(name))
        del self._capsItems[name]
        self._capsObj.deleteAttr(name)

    def __delitem__(self, key):
        if key in self._capsItems:
            del self._capsItems[key]
            self._capsObj.deleteAttr(key)
        else:
            raise KeyError("No CAPS Attribute: {!r}".format(key))

#==============================================================================
## Defines a Sequence of CAPS Parameter Value Objects
class ParamSequence(Sequence):

    __slots__ = ["_problemObj"]

    def __init__(self, problemObj):
        super(self.__class__, self).__init__('Parameter')

        super(self.__class__, self).__setattr__('_problemObj', problemObj)

        nParam = problemObj.size(caps.oType.VALUE, caps.sType.PARAMETER)

        for i in range(nParam):
            valObj = problemObj.childByIndex(caps.oType.VALUE, caps.sType.PARAMETER, i+1)
            name, otype, stype, link, parent, last = valObj.info()
            self._capsItems[name] = ValueIn(valObj)

    ## Create an parameter CAPS Value Object.
    #
    # \param name Name used to define the parameter.
    #
    # \param data Initial data value(s) for the parameter. Note that type casting in done
    # automatically based on the determined type of the Python object.
    #
    # \param limits Limits on the parameter values
    #
    # \param fixedLength Boolean if the value is fixed length
    #
    # \param fixedShape Boolean if the value is fixed shape
    #
    # \return The new Value Object is added to the sequence and returned
    def create(self, name, data, limits=None, fixedLength=True, fixedShape=True):

        # Create the new Value Object
        valObj = self._problemObj.makeValue(name, caps.sType.PARAMETER, data)

        # Get shape parameters
        dim, pmtr, lfixed, sfixed, ntype = valObj.getValueProps()

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
        valObj.setValueProps(dim, lfixed, sfixed, ntype)

        # Set limits
        if limits:
            valObj.setLimits(limits)

        self._capsItems[name] = ValueIn(valObj)
        return self._capsItems[name]

    def __getattr__(self, name):
        if name not in self._capsItems:
            raise AttributeError("No Parameter: {!r}".format(name))
        return self._capsItems[name].value

    def __setattr__(self, name, value):
        if name not in self._capsItems:
            raise AttributeError("No Parameter: {!r}".format(name))
        self._capsItems[name].value = value

    def __delattr__(self, name):
        if name in self._capsItems:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot del Parameter: {!r}".format(name))
        else:
            raise AttributeError("No Parameter: {!r}".format(name))

#==============================================================================
## Defines a CAPS output Value Object
# Not a standalone class.
class ValueOut(object):

    __slots__ = ["_valObj", "_name"]

    def __init__(self, valObj):
        self._valObj = valObj
        self._name, otype, stype, link, parent, last = valObj.info()

    ## Property getter returns a copy the values stored in the CAPS Value Object
    @property
    def value(self):
        return self._valObj.getValue()

    @value.setter
    def value(self, val):
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot assign Out Value")

    ## Property returns the name of the CAPS Value Object
    @property
    def name(self):
        return self._name

    ## Property getter returns a copy the values stored in the CAPS Value Object
    def props(self):
        dim, pmtr, lfix, sfix, ntype = self._valObj.getValueProps()
        return {"dim":dim, "pmtr":pmtr, "lfix":lfix, "sfix":sfix, "ntype":ntype}

    ## Returns a string list of of the input Value Object names that can be used in \ref deriv
    def hasDeriv(self):
        return self._valObj.hasDeriv()

    ## Returns derivatives of the output Value Object
    #
    # \param name Name of the input Value Object to take derivative w.r.t.
    #        if name is None then a dictionary with all dervatives from \ref hasDeriv are returned
    def deriv(self, name=None):
        if name is not None:
            return self._valObj.getDeriv(name)
        else:
            names = self._valObj.hasDeriv()
            derivs = {}
            for name in names:
                derivs[name] = self._valObj.getDeriv(name)

            return derivs

#==============================================================================
## Defines a Sequence of CAPS output Value Objects
class ValueOutSequence(Sequence):

    __slots__ = ["_subVals"]

    def __init__(self, valObjs, level=0):
        super(self.__class__, self).__init__('Out Value')

        subVals = {}

        for valObj in valObjs:
            name, otype, stype, link, parent, last = valObj.info()
            if ":" in name:
                names = name.split(":")
                self._capsItems[":".join(names[level:])] = ValueOut(valObj)
                if level+1 == len(names):
                    continue

                if names[level] not in subVals:
                    subVals[names[level]] = []
                subVals[names[level]].append(valObj)
            else:
                self._capsItems[name] = ValueOut(valObj)

        super(self.__class__, self).__setattr__('_subVals', dict())
        for name, vals in subVals.items():
            self._subVals[name] = ValueOutSequence(vals, level+1)

    def __getattr__(self, name):
        if name in self._capsItems:
            return self._capsItems[name].value
        elif name in self._subVals:
            return self._subVals[name]

        raise AttributeError("No Out Value: {!r}".format(name))

    def __setattr__(self, name, value):
        if name in self._capsItems or name in self._subVals:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot set Out Value: {!r}".format(name))
        else:
            raise AttributeError("No Out Value: {!r}".format(name))

    def __delattr__(self, name):
        if name in self._capsItems or name in self._subVals:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot del Out Value: {!r}".format(name))
        else:
            raise AttributeError("No Out Value: {!r}".format(name))

    def __getitem__(self, name):
        if name in self._capsItems:
            return self._capsItems[name]
        elif name in self._subVals:
            return self._subVals[name]

        raise KeyError("No Out Value: {!r}".format(name))


#==============================================================================
## Defines a CAPS dynamic output Value Object
# Not a standalone class.
class ValueDynOut(object):

    __slots__ = ["_analysisObj", "_name"]

    def __init__(self, analysisObj, name):
        self._analysisObj = analysisObj
        self._name = name

    ## Property getter returns a copy the values stored in the CAPS Value Object
    @property
    def value(self):
        valObj = self._analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISDYNO, self._name)
        return valObj.getValue()

    @value.setter
    def value(self, val):
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot assign Out Value")

    ## Property returns the name of the CAPS Value Object
    @property
    def name(self):
        return self._name

    ## Property getter returns a copy the values stored in the CAPS Value Object
    def props(self):
        valObj = self._analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISDYNO, self._name)
        dim, pmtr, lfix, sfix, ntype = valObj.getValueProps()
        return {"dim":dim, "pmtr":pmtr, "lfix":lfix, "sfix":sfix, "ntype":ntype}

    ## Returns a string list of of the input Value Object names that can be used in \ref deriv
    def hasDeriv(self):
        valObj = self._analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISDYNO, self._name)
        return valObj.hasDeriv()

    ## Returns derivatives of the output Value Object
    #
    # \param name Name of the input Value Object to take derivative w.r.t.
    #        if name is None then a dictionary with all dervatives from \ref hasDeriv are returned
    def deriv(self, name=None):
        valObj = self._analysisObj.childByName(caps.oType.VALUE, caps.sType.ANALYSISDYNO, self._name)
        if name is not None:
            return valObj.getDeriv(name)
        else:
            names = valObj.hasDeriv()
            if names is None: return {}
            derivs = {}
            for name in names:
                derivs[name] = valObj.getDeriv(name)

            return derivs

#==============================================================================
## Defines a Sequence of CAPS dynamic output Value Objects
class ValueDynOutSequence(object):

    __slots__ = ["_analysisObj"]

    def __init__(self, analysisObj):
        super(self.__class__, self).__setattr__('_analysisObj', analysisObj)

    def __getattr__(self, name):
        return ValueDynOut(self._analysisObj, name).value

    def __setattr__(self, name, value):
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot set DynOut Value: {!r}".format(self._typeName, key))

    def __delattr__(self, name):
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot del DynOut Value: {!r}".format(self._typeName, key))

    def __getitem__(self, name):
        return ValueDynOut(self._analysisObj, name)

    def __setitem__(self, key, value):
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot set DynOut Value: {!r}".format(self._typeName, key))

    def __delitem__(self, key):
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Cannot del DynOut Value: {!r}".format(self._typeName, key))

    def _keys(self):
        size = self._analysisObj.size(caps.oType.VALUE, caps.sType.ANALYSISDYNO)
        keys = [None]*size
        for i in range(size):
            valObj = self._analysisObj.childByIndex(caps.oType.VALUE, caps.sType.ANALYSISDYNO, i+1)
            keys[i], otype, stype, link, parent, last = valObj.info()
        return keys

    def _capsItems(self):
        keys = self._keys()
        capsItems = {}
        for key in keys:
            capsItems[key] = ValueDynOut(self._analysisObj, key)
        return capsItems

    def __contains__(self, item):
        return item in self._keys()

    def __len__(self):
        return self._analysisObj.size(caps.oType.VALUE, caps.sType.ANALYSISDYNO)

    def __iter__(self):
        return self._capsItems().__iter__()

    def __reversed__(self):
        return self._capsItems().__reversed__()

    def __bool__(self):
        return bool(self._capsItems())

    ## Returns the keys of the ValueDynOutSequence
    def keys(self):
        return self._keys()

    ## Returns the values of the ValueDynOutSequence
    def values(self):
        keys = self._keys()
        for key in keys:
            values[i] = ValueDynOut(self._analysisObj, key)
        return values

    ## Returns the items of the ValueDynOutSequence
    def items(self):
        return self._capsItems().items()

#==============================================================================
## Defines a CAPS Analysis Object.
# Created via Problem.analysis.create().
#
# \param Analysis.geometry  \ref AnalysisGeometry instances representing the bodies associated with the analysis
# \param Analysis.input  \ref ValueInSequence of \ref ValueIn inputs
# \param Analysis.output \ref ValueOutSequence of \ref ValueOut outputs
# \param Analysis.attr   \ref AttrSequence of \ref ValueIn attributes
class Analysis(object):

    __slots__ = ["_analysisObj", "input", "output", "dynout", "geometry", "attr", "_analysisDir", "_name"]

    def __init__(self, analysisObj):
        self._analysisObj = analysisObj

        self.attr     = AttrSequence(self._analysisObj)
        self.input    = ValueInSequence (_values(self._analysisObj, caps.sType.ANALYSISIN ))
        self.output   = ValueOutSequence(_values(self._analysisObj, caps.sType.ANALYSISOUT))
        self.dynout   = ValueDynOutSequence(self._analysisObj)
        self.geometry = AnalysisGeometry(self._analysisObj)

        self._name, otype, stype, link, parent, last = self._analysisObj.info()

        self._analysisDir, unitSystem, major, minor, capsIntent, fnames, franks, fInOut, execution, cleanliness = self._analysisObj.analysisInfo()

#     def __setattr__(self, name, data):
#         raise AttributeError("Cannot add attribute: {!r}".format(name))
#
#     def __delattr__(self, name):
#         raise AttributeError("Cannot del attribute: {!r}".format(name)) #from None

    ## Run the pre-analysis function for the AIM.
    def preAnalysis(self):
        self._analysisObj.preAnalysis()

    ## Run the pre/exec/post functions for the AIM (if AIM execution is available).
    def runAnalysis(self):
        self._analysisObj.execute()

    ## Execute the Command Line String
    #    Notes:
    #    1. only needed when explicitly executing the appropriate analysis solver (i.e., not using the AIM)
    #    2. should be invoked after caps_preAnalysis and before caps_postAnalysis
    #    3. this must be used instead of the OS system call to ensure that journaling properly functions
    #
    # \param cmd  the command line string to execute
    # \param rpath  the relative path from the Analysis' directory or None (in the Analysis path)
    def system(self, cmd, rpath=None):
        self._analysisObj.system(cmd, rpath)

    ## Run post-analysis function for the AIM.
    def postAnalysis(self):
        self._analysisObj.postAnalysis()

    ## Property returns the path to the analysis directory
    @property
    def analysisDir(self):
        return self._analysisDir

    ## Property returns the name of the CAPS Analysis Object
    @property
    def name(self):
        return self._name

    ## Returns linked analyses that are dirty.
    #
    # \return A list of dirty analyses that need to be exeuted before executing this analysis. An empty list is returned if no
    # linked analyses are dirty.
    def dirty(self):

        dirtyAnalysisObjs = self._analysisObj.dirtyAnalysis()

        dirtyAnalysis = [Analysis(dirty) for dirty in dirtyAnalysisObjs]

        return dirtyAnalysis

    ## Gets analysis information for the analysis object.
    #
    # \param printInfo Print information to sceen (default - False).
    # \param **kwargs See below.
    #
    # \return Cleanliness state of analysis object or a dictionary containing analysis
    # information (infoDict must be set to True). For cleanliness state:
    # 0 = "Up to date",
    # 1 = "Dirty analysis inputs",
    # 2 = "Dirty geometry inputs",
    # 3 = "Both analysis and geometry inputs are dirty",
    # 4 = "New geometry",
    # 5 = "Post analysis required",
    # 6 = "Execution and Post analysis required"
    #
    # Valid keywords:
    # \param infoDict Return a dictionary containing analysis information instead
    # of just the cleanliness state (default - False).
    def info(self, printInfo=False, **kwargs):

        dirtyInfo = {0 : "Up to date",
                     1 : "Dirty analysis inputs",
                     2 : "Dirty geometry inputs",
                     3 : "Both analysis and geometry inputs are dirty",
                     4 : "New geometry",
                     5 : "Post analysis required",
                     6 : "Execution and Post analysis required"}

        name, otype, stype, link, parent, last = self._analysisObj.info()
        analysisPath, unitSystem, major, minor, capsIntent, fnames, franks, fInOut, execution, cleanliness = self._analysisObj.analysisInfo()

        if printInfo:
            numField = len(fnames)
            print("Analysis Info:")
            print("\tName           = ", name)
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

        if execution == 0:
            executionFlag = False
        else:
            executionFlag = True

        if kwargs.pop("infoDict", False):
            return {"name"          : name,
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

    ## Create a HTML dendrogram/tree of the current state of the analysis.
    # The HTML file relies on the open-source JavaScript library, D3, to visualize the data.
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
    # \param reverseMap Reverse the attribute map (default - False). See \ref attrMap for details.
    def createTree(self, filename = "name", **kwargs):

        if filename == "name":
            filename = self.name

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
    # \return Returns the reference to the OpenMDAO component object created.
    def createOpenMDAOComponent(self, inputVariable, outputVariable, **kwargs):

         print("Creating an OpenMDAO component for analysis, " + str(self.name))

         openMDAOComponent = createOpenMDAOComponent(self,
                                                     inputVariable,
                                                     outputVariable,
                                                     kwargs.pop("executeCommand", None),
                                                     kwargs.pop("inputFile", None),
                                                     kwargs.pop("outputFile", None),
                                                     kwargs.pop("setSensitivity", {}),
                                                     kwargs.pop("changeDir", True),
                                                     kwargs.pop("saveIteration", False))

         if kwargs.get("stdin", None) is not None:
             openMDAOComponent.stdin = kwargs.get("stdin", None)

         if kwargs.get("stdout", None) is not None:
             openMDAOComponent.stdout = kwargs.get("stdout", None)

         return openMDAOComponent

#=============================================================================
## Defines Analysis Geometry Object
#
# \param AnalysisGeometry.despmtr \ref ValueInSequence of \ref ValueIn CSM design parameters
# \param AnalysisGeometry.cfgpmtr \ref ValueInSequence of \ref ValueIn CSM configuration parameters
# \param AnalysisGeometry.conpmtr \ref ValueInSequence of \ref ValueIn CSM constant parameters
# \param AnalysisGeometry.outpmtr \ref ValueOutSequence of \ref ValueOut CSM outputs

class AnalysisGeometry(object):

    __slots__ = ['_analysisObj', 'despmtr', 'cfgpmtr', 'conpmtr', 'outpmtr']

    def __init__(self, analysisObj):
        super(AnalysisGeometry, self).__setattr__('_analysisObj', analysisObj)

        super(AnalysisGeometry, self).__setattr__("despmtr", ValueInSequence (_values(analysisObj.problemObj(), caps.sType.GEOMETRYIN)))
        super(AnalysisGeometry, self).__setattr__("cfgpmtr", ValueInSequence (_values(analysisObj.problemObj(), caps.sType.GEOMETRYIN, 1)))
        super(AnalysisGeometry, self).__setattr__("conpmtr", ValueInSequence (_values(analysisObj.problemObj(), caps.sType.GEOMETRYIN, 2)))
        super(AnalysisGeometry, self).__setattr__('outpmtr', ValueOutSequence(_values(analysisObj.problemObj(), caps.sType.GEOMETRYOUT)))

    def __setattr__(self, name, data):
        raise AttributeError("Cannot add attribute: {!r}".format(name))

    def __delattr__(self, name):
        raise AttributeError("Cannot del attribute: {!r}".format(name)) #from None

    ## Get dict of geometric bodies.
    #
    # \return Returns a dictionary of the bodies in the Analysis Object, as well a the capsLength unit.
    # Keys use the body "_name" attribute or "Body_#".
    def bodies(self):
        bodies = self._analysisObj.getBodies()
        nameList = []

        if len(bodies) < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "The number of bodies in the AIM is less than 1!!!" )

        body, units = self._analysisObj.problemObj().bodyByIndex(1)

        bodyDict = {}
        for i in range(len(bodies)):
            body = bodies[i]
            name = _getBodyName(i, nameList, body)

            bodyDict[name] = body

        return bodyDict, units

    ## Save the current geometry used by the AIM to a file.
    # \param filename File name to use when saving geometry file.
    # \param directory Directory where to save file. Default current working directory.
    # \param extension Extension type for file if filename does not contain an extension.
    # \param writeTess Write tessellations to the EGADS file (only applies to .egads extension)
    def save(self, filename, directory=os.getcwd(), extension=".egads", writeTess=True):

        filename = _createFilename(filename, directory, extension)

        flag = 1 if writeTess else 0
        self._analysisObj.writeGeometry(filename, flag)

    ## View the geometry associated with the analysis.
    #  If the analysis produces a surface tessellation, then that is shown.
    #  Otherwise the bodies are shown with default tessellation parameters.
    #  Note that the geometry must be built and will not autoamtically be built by this function.
    #
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param portNumber Port number to start the server listening on (default - 7681).
    def view(self, **kwargs):

        analysisDir, unitSystem, major, minor, capsIntent, fnames, franks, fInOut, execution, cleanliness = self._analysisObj.analysisInfo()

        egadsFile = os.path.join(analysisDir, "viewAnalysisGeometry.egads")
        relFile = os.path.relpath(egadsFile, os.getcwd())

        # use the shortest filename
        if len(relFile) < len(egadsFile):
            egadsFile = relFile

        self.save( egadsFile )

        # dissable the viewer (used for testing)
        if os.getenv("CAPS_BATCH"):
            os.remove(egadsFile)
            return

        egadsFile = _ocsmPath(egadsFile)
        os.system("serveCSM -port " + str(kwargs.pop("portNumber", 7681)) + " -outLevel 0 " + egadsFile)
        os.remove(egadsFile)
        os.remove("autoEgads.csm")

        return

        # Initialize the viewer class
        viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))

        # Get the bodies from the analysis
        bodies, units = self.bodies()

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
    def attrList(self, attributeName, **kwargs):
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
            raise caps.CAPSError(CAPS_BADVALUE, "Error: Invalid attribute level!")

        # How many bodies do we have
        numBody = self._analysisObj.size(caps.oType.BODIES, caps.sType.NONE)

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

            body, units = self._analysisObj.bodyByIndex(i+1)

            attributeList += _createAttributeValList(body,
                                                     attrLevel,
                                                     attributeName)

        try:
            return list(set(attributeList))
        except TypeError:
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
    def attrMap(self, getInternal = False, **kwargs):

        attributeMap = {}

        # Get the bodies
        bodies, units = self.bodies()

        if len(bodies) < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating attribute map"
                                                         + " for aim " + str(self.aimName)
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        # Loop through bodies
        for name, body in bodies.items():
            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            attributeMap[name] = _createAttributeMap(body, getInternal)

        if kwargs.pop("reverseMap", False):
            return _reverseAttributeMap(attributeMap)
        else:
            return attributeMap



#==============================================================================
# Legacy depricated functions
class capsAnalysis(Analysis):

    __slots__ = ['officialName']

    def __init__(self, problem, aim, altName=None, unitSystem=None, capsIntent=None, analysisDir=None, parents=None, copyAIM=None, analysisObj=None):

        if analysisObj is not None:
            super(capsAnalysis,self).__init__(analysisObj)
            return

        deprecate("capsAnalysis.__init__","'Problem.analysis.create'")

        # Check problem object
        if not isinstance(problem, capsProblem) and not isinstance(problem, CapsProblem):
            raise TypeError("Provided problem is not a capsProblem or CapsProblem class")

        # Check the analysis dictionary for redundant names - need unique names
        def checkAnalysisDict(name, printInfo=True):
            if name in problem.analysis.keys():
                if printInfo:
                    print("\nAnalysis object with the name, " + str(name) +  ", has already been initilized")
                return True
            else:
                return False

        currentAIM = aim
        if altName is None: altName = currentAIM
        if currentAIM is None: currentAIM = ""

        if not analysisDir or analysisDir == ".":
            print("No analysis directory provided - defaulting to altName")
            analysisDir = altName

        if analysisDir is not None:
            ## Adjust legacy directory names
            analysisDir = analysisDir.replace(".","_")
            analysisDir = analysisDir.replace("/","_")
            analysisDir = analysisDir.replace("\\","_")

            ## Analysis directory of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$analysis).
            # If the directory does not exist it will be made automatically.
            analysisDir = analysisDir + "_" + currentAIM

        if problem.geometry == None:
            raise caps.CAPSError(caps.CAPS_NULLOBJ, msg = "while loading aim"
                                                        + "\nError: A *.csm, *.caps, or *.egads file has not been loaded!")

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


        if analysisDir:
            root = problem._problemObj.getRootPath()
            fullPath = os.path.join(root,analysisDir)

            # Check to see if directory exists
            if not os.path.isdir(fullPath):
                #print("Analysis directory does not currently exist - it will be made automatically")
                os.makedirs(fullPath)

            # Remove the top directory as cCAPS will create it
            if os.path.isdir(fullPath):
                shutil.rmtree(fullPath)


        # Create the new analysisObj
        if copyAIM is not None:
            analysisObj = problem.analysis._capsItems[copyAIM]._analysisObj.dupAnalysis(name=analysisDir)
        else:
            analysisObj = problem._problemObj.makeAnalysis(aim,
                                                           name=analysisDir,
                                                           unitSys=unitSystem,
                                                           intent=capsIntent,
                                                           execute=0)

        # Store the analysisObj
        super(capsAnalysis,self).__init__(analysisObj)

        # Store capsAnalysis in the AnalysisSequence
        problem.analysis._capsItems[altName] = self

        self.officialName = currentAIM
        self._name = altName

        if parents:
            if not isinstance(parents,list):
                parents = [parents]

            inNames = self.input.keys()

            inVal = None

            if "Mesh" in inNames:
                inVal = self.input["Mesh"]
            elif "Surface_Mesh" in inNames:
                inVal = self.input["Surface_Mesh"]
            elif "Volume_Mesh" in inNames:
                inVal = self.input["Volume_Mesh"]

            if inVal is not None:
                for parent in parents:
                    if isinstance(parent, capsAnalysis):
                        parent = parent.name
                    elif not isinstance(parent,str):
                        raise TypeError("Provided parent is not a capsAnalysis or CapsAnalysis class or a string name")

                    if parent not in problem.analysis:
                        raise caps.CAPSError(caps.CAPS_NULLOBJ, msg = "while loading aim: "
                                                                    + "Unable to find parents")

                    outNames = problem.analysis[parent].output.keys()

                    if "Surface_Mesh" in outNames:
                        outVal = problem.analysis[parent].output["Surface_Mesh"]
                    elif "Area_Mesh" in outNames:
                        outVal = problem.analysis[parent].output["Area_Mesh"]
                    elif "Volume_Mesh" in outNames:
                        outVal = problem.analysis[parent].output["Volume_Mesh"]
                    else:
                        continue

                    inVal.link(outVal)


    @deprecated("'Analysis.input[\"varname\"].value'")
    def setAnalysisVal(self, varname, value, units=None):

        # This is a complete HACK, but consistent with previous lack of using units...
        if units is None:
            units = self.input[varname]._valObj.getValueUnit()

        if isinstance(value,list) and \
           len(value) > 0 and \
           isinstance(value[0], tuple) and \
           len(value[0]) == 2 and \
           isinstance(value[0][0], str):
            data = {}
            for i in range(len(value)):
                data[value[i][0]] = value[i][1]
            value = data

        if isinstance(value,tuple) and \
           len(value) == 2 and \
           isinstance(value[0], str):
            value = {value[0]: value[1]}

        self.input[varname].value = caps._withUnits(value, units)

    @deprecated("'Analysis.input[\"varname\"].value'")
    def getAnalysisVal(self, varname = None, **kwargs):
        namesOnly = kwargs.pop("namesOnly", False)
        units = kwargs.pop("units", None)

        if namesOnly:
            return self.input.keys()

        if varname is None:
            values = {}
            for i in self.input.keys():
                values[i] = self.input[i].value
                if isinstance(values[i], caps.Quantity):
                    values[i] = values[i]._value

            return values

        value = self.input[varname].value

        if isinstance(value,dict):
            data = []
            for key in sorted(value.keys()):
                data.append( (key, value[key]) )
            value = data
            if len(value) == 1: value = value[0]

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

    @deprecated("'Analysis.output[\"varname\"].value'")
    def getAnalysisOutVal(self, varname = None, **kwargs):
        namesOnly = kwargs.pop("namesOnly", False)
        units = kwargs.pop("units", None)

        if namesOnly:
            return self.output.keys()

        stateOfAnalysis = self.info(printInfo=False,infoDict=True)

        if stateOfAnalysis["status"] != 0 and stateOfAnalysis["executionFlag"] == False:
            raise caps.CAPSError(caps.CAPS_DIRTY, msg = "while getting analysis out variables"
                                                      + "Analysis state isn't clean (state = " + str(stateOfAnalysis) + ") can only return names of ouput variables; set namesOnly keyword to True")

        if varname is None:
            values = {}
            for i in self.output.keys():
                values[i] = self.output[i].value
                if isinstance(values[i], caps.Quantity):
                    values[i] = values[i]._value

            return values

        value = self.output[varname].value

        if isinstance(value,dict):
            data = []
            for key in sorted(value.keys()):
                data.append( (key, value[key]) )
            value = data
            if len(value) == 1: value = value[0]

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

    @deprecated("'Analysis.preAnalysis'")
    def aimPreAnalysis(self):
        self.preAnalysis()

    @deprecated("'Analysis.postAnalysis'")
    def aimPostAnalysis(self):
        self.postAnalysis()

    @deprecated("'Analysis.info'")
    def getAnalysisInfo(self, printInfo=True, **kwargs):
        return self.info(printInfo, **kwargs)

    @property
    def aimName(self):
        deprecate("capsAnalysis.aimName","'name'")
        return self.name

    @deprecated("'Analysis.attr.create'")
    def addAttribute(self, name, data):
        self.attr.create(name,data,True)

    @deprecated("'Analysis.attr[\"name\"].value'")
    def getAttribute(self, name):
        value = self.attr[name].value

        if isinstance(value,caps.Quantity):
            value = value._value

        return value

    @deprecated("'Analysis.geometry.attrList'")
    def getAttributeVal(self, attributeName, **kwargs):

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
            kwargs["attrLevel"] = attrLevel

        return self.geometry.attrList(attributeName, **kwargs)

    @deprecated("'Analysis.geometry.attrMap'")
    def getAttributeMap(self, getInternal = False, **kwargs):

        attributeMap = {}

        # How many bodies do we have
        numBody = self._analysisObj.size(caps.oType.BODIES, caps.sType.NONE)

        if numBody < 1:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while creating attribute map"
                                                         + " for aim " + str(self.aimName)
                                                         + "\nThe number of bodies in the problem is less than 1!!!" )

        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i

            body, units = self._analysisObj.bodyByIndex(i+1)

            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            attributeMap["Body " + str(i + 1)] = _createAttributeMap(body, getInternal)

        if kwargs.pop("reverseMap", False):
            return _reverseAttributeMap(attributeMap)
        else:
            return attributeMap


    @deprecated("'Analysis.geometry.save'")
    def saveGeometry(self, filename, directory=os.getcwd(), extension=".egads"):
        return self.geometry.save(filename, directory, extension)

    @deprecated("'Analysis.geometry.view'")
    def viewGeometry(self, **kwargs):
        return self.geometry.view(**kwargs)


    @deprecated("'Analysis.output[\"varname\"].dot()'")
    def getSensitivity(self, inputVar, outputVar):
        JSONin = jsonDumps({"mode": "Sensitivity",
                            "outputVar": str(outputVar),
                            "inputVar" : str(inputVar)})

        return self.backDoor(JSONin)

    @deprecated("'None'")
    def aimBackDoor(self, JSONin):
        self.backDoor(JSONin)

    @deprecated("'None'")
    def backDoor(self, JSONin):

        if not isinstance(JSONin, str):
            JSONin = jsonDumps(JSONin)

        JSONout = self._analysisObj.AIMbackdoor(JSONin)

        if JSONout:
            temp =  jsonLoads( JSONout )
            return float(temp["sensitivity"])

        return None

    @deprecated("'Analysis.bodies'")
    def getBoundingBox(self):

        boundingBox = {}

        bodies, units = self.geometry.bodies()

        if bodies is None:
            raise caps.CAPSError(caps.CAPS_SOURCEERR, msg = "No bodies available to retrieve the bounding box!")

        for name, body in bodies.items():
            box  = body.getBoundingBox()

            boundingBox[name] = [box[0][0], box[0][1], box[0][2],
                                 box[1][0], box[1][1], box[1][2]]

        return boundingBox

## Pythonic alias to capsAnalysis (deprecated)
CapsAnalysis = capsAnalysis

#==============================================================================
## Defines a Sequence of CAPS Analysis Objects
class AnalysisSequence(Sequence):

    __slots__ = ['_problemObj', "_deprecated"]

    def __init__(self, problem):
        super(self.__class__, self).__init__('Analysis')

        super(Sequence, self).__setattr__('_problemObj', problem._problemObj)
        super(Sequence, self).__setattr__('_deprecated', isinstance(problem, capsProblem))

        nAnalyses = self._problemObj.size(caps.oType.ANALYSIS, caps.sType.NONE)

        for i in range(nAnalyses):
            analysisObj = self._problemObj.childByIndex(caps.oType.ANALYSIS, caps.sType.NONE, i+1)
            name, otype, stype, link, parent, last = analysisObj.info()
            if self._deprecated:
                self._capsItems[name] = capsAnalysis(problem=None,aim=None,analysisObj=analysisObj)
            else:
                self._capsItems[name] = Analysis(analysisObj)

    ## Create a CAPS Analysis Object.
    #
    # \param aim Name of the AIM module
    #
    # \param name Name (e.g. key) of the Analysis Object.
    # Must be unique if specified.
    # If None, the defalt is aim+str(instanceCount) where
    # instanceCount is the count of the existing 'aim' instances.
    #
    # \param capsIntent Analysis intention in which to invoke the AIM.
    #
    # \param unitSystem See AIM documentation for usage.
    #
    # \param autoExec If false dissable any automatic execution of the AIM.
    #
    # \return The new Analysis Object is added to the sequence and returned
    def create(self, aim, name=None, capsIntent=None, unitSystem=None, autoExec=True):

        analysisObj = self._problemObj.makeAnalysis(aim,
                                                    name,
                                                    unitSys=unitSystem,
                                                    intent=capsIntent,
                                                    execute = 1 if autoExec else 0)

        if name is None:
            name, otype, stype, link, parent, last = analysisObj.info()

        if self._deprecated:
            self._capsItems[name] = capsAnalysis(problem=None,aim=None,analysisObj=analysisObj)
        else:
            self._capsItems[name] = Analysis(analysisObj)

        return self._capsItems[name]

    ## Create a copy of an CAPS Analysis Object.
    #
    # \param src Name of the source Analysis Object or an Analysis Object
    # \param name Name of the new Analysis Object copy
    def copy(self, src, name=None):

        if isinstance(src, Analysis):
            src = src.name

        analysisObj = self._capsItems[src]._analysisObj.dupAnalysis(name)

        if self._deprecated:
            self._capsItems[name] = capsAnalysis(problem=None,aim=None,analysisObj=analysisObj)
        else:
            self._capsItems[name] = Analysis(analysisObj)

        return self._capsItems[name]

    ## Returns analyses that are dirty.
    #
    # \return A list of dirty analyses. An empty list is returned if no
    # analyses are dirty.
    def dirty(self):

        dirtyAnalysisObjs = self._problemObj.dirtyAnalysis()

        dirtyAnalysis = [Analysis(dirty) for dirty in dirtyAnalysisObjs]

        return dirtyAnalysis

#==============================================================================
## Defines a CAPS Bound Object.
# Created via Problem.bound.create().
#
# \param Bound.vertexSet \ref VertexSetSequence of \ref VertexSet instances
# \param Bound.attr      \ref AttrSequence of \ref ValueIn attributes
class Bound(object):

    __slots__ = ["_boundObj", "_name", "vertexSet", "attr"]

    def __init__(self, boundObj):
        self._boundObj = boundObj
        self.attr      = AttrSequence(self._boundObj)
        self.vertexSet = VertexSetSequence(self._boundObj)

        self._name, otype, stype, link, parent, last = self._boundObj.info()

    ## Property returns the name of the CAPS Bound Object
    @property
    def name(self):
        return self._name

    ## Closes the bound indicating it's complete
    def close(self):
        self._boundObj.closeBound()

    def complete(self):
        self._boundObj.closeBound()

    ## Gets information for the bound object.
    #
    # \param printInfo Print information to sceen if True.
    # \param **kwargs See below.
    #
    # \return State of bound object.
    #
    # Valid keywords:
    # \param infoDict Return a dictionary containing bound information instead
    # of just the state (default - False)
    def info(self, printInfo=False, **kwargs):

        boundInfo = {-2 : "multiple BRep entities - Error in reparameterization!y",
                     -1 : "Open",
                      0 : "Empty and Closed",
                      1 : "Single BRep entity",
                      2 : "Multiple BRep entities"}

        state, dimension, limits = self._boundObj.boundInfo()

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


#==============================================================================
class capsBound(Bound):

    __slots__ = ['boundObj','boundName','variables','dataSetSrc','dataSetDest']

    def __init__(self, problem, **kwargs):

        # Check problem object
        if not isinstance(problem, capsProblem):
            raise TypeError("Provided problem is not a capsProblem class")

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
                if i not in problem.analysis.keys():
                    raise caps.CAPSError(caps.CAPS_BADVALUE, msg="AIM name = " + i + " not found in analysis dictionary")

        if "capsTransfer" in kwargs.keys():
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="The keyword, capsTransfer, should be changed to capsBound")

        if "variableRank" in kwargs.keys():
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Keyword variableRank is not necessary. Please remove it from your script.")

        boundName = kwargs.pop("capsBound", None)
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

            fieldName = problem.analysis[src].getAnalysisInfo(printInfo = False, infoDict = True)["fieldName"]

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

        # Set the bound in the problem
        problem.dataBound[boundName] = self

        ## Bound/transfer name used to set up the data bound.
        self.boundName = boundName

        ## List of variables in the bound.
        self.variables = varName

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
        self.boundObj = problem._problemObj.makeBound(2, boundName)
        super(capsBound,self).__init__(self.boundObj)

        problem.bound._capsItems[boundName] = self

        # Make all of our vertex sets
        for i in set(aimSrc+aimDest):
            #print"Vertex Set - ", i
            self.vertexSet.create(problem.analysis[i])

        if len(varName) != len(aimSrc):
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg="ERROR: varName and aimSrc must all be the same length")

        # Make all of our source data sets
        for var, src in zip(varName, aimSrc):
            #print"DataSetSrc - ", var, src
            self.dataSetSrc[var] = self.vertexSet[src].dataSet.create(var, caps.fType.FieldOut)

        # Make destination data set
        if len(aimDest) > 0:

            if len(varName) != len(aimDest) or \
               len(varName) != len(transferMethod) or \
               len(varName) != len(aimSrc) or \
               len(varName) != len(initValueDest):
                raise caps.CAPSError(caps.CAPS_BADVALUE, msg="ERROR: varName, aimDest, transferMethod, aimSrc, and initValueDest must all be the same length")

            for var, dest, src, initVal in zip(varName, aimDest, aimSrc, initValueDest):
                #print"DataSetDest - ", var, dest, method

                self.dataSetDest[var] = self.vertexSet[dest].dataSet.create(var, caps.fType.FieldIn, initVal)

        for var, method in zip(varName, transferMethod):
            self.dataSetDest[var].link(self.dataSetSrc[var], method)

        # Close bound
        self.boundObj.closeBound()


    @deprecated("'Bound.info'")
    def getBoundInfo(self, printInfo=True, **kwargs):

        boundInfo = {-2 : "multiple BRep entities - Error in reparameterization!y",
                     -1 : "Open",
                      0 : "Empty and Closed",
                      1 : "Single BRep entity",
                      2 : "Multiple BRep entities"}

        state, dimension, limits = self._boundObj.boundInfo()

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

        return self.info(printInfo, **kwargs)

    @deprecated("'None'")
    def fillVertexSets(self):
        pass

    @deprecated("'None'")
    def executeTransfer(self, variableName = None):

        if not self.dataSetDest:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "No destination data exists in bound!")

        self.getDestData(variableName)

    @deprecated("'vertexSet[\"name\"].dataSet[\"variableName\"].getData'")
    def getDestData(self, variableName = None):

        if not self.dataSetDest:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "No destination data exists in bound!")

        if variableName is None:
            print("No variable name provided - defaulting to variable - " + self.variables[0])
            variableName = self.variables[0]

        elif variableName not in self.variables:
            print("Unreconized varible name " + str(variableName))
            return

        self.dataSetDest[variableName].getData()

    @deprecated("'vertexSet[\"name\"].dataSet[\"variableName\"].view'")
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
                raise ImportError("Error: Unable to import matplotlib - viewing the data is not possible")

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

    @deprecated("'vertexSet[\"name\"].dataSet[\"variableName\"].writeTecplot'")
    def writeTecplot(self, filename, variableName=None):

        if filename is None:
            print("No filename provided")
            return

        if "." not in filename:
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

            self.dataSetSrc[var].writeTecplot(file=file)

            if self.dataSetDest:

                self.dataSetDest[var].writeTecplot(file=file)

        # Close the file
        file.close()

## Pythonic alias to capsBound (deprecated)
CapsBound = capsBound


#==============================================================================
## Defines a Sequence of CAPS Bound Objects
class BoundSequence(Sequence):

    __slots__ = ['_problemObj']

    def __init__(self, problemObj):
        super(self.__class__, self).__init__('Bound')

        super(Sequence, self).__setattr__('_problemObj', problemObj)

        nBound = problemObj.size(caps.oType.BOUND, caps.sType.NONE)

        for i in range(nBound):
            boundObj = problemObj.childByIndex(caps.oType.BOUND, caps.sType.NONE, i+1)
            name, otype, stype, link, parent, last = boundObj.info()
            self._capsItems[name] = Bound(boundObj)

    ## Create a CAPS Bound Object.
    #
    # \param capsBound The string value of the capsBound geometry attributes
    #
    # \param dim The dimension of the bound
    #
    # \return The new Bound Object is added to the sequence and returned
    def create(self, capsBound, dim = 2):

        boundObj = self._problemObj.makeBound(dim, capsBound)

        self._capsItems[capsBound] = Bound(boundObj)
        return self._capsItems[capsBound]

#==============================================================================
## Defines a CAPS VertexSet Object.
# Created via Bound.vertexSet.create().
#
# \param VertexSet.dataSet \ref DataSetSequence of \ref DataSet instances
# \param VertexSet.attr    \ref AttrSequence of \ref ValueIn attributes
class VertexSet(object):

    __slots__ = ["_vertexSetObj", "_name", "attr", "dataSet"]

    def __init__(self, vertexSetObj, fieldRank):
        self._vertexSetObj = vertexSetObj
        self.attr          = AttrSequence(self._vertexSetObj)
        self.dataSet       = DataSetSequence(self._vertexSetObj, fieldRank)

        self._name, otype, stype, link, parent, last = self._vertexSetObj.info()

    ## Property returns the name of the CAPS VertexSet Object
    @property
    def name(self):
        return self._name

    ## Executes caps_getTriangles on data set's vertex set to retrieve the connectivity (triangles only) information
    # for the data set.
    # \return Optionally returns a list of lists of connectivity values
    # (e.g. [ [node1, node2, node3], [node2, node3, node7], etc. ] ) and a list of lists of data connectivity (not this is
    # an empty list if the data is node-based) (eg. [ [node1, node2, node3], [node2, node3, node7], etc. ]
    def getDataConnect(self):
        Gtris, Gsegs, Dtris, Dsegs = self._vertexSetObj.getTriangles()
        return Gtris, Gsegs, Dtris, Dsegs

#==============================================================================
## Defines a Sequence of CAPS Bound Objects
class VertexSetSequence(Sequence):

    __slots__ = ['_boundObj']

    def __init__(self, boundObj):
        super(self.__class__, self).__init__('VertexSet')

        super(Sequence, self).__setattr__('_boundObj', boundObj)

        nVertexSet = self._boundObj.size(caps.oType.VERTEXSET, caps.sType.NONE)

        for i in range(nVertexSet):
            vertexSetObj = boundObj.childByIndex(caps.oType.VERTEXSET, caps.sType.NONE, i+1)
            name, otype, stype, link, parent, last = vertexSetObj.info()
            self._capsItems[name] = VertexSet(vertexSetObj)

    ## Create a CAPS VertexSet Object.
    #
    # \param analysis A CAPS Analysis Object or the string name of an Analsysis Object instance
    #
    # \param vname Name of the VertexSet (same as the Analysis Object if None)
    #
    # \return The new VertexSet Object is added to the sequence and returned
    def create(self, analysis, vname=None):

        if isinstance(analysis,Analysis):
            analysisObj = analysis._analysisObj
        elif isinstance(analysis,str):
            analysisObj = self._boundObj.problemObj().childByName(caps.oType.ANALYSIS, caps.sType.NONE, analysis)
        else:
            raise TypeError("'analysis' must be an Analsysis Object or name of an Analysis Object! type(analysis)={!r}".format(type(analysis)))

        vertexSetObj = self._boundObj.makeVertexSet(analysisObj, vname)

        if isinstance(analysis,capsAnalysis):
            name = analysis.name
        else:
            name, otype, stype, link, parent, last = vertexSetObj.info()

        analysisPath, unitSystem, major, minor, capsIntent, fnames, franks, fInOut, execution, cleanliness = analysis._analysisObj.analysisInfo()

        self._capsItems[name] = VertexSet(vertexSetObj, dict(zip(fnames, franks)))
        return self._capsItems[name]

#==============================================================================
## Defines a CAPS DataSet Object.
# Created via VertexSet.dataSet.create().
#
# \param DataSet.attr \ref AttrSequence of \ref ValueIn attributes
class DataSet(object):

    __slots__ = ["_dataSetObj", "_vertexSetObj", "_name", "attr"]

    def __init__(self, dataSetObj):
        self._dataSetObj = dataSetObj
        self.attr        = AttrSequence(self._dataSetObj)

        self._name, otype, stype, link, self._vertexSetObj, last = self._dataSetObj.info()

    ## Property returns the name of the CAPS DataSet Object
    @property
    def name(self):
        return self._name

    ## Executes caps_getData on data set object to retrieve data set variable.
    # \return Optionally returns a list of data values. Data with a rank greater than 1 returns a list of lists (e.g.
    # data representing a displacement would return [ [Node1_xDisplacement, Node1_yDisplacement, Node1_zDisplacement],
    # [Node2_xDisplacement, Node2_yDisplacement, Node2_zDisplacement], etc. ]
    def data(self):
        return self._dataSetObj.getData()

    @deprecated("'DataSet.data'")
    def getData(self):
        return self.data()

    ## Executes caps_getData on data set object to retrieve XYZ coordinates of the data set.
    # \return Optionally returns a list of lists of x,y, z values (e.g. [ [x2, y2, z2], [x2, y2, z2],
    # [x3, y3, z3], etc. ] )
    def xyz(self):
        xyzDataObj = self._vertexSetObj.childByName(caps.oType.DATASET, caps.sType.NONE, "xyz")
        return xyzDataObj.getData()

    @deprecated("'DataSet.xyz'")
    def getDataXYZ(self):
        return self.xyz()

    ## Executes caps_getTriangles on data set's vertex set to retrieve the connectivity (triangles only) information
    # for the data set.
    # \return Optionally returns a list of lists of connectivity values
    # (e.g. [ [node1, node2, node3], [node2, node3, node7], etc. ] ) and a list of lists of data connectivity (not this is
    # an empty list if the data is node-based) (eg. [ [node1, node2, node3], [node2, node3, node7], etc. ]
    def connectivity(self):
        Gtris, Gsegs, Dtris, Dsegs = self._vertexSetObj.getTriangles()
        return Gtris, Gsegs, Dtris, Dsegs

    @deprecated("'DataSet.connectivity'")
    def getDataConnect(self):
        Gtris, Gsegs, Dtris, Dsegs = self._vertexSetObj.getTriangles()
        return Gtris, Gsegs, Dtris, Dsegs

    ## Link this DataSet to an other CAPS DataSet Object
    # \param source The source DataSEt Object
    # \param dmethod Transfter method: dMethod.Interpolate or "Interpolate", tMethod.Conserve or "Conserve"
    def link(self, source, dmethod = caps.dMethod.Interpolate):

        dMethods = {"interpolate" : caps.dMethod.Interpolate,
                    "conserve"    : caps.dMethod.Conserve}

        if not isinstance(dmethod, (int,float)):
            if dmethod.lower() in dMethods.keys():
                dmethod = dMethods[dmethod.lower()]
            else:
                raise caps.CAPSError(caps.CAPS_BADVALUE, "Unknown dmethod = {!r}".format(dmethod))

        if not (isinstance(source, DataSet)):
            raise caps.CAPSError(caps.CAPS_BADVALUE, "source must be a DataSet Object")

        self._dataSetObj.linkDataSet(source._dataSetObj, dmethod)

    ## Visualize data set.
    # The function currently relies on matplotlib to plot the data.
    # \param fig Figure object (matplotlib::figure) to append image to.
    # \param numDataSet Number of data sets in \ref DataSet.view$fig.
    # \param dataSetIndex Index of data set being added to \ref DataSet.view$fig.
    # \param **kwargs See below.
    #
    # Valid keywords:
    # \param filename Save image(s) to file specified (default - None).
    # \param colorMap  Valid string for a, matplotlib::cm,  colormap (default - 'Blues').
    # \param showImage Show image(s) (default - True).
    # \param title Set a custom title on the plot (default - VertexSet= 'name', DataSet = 'name', (Var. '#') ).
    def view(self, fig=None, numDataSet=1, dataSetIndex=0, **kwargs):

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
            cMap = colorMap.get_cmap(cMap)
        except ValueError:
            print("Colormap ",  cMap, "is not recognized. Defaulting to 'Blues'")
            cMap = colorMap.get_cmap("Blues")

        # Initialize values
        colorArray = []

        data = self.data()
        xyz = self.xyz()

        triConn, Gsegs, Dtris, dataConn = self._vertexSetObj.getTriangles()

        numData = len(data)
        numTri = len(triConn)

        try:
            dataRank = len(data[0])
        except TypeError:
            dataRank = 1

        if not ((numData == len(xyz)) or (numData == numTri)):
            raise caps.CAPSError(caps.CAPS_MISMATCH, msg = "viewData only supports node-centered or cell-centered data!")

        if dataConn is not None:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "viewData does NOT currently support cell-centered based data sets!")

        if fig is None:
            fig = plt.figure() #figsize=(9,9)

        cellCenter = (numData == numTri)

        ax = []
        for j in range(dataRank):

            if dataRank > 1 or numDataSet > 1:
                ax.append(fig.add_subplot(dataRank, numDataSet, numDataSet*j + 1 + dataSetIndex, projection='3d'))
            else:
                ax.append(fig.gca(projection='3d'))

            if dataRank == 1:
                colorArray = data
            else:
                colorArray = [None]*numData
                for i in range(numData):
                    colorArray[i] = data[i][j]

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

            for i in range(numTri):
                tri = []

                index = triConn[i][0]-1
                tri.append(xyz[index])

                if zeroP == False and not cellCenter:
                    p = (colorArray[index] - minColor)/dColor

                index = triConn[i][1]-1
                tri.append(xyz[index])

                if zeroP == False and not cellCenter:
                    p += (colorArray[index] - minColor)/dColor

                index = triConn[i][2]-1
                tri.append(xyz[index])

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

            vertexSetName, otype, stype, link, parent, last = self._vertexSetObj.info()

            #ax[-1].autoscale_view()
            title = "DataSet = " + str(self.name)
            if dataRank > 1:
                title += " (Var. " + str(j+1) + ")"

            if j == 0:
                title = "VertexSet = " + str(vertexSetName) + \
                        "\n" + title

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

        #fig.subplots_adjust(hspace=0.1) #left=0.0, bottom=0.0, right=1.0, top=.93, wspace=0.0 )

        if showImage:
            plt.show()

        if filename is not None:
            print("Saving figure - ", filename)
            fig.savefig(filename)

    ## Write data set to a Tecplot compatible data file. A triagulation of the data set will be used
    # for the connectivity.
    # \param file Optional open file object to append data to. If not provided a
    # filename must be given via the keyword arguement \ref DataSet.writeTecplot$filename.
    # \param filename Write Tecplot file with the specified name.
    def writeTecplot(self, filename=None, file=None):

        # If no file
        if file == None:

            if filename == None:
                raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while writing Tecplot file for data set - " + str(self.dataSetName)
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

        data = self.data()
        if not data: # Do we have empty data
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "the data set is empty for - "
                                                         + str(self.name))
        else:
            if not isinstance(data[0], tuple): # Only returns a list if rank 1
                temp = []
                for d in data:
                    temp.append([d])
                data = temp

        xyz = self.xyz()

        connectivity, Gsegs, dconnectivity, Dsegs = self.connectivity()

        #TODO: dconnectivity does not appear to indicate a cell centered data set?
        if dconnectivity:
            raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "writeTecplot does NOT currently support dconnectivity based data sets!")

        if not ((len(data) == len(xyz)) or (len(data) == len(connectivity))):
            raise caps.CAPSError(caps.CAPS_MISMATCH, msg = "writeTecplot only supports node-centered or cell-centered data!")

        vertexSetName, otype, stype, link, parent, last = self._vertexSetObj.info()

        title = 'TITLE = "VertexSet: ' + vertexSetName + '"\n'
        fp.write(title)

        title = 'VARIABLES = "x", "y", "z"'
        if isinstance(data[0], tuple):
            for i in range(len(data[0])):
                title += ', "' + self.name + '_' + str(i) + '"'
        else:
            title += ', "' + self.name + '"'
        title += "\n"
        fp.write(title)

        #numNode = self.numNode if len(data) == len(xyz) else 3*len(connectivity)
        numNode = len(xyz)

        title = "ZONE T = " + \
                '"VertexSet: ' + vertexSetName + '"' +\
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

#==============================================================================
## Defines a Sequence of CAPS DataSet Objects
class DataSetSequence(Sequence):

    __slots__ = ['_vertexSetObj', '_fieldRank']

    def __init__(self, vertexSetObj, fieldRank):
        super(self.__class__, self).__init__('DataSet')

        super(Sequence, self).__setattr__('_vertexSetObj', vertexSetObj)
        super(Sequence, self).__setattr__('_fieldRank', fieldRank)

        nDataSet = vertexSetObj.size(caps.oType.DATASET, caps.sType.NONE)

        for i in range(nDataSet):
            dataSetObj = vertexSetObj.childByIndex(caps.oType.DATASET, caps.sType.NONE, i+1)
            name, otype, stype, link, parent, last = dataSetObj.info()
            self._capsItems[name] = DataSet(dataSetObj)

    ## Create a CAPS DataSet Object.
    #
    # \param dname The name of the data set
    #
    # \param ftype The field type (FieldIn, FieldOut, GeomSens, TessSens, User)
    #
    # \param init Inital value assiged to the DataSet. Length must be consistent with the rank.
    #
    # \param rank The rank of the data set (only needed for un-connected data set)
    #
    # \return The new DataSet Object is added to the sequence and returned
    def create(self, dname, ftype, init=None, rank=None):

        if dname not in self._fieldRank:
            for key in self._fieldRank:
                if '#' not in key: continue
                pattern = key.replace('#',"[0-9]+")
                if re.search(pattern, dname):
                    rank = self._fieldRank[key]
                    break

            if rank is None:
                raise caps.CAPSError(caps.CAPS_BADVALUE, "No DataField: {!r}\n Available fields: {!r}\n A rank must be specified.".format(dname, self.fields()))
        else:
            rank = self._fieldRank[dname]

        dataSetObj = self._vertexSetObj.makeDataSet(dname, ftype, rank)

        if init:
            dataSetObj.initDataSet(init)

        self._capsItems[dname] = DataSet(dataSetObj)
        return self._capsItems[dname]

    ## Returns a list of the fields in the Analysis Object associated with this DataSet
    def fields(self):
        return list(self._fieldRank.keys())
