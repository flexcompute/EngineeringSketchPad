#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

from struct import unpack as structUnpack
from struct import pack as structPack
from math import cos, sin, sqrt
from copy import deepcopy

from json import dumps as jsonDumps
from json import loads as jsonLoads

cdef class capsTecplotZone:
    
    cdef:
        
        readonly object zoneName
        
        public numElement
        public numNode
        
        public i
        public j
        public k
       
        readonly object variable # Needed to keep track of the order in which the data was stored
        readonly object data
        readonly object connectivity
        
        public object dataPacking
        public object zoneType
        public fileName
        
        readonly _lineOffset
    
    def __cinit__(self, name,  variables, fileName=None):
        
        self.zoneName = name
        
        self.numElement = None
        self.numNode = None
        
        self.i = None
        self.j = None
        self.k = None
        
        self.fileName = fileName
        
        self._lineOffset = 0
        
        if not isinstance(variables, list):
             self.variable = [variables]
        else:
             self.variable = variables
                
        self.data = {}
        for i in self.variable:
            if i: # Deal with empty characters in list 
                self.data[i] = [] # Dictionary [variable] = list of data 
                
        self.connectivity = []# List of list element connectivity 
            
        self.zoneType = "Ordered".upper()
        
        self.dataPacking = "Block".upper()
        
    def _switchIJK(self): # Help to handle archaic version of ASCII formated file 
        
        if "FE" in self.zoneType:
            self.i = None
            self.j = None
            self.k = None
       
        else:
            self.i = self.numNode
            self.j = self.numElement
        
            self.numNode = None
            self.numElement = None
            
    def _checkZoneType(self, zoneType):
        
        if "ORDERED" in zoneType:
            self.zoneType = "ORDERED"
            
        elif "FELINESEG" in zoneType:
            self.zoneType = "FELINESEG"
            
        elif "FETRIANGLE" in zoneType:
            self.zoneType = "FETRIANGLE" 
            
        elif "FEQUADRILATERAL" in zoneType:
            self.zoneType = "FEQUADRILATERAL"
            
        elif "FETETRAHEDRON" in zoneType:
            self.zoneType = "FETETRAHEDRON"
            
        elif "FEBRICK" in zoneType:
            self.zoneType = "FEBRICK"
            
        elif "FEPOLYGON" in zoneType:
            self.zoneType = "FEPOLYGON" 
            
        elif "FEPOLYHEDRAL" in zoneType:
            self.zoneType = "FEPOLYHEDRAL" 
            
        else:
            print("Invalid ZONETYPE - ", zoneType)
            raise TypeError
        
    def _checkDataPacking(self, dataPacking):
        
        if "BLOCK" in dataPacking:
            self.dataPacking = "BLOCK"
            
        elif "POINT" in dataPacking:
            self.dataPacking = "POINT"
            
        else:
            print("Invalid DATAPACKING - ", dataPacking)
            raise TypeError
    
    # Return true for correct zone on a given line, false otherwise
    def _checkZoneTitle(self, line, zoneCount):
        
        if "Zone_" in self.zoneName:
            
            if self.zoneName != "Zone_" + str(zoneCount):
                return False 
        else:
      
            # Check for zone title 
            if not regSearch(r'[,\s]T\s*=' , line, regMultiLine | regIgnoreCase):
                return False
            
            # Check for zone name 
            if not regSearch(self.zoneName, line, regMultiLine | regIgnoreCase):
                return False
            
        return True
        
    # Parse a zone header
    def readZoneHeader(self):
        
        reZone = regCompile(r'Zone',  regMultiLine | regIgnoreCase)
        
        reNode  = regCompile(r'Nodes\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        reNode1 = regCompile(r'[,\s]N\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        reNode2 = regCompile(r'[,\s]I\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        
        reElement  = regCompile(r'Elements\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        reElement1 = regCompile(r'[,\s]E\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        reElement2 = regCompile(r'[,\s]J\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        
        reK = regCompile(r'K\s*=\s*(\d*)', regMultiLine | regIgnoreCase)
        
        rePack = regCompile(r'DataPacking\s*=\s*(.*)', regMultiLine | regIgnoreCase)
        rePack1 = regCompile(r'[,\s]F\s*=\s*(.*)', regMultiLine | regIgnoreCase)
        
        reZoneType = regCompile(r'ZoneType\s*=\s*(.*)', regMultiLine | regIgnoreCase)
        
        zoneCount = 0
        self._lineOffset = 0
        with open(self.fileName, 'r') as fp:
            
            # Loop through the file and find the zone
            for line in fp.readlines():
                self._lineOffset += len(line)
                
                if not reZone.search(line):
#                 if not regSearch(r'Zone', line, regMultiLine | regIgnoreCase):    
                    continue
                
                zoneCount += 1
                if not self._checkZoneTitle(line, zoneCount-1):
                    continue
                
                # Nodes or i
#                 numNode = regSearch(r'Nodes\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                numNode = reNode.search(line)
                if not numNode:
#                     numNode = regSearch(r'[,\s]N\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                    numNode = reNode1.search(line)
                    if not numNode:
#                         numNode = regSearch(r'[,\s]I\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                        numNode = reNode2.search(line)
                        
                        if not numNode:
                            numNode = None
            
                if numNode:
                    numNode = int(numNode.group(1).replace(",",""))
            
                self.numNode = numNode
            
                # Elements or j
#                 numElement = regSearch(r'Elements\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                numElement = reElement.search(line)
                if not numElement:
#                     numElement = regSearch(r'[,\s]E\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                    numElement = reElement1.search(line)
                    
                    if not numElement:
#                         numElement = regSearch(r'[,\s]J\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                        numElement = reElement2.search(line)
                        if not numElement:
                            numElement = None
            
                if numElement:
                    numElement = int(numElement.group(1).replace(",",""))
                    
                self.numElement = numElement
                
                # k
#                 k = regSearch(r'K\s*=\s*(\d*)', line, regMultiLine | regIgnoreCase)
                k = reK.search(line)
                if k:
                    k = int(k.group(1).replace(",",""))
                else:
                    k = 1
                    
                self.k = k
                
#                 dataPacking = regSearch(r'DataPacking\s*=\s*(.*)', line, regMultiLine | regIgnoreCase)
                dataPacking = rePack.search(line)
                if not dataPacking:
#                     dataPacking = regSearch(r'[,\s]F\s*=\s*(.*)', line, regMultiLine | regIgnoreCase)
                    dataPacking = rePack.search(line)
            
                if dataPacking:
                    self._checkDataPacking(dataPacking.group(1).upper())
                
                else:
                    print("No DATAPACKING specified - defaulting to", self.dataPacking)
                
#                 zoneType = regSearch(r'ZoneType\s*=\s*(.*)', line, regMultiLine | regIgnoreCase)
                zoneType = reZoneType.search(line)
                if zoneType:
                    
                    self._checkZoneType(zoneType.group(1).upper())
                
                else:
                    if "FE" in dataPacking.group(1).upper(): # "FEPOINT and FEBLOCK is archaic
                        print("No ZONETYPE specified, but 'FE' found for DataPacking, defaulting to 'FEQUADRILATERAL'")
                        self._checkZoneType("FEQUADRILATERAL")
                    
                    else:
                        print("No ZONETYPE specified - defaulting to", self.zoneType)
                
                break # Quit after we have read the zone we were look for
            
        self._switchIJK()
    
    def _checkConnectivity(self, length, line):
        print("ZoneType set to", self.zoneType, "but connectivity length is", length)
        print("Line -")
        print(line)
                           
    # Given a open file parse it to populate the data and connectivity information
    # Assumed the file pointer is on the correct line!
    def _readZoneData(self, repeatConnect = False):
        
        if not self.data:
            print("Data dictionary hasn't been setup yet!")
            raise ValueError
        
        if self.zoneType.upper() == "ORDERED":
            if not self.i:
                print("'I' hasn't been setup yet for ordered data!")
                raise ValueError
            
            if not self.j:
                print("'J' hasn't been setup yet for ordered data!")
                raise ValueError
            
            if not self.k:
                print("'K' hasn't been setup yet for ordered data!")
                raise ValueError
        else:
            if not self.numElement:
                print("Number of elements hasn't been setup yet!")
                raise ValueError
        
            if not self.numNode:
                print("Number of nodes hasn't been setup yet!")
                raise ValueError
    
        print("Reading zone -", self.zoneName)
#         print "Line offset-  ", self._lineOffset
        
        # Open the file and read all the lines
        with open(self.fileName, 'r') as fp:
            fp.seek(self._lineOffset)
            allLine = fp.readlines()

# Obsolete with the introduction of _lineOffset
#         # Loop through the file and find the zone
#         zoneCount = 0
#         for i, line in enumerate(allLine):
#
#             if not regSearch(r'Zone', line, regMultiLine | regIgnoreCase):
#                 continue
#                 
#             zoneCount += 1
#             if not self._checkZoneTitle(line, zoneCount-1):
#                 continue
#         
#             break
#         
#         # Remove proceeding lines - don't care about anything prior to our zone 
#         for j in range(i,-1,-1):
#             allLine.pop(j)
                 
        if self.dataPacking.upper() == "POINT":
            
            if self.zoneType.upper() == "ORDERED":
                size = self.i*self.j*self.k
            else:
                size = self.numNode
            
            # Pre-allocate data array
            for j, key in enumerate(self.variable):
                self.data[key] = [None]*size
            
            # Read data 
            for i, line in enumerate(allLine):
                line = line.lstrip().rstrip().split()
                
                for j, key in enumerate(self.variable):
                    self.data[key][i] = float(line[j])
#                     self.data[key].append(float(line.lstrip().rstrip().split(" ")[j]))

                if self.zoneType.upper() == "ORDERED":
                    if i == self.i*self.j*self.k -1:
                        break
                else:
                    if i == self.numNode-1:
                        break

            if self.zoneType.upper() == "ORDERED": # Don't have connectivity data for "ordered" zones
                return  # We are done
            
#             print len(self.data["x"]), self.numNode, self.data["x"][0], self.data["x"][-1]
            
            # Remove all the lines we just read
            allLine = allLine[i+1:]
            
            # Pre-allocate connectivity
            self.connectivity = [None]*self.numElement
            
            # Read connectivity
            if (self.zoneType.upper() == "FELINESEG" or
                self.zoneType.upper() == "FETRIANGLE" or 
                self.zoneType.upper() == "FEQUADRILATERAL" or 
                self.zoneType.upper() == "FETETRAHEDRON" or 
                self.zoneType.upper() == "FEBRICK"):
                
                for i, line in enumerate(allLine):
                    connect = [None]*8
                
                    temp = line.lstrip().rstrip().replace("\n","").split()
                    length = 0
                    for j in temp:
                        
                        val = int(j)
                            
                        if not repeatConnect: # Should we allow for repeated connectivity indices?
                            
                            if val not in connect:
                                connect[length] = val
                                length += 1
                        else:
                            connect[length] = val
                            length += 1
                    
                    if self.zoneType.upper() == "FELINESEG" and length != 2:
                        self._checkConnectivity(length, line)
                        
                    elif self.zoneType.upper() == "FETRIANGLE" and length != 3:
                        self._checkConnectivity(length, line)
                        
                    elif self.zoneType.upper() == "FEQUADRILATERAL" and (length != 4 and length != 3):
                        self._checkConnectivity(length, line)
                        
                    elif self.zoneType.upper() == "FETETRAHEDRON" and length != 4:
                        self._checkConnectivity(length, line)
                        
                    elif self.zoneType.upper() == "FEBRICK" and length >= 8:
                        self._checkConnectivity(length, line)
                         
                    self.connectivity[i] = connect[0:length]
                
                    if i == self.numElement-1:
                        break
            else:
                print("Tecplot zone reader does not fully support", self.zoneType, "zones yet!")
                raise TypeError
                
        elif self.dataPacking.upper() == "BLOCK":
            print("Tecplot zone reader does not fully support ",  self.dataPacking, "packing yet!")
            raise TypeError
        
        else:
            print("Unsupport packing type ",  self.dataPacking, "packing yet!")
            raise TypeError
                
        return 
    
cdef class capsTecplot:
    
    cdef: 
        readonly object title 
        readonly object variable 
        
        readonly object zone 
        
        public object fileName 
    
    # \param fileName Name of the file that 
    def __cinit__(self, fileName=None):
      
        self.title = None
        self.variable = []
         
        # Dictionary con
        self.zone = {} 
         
        self.fileName = fileName

    # Read the header information for the files 
    def loadHeader(self):
        
        if not self.fileName: 
            print("No fileName as been specified")
            raise IOError
        
        with open(self.fileName, 'r') as fp:
            allLine = fp.readlines() # Read all the lines and up them in one big string
            
            # Find title
            for line in allLine:
                title= regSearch(r'Title\s*=\s*(.*)', line, regMultiLine | regIgnoreCase)
                if not title:
                    continue
                else:
                    self.title = title.group(1).replace('"', "").replace("'", "")
                    break
            
            if not title:
                print("No title found - defaulting to", self.title)
                
            # Find variables 
            for line in allLine:
                variables = regSearch(r'Variables\s*=\s*(.*)', line, regMultiLine | regIgnoreCase)
                if not variables:
                    continue
                else:
                    break
            
            if not variables:
                print("No variables found!")
                raise ValueError
            else:
                variables = variables.group(1)
                
            if "," in variables:
                variables = variables.split(",")
            elif " " in variables:
                variables = variables.split(" ")
            else:
                print("Unable to parse variable list - ", variables)
                raise TypeError
            
            for i in variables:
                if i: # Deal with empty characters in list
                    self.variable.append(i.replace('"', "").replace("'", "").replace(" ", "").replace("\n", ""))
            
            # Figure out how many zones there are 
            for line in allLine:
                
                # Look for Zone line 
                zone = regSearch(r"Zone", line, regMultiLine | regIgnoreCase) 
                if zone:
                    
                    # Get title
                    zone = regSearch(r"\s*T\s*=\s*(\".*\")", line, regMultiLine | regIgnoreCase) # First look for "
                    if not zone:
                        
                        zone = regSearch(r"\s*T\s*=\s*(\'.*\')", line, regMultiLine | regIgnoreCase) # Next try '
                        if not zone:
                            
                            title = "Zone_" + str(len(self.zone.keys()))
                            
                            print("No zone title found!, Defaulting to", title)
                        
                        else:
                            title = zone.group(1).replace('"', "").replace("'", "")
                    else:
                        title = zone.group(1).replace('"', "").replace("'", "")
                  
                else:
                    continue
            
                # Setup zone dictionary 
                self.zone[title] = capsTecplotZone(title, self.variable, self.fileName) # Set the zone name and variables in the zone
            
            if not self.zone:
                print("No zones found!")
                raise ValueError
            
    # Load all zones in the file  
    def loadZone(self):
        
        if not self.zone:
            self.loadHeader()
        
        for i in self.zone.keys():
            self.zone[i].readZoneHeader()
            self.zone[i]._readZoneData()


cdef class capsUnstructMesh:
    
    cdef:
        
        object filename
                
        public name 
        
        readonly numElement
        readonly numNode
       
        readonly variable
        readonly xyz 
        readonly object data
        readonly object connectivity
        readonly object oneBias
        
        readonly object meshType
    
    def __cinit__(self, filename):
        
        self.name = "UnstructuredMesh"
        self.numElement = 0
        self.numNode = 0
        
        self.xyz = []
        self.connectivity = []

        self.variable = []
        self.data = {}
        
        self.filename = filename
        
        self.meshType = "unknown"
        
        self.oneBias = True
        
    def loadFile(self, **kwargs):
        
        if ".origami" in self.filename:
            self._loadOrigami()
        elif ".ugrid" in self.filename:
            self._loadAFLR3()
        elif ".stl" in self.filename:
            self._loadSTL(**kwargs)
        elif ".igs" in self.filename or \
             ".iges" in self.filename or \
             ".stp" in self.filename or \
             ".step" in self.filename or \
             ".brep" in self.filename or \
             ".egads" in self.filename or \
             ".csm"in self.filename:
            self._loadGeometry(self.filename, **kwargs)
        else:
            print("Unable to determine mesh format!" )
            raise ValueError
    
    def _printLoading(self):
        print("Loading ( type =", self.meshType, ") file", self.filename)
        
    def _parseFileName(self):
        self.name = os.path.basename(self.filename)
        
    def _loadAFLR3(self):
        
        sint = 4 # bytes
        sdouble = 8 # bytes 
        
        self._parseFileName()
        
        self.meshType = "aflr3"
                
        self.oneBias = True
        
        self.variable = ["BCID"]
        self.data["BCID"] = []
        
        numTriangle = 0
        numQuadrilateral = 0
        numTetrahedral = 0
        numPyramid = 0
        numPrism = 0
        numHexahedral = 0
        
        self._printLoading()
        with open(self.filename, 'r') as fp:
            
            # Little endian  or big endian 
            if "lb8" in self.name or "b8" in self.name:
            
                self.numNode     = structUnpack("i", fp.read(sint))[0]
                
                numTriangle      = structUnpack("i", fp.read(sint))[0]
                numQuadrilateral = structUnpack("i", fp.read(sint))[0]
                numTetrahedral   = structUnpack("i", fp.read(sint))[0]
                numPyramid       = structUnpack("i", fp.read(sint))[0]
                numPrism         = structUnpack("i", fp.read(sint))[0]
                numHexahedral    = structUnpack("i", fp.read(sint))[0]
                
                self.numElement = numTriangle + \
                                  numQuadrilateral + \
                                  numTetrahedral + \
                                  numPyramid + \
                                  numPrism + \
                                  numHexahedral
                
                print "AFLR3 header = ", self.numNode, numTriangle, numQuadrilateral, numTetrahedral, numPyramid, numPrism, numHexahedral
                
                
                # Nodes
                i = 0 
                while i < self.numNode:
                    self.xyz.append( [ structUnpack("d", fp.read(sdouble))[0], 
                                       structUnpack("d", fp.read(sdouble))[0], 
                                       structUnpack("d", fp.read(sdouble))[0] ])
                    i += 1
                    
                # Triangles
                i = 0 
                while i < numTriangle:
                    
                    self.connectivity.append([ structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0]])
                    i += 1
                    
                # Quadrilaterals 
                i = 0
                while i < numQuadrilateral:
                    
                    self.connectivity.append([ structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0]])
                    i += 1
                    
                # Boundary ids
                self.data["BCID"] = [0]*(numTriangle + numQuadrilateral) #[0 for x in range(numTriangle + numQuadrilateral)]
                
                
                #for i in range(numTriangle + numQuadrilateral):
                #    self.data["BCID"][i] = structUnpack("i", fp.read(sint))[0]
                
                self.data["BCID"] = list(structUnpack("i"*(numTriangle + numQuadrilateral), fp.read(sint*(numTriangle + numQuadrilateral))))
                
                # Tetrahedrals
                i = 0
                while i < numTetrahedral:
                    
                    self.connectivity.append([ structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0]])
                    i += 1
                    
                # Pyramids
                i = 0
                while i < numPyramid:
                    
                    self.connectivity.append([ structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0]])
                    i += 1
                    
                # Prisms
                i = 0
                while i < numPrism:
                    
                    self.connectivity.append([ structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0],
                                               structUnpack("i", fp.read(sint))[0]])
                    i +=1 
                    
                # Hexahedrals
                i = 0
                while i < numHexahedral:
                    
                    self.connectivity.append([ structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0], 
                                               structUnpack("i", fp.read(sint))[0],
                                               structUnpack("i", fp.read(sint))[0],
                                               structUnpack("i", fp.read(sint))[0],
                                               structUnpack("i", fp.read(sint))[0]])
                    i += 1
                    
            # ASCII 
            elif ".ugrid" in self.filename:
                
                line = fp.readline()
                if "," in line:
                    temp = line.split(",")
                else:
                    temp = line.split(" ")
     
                self.numNode     = int(temp[0]) 
                
                numTriangle      = int(temp[1]) 
                numQuadrilateral = int(temp[2]) 
                numTetrahedral   = int(temp[3]) 
                numPyramid       = int(temp[4]) 
                numPrism         = int(temp[5]) 
                numHexahedral    = int(temp[6]) 
                
                self.numElement = numTriangle + \
                                  numQuadrilateral + \
                                  numTetrahedral + \
                                  numPyramid + \
                                  numPrism + \
                                  numHexahedral
                
                # Nodes 
                for i in range(self.numNode):
                    temp = fp.readline().split(" ")
                    
                    self.xyz.append( [ float(temp[0]), 
                                       float(temp[1]), 
                                       float(temp[2]) ])
                    
                # Triangles 
                for i in range(numTriangle):
                    temp = fp.readline().split(" ")
                    
                    self.connectivity.append([ int(temp[0]), 
                                               int(temp[1]), 
                                               int(temp[2])])
                    
                # Quadrilaterals 
                for i in range(numQuadrilateral):
                    temp = fp.readline().split(" ")
                    
                    self.connectivity.append([ int(temp[0]), 
                                               int(temp[1]), 
                                               int(temp[2]),
                                               int(temp[3])])
                
                # Boundary ids
                self.data["BCID"] = [0 for x in range(numTriangle + numQuadrilateral)]
                
                for i in range(numTriangle + numQuadrilateral):
                    temp = fp.readline().split(" ")
                    
                    self.data["BCID"][i] = int(temp[0])
                    
                # Tetrahedrals
                for i in range(numTetrahedral):
                    temp = fp.readline().split(" ")
                     
                    self.connectivity.append([ int(temp[0]), 
                                               int(temp[1]), 
                                               int(temp[2]),
                                               int(temp[3])])
                    
                # Pyramids
                for i in range(numPyramid):
                    temp = fp.readline().split(" ")
                    
                    self.connectivity.append([ int(temp[0]), 
                                               int(temp[1]), 
                                               int(temp[2]),
                                               int(temp[3]),
                                               int(temp[4])])
                    
                # Prisms
                for i in range(numPrism):
                    temp = fp.readline().split(" ")
                    
                    self.connectivity.append([ int(temp[0]), 
                                               int(temp[1]), 
                                               int(temp[2]),
                                               int(temp[3]),
                                               int(temp[4]),
                                               int(temp[5])])
                # Hexahedrals
                for i in range(numHexahedral):
                    temp = fp.readline().split(" ")
                    
                    self.connectivity.append([ int(temp[0]), 
                                               int(temp[1]), 
                                               int(temp[2]),
                                               int(temp[3]),
                                               int(temp[4]),
                                               int(temp[5]),
                                               int(temp[6]),
                                               int(temp[7])])
            else:
                raise TypeError
        
    def _loadOrigami(self):
        
        self._parseFileName()
        
        self.meshType = "origami"
        
        self.oneBias = True
        
        nodeCount = 0
        elementCount = 0
        
        self._printLoading()
        with open(self.filename, 'r') as fp:
            
            for line in fp:
                
                if line.startswith("#"):
                    continue
                    
                if self.numNode == 0 and self.numElement == 0:
                    
                    self.numNode = int(line.split(" ")[0])
                    self.numElement = int(line.split(" ")[1])
                    
                    self.xyz = [[] for x in range(self.numNode)]
                    self.connectivity = [[] for x in range(self.numElement)]
                    
                    self.variable = ["ClosetNeighbor_1", "ClosetNeighbor_2", "FoldLine"]
                    
                    self.data["ClosetNeighbor_1"] = [0 for x in range(self.numElement)]
                    self.data["ClosetNeighbor_2"] = [0 for x in range(self.numElement)]
                    self.data["FoldLine"] = [0 for x in range(self.numElement)]
                    
                    continue
                    
                if nodeCount < self.numNode:
                    xyz = line.split(" ")
                    
                    index = int(xyz[0])-1
                    
                    self.xyz[index].append(float(xyz[1]))
                    self.xyz[index].append(float(xyz[2]))
                    self.xyz[index].append(float(xyz[3]))
                    
                    nodeCount += 1
                    
                    continue
                
                if elementCount < self.numElement:
                    
                    element = line.split(" ")
                    
                    index = int(element[0])-1
                    
                    self.connectivity[index].append(int(element[1]))
                    self.connectivity[index].append(int(element[2]))
                    
                    self.data["ClosetNeighbor_1"][index] = int(element[3])
                    self.data["ClosetNeighbor_2"][index] = int(element[4])
                    self.data["FoldLine"][index]         = int(element[5])
                    
                    elementCount += 1
                    
                    continue
                
    def _loadSTL(self, stitch=True):
        
        def sortNode(node, connect):
            
            if stitch == True:
                # Does the vertex already exist?           
                try:
                    index = self.xyz.index(node)
                    connect.append(index+1) # Yes, get the index for the connectivity
                                
                except ValueError: # No! Append it 
                    self.xyz.append(node)
                    
                    self.numNode += 1
                    connect.append(self.numNode)
                except:
                    raise
            else: 
#                 print "Not stitching stl file!"
                self.xyz.append(node)
                    
                self.numNode += 1
                connect.append(self.numNode)
            
        sint  = 4 # bytes
        schar = 1 # bytes
        sfloat = 4 # bytes 
        sshort = 2 # bytes
                    
        self._parseFileName()
        
        self.meshType = "stl"
        
        self.oneBias = True
        
        self._printLoading()
        
        ascii = True
        
        with open(self.filename, 'r') as fp:
            
            # ASCII or binary
            line = fp.read(schar*80)
            
            if "solid" in line: # ASCII
                
                print("'solid' tag found! Assuming STL is in ASCII format....")
                
                fp.seek(0) # Return to the beginning and get the whole line
                
                line = fp.readline()
                
                line = line.split()
                if len(line) > 1:
                    self.name = str(line[1])
                else:
                    self.name = self.filename
         
            else: # Binary
                    
                print("No 'solid' tag found! Assuming STL is in binary format....")
            
                ascii = False
                    
                fp.seek(0)
                
                self.name = str(fp.read(schar*80)) #structUnpack("s", fp.read(schar*80))[0]
                    
                self.numElement = structUnpack("I", fp.read(sint))[0]
                print "STL name - ", self.name
                print "Number of elements - ", self.numElement
            
            if ascii: # STL appears to be ASCII       

                for line in fp:
                    if "outer loop" in line:
                        
                        connect = []
                        for line in fp:
                            
                            if "vertex" in line:
                                data = line.split()
                                # Read vertex
                                node = [float(data[1]), 
                                        float(data[2]),
                                        float(data[3]) ]
                            
                                # Determine if node is new or not 
                                sortNode(node, connect)
                                
                            if "endloop" in line:
                                self.connectivity.append(connect)
                                self.numElement += 1
                                
                                break 
                    
                    if "endsolid" in line:
                        break
                    
            else: # Binary 
                
                # Pre-allocate connectivity array
                self.connectivity = [None]*self.numElement
                
                if stitch == False:
                    self.xyz = [None]*self.numElement*3
                    
                i = 0
                while i < self.numElement:
#                     start = time.time()
                    # Read the normals
                    fp.read(sfloat*3)
#                     fp.read(sfloat)
#                     fp.read(sfloat)
#                     fp.read(sfloat)
#                     
                    connect = []
                    j = 0
                    while j < 3:
                        
                        # Read vertex
                        node = structUnpack("f"*3, fp.read(sfloat*3))
                        
                        if stitch == True:
                            # Determine if node is new or not 
                            sortNode([node[0], node[1], node[2]], connect)
                        else:
                            self.xyz[self.numNode] = [node[0], node[1], node[2]]
                            self.numNode +=1
                            connect.append(self.numNode)
                         
                        j += 1
                    
                    # Add element to connectivity 
                    #self.connectivity.append(connect)
                    self.connectivity[i] = connect[:]
                         
                    # Read attribute bytes - currently not doing anything with it - could add it to the data 
                    attrByte = structUnpack("H", fp.read(sshort))[0]
                    
                    i += 1
#                     print "Time - ", time.time() - start
                     
       # print self.numNode, self.numElement 
        
        if self.numNode == 0 or self.numElement == 0:
            print "No nodes or elements created.... something went wrong"
            raise ValueError
    
    
    def _loadGeometry(self, geomFile, flags=0, dictFile = None):
        
        cdef: 
            cEGADS.ego context
            
            cEGADS.ego model
            cEGADS.ego ref
            int oclass, mtype
            double[4] data
            int *senses
            
            int numBody
            cEGADS.ego *body = NULL
            object freeBody = False
            
            # OpenCSM
            void *modl = NULL
            cOCSM.modl_T *MODL=NULL
            
            int builtTo
     
        def cleanup():
            
            cEGADS.EG_close(context)
     
            MODL = NULL
                
            if modl: 
                cOCSM.ocsmFree(modl)
                
            if freeBody:
                if body:
                    cEGADS.EG_free(body)
            
        def raiseError(msg):
            cleanup() 
               
            raise ValueError(msg=msg)
                        
        if not isinstance(geomFile, str):
            raise TypeError
        
        fileName = _byteify(geomFile)
        
        if ".csm" in geomFile: # If a ESP file - use OpenCSM loader
        
            # Set OCSMs output level
            cOCSM.ocsmSetOutLevel(<int> 1)
            status = cOCSM.ocsmLoad(fileName, &modl)
            if status < 0: 
                raiseError("while getting bodies from model (during a call to ocsmLoad)")
            
            # Cast void to model 
            MODL = <cOCSM.modl_T *> modl
            
            # If there was a dictionary, load the definitions from it now */
            if dictFile:
                dictName = _byteify(dictFile)
                status = cOCSM.ocsmLoadDict(modl, dictName)
                          
                if status < 0: 
                    raiseError("while loading dictionary (during a call to ocsmLoadDict)")
                
        
            status = cOCSM.ocsmBuild( modl,     # (in)  pointer to MODL
                                           <int> 0,  # (in)  last Branch to execute (or 0 for all, or -1 for no recycling)
                                           &builtTo, # (out) last Branch executed successfully
                                           &numBody, # (in)  number of entries allocated in body[], (out) number of Bodys on the stack
                                           NULL)     # (out) array  of Bodys on the stack (LIFO) (at least *nbody long)
            if status < 0: 
                raiseError("while building model (during a call to ocsmBuild)")
            
            numBody = 0
            for i  in range(MODL.nbody):
           
                if MODL.body[i+1].onstack != 1:
                    continue

                numBody += 1
                    
                body = <cEGADS.ego *> cEGADS.EG_reall(body, numBody*sizeof(cEGADS.ego))
                     
                if body == NULL:
                    status = cEGADS.EGADS_MALLOC
                    raiseError("while allocating bodies  (during a call to EG_reall)")
                     
                freeBody = True
                
                status = cEGADS.EG_copyObject(MODL.body[i+1].ebody, NULL, &body[numBody -1])
                if status != cEGADS.EGADS_SUCCESS:
                    raiseError("while getting bodies from model (during a call to EG_copyObject)") 
              
        else: # else use a EGADS model loader
        
            status = cEGADS.EG_open(&context)
            if status != cEGADS.EGADS_SUCCESS:
                raise CAPSError(status, msg="while opening context  (during a call to EG_open)")
            
            status = cEGADS.EG_loadModel(context, <int> flags, <char *>fileName, &model)
            if status != cEGADS.EGADS_SUCCESS: 
                raiseError("while loading model from geometry file (during a call to EG_loadModel)")
            
            status = cEGADS.EG_getTopology(model, &ref, &oclass, &mtype, data, &numBody, &body, &senses)
            if status != cEGADS.EGADS_SUCCESS or oclass != cEGADS.MODEL:
                raiseError("while getting bodies from model (during a call to EG_getTopology)")
        
        try:
            self.__loadEGADS(numBody, body)
        except CAPSError as e:
            status = e.errorCode
            raiseError(e.args)
            
        cleanup()
    
    cdef __loadEGADS(self, int numBody, cEGADS.ego *body): 
        cdef: 
            cEGADS.ego ref
            int oclass, mtype
            double[4] data
            int *senses
            
            int numFace, numChildren
            cEGADS.ego tempBody
            cEGADS.ego *face = NULL
            cEGADS.ego *children
            
            cEGADS.ego tess
            
            double box[6]
            double size = 0
            double params[3]
            
            #EGADS tess return variables
            int numPoint = 0, numTri = 0
            const int *tris = NULL
            const int *triNeighbor = NULL
            const int *ptype = NULL
            const int *pindex = NULL
            
            const double *points = NULL
            const double *uv = NULL

            int     *triConn = NULL
            double[3] xyzs 
            
            int dummyInt
            int globalID = 0
            int pointType
            int pointIndex
            
            object xyz
            object connectivity
            
            # Attributes
            int atype
            int alen
            const int *ints
            const double *reals 
            const char *string
            
        def cleanup():
            
            if face:
                cEGADS.EG_free(face)
            
        def raiseError(msg):
            cleanup() 
            raise CAPSError(status, msg=msg)
        
        # Loop through bodies - tesselate and create gPrims 
        for i in range(numBody):
            
            #print "Tessellating Body", i+1
            
            xyz = []
            connectivity = []
            
            status = cEGADS.EG_getTopology(body[i], &ref, &oclass, &mtype, data, &numChildren, &children, &senses)
            if status != cEGADS.EGADS_SUCCESS:
                raiseError("while getting bodies from model (during a call to EG_getTopology)")
            
            # TODO: Add support for line bodies 
            if mtype != cEGADS.FACEBODY and mtype != cEGADS.SHEETBODY and mtype != cEGADS.SOLIDBODY:
                print "Unsupported body type... skipping tessellation" 
                continue
            
            # Get tessellation             
            status = cEGADS.EG_getBoundingBox(body[i], box)
            if status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting bounding box of body (during a call to (during a call to EG_makeTessBody)")
            
            size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
                        (box[1]-box[4])*(box[1]-box[4]) +
                        (box[2]-box[5])*(box[2]-box[5]))
            
            params[0] = 0.00250*size
            params[1] = 0.00010*size
            params[2] = 15.0

            status = cEGADS.EG_makeTessBody(body[i], params, &tess)
            if status != cEGADS.EGADS_SUCCESS: 
                 raiseError("while making tesselation (during a call to EG_makeTessBody)")
        
            # How many faces do we have? 
            status = cEGADS.EG_getBodyTopos(body[i], 
                                                 <cEGADS.ego> NULL, 
                                                 cEGADS.FACE, 
                                                 &numFace, 
                                                 &face)
            if status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting number of faces (during a call to EG_getBodyTopos)")
            
#             if face:
#                 cEGADS.EG_free(face)
#                 face = NULL
                
            status = cEGADS.EG_statusTessBody(tess, &tempBody, &dummyInt, &numPoint)
            if status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting the status of the tessellation (during a call to EG_statusTessBody)")
            
            # Get the points 
            for j in range(numPoint):
      
                status = cEGADS.EG_getGlobal(tess, j+1, &pointType, &pointIndex, xyzs)
                if status != cEGADS.EGADS_SUCCESS: 
                    raiseError("while getting coordinates for the tessellation (during a call to EG_getGlobal)")
                
                xyz.append([xyzs[0], xyzs[1], xyzs[2]])
            
            color = '#ffffa3' # yellow - Default color for the body/face
            
            # Loop through face by face and get tessellation information 
            for faceIndex in range(numFace):
                
                status = cEGADS.EG_getTessFace(tess, faceIndex+1, &numPoint, &points, &uv, &ptype, 
                                                    &pindex, &numTri, &tris, &triNeighbor)
                if status != cEGADS.EGADS_SUCCESS: 
                    raiseError("while getting the face tessellation (during a call to EG_getTessFace)")

                for j in range(numTri):
                    tri = []
                    status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[3*j + 0], &globalID)
                    if status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting local-to-global indexing [connectivity - 0] (during a call to EG_localToGlobal)")
                
                    tri.append(globalID)
                
                    status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[3*j + 1], &globalID)
                    if status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting local-to-global indexing [connectivity - 1] (during a call to EG_localToGlobal)")
                
                    tri.append(globalID)
                
                    status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[3*j + 2], &globalID)
                    if status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting local-to-global indexing [connectivity - 2] (during a call to EG_localToGlobal)")                 
                    
                    tri.append(globalID)
                
                    connectivity.append(tri)
                
                try:
                    temp = __getColor(face[faceIndex])
                except CAPSError as e: 
                    status = e.errorCode
                    raiseError("while getting color")
                except Exception as e:
                    raise(e)
             
                if temp:
                    color = temp
            
            # Get color
            try: 
                 color = __getColor(body[i]) or color or '#ffffa3' # yellow - Default color for the body/face
            except CAPSError as e: 
                    status = e.errorCode
                    raiseError("while getting color")
            except Exception as e:
                    raise(e)
                    
            #print xyz
            #print connectivity
            self.xyz = xyz
            self.numNode = len(self.xyz)
            self.connectivity = connectivity
            self.numElement = len(self.connectivity)
            
            #self.addColor(color)
            #self.addVertex(xyz)
            #self.addIndex(connectivity)
            #self.addLineIndex(connectivity)
            #self.addLineColor(1, colorMap = "grey")
            
            attr = _byteify("_name")
            # Get body name - if available 
            status = cEGADS.EG_attributeRet(body[i], <char *>attr, &atype, &numPoint, &ints, &reals, &string)
            if status == cEGADS.EGADS_NOTFOUND:
                name = None
            elif status != cEGADS.EGADS_SUCCESS:
                raiseError("while getting trying to retrieve body name (during a call to EG_attributeRet)")    
            else: 
                if atype == cEGADS.ATTRSTRING:
                    name = _strify(string)
                else:
                    name = None
            #self.createTriangle(name=name)
                
        cleanup()
        return 
        
cdef class capsUnstructData:
    
    cdef:
        public name 
        
        readonly numElement
        readonly numNode
        readonly object variable # List of variables found in _data to keep track of variable order
        
        object _xyz 
        object _data
        object _connectivity
        
        #readonly blankParam
        object _blank #
        
        object _metaData
        
        object fileExt
    
    def __cinit__(self, file=None):
        
        self.name = "capsUnstructuredData"
        self.numElement = 0
        self.numNode = 0
        
        self.xyz = []
        self.connectivity = []
        
        self._blank = []
        
        self.metaData = None
        
        self.variable = []
        self.data = {}
        
        self.fileExt = ".cud"
        
        if file:
            self.loadFile(file)
    
    @property
    def xyz(self):
        return self._xyz
    
    @xyz.setter
    def xyz(self, data):
        
        if not data:
            self._xyz = []
            self.numNode = 0
            return 
        
        if not isinstance(data, list):
            data = [data]
        
        if not isinstance(data[0], list):
            data = [data]
            
        self._xyz = data[:]
        self.numNode = len(data)
        
#     @xyz.append
#     def xyz(self, data):
#         
#         if not data:
#             return
#         
#        
#         self._xyz += data[:]
#         self.numNode += len(data)
        
    
    @property
    def connectivity(self):
        return self._connectivity
    
    @connectivity.setter
    def connectivity(self, data):
        
        if not data:
            self._connectivity = []
            self.numElement = 0
            return 
        
        if not isinstance(data, list):
            data = [data]
        
        if not isinstance(data[0], list):
            data = [data]
            
        self._connectivity = data[:]
        self.numElement = len(data)

    @property
    def data(self): 
        return self._data
    
    @data.setter
    def data(self, dataDict):
        
        if not isinstance(dataDict, dict):
            raise TypeError("Data should be a dictionary!")
        
        self._data = dataDict.copy()
        
        for i in self._data.keys():
            if len(self._data[i]) != self.numNode:
                raise IndexError("Data array, " + str(i) +", not equal to the number of nodes!")
            
        self.variable = self._data.keys()
        
    @property
    def metaData(self):
        return self._metaData
    
    @metaData.setter
    def metaData(self, data):
        self._metaData = deepcopy(data)
        
#     @property
#     def blank(self):
#         return self._blank
#     
#     @blank.setter
#     def blank(self, value, data):
#         if len(data)
    
    def boundingBox(self):
        box = [1.0E100, 1.0E100, 1.0E100, -1.0E100, -1.0E100, -1.0E100]
        i = 0
        while i < self.numNode:
    
            if (box[0] > self._xyz[i][0]): box[0] = self._xyz[i][0]
            if (box[3] < self._xyz[i][0]): box[3] = self._xyz[i][0]
            
            if (box[1] > self._xyz[i][1]): box[1] = self._xyz[i][1]
            if (box[4] < self._xyz[i][1]): box[4] = self._xyz[i][1]
            
            if (box[2] > self._xyz[i][2]): box[2] = self._xyz[i][2]
            if (box[5] < self._xyz[i][2]): box[5] = self._xyz[i][2]
            
            i += 1
        
        return box
    
    def blankData(self, param, val, operator="<="):
        
        import operator as op

        opMap = {'<' : op.lt,
                 '<=' : op.le,
                 '==' : op.eq,
                 '!=' : op.ne,
                 '>=' : op.ge,
                 '>' : op.gt}
        
        param = param.lower()
        
        if param != "x" and param != "y" and param != "z" :
            raise ValueError("Unsupported param - " + str(param))
        
        if operator not in opMap.keys():
            raise ValueError("Unsupported operator - " + str(operator))
        
        if not self._blank: 
            self._blank = [0]*self.numNode
       
        opFn = opMap[operator] 
        
        numBlank = 0
        i = 0
        while i < self.numNode:
            
            if param == "x":
                if opFn(self._xyz[i][0], val): 
                    self._blank[i] = 1
                    numBlank += 1
    
            if param == "y":
                if opFn(self._xyz[i][1], val):
                    self._blank[i] = 1
                    numBlank += 1
    
            if param == "z":
                if opFn(self._xyz[i][2], val):
                    self._blank[i] = 1
                    numBlank += 1
                     
            i += 1
        print "Number of blanked nodes = ", numBlank
        return numBlank
   
    def writeBlankFile(self, filename):
 
        if not self._blank:
            print "No blanking done - just writing file"
            self.writeFile(filename)
            return 
               
        mesh = capsUnstructData()
        
        mesh.numNode = 0
        mesh.numElement = 0
        
        nodeMap = [0]*self.numNode
        
        i = 0
        while i < self.numNode:
            
            if self._blank[i] == 0: 
                mesh.numNode += 1
                nodeMap[i] = mesh.numNode # Old node -> new node
                
            i += 1
                
        # Pre-allocate arrays
        mesh._xyz = [[0.0,0.0,0.0]]*mesh.numNode
        mesh._connectivity = [None]*self.numElement

        i = 0
        while i < self.numNode:
            if nodeMap[i] != 0:
                mesh._xyz[nodeMap[i]-1] = self._xyz[i][:]
            i += 1
            
        connect = [0]*30
        
        i = 0
        while i < self.numElement:
            blanked = False
            
            for k, j in enumerate(self.connectivity[i]): 
            
                if self._blank[j-1] == 1: #Blanked
                    blanked = True
                    break
            
                connect[k] = nodeMap[j-1]
            
            if blanked == False:
                mesh._connectivity[mesh.numElement] = connect[:len(self.connectivity[i])] 
                mesh.numElement += 1
                
            i += 1
        
        mesh._connectivity = mesh._connectivity[:mesh.numElement]        
                
        print "NOT CURRENTLY COPY DATA AS WELL!"
        
        print "Blank file:"
        print "Number of nodes,", mesh.numNode
        print "Number of elements",mesh.numElement
        
        mesh.writeFile(filename)
    
    def scale(self, scaleFactor):
        i = 0
        while i < self.numNode:
            self._xyz[i][0] *= scaleFactor
            self._xyz[i][1] *= scaleFactor
            self._xyz[i][2] *= scaleFactor
            i += 1
    
    def translate(self, translationVector):
        
        if len(translationVector) != 3 or not isinstance(translationVector, list):
            raise ValueError("Variable, transVector, should be a list of len 3.")
        
        i = 0
        while i < self.numNode:
            self._xyz[i][0] += translationVector[0]
            self._xyz[i][1] += translationVector[1]
            self._xyz[i][2] += translationVector[2]
            i += 1
            
    def flipCoordinate(self, coordinate="x"):

        if coordinate.lower() == "x":
            self._flipXCoord()
        elif coordinate.lower() == "y":
            self._flipYCoord()
        elif coordinate.lower() == "z":
            self._flipZCoord()
        else: 
            raise ValueError("Unreconized value for 'coordinate'. Valid options = 'x', 'y', 'z'")
        
    def _flipXCoord(self):
        i = 0
        while i < self.numNode:
            self._xyz[i][0] *= -1
            i += 1
            
    def _flipYCoord(self):
        i = 0
        while i < self.numNode:
            self._xyz[i][1] *= -1
            i += 1

    def _flipZCoord(self):
        i = 0
        while i < self.numNode:
            self._xyz[i][2] *= -1
            i += 1
    
    def flipNormal(self):
        i = 0
        while i < self.numElement:
            self._connectivity[i].reverse()
            i += 1
        
    def mirror(self, plane):

        if plane.lower() == "xy":
            self._mirrorXYPlane()
        elif plane.lower() == "xz":
            self._mirrorXZPlane()
        elif plane.lower() == "yz":
            self._mirrorYZPlane()
        else: 
            raise ValueError("Unreconized value for 'plane'. Valid options = 'xy', 'xz', 'yz'")
    
    def __mirrorData(self):
        for i in self._data.keys():
            self._data[i] += self._data[i][:]
            
    def __mirrorConnect(self):

        connect = deepcopy(self._connectivity)
        
        i = 0
        while i < self.numElement:
            j = 0
            while j < len(connect[i]):
                connect[i][j] += self.numNode
                j += 1
            connect[i].reverse()
            i += 1
        self._connectivity += connect
        self.numElement *= 2
         
    # No stiching at mirror plane
    def _mirrorXYPlane(self):
        
        mirror = deepcopy(self._xyz)
        
        i = 0
        while i < self.numNode:
            mirror[i][2] *= -1
            i += 1
        self._xyz += mirror
        
        self.__mirrorConnect()
            
        self.numNode *= 2
        
        self.__mirrorData()
        
    # No stiching at mirror plane   
    def _mirrorXZPlane(self):
        
        mirror = deepcopy(self._xyz)
        
        i = 0
        while i < self.numNode:
            mirror[i][1] *= -1
            i += 1
        self._xyz += mirror
        
        self.__mirrorConnect()
            
        self.numNode *= 2
        
        self.__mirrorData()
        
    # No stiching at mirror plane
    def _mirrorYZPlane(self):
        
        mirror = deepcopy(self._xyz)
        
        i = 0
        while i < self.numNode:
            mirror[i][0] *= -1
            i += 1
        self._xyz += mirror
        
        self.__mirrorConnect()
            
        self.numNode *= 2
        
        self.__mirrorData()

    def rotate(self, theta, axis="x"):

        if axis.lower() == "x":
            self._rotateXAxis(theta)
        elif axis.lower() == "y":
            self._rotateYAxis(theta)
        elif axis.lower() == "z":
            self._rotateZAxis(theta)
        else: 
            raise ValueError("Unreconized value for 'axis'. Valid options = 'x', 'y', 'z'")
        
    # theta in radions
    def _rotateXAxis(self,theta):
        i = 0
        while i < self.numNode:
            x = self._xyz[i][:]
            #self._xyz[i][0] = self._xyz[i][0]
            self._xyz[i][1] = cos(theta)*x[1] - sin(theta)*x[2]
            self._xyz[i][2] = sin(theta)*x[1] + cos(theta)*x[2]     
            i += 1

    # theta in radions
    def _rotateYAxis(self,theta):
        i = 0
        while i < self.numNode:
            x = self._xyz[i][:]
            self._xyz[i][0] =  cos(theta)*x[0] + sin(theta)*x[2]
            #self._xyz[i][1] = self._xyz[i][1]
            self._xyz[i][2] = -sin(theta)*x[0] + cos(theta)*x[2]
            i += 1
    
    # theta in radions
    def _rotateZAxis(self,theta):
        i = 0
        while i < self.numNode:
            x = self._xyz[i][:]
            self._xyz[i][0] = cos(theta)*x[0] - sin(theta)*x[1]
            self._xyz[i][1] = sin(theta)*x[0] + cos(theta)*x[1]
            #self._xyz[i][2] = self._xyz[i][2]
            i += 1
   
    # Determine if a ray intersects a given triangle 
    # 
    # \return Returns a boolean indicating intersection and the intersection 
    # point (intersectBool, [X_intersect, Y_intersect, Z_intersect])
    def _rayIntersectTriangle(self, origin, ray, triangle, epsilon=0.0000001):
     
        def diff(a, b):
            return [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
     
        def cross(a, b):
            return [a[1]*b[2] - a[2]*b[1],
                    a[2]*b[0] - a[0]*b[2],
                    a[0]*b[1] - a[1]*b[0]]

        def dot(a, b):
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
     
        vertex0 = triangle[0]
        vertex1 = triangle[1]  
        vertex2 = triangle[2]
        
        edge1 = diff(vertex0, vertex1)
        edge2 = diff(vertex0, vertex2)
        
        h = cross(ray, edge2)
        
        g = dot(edge1,h)
        
        if (g > -epsilon and g < epsilon):
            #print "This ray is parallel to this triangle."
            return False, []   # This ray is parallel to this triangle.
        
        f = 1.0/g;
        s = diff(vertex0, origin)
        
        u = f * dot(s, h)
        if (u < 0.0 or u > 1.0):
            return False, []
        
        q = cross(s, edge1)
        
        v = f * dot(ray, q)
        
        if (v < 0.0 or u + v > 1.0):
            return False, []
        # At this stage we can compute t to find out where the intersection point is on the line.
        t = f * dot(edge2, q)
        
        if t > epsilon: # ray intersection
            intersectionPoint = [origin[0] + ray[0]*t,
                                 origin[1] + ray[1]*t,
                                 origin[2] + ray[2]*t]
            return True, intersectionPoint
        else: #This means that there is a line intersection but not a ray intersection.
            return False, []
    
    # Determine if a ray intersects the mesh using the Möller–Trumbore intersection algorithm.
    # \param origin Ray origin [x_origin, y_origin, z_origin].
    # \param ray Ray vector [x_direction, y_direction, z_direction]. 
    # \param all Boolean indictating whether or not all intersections are evaluated or if the search 
    # stops after finding a single hit.
    # \return Returns the intersection 
    # point(s)  [[X1_intersect, Y1_intersect, Z1_intersect], [X2_intersect, Y2_intersect, Z2_intersect], ...], an
    # empty list is returned if no intersections are found.
    def rayIntersect(self, origin, ray, all=False):
        
        i = 0
        intersect = []
        
        while i < self.numElement:
            
            val = self._connectivity[i]
                
            j = 0
            while j < len(val)-2: # 4 goes to 2 5 -> 3, 6 -> 4 
                connect = [val[0], val[j+1], val[j+2]]
                
                triangle = [self._xyz[connect[0]-1], # Assumed one bias
                            self._xyz[connect[1]-1],
                            self._xyz[connect[2]-1]]
                
                j += 1
               
                hit, point = self._rayIntersectTriangle(origin, ray, triangle)
                if not hit: continue
                break
            
            i += 1
            
            if not hit: continue
            
            intersect.append(point[:])
            
            if not all: break
    
        return intersect
        
    def writeFile(self, filename):
        
        if not isinstance(filename, str):
            raise TypeError("Filename should be of type str!")
            
        if self.fileExt not in filename:
            filename = filename+self.fileExt
        
        print "Writing file - ", filename
        
        with open(filename, 'wb') as fp:
            fp.write("Type=2\n")
            
            fp.write("Zones=1\n") # Number of zones
            fp.write("MetaData=")
            if (self._metaData): fp.write(jsonDumps(self._metaData))
            fp.write("\n")
            
            fp.write("Variables=")
            for i in self.variable:
                if i == self.variable[-1]:
                    fp.write("{}".format(i))
                else:
                    fp.write("{},".format(i))
            fp.write("\n")
            
            fp.write("{:d} {:d}\n".format(self.numNode, self.numElement))
            
            # Nodes
            i = 0
            while i < self.numNode:
                s = structPack("d"*3, *self._xyz[i])
                fp.write(s)
                i += 1
            
            # Connectivity
            i = 0
            while i < self.numElement:
                lenCon = len(self._connectivity[i])
                s = structPack("i"*(1+lenCon), lenCon, *self._connectivity[i])
                fp.write(s)
                i += 1
            
            # Data - Only support node-based data
            for i in self.variable:
                s = structPack("d"*len(self._data[i]), *self._data[i])
                fp.write(s)
    
    def _writeSTL(self, filename, ascii = False):
        
        def cross(a, b):
            c = [a[1]*b[2] - a[2]*b[1],
                 a[2]*b[0] - a[0]*b[2],
                 a[0]*b[1] - a[1]*b[0]]

            return c
        
        def dot(a, b):
             return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        
        
        file = filename + ".stl"
        
        numTri = 0
        i = 0
        while i < self.numElement:
            numTri += len(self._connectivity[i]) - 2
            i += 1
        
        print "Write STL file"
        print "Number of triangles = ", numTri
        
        with open(file, 'wb') as fp:
                
            fp.write(self.name)
            fp.write(' '*(80-len(self.name)))
                
            if ascii: 
                raise ("ASCII format not supported yet!")
            else:
                
                s = structPack("I", numTri)
                fp.write(s)
                
                i = 0
                while i < self.numElement:
   
                    val = self._connectivity[i]
                    
                    j = 0
                    while j < len(val)-2: # 4 goes to 2 5 -> 3, 6 -> 4 
                        connect = [val[0], val[j+1], val[j+2]]
                    
                        a = self._xyz[connect[0]-1]
                        b = self._xyz[connect[1]-1]
                        c = self._xyz[connect[2]-1]
                        
                        norm = cross([b[0]-a[0], b[1]-a[1], b[2]-a[2]],
                                     [c[0]-a[0], c[1]-a[1], c[2]-a[2]])
                        
                        mag = dot(norm, norm)
                        norm[0] /= sqrt(mag)
                        norm[1] /= sqrt(mag)
                        norm[2] /= sqrt(mag)
                       
                        s = structPack("3f", *norm)
                        fp.write(s)
                        
                        s = structPack("3f", *a)
                        fp.write(s)
                        
                        s = structPack("3f", *b)
                        fp.write(s)
                        
                        s = structPack("3f", *c)
                        fp.write(s)
                        
                        s = structPack("H", 0)
                        fp.write(s)
                        
                        j += 1
                    i += 1
                   
    def loadFile(self, filename, zoneIndex=0, timeFunction=False):
        
        if timeFunction: import time
        
        if not isinstance(filename, str):
            raise TypeError("Filename should be of type str!")
            
        if self.fileExt not in filename:
            filename = filename+self.fileExt
        
        print "Reading file - ", filename
        
        sint = 4 # bytes
        sdouble = 8 # bytes 
        
        with open(filename, 'rb') as fp:
            
            type = int(fp.readline().replace("Type=","").replace("\n","")) # File type
           
            fp.readline() # Zone number
            if type == 2:
                meta = fp.readline().replace("MetaData=", "").replace("\n","")  # Meta data
                
                if meta: 
                    self._metaData = jsonLoads(meta, object_hook=_byteify)
            else :
                pass
            
            line = fp.readline().replace("Variables=", "").replace("\n","") # Variables
            
            if not line:
                self.variable = []
            else:
                self.variable = line.split(",")
     
            line = fp.readline()
            
            self.numNode = int(line.split()[0])
            self.numElement = int(line.split()[1])
            
            print "Number of nodes", self.numNode
            print "Number of elements", self.numElement
            print "Variables = ", self.variable
            
            # Pre-allocate arrays
            self._xyz = [[0.0,0.0,0.0]]*self.numNode
            self._connectivity = [None]*self.numElement
            for i in self.variable:
                self._data[i] = [] #[0.0]*self.numNode
            
            # Nodes
            print "Reading nodes"
            if timeFunction: start = time.time()
           
            i = 0
            while i < self.numNode:
                self._xyz[i] = list(structUnpack("d"*3, fp.read(sdouble*3)))
                i += 1

            if timeFunction: 
                end = time.time()
                print "Node read - ", end-start, "(s)"
            
            # Connectivity
            print "Reading elements"   
            if timeFunction: start = time.time() 
    
            i = 0
            while i < self.numElement:
                lenCon = structUnpack("i", fp.read(sint))[0]
                self._connectivity[i] = list(structUnpack("i"*lenCon, fp.read(sint*lenCon)))
                i += 1
            if timeFunction:
                end = time.time()
                print "Element read - ", end-start, "(s)"
            
            # Data
            print "Reading data"
            if timeFunction: start = time.time()
            
            for i in self.variable:
                self._data[i] = list(structUnpack("d"*self.numNode, fp.read(sdouble*self.numNode)))
            if timeFunction:
                end = time.time()
                print "Data read - ", end-start, "(s)"
    
    def _loadSTL(self, filename, stitch=True):
        
        def sortNode(node, connect):
            
            if stitch == True:
                # Does the vertex already exist?           
                try:
                    index = self._xyz.index(node)
                    connect.append(index+1) # Yes, get the index for the connectivity
                                
                except ValueError: # No! Append it 
                    self._xyz.append(node)
                    
                    self.numNode += 1
                    connect.append(self.numNode)
                except:
                    raise
            else: 
#                 print "Not stitching stl file!"
                self._xyz.append(node)
                    
                self.numNode += 1
                connect.append(self.numNode)
            
        sint  = 4 # bytes
        schar = 1 # bytes
        sfloat = 4 # bytes 
        sshort = 2 # bytes
                            
        ascii = True
        
        with open(filename, 'r') as fp:
            
            # ASCII or binary
            line = fp.read(schar*80)
            
            if "solid" in line: # ASCII
                
                print "'solid' tag found! Assuming STL is in ASCII format...."
                
                fp.seek(0) # Return to the beginning and get the whole line
                
                line = fp.readline()
                
                line = line.split()
                if len(line) > 1:
                    self.name = str(line[1])
                else:
                    self.name = self.filename
         
            else: # Binary
                    
                print "No 'solid' tag found! Assuming STL is in binary format...."
            
                ascii = False
                    
                fp.seek(0)
                
                self.name = str(fp.read(schar*80)) #structUnpack("s", fp.read(schar*80))[0]
                    
                self.numElement = structUnpack("I", fp.read(sint))[0]
                print "STL name - ", self.name
                print "Number of elements - ", self.numElement
            
            print "ASCII STL", ascii
            if ascii: # STL appears to be ASCII       
               
                for line in fp:
                    if "outer loop" in line:
                        
                        connect = []
                        for line in fp:
                            
                            if "vertex" in line:
                                data = line.split()
                                # Read vertex
                                node = [float(data[1]), 
                                        float(data[2]),
                                        float(data[3]) ]
                            
                                # Determine if node is new or not 
                                sortNode(node, connect)
                                
                            if "endloop" in line:
                                self._connectivity.append(connect)
                                self.numElement += 1
                                
                                break 
                    
                    if "endsolid" in line:
                        break
                    
            else: # Binary 
                
                # Pre-allocate connectivity array
                self._connectivity = [None]*self.numElement
                
                if stitch == False:
                    self._xyz = [None]*self.numElement*3
                    
                i = 0
                while i < self.numElement:
#                     start = time.time()
                    # Read the normals
                    fp.read(sfloat*3)
#                     fp.read(sfloat)
#                     fp.read(sfloat)
#                     fp.read(sfloat)
#                     
                    connect = []
                    j = 0
                    while j < 3:
                        
                        # Read vertex
                        node = structUnpack("f"*3, fp.read(sfloat*3))
                        
                        if stitch == True:
                            # Determine if node is new or not 
                            sortNode([node[0], node[1], node[2]], connect)
                        else:
                            self._xyz[self.numNode] = [node[0], node[1], node[2]]
                            self.numNode +=1
                            connect.append(self.numNode)
                         
                        j += 1
                    
                    # Add element to connectivity 
                    #self.connectivity.append(connect)
                    self._connectivity[i] = connect[:]
                         
                    # Read attribute bytes - currently not doing anything with it - could add it to the data 
                    attrByte = structUnpack("H", fp.read(sshort))[0]
                    
                    i += 1
#                     print "Time - ", time.time() - start
                     
       # print self.numNode, self.numElement 
        
        if self.numNode == 0 or self.numElement == 0:
            print("No nodes or elements created.... something went wrong")
            raise ValueError
        
cdef class capsFieldView:
    
    cdef:
        
        object filename
                
        public name 
        
        readonly numElement
        readonly numNode
       
        readonly variable
        readonly xyz 
        readonly object data
        readonly object connectivity
        
    
    def __cinit__(self, filename):
        
        self.name = "FieldView"
        self.numElement = 0
        self.numNode = 0
        
        self.xyz = []
        self.connectivity = []

        self.variable = []
        self.data = {}
        
        self.filename = filename
                
    def loadFile(self):
        
        if ".txt" in self.filename:
            self._loadText()
        else:
            print "Unable to determine FieldView  format!" 
            raise ValueError
    
    def _printLoading(self):
        print "Loading file", self.filename
        
    def _parseFileName(self):
        self.name = os.path.basename(self.filename)
    
    def _loadText(self):
        
        self._parseFileName()
        
        self.variable = ["BCID"]
        self.data["BCID"] = []
        
        self._printLoading()
        
        with open(self.filename, 'r') as fp:
        
            fp.readline() # Header  
            self.numNode = int(fp.readline()) # Number of nodes
            print "Number of nodes", self.numNode
            
            self.variable = fp.readline().lstrip().rstrip().split()[3:]
          
#             print self.variable
          
            # Pre-allocate xyz array 
            self.xyz = [None]*self.numNode
            
            # Pre-allocate data array
            for j, key in enumerate(self.variable):
                self.data[key] = [None]*self.numNode
            
            # Nodes xyz and data 
            for i in range(self.numNode):
                
                line = fp.readline().lstrip().rstrip().split()
                
                self.xyz[i] = [float(line[0]),
                               float(line[1]),
                               float(line[2])]
                
                for j, key in enumerate(self.variable):
                    self.data[key][i] = float(line[j+3])
            
            fp.readline() # Geometry line 
            self.numElement = int(fp.readline()) # Number of elements
            print "Number of elements", self.numElement
            
            # Pre-allocate connectivity
            self.connectivity = [None]*self.numElement
            
            for i in range(self.numElement):
                line = fp.readline().lstrip().rstrip().split()
                
                self.connectivity[i] = [int(x) for x in line[1:]] 
               