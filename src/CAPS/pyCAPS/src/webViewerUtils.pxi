#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

# Python version
try:
    from cpython.version cimport PY_MAJOR_VERSION
except: 
    print("Error: Unable to import cpython.version")
    raise ImportError

# C-Imports 
cimport cWV
cimport cCAPS

import webbrowser
import time 
import struct

from libc.stdio cimport FILE, fopen, fclose, fprintf, fscanf, printf


#cdef int numInstance

ctypedef struct globalContainer:
    void ** gprim
    int numPrim

cdef globalContainer globalInstance

#cdef cWV.wvContext *globalInstanceontext

cdef clearGlobal():
    global globalInstance
    
    if globalInstance.gprim:
        free(globalInstance.gprim)
        
    globalInstance.gprim = NULL 
    
    globalInstance.numPrim = 0
    
    
# Set the scalar color map associated with data. It is assumed data is in a list. 
# Return a void * (must be freed) for the scalar color map (which if of float type) 
cdef void * setScalarColor(object data, double dataMin, double dataMax, int numContour, int reverseMap, char * cMapName):

    cdef: 
        float *color = NULL
        double scalar, frac
        int index, numColor
    
    if isinstance(data, str):
        numData = 1
    else:
        numData = len(data)

    color = <float *> malloc(3*numData*sizeof(float))
    if not color:
        raise MemoryError()
    
    if isinstance(data, str):
        
        h = data.lstrip('#')
      
        color[0] = <float> int(h[0:2], 16)/255.0
        color[1] = <float> int(h[2:4], 16)/255.0
        color[2] = <float> int(h[4:6], 16)/255.0
        
    else:
        if isinstance(data[0], str):

            i = 0
            while i < numData:
                h = data[i].lstrip('#')
                color[3*i + 0] = <float> int(h[0:2], 16)/255.0
                color[3*i + 1] = <float> int(h[2:4], 16)/255.0
                color[3*i + 2] = <float> int(h[4:6], 16)/255.0
                
                i += 1
        else:
     
            colorMap = getColorMap(<object> cMapName, numContour, reverseMap)
            numColor = int(len(colorMap)/3) - 1
            
            if numColor == 0:
                dataMin = dataMax
        
            i = 0
            for scalar in data:
                if (dataMin == dataMax):
                    color[3*i + 0] = <float> colorMap[-3]
                    color[3*i + 1] = <float> colorMap[-2]
                    color[3*i + 2] = <float> colorMap[-1]
          
                elif scalar <= dataMin:
                    color[3*i + 0] = <float> colorMap[0]
                    color[3*i + 1] = <float> colorMap[1]
                    color[3*i + 2] = <float> colorMap[2]
          
                elif scalar >= dataMax:
                    color[3*i + 0] = <float> colorMap[-3]
                    color[3*i + 1] = <float> colorMap[-2]
                    color[3*i + 2] = <float> colorMap[-1]
                    
                else:
            
                    frac  = numColor * (scalar- dataMin) / (dataMax- dataMin)
                    
                    if (frac < 0):
                        frac = 0
                    if (frac > numColor):
                        frac = numColor
            
                    index  = int(frac)
                    
                    #frac -= index
                    
                    #if index == numColor:
                    #    index = index -1
                    #   frac += 1.0
            
                    #color[3*i + 0] = <float > frac * colorMap[3*(index+1)+0] + (1.0-frac) * colorMap[3*index+0]
                    #color[3*i + 1] = <float > frac * colorMap[3*(index+1)+1] + (1.0-frac) * colorMap[3*index+1]
                    #color[3*i + 2] = <float > frac * colorMap[3*(index+1)+2] + (1.0-frac) * colorMap[3*index+2]
                    
                    color[3*i + 0] = <float> colorMap[3*(index+1)+0]
                    color[3*i + 1] = <float> colorMap[3*(index+1)+1]
                    color[3*i + 2] = <float> colorMap[3*(index+1)+2]
                
                #print color[3*i + 0], color[3*i + 1], color[3*i + 2]
                i = i + 1
        
    return <void *> color


# Look for color attribute on the object - if available; return None if not found
cdef object __getColor(cEGADS.ego egoObj):
    cdef:
        int status
        int atype
        int alen
        const int *ints
        const double *reals
        const char *string
    
    color = None
  
    attr = _byteify("_color")
    status = cEGADS.EG_attributeRet(egoObj, <char *>attr, &atype, &alen, &ints, &reals, &string)
    if status != cEGADS.EGADS_SUCCESS and \
       status != cEGADS.EGADS_NOTFOUND: 
        raise CAPSError(status, "while getting color (during a call to EG_attributeRet)");
    
    if status == cEGADS.EGADS_SUCCESS:
        if atype == cEGADS.ATTRREAL and alen == 3:
            
            red   = min(max(0, int(255*reals[0])), 255);
            green = min(max(0, int(255*reals[1])), 255);
            blue  = min(max(0, int(255*reals[2])), 255);

            color ='#%02x%02x%02x' % (red, green, blue)
            
        elif atype == cEGADS.ATTRSTRING:
            colorString = _strify(string)
    
            if "red" in colorString:
                color = '#ff0000'
            elif "green" in colorString:
                color = '#00ff00'
            elif "blue" in colorString:
                color = '#0000ff'
            elif "yellow" in colorString:
                color = '#ffff00'
            elif "magenta" in colorString:
                color = '#ff00ff'
            elif "cyan" in colorString:
                color = '#00ffff'
            elif "white"  in colorString:
                color = '#ffffff'
            elif "black" in colorString:
                color = '#000000'
                
    return color


# Required call-back to catch browser messages
# 'with gil' enables calling CPython functions from within browserMessage
# https://cython.readthedocs.io/en/latest/src/userguide/external_C_code.html#acquiring-the-gil
cdef public void browserMessage(void *wsi, char *inText, int lenText) with gil:
    
    cdef:
        void *temp
        char *msg
        double upper, lower
        
        int i=0, istart=0, iend=0
        
        double *colorRange 
        
        int numContour
        
        FILE  *fp=NULL
        
        char msgFinite[1024]
        char cscale[512]
        char cmatrix[512]

    # convert text to a python string
    if PY_MAJOR_VERSION >= 3:
        Text = inText[:lenText].decode('ascii', 'strict')
    else:
        Text = <bytes> inText[:lenText]

    cscale[0] = b'\0'
    cmatrix[0] = b'\0'
    msgFinite[0] = b'\0'
    msg = msgFinite

    #bText = bytes(Text, "ascii")
    #printf("Server message: %s\n", <char*> bText)

    if 'nextScalar' == Text:
        
        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
            
            # build the response 
            msg = prim.nextColor()
            
            if (msg == NULL):
                msg = "No colors!"
            
            prim.raiseException = True
            
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break

    elif 'previousScalar' == Text:
        
        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
            
            # build the response 
            msg = prim.previousColor()
            if (msg == NULL):
                msg = "No colors!"
            
            prim.raiseException = True
            
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break
                 
    elif "limits" in Text:
        
        for i in range(lenText):
                
            if b',' == inText[i] and istart == 0:
                istart = i
            elif b',' == inText[i]:
                iend = i
        
        # Need special treatment for negatives to a avoid a seg-fault in PyFloat_FromString
        if inText[istart+1] == b"-": 
            lower = <double > float(inText[istart+2:iend])
            lower = -lower
        else: 
            lower = <double > float(inText[istart+1:iend])
        
        if inText[iend+1] == b"-":
            upper = <double > float(inText[iend+2:])
            upper = -upper 
        else:
            upper = <double > float(inText[iend+1:])

        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
              
            prim.changeLimit(lower, upper)
            
            prim.raiseException = True
            
            # build the response 
            msg = "Limits updated"
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break
        
    elif "colorMap," in Text:
        
        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
              
            prim.changeColorMap(inText[9:])
            
            prim.raiseException = True
            
            # build the response 
            msg = "ColorMap updated"
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break
        
    elif "reverseMap" in Text:
        
        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
              
            prim.reverseColorMap()
            
            prim.raiseException = True
            
            # build the response 
            msg = "ColorMap reversed/inverted"
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break
        
    elif "getRange" in Text:
    
        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
            
            colorRange = prim.getRange()
            
            # build the response 
            sprintf(msgFinite, "Range:%10.6g,%10.6g", colorRange[0], colorRange[1])
            
            prim.raiseException = True
            
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break
        
    elif "contourLevel," in Text:
        
        istart = 12
        
        # Need special treatment for negatives to a avoid a seg-fault in PyFloat_FromString
        if inText[istart+1] == b"-": 
            numContour = < int > float(inText[istart+2:])
            numContour = -numContour
        else: 
            numContour = < int > float(inText[istart+1:])
            
        for i in range(globalInstance.numPrim):
            temp = globalInstance.gprim[i]
        
            prim = <graphicPrimitive > temp
            
            prim.raiseException = False
              
            prim.changeContourLevel(numContour)
            
            prim.raiseException = True
            
            # build the response 
            msg = "Contour level updated"
            if prim.status < 0:
                msg = msgFinite
                sprintf(msg, "Error:%+d", prim.status)
                break
 
    # "saveView|viewfile|scale|matrix" 
    elif "saveView|" in Text:

        # extract arguments
        viewfile, scale, matrix = Text.split("|")[1:4]

        viewfile = bytes(viewfile, 'ascii')
        scale    = bytes(scale   , 'ascii')
        matrix   = bytes(matrix  , 'ascii')

        # save the view in viewfile 
        fp = fopen(<char*> viewfile, "w")
        if fp != <FILE*>NULL:
            fprintf(fp, "%s\n", <char*>scale);
            fprintf(fp, "%s\n", <char*>matrix);
            fclose(fp);
        
        # build the response 
        msg = "saveView|";

    # "readView|viewfile" 
    elif "readView|" in Text:
        # extract arguments 
        viewfile = Text.split("|")[1]

        # read the view from viewfile 
        viewfile = bytes(viewfile, 'ascii')
        fp = fopen(<char*> viewfile, "r");
        if fp != <FILE*>NULL:
            fscanf(fp, "%s", cscale);
            fscanf(fp, "%s", cmatrix);
            fclose(fp);

            # build the response 
            sprintf(msg, "readView|%s|%s|", cscale, cmatrix);
        else:
            # build the response that the file did not exist
            sprintf(msg, "readView|%s", <char*>viewfile);

    else:
        cWV.wv_sendText(wsi, "Command not recognized!")

    #printf("Sending message: %s\n", msg);
    cWV.wv_sendText(wsi, msg)
    
    return 

# Copy an item. Note that isn't a deep copy - the data in the data pointer isn't being copied 
cdef void copyItem(cWV.wvData *inVal, cWV.wvData *outVal):
    
    outVal.dataType = inVal.dataType   #     /* VBO type */
    outVal.dataLen  = inVal.dataLen    #    /* length of data */
    outVal.dataPtr  = inVal.dataPtr    #     /* pointer to data */

    outVal.data[0] = inVal.data[0] #    /* size = 1 data (not in ptr) *
    outVal.data[1] = inVal.data[1] #    /* size = 1 data (not in ptr) *
    outVal.data[2] = inVal.data[2] #    /* size = 1 data (not in ptr) *
        
cdef class wvData:
    cdef:
        int status
        
        cWV.wvData item
        cWV.wvData *itemPointer
        
        public object name # Name of the item. Keep as an object to avoid string memory complications
        
        public object colorData # Need to store way the color data in case we are going to change the limits
        public object colorMap # Name of the color map (if any). Keep as an object to avoid string memory complications
        
        public double minColor 
        public double maxColor 
        
        public double minColorRange  # Initial max and min color value 
        public double maxColorRange   
        
        public int numContour
        public int reverseMap
    
        public int raiseException
        
        #public float focus[4]
        public float boundingBox[6] # Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
        
    def __cinit__(self):
        
        cdef: 
            int i
            char *nameTemp = "UnknownVariable"
            
        self.status = 0
        
        self.itemPointer = NULL
        
        self.name = nameTemp
               
        self.colorData
        self.colorMap = None
        
        self.minColor = -1.0
        self.maxColor = 1.0
        
        self.minColorRange = cCAPS.CAPSMAGIC
        self.maxColorRange = cCAPS.CAPSMAGIC
        
        self.numContour = 0
        
        self.reverseMap = <int> False
        
        self.raiseException = True
        
        #self.focus[0] = 0.0
        #self.focus[1] = 0.0
        #elf.focus[2] = 0.0
        #self.focus[3] = 0.0
        
        for i in range(6):
            self.boundingBox[i] = 0.0
        
    def __del__(self):
        cWV.wv_freeItem(&self.item)
        
        self.itemPointer = NULL
        
    def __dealloc__(self):
        cWV.wv_freeItem(&self.item)
        
        self.itemPointer = NULL
    
    cdef cWV.wvData* __wvDataItem(self):
        
        return self.itemPointer
    
    cdef freeItem(self):
        cWV.wv_freeItem(&self.item)
        
    cdef createVertex(self, object xyz, char *name = NULL):
        
        if name:
            self.name = name
            
        if not isinstance(xyz, list):
            xyz = [xyz]
        
        length = len(xyz)
        
        data = castValue2VoidP(xyz, None, cCAPS.Double)
                            
        self.status = cWV.wv_setData(cWV.WV_REAL64, <int> length, data,  cWV.WV_VERTICES, &self.item)
        if data:
            free(data)
        
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createVertex), code error "+str(self.status))
        
        return self 
    
    cdef createIndex(self, object connectivity, char *name = NULL):
        
        cdef: 
            int length=0
            object removeList = []

        if name:
            self.name = name
        
        guess = len(connectivity)
        connect = [None]*guess
        
        numElement = 0
        
        i = 0
        while i < len(connectivity):
        #for i, val in enumerate(connectivity):
            val = connectivity[i]
            i += 1
            
            if len(val) > 3:
                
                j = 0
                while j < len(val)-2: # 4 goes to 2 5 -> 3, 6 -> 4 
                #for j in range(len(val)-2): # 4 goes to 2 5 -> 3, 6 -> 4 
                    
                    if numElement < guess:
                        connect[numElement] = [val[0], val[j+1], val[j+2]]
                        # length += len(connect[numElement])
                    else:
                        connect.append([val[0], val[j+1], val[j+2]])
                        # length += len(connect[-1])
                    
                    length += 3
                    numElement += 1
                    
                    j += 1
            else:
                length += len(val)
                
                if numElement < guess:
                    connect[numElement] = val[:]
                else:
                    connect.append(val)
                
                numElement += 1
            
        data = castValue2VoidP(connect, None, cCAPS.Integer)

        self.status = cWV.wv_setData(cWV.WV_INT32, length, data, cWV.WV_INDICES, &self.item)
        if data:
            free(data)
            
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createIndex), code error "+str(self.status))
        
        return self
    
    cdef createLineIndex(self, object connectivity, char *name = NULL):
        
        cdef: 
            int length=0 

        if name:
            self.name = name
            
        segment = []
        for i in connectivity:
            
            if len(i) == 2: # If just a line - don't loop back on self
                for j in i:
                    segment.append(j)
                    length += 1
            else:
                j = 0
                while j < len(i):
                #for j in range(len(i)):
                    
                    segment.append(i[j])
                    length += 1
                    
                    if j == len(i)-1:
                        segment.append(i[0])
                    
                    else:
                        segment.append(i[j+1])
                    
                    length += 1
                    
                    j += 1
    
        data = castValue2VoidP(segment, None, cCAPS.Integer)

        self.status = cWV.wv_setData(cWV.WV_INT32, length, data, cWV.WV_LINDICES, &self.item)
        if data:
            free(data)
    
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createLineIndex), code error "+str(self.status))
        
        return self
    
    cdef createPointIndex(self, object connectivity, char *name = NULL):
        
        cdef: 
            int length=0 

        if name:
            self.name = name

        length = len(connectivity)
        data = castValue2VoidP(connectivity, None, cCAPS.Integer)

        self.status = cWV.wv_setData(cWV.WV_INT32, length, data, cWV.WV_PINDICES, &self.item)
        if data:
            free(data)
    
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createPointIndex), code error "+str(self.status))
        
        return self
        
    cpdef updateColor(self, double minColor, double  maxColor, int numContour, int reverseMap, char *colorMap):
    
        cdef: 
            int length
            void *data = NULL
        
        self.freeItem()
        
        if isinstance(self.colorData, str):
            length = 1
        else:
            length = len(self.colorData)
        
        self.minColor = minColor 
        self.maxColor = maxColor
        self.numContour = numContour
        self.reverseMap = reverseMap
        
        if colorMap:
            self.colorMap = colorMap    
        
        # Set color returns a void * to a float *
        data = setScalarColor(self.colorData, 
                              self.minColor, self.maxColor, 
                              self.numContour, self.reverseMap, 
                              <char *> self.colorMap)

        self.status = cWV.wv_setData(cWV.WV_REAL32, length, data, cWV.WV_COLORS, &self.item)
        if data:
            free(data)
        
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (updateColor), code error "+str(self.status))
        
        return self
    
    cdef createColor(self, object colorData, double minColor, double  maxColor, int numContour, int reverseMap, char * colorMap, char * name=NULL): 
        cdef: 
            int length
            void *data = NULL
        
        self.freeItem()
        
        if colorMap:
            self.colorMap = colorMap 
                 
        if self.minColorRange == cCAPS.CAPSMAGIC:
            self.minColorRange = minColor
        
        if self.maxColorRange == cCAPS.CAPSMAGIC:
            self.maxColorRange = maxColor 
        
        self.minColor = minColor 
        self.maxColor = maxColor
        
        self.numContour = numContour
        
        self.reverseMap = reverseMap
        
        if name:
            self.name = name
        
        if isinstance(colorData, str):
           
            length = 1
            
        else:
            if not isinstance(colorData, list):
                colorData = [colorData]
        
            length = len(colorData)
        
        self.colorData = colorData
            
        # Set color returns a void * to a float *
        data = setScalarColor(colorData, self.minColor, self.maxColor, self.numContour, self.reverseMap, self.colorMap)

        self.status = cWV.wv_setData(cWV.WV_REAL32, length, data, cWV.WV_COLORS, &self.item)
        if data:
            free(data)
        
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createColor), code error "+str(self.status))
        
        return self
    
    cdef createLineColor(self, object colorData, double minColor, double maxColor, int numContour, int reverseMap, char *colorMap=NULL, char *name = NULL):
        
        cdef: 
            int length
            void *data = NULL
        
        if colorMap:
            self.colorMap = colorMap
            
        if name:
            self.name = name
        
        if isinstance(colorData, str):
           
            length = 1
        else: 
            if not isinstance(colorData, list):
                colorData = [colorData]
            
            length = len(colorData)
        
        # Set color returns a void * to a float *
        data = setScalarColor(colorData, minColor, maxColor, numContour, reverseMap, colorMap)

        self.status = cWV.wv_setData(cWV.WV_REAL32, length, data, cWV.WV_LCOLOR, &self.item)
        if data:
            free(data)
        
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createLineColor), code error "+str(self.status))
        
        return self
    
    cdef createPointColor(self, object colorData, double minColor, double maxColor, int numContour, int reverseMap, char *colorMap=NULL, char *name = NULL):
        
        cdef: 
            int length
            void *data = NULL
        
        if colorMap:
            self.colorMap = colorMap
            
        if name:
            self.name = name
        
        if isinstance(colorData, str):
           
            length = 1
        else: 
            if not isinstance(colorData, list):
                colorData = [colorData]
            
            length = len(colorData)
        
        # Set color returns a void * to a float *
        data = setScalarColor(colorData, minColor, maxColor, numContour, reverseMap, colorMap)

        self.status = cWV.wv_setData(cWV.WV_REAL32, length, data, cWV.WV_PCOLOR, &self.item)
        if data:
            free(data)
        
        self.itemPointer = &self.item
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting data (createLineColor), code error "+str(self.status))
        
        return self
    
    def setBoundingBox(self, object box):
        if len(box) < 6:
            self.status = cCAPS.CAPS_RANGEERR
            
            if self.raiseException:
                raise CAPSError(self.status, msg="while setting bounding box data")
        
        for i in range(6):
            self.boundingBox[i] = box[i]
            
    def adjustVertex(self, object focus):
        cdef: 
            float focal[4]
            
        if len(focus) < 4:
            self.status = cCAPS.CAPS_RANGEERR
            
            if self.raiseException:
                raise CAPSError(self.status, msg="while setting adjust vertex data")
            else:
                return 
        
        for i in range(4):
            focal[i] = focus[i]   
            
        cWV.wv_adjustVerts(&self.item, focal)


cdef class graphicPrimitive:
    cdef:
        _wvContext context # Context class
        
        int status 
        
        readonly object colors # List of color wvData object 
        
        int currentColor
        
        int primIndex
        
        int raiseException
        
        readonly object name # Name of the primitive
        
        double range[2]
        
        readonly object vertex # Copy of vertices wvData object
        
    def __cinit__(self, context):
        
        cdef: 
            char *tempName = "UnknownPrimitive"
            
        if not isinstance(context, _wvContext):
            raise TypeError
        
        self.context = context
        
        self.colors = []
        self.currentColor = 0
        
        self.status = 0
        
        self.primIndex = 0
        
        self.raiseException = True
        
        self.name = tempName
        
        self.range[0] = -1.0
        self.range[1] = 1.0
        
        self.vertex = None
        
    def __del__(self):
        
        self.primIndex = 0
        
        self.colors =[]
        self.currentColor = 0
        
        self.context = None
        
    def __dealloc__(self):
        
        self.primIndex = 0
        
        self.colors =[]
        self.currentColor = 0 
        
        self.context = None
 
    # Append a or multiple colors to the color list after even after a primitive has been created 
    cpdef appendColor(self, object items):
        cdef: 
           cWV.wvData *temp =NULL
    
        if not isinstance(items, list):
            items = [items]
            
        if not isinstance(items[0], wvData):
            if self.raiseException: 
                raise TypeError
            else:
                self.status = -99
                return self.status
            
        for i in items:
            temp = (<wvData>i).__wvDataItem()
           
            if temp.dataType == cWV.WV_COLORS:
                self.colors.append(i)
                
    # Update vertex - assumes the item the "vertex" corresponds to has been changed (i.e. a adjustVertex call was made)
    cpdef updateVertex(self):
        
        if self.colors:
            
            self.modifyPrimitive([self.vertex, self.colors[self.currentColor-1]])
        else : 
            
            self.modifyPrimitive(self.vertex)
    
    # \param items List of wvData objects used to create the triangle 
    # \param name Name of graphic primitive. Important: Name must be unique 
    # \return Instance of  graphicPrimitive class 
    def createPoint(self, items, name):
        cdef:
           int attrs = 0
        
        if not isinstance(items, list):
            items = [items]
            
        if not isinstance(items[0], wvData):
            raise TypeError
            
        return self.createPrimative(items, name, "point", attrs)
    
    # \param items List of wvData objects used to create the triangle 
    # \param name Name of graphic primitive. Important: Name must be unique 
    # \param lines Turn lines in the primitive on or off.
    # \return Instance of  graphicPrimitive class 
    def createLine(self, items, name, lines=True):
        cdef:
           int attrs = 0
           
        if not isinstance(items, list):
            items = [items]
            
        if not isinstance(items[0], wvData):
            raise TypeError
        
        if lines: attrs = cWV.WV_LINES
        
        return self.createPrimative(items, name, "line", attrs)
        
    # \param items List of wvData objects used to create the triangle 
    # \param name Name of graphic primitive. Important: Name must be unique 
    # \param lines Turn lines in the primitive on or off.
    # \return Instance of  graphicPrimitive class 
    def createTriangle(self, items, name, lines=True):
        cdef:
           int attrs = 0
        
        if not isinstance(items, list):
            items = [items]
            
        if not isinstance(items[0], wvData):
            raise TypeError
            
        if lines: attrs = cWV.WV_LINES
        
        return self.createPrimative(items, name, "triangle", attrs)
    
    # \return The item name for the color being set 
    cdef char *nextColor(self):
        
        cdef:
            object item 
            
        if len(self.colors) == 0:
            return <char *> NULL 
        
        if self.currentColor == len(self.colors):
        
            self.currentColor = 1
        
        else:
            self.currentColor += 1 # if not the last color
            
        item = self.colors[self.currentColor -1]
        
        # Update primitive with the new color 
        self.modifyPrimitive(item)
        
        if not self.raiseException:
            if self.status < 0:
                return NULL
            
        # Update the key
        self.setKey()
        
        if not self.raiseException:
            if self.status < 0:
                return NULL
        
        # Return the item's name 
        return <char *> item.name
    
    # \return The item name for the color being set 
    cdef char * previousColor(self):
        
        cdef:
            object item 
        
        if len(self.colors) == 0:
            return <char *> NULL
        
        if self.currentColor == 1:
        
            self.currentColor = len(self.colors)
        
        else:
          
            self.currentColor -= 1 # if not the first color 
            
        item = self.colors[self.currentColor -1]
        
        # Update primitive with the new color 
        self.modifyPrimitive(item)
        
        if not self.raiseException:
            if self.status < 0:
                return NULL
        
        # Update the key
        self.setKey()
        
        if not self.raiseException:
            if self.status < 0:
                return NULL
    
        # Return the item's name 
        return <char *> item.name
    
    # Change the limits of the current data set 
    cdef void changeLimit(self, double minColor, double maxColor):
        cdef:
            wvData item
        
        if len(self.colors) == 0:
            return
        
        # Create a new color based on the current colors data, but change the limits
        item = self.colors[self.currentColor -1]
       
        item.raiseException = False
        
        item.updateColor(minColor, maxColor, item.numContour, item.reverseMap, item.colorMap)
         
        if not item.raiseException:
            if item.status < 0:
                
                item.raiseException = True
               
                self.status = item.status
                
                return 
        
        item.raiseException = True
        
        # Update primitive with the new color 
        self.modifyPrimitive(item)
        
        if not self.raiseException:
            if self.status < 0:
                return 
         
        # Update the key
        self.setKey()
        
        if not self.raiseException:
            if self.status < 0:
                return 
    
    # Change all the color maps 
    cdef void changeColorMap(self, char *colorMap):
        cdef:
            wvData item
            int i 
        
        if len(self.colors) == 0:
            return
        
        for i in range(len(self.colors)):
             
            # Create a new color based on the current colors data, but change the colorMap
            item = self.colors[i]
            
            item.raiseException = False
            
            item.updateColor(item.minColor, item.maxColor, item.numContour, item.reverseMap, colorMap)
            
            if not item.raiseException:
                if item.status < 0:
                
                    item.raiseException = True
               
                    self.status = item.status
                
                    return 
        
            item.raiseException = True
            
            # Only update the primitive if it is our current color 
            if i == self.currentColor-1:
                self.modifyPrimitive(item)
                
                if not self.raiseException:
                    if self.status < 0:
                        return 
                
        # Update the key
        self.setKey()
        
        if not self.raiseException:
            if self.status < 0:
                return 
    
    # Change all the contour levels  
    cdef void changeContourLevel(self, int numContour):
        cdef:
            wvData item
            int i 
        
        if len(self.colors) == 0:
            return
        
        for i in range(len(self.colors)):
             
            # Create a new color based on the current colors data, but change the colorMap
            item = self.colors[i]
            
            item.raiseException = False
            
            item.updateColor(item.minColor, item.maxColor, numContour, item.reverseMap, item.colorMap)
            
            if not item.raiseException:
                if item.status < 0:
                
                    item.raiseException = True
               
                    self.status = item.status
                
                    return 
        
            item.raiseException = True
            
            # Only update the primitive if it is our current color 
            if i == self.currentColor-1:
                self.modifyPrimitive(item)
                
                if not self.raiseException:
                    if self.status < 0:
                        return 
    
        # Update the key
        self.setKey()
        
        if not self.raiseException:
            if self.status < 0:
                return 
    
    # Reverse/invert all the color maps  
    cdef void reverseColorMap(self):
        cdef:
            wvData item
            int i 
        
        if len(self.colors) == 0:
            return
        
        for i in range(len(self.colors)):
             
            # Create a new color based on the current colors data, but change the colorMap
            item = self.colors[i]
            
            item.raiseException = False
            
            if item.reverseMap == True:
                item.updateColor(item.minColor, item.maxColor, item.numContour, False, item.colorMap)
            else:
                item.updateColor(item.minColor, item.maxColor, item.numContour, True, item.colorMap)
            
            if not item.raiseException:
                if item.status < 0:
                
                    item.raiseException = True
               
                    self.status = item.status
                
                    return 
        
            item.raiseException = True
            
            # Only update the primitive if it is our current color 
            if i == self.currentColor-1:
                self.modifyPrimitive(item)
                
                if not self.raiseException:
                    if self.status < 0:
                        return 
    
        # Update the key
        self.setKey()
        
        if not self.raiseException:
            if self.status < 0:
                return

    # Return the max and min initial color values (range)
    cdef double *getRange(self):
        cdef:
            double * range = NULL
        
        if len(self.colors) == 0:
            return self.range
        
        # Create a new color based on the current colors data, but change the limits
        item = self.colors[self.currentColor -1]
        
        self.range[0] = item.minColorRange
        self.range[1] = item.maxColorRange
        
        return self.range
    
    # Set the key for the context
    #cdef void setKey(self):
    cpdef setKey(self):
        
        cdef:
            wvData item
            
            float *keyData = NULL
            int i, numCol
            
            char *title
            
            object cMapName
            
            float minColor 
            float maxColor 
            
            object colorMap
            
            int numContour 
            int reverseMap
            
        if len(self.colors) == 0:
            return
        
        # Get current color and it's information 
        item = self.colors[self.currentColor -1]
        
        
        title    = <char*>  item.name
        minColor = <float>  item.minColor
        maxColor = <float>  item.maxColor
        cMapName = <object> item.colorMap
        numContour = <int>  item.numContour
        reverseMap = <int>  item.reverseMap
        
#         print cMapName, numContour, minColor, maxColor
        # Get color map 
        colorMap = getColorMap(cMapName, numContour, reverseMap)
       
        # Create key data 
        if minColor == maxColor:
            numCol = 1
        else:
            numCol = int(len(colorMap)/3)
        
        keyData = <float *> malloc(3*numCol*sizeof(float))
        if not keyData:
            raise MemoryError

        if minColor == maxColor:
            keyData[0] = <float> colorMap[-3]
            keyData[1] = <float> colorMap[-2]
            keyData[2] = <float> colorMap[-1]
            
        else:
            i = 0
            while i < 3*numCol:
                keyData[i] = <float> colorMap[i]
                i += 1
#                 print keyData[i-1],colorMap[i-1]
                
        self.status = cWV.wv_setKey(self.context.context, numCol,  
                                                          keyData, 
                                                          minColor, maxColor, 
                                                          title)
        
        if keyData:
           free(keyData)
        
        if self.raiseException:
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while setting key (setKey), code error "+str(self.status))
    
    cdef int createPrimative(self, list items, char * name, char *type, int attrs) except -1:
        cdef:
            
            int numItem = 0
            cWV.wvData *tempItem =NULL
            
            cWV.wvData *temp =NULL
            
            #int colorCount = 0
        
        # Store away the name
        self.name = name
                
        for i in items:
            
            temp = (<wvData>i).__wvDataItem()
            
            if temp.dataType == cWV.WV_COLORS:
                self.appendColor(i)
            
                if len(self.colors) > 1:
                    continue
                elif len(self.colors) == 1:
                    self.currentColor = 1
            
            #   self.currentColor = 1
            #   self.colors.append(i) # Store the wvData instance 
                
            #    colorCount += 1
                
            #if  temp.dataType == cWV.WV_COLORS and colorCount > 1: 
            # 
            #     continue
            
            if not self.vertex:
                if temp.dataType == cWV.WV_VERTICES:
                    self.vertex = i
        
            # Copy the item data temporarily 
            numItem += 1
            
            tempItem = <cWV.wvData *> realloc(tempItem, <int>numItem*sizeof(cWV.wvData))
            if not tempItem:
                if self.raiseException: 
                    raise MemoryError
                else:
                    self.status = -99
                    return self.status 
                
            # Note that isn't a deep copy - the data in the data pointer isn't being copied 
            copyItem(temp, &tempItem[numItem-1])
        
        if _strify(type) == "point":
            self.status = cWV.wv_addGPrim(self.context.context,
                                          <char *> name, 
                                          cWV.WV_POINT, 
                                          cWV.WV_ON+cWV.WV_ORIENTATION+attrs,
                                          numItem, 
                                          tempItem)
        elif _strify(type) == "line":
            self.status = cWV.wv_addGPrim(self.context.context,
                                          <char *> name, 
                                          cWV.WV_LINE, 
                                          cWV.WV_ON+cWV.WV_ORIENTATION+attrs,
                                          numItem, 
                                          tempItem)
        elif _strify(type) == "triangle":
            self.status = cWV.wv_addGPrim(self.context.context,
                                          <char *> name, 
                                          cWV.WV_TRIANGLE, 
                                          cWV.WV_ON+cWV.WV_ORIENTATION+cWV.WV_SHADING+attrs,
                                          numItem, 
                                          tempItem)
        else: 
            
            if tempItem:
                free(tempItem)
            
            tempItem = NULL
            
            if self.raiseException: 
                raise TypeError
            else:
                self.status = -99
                return self.status 

        #################################################
        #################################################
        # TODO: Ryan we need to expose pSize somehow...
        #################################################
        #################################################
        if self.status >= 0:
            self.context.context.gPrims[self.status].pSize = 6.0 

        if tempItem:
                free(tempItem)
        
        if self.raiseException:      
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while adding graphic primitive (wv_addGPrim), code error "+str(self.status))

        self.primIndex = self.status
        
        return <int> self.primIndex
              
    cdef void modifyPrimitive(self, items):
        
        cdef: 
           cWV.wvData *tempItem =NULL
           cWV.wvData *temp =NULL
           
           
        if not isinstance(items, list):
            items = [items]
            
        if not isinstance(items[0], wvData):
            raise TypeError
        
        tempItem = <cWV.wvData *> malloc(len(items)*sizeof(cWV.wvData))
        if not tempItem:
            raise MemoryError
        
        for i, item in enumerate(items):
            
            temp = (<wvData>item).__wvDataItem()
            
            # Note that isn't a deep copy - the data in the data pointer isn't being copied 
            copyItem(temp, &tempItem[i])
            
        self.status = cWV.wv_modGPrim(self.context.context, 
                                      self.primIndex,
                                      <int> len(items), 
                                      tempItem)
        if tempItem:
            free(tempItem)
        
        if self.raiseException:
            
            if self.status < 0:
                raise CAPSError(self.status-1000, msg="while modifying graphic primitive (wv_modGPrim), code error"+str(self.status))


cdef class _wvContext:
    cdef:
        int status 
        cWV.wvContext *context
        
        readonly object portNumber
        int serverInstance
                
        readonly int numGPrim
        object graphicPrim 
        object graphicPrimName
        
        # KeyData
        int keyNumCol
        float *keyData
        float keyMax
        float keyMin
        object keyTitle
            
        readonly object oneBias
        
        object buffer
        object bufferIndex
        
        int bufferLen
        
    def __cinit__(self, portNumber=7681,  oneBias=True, bufferLen = 3205696):
        cdef:
            float eye[3]
            float center[3]
            float up[3]
        
        eye[:] = [0.0, 0.0, 7.0]
        center[:] = [0.0, 0.0, 0.0]
        up[:] = [0.0, 1.0, 0.0]
        
        self.graphicPrim = []
        self.graphicPrimName = []
        
        self.portNumber = portNumber
        self.oneBias = oneBias 
        
        self.context = NULL
        
        self.context = cWV.wv_createContext(<int> int(self.oneBias), 30.0, 1.0, 10.0, eye, center, up)
        if self.context is NULL: 
            self.status = cCAPS.CAPS_NULLOBJ
            raise CAPSError(self.status, msg="while creating wvContext")
        
        self.context.keepItems = 1
        
        self.bufferLen = bufferLen
        self.buffer = bytearray(self.bufferLen)
        self.bufferIndex = 0
            
    def setKey(self, primitiveIndex=0):
        self.graphicPrim[primitiveIndex].setKey()
        
    def addPrimative(self, type, items, name=None, lines=True ):
         
        if name is None:
            name = "GPrim_" + str(len(self.graphicPrim)+1)
              
        gPrim = graphicPrimitive(self)
         
        if not isinstance(name, bytes) or isinstance(name, str):
            name = str(name)
     
        if type.lower() == "point":
            gPrim.createPoint(items, _byteify(name))
        elif type.lower() == "line":
            gPrim.createLine(items, _byteify(name), lines)
        elif type.lower() == "triangle":
            gPrim.createTriangle(items, _byteify(name), lines)
        else:
            raise TypeError
        
        self.graphicPrim.append(gPrim)
  
        return gPrim
    
    def sendInit(self, client):
        self.sendPrimative(client, 1)
        self.sendPrimative(client, -1)
        
    def sendUpdate(self, client):
        self.sendPrimative(client, 0)

    def sendPrimative(self, client, flag): 
        
        cdef :
            cWV.wvGPrim *gp
            int i4
            unsigned short *s2 = <unsigned short *> &i4
            unsigned char *c1 = <unsigned char *> &i4
        
        # TODO: Check client type for a send function
        
        self.bufferIndex = 0 
        
        print("Flag - ", flag)
        if flag == 1: # Intiatial message
            self.__sendInitMessage(client)
            return 0
        
        self.__sendThumbnail(client)
        
        self.__sendKey(client, flag)
        
        if self.context.gPrims == NULL: 
            print("No GPrims")
            return -1
        
        # Any changes
        if (flag == 0) and (self.context.cleanAll == 0):
            
            for i in range(self.context.nGPrim):
                if self.context.gPrims[i].updateFlg != 0:
                    break

            if i == self.context.nGPrim:
                return 0
  
        # Put out the new data
        self.bufferIndex = 0
 
        if self.context.cleanAll != 0:
        
            npack = 8
            self.__sendBuffer(client, npack)
            
            struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+3, 2) # delete opcode for all
     
            struct.pack_into('B', self.buffer, self.bufferIndex+4, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+5, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+6, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+7, 0) 
            
            self.bufferIndex += npack
            self.context.cleanAll = 0
        
        for i in range(self.context.nGPrim):
            print("Primative - ", i)
            
            gp = &self.context.gPrims[i]
            
            if (gp.updateFlg == 0) and (flag != -1):
                continue  
            if (gp.updateFlg == cWV.WV_DELETE) and (flag == -1):
                continue
            
            if gp.updateFlg == cWV.WV_DELETE:
                print("Deleting all ")
                # Delete the gPrim
                npack = 8 + gp.nameLen
                self.__sendBuffer(client, npack)
                
                struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 2) # delete opcode for all
                
                self.bufferIndex += 4
                
#                 s2[0] = gp.nameLen # TODO: FIX ME 
#                 s2[1] = 0
#                 print "Addeding s"
#                 struct.pack_into('4B', self.buffer, self.bufferIndex+4, s2[0], s2[1], s2[2], s2[3])
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 0) # delete opcode for all
                self.bufferIndex += 4
                
#                 print "Done adding s"
                #struct.pack_into('i', self.buffer, self.bufferIndex+4, gp.nameLen)
                
                self.__addStringToBuffer(gp.name, gp.nameLen)
                
                gp.updateFlg |= cWV.WV_DONE
                
            elif (gp.updateFlg == cWV.WV_PCOLOR) or (flag == -1):
                #print "New GPrim"
                print("New Primative (# " , i, ") - ", gp.name)
                
                # New gPrim
                npack = 8 + gp.nameLen + 16
                if gp.gtype > 0:
                    npack += 16
                if gp.gtype > 1: 
                    npack += 36
                
                self.__sendBuffer(client, npack)    
                 
#                 i4    = gp.nStripe
#                 c1[3] = 1 # New opcode
#                 print "Adding N-Stripe",  gp.nStripe, self.bufferIndex
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nStripe)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 1) # New opcode
                
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nStripe)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 1) # New opcode
#                 
                self.bufferIndex += 4
#                 i4    = gp.nameLen
#                 c1[2] = 0
#                 c1[3] = gp.gtype
                 
#                 print "Adding name len", gp.nameLen, gp.gtype, self.bufferIndex
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype)
                
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype)                  
#                 
                self.bufferIndex += 4
                
                self.__addStringToBuffer(gp.name, gp.nameLen)
                
                struct.pack_into('!i', self.buffer, self.bufferIndex, gp.attrs)
                self.bufferIndex +=4
                
                struct.pack_into('f', self.buffer, self.bufferIndex, gp.pSize)
                self.bufferIndex +=4
                
                struct.pack_into('3f', self.buffer, self.bufferIndex, gp.pColor[0], gp.pColor[1], gp.pColor[2])
                self.bufferIndex +=12
            
                if gp.gtype > 0:
                    struct.pack_into('f', self.buffer, self.bufferIndex+0, gp.lWidth)
                    struct.pack_into('3f', self.buffer, self.bufferIndex+4 , gp.lColor[0], gp.lColor[1], gp.lColor[2])
                    struct.pack_into('3f', self.buffer, self.bufferIndex+16, gp.fColor[0], gp.fColor[1], gp.fColor[2])
                    struct.pack_into('3f', self.buffer, self.bufferIndex+28, gp.bColor[0], gp.bColor[1], gp.bColor[2])
                    self.bufferIndex += 40
            
                if gp.gtype > 1:
                    struct.pack_into('3f', self.buffer, self.bufferIndex, gp.normal[0], gp.normal[1], gp.normal[2])
                    self.bufferIndex += 12
                
                self.__writeGPrim(client, gp) 
            else:
                print("Update GPrim")
                self.__updateGPrim(client, gp)
                
                
        # End of Franme 
        struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
        struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
        struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
        struct.pack_into('B', self.buffer, self.bufferIndex+3, 7)
        
        self.bufferIndex += 4
        
        client.sendMessage(self.buffer[:self.bufferIndex])
        self.bufferIndex = 0
        
        return
        
    
    cdef __writeGPrim(self, client, cWV.wvGPrim *gp):
        
#         print "Write Gprim"
            
        i = -1
        while i+1 < gp.nStripe:
        #for i in range(gp.nStripe):
            i += 1
        
            npack = 12+gp.nameLen
            vflag = cWV.WV_VERTICES
            
            if (gp.stripes[i].nsVerts  == 0) or \
               (gp.stripes[i].vertices == NULL):
                continue
            
            npack += 3*4*gp.stripes[i].nsVerts
            
            if (gp.stripes[i].nsIndices != 0) and \
               (gp.stripes[i].sIndice2 != NULL):
                
                npack += 2*gp.stripes[i].nsIndices + 4
                
                if (gp.stripes[i].nsIndices%2) != 0:
                    npack += 2
                    
                vflag |= cWV.WV_INDICES
        
            if gp.stripes[i].colors != NULL:
                
                npack += 3*gp.stripes[i].nsVerts + 4
                
                if ((3*gp.stripes[i].nsVerts)%4) != 0:
                    npack += 4 - (3*gp.stripes[i].nsVerts)%4
                    
                vflag |= cWV.WV_COLORS
       
            if gp.stripes[i].normals != NULL:
                npack += 3*4*gp.stripes[i].nsVerts + 4
                vflag |= cWV.WV_NORMALS
    
            if (gp.gtype == cWV.WV_LINE) and (gp.normals != NULL) and (i == 0):
                npack += 3*4*gp.nlIndex + 4
                vflag |= cWV.WV_NORMALS
        
            self.__sendBuffer(client, npack)
            
            if npack > self.bufferLen:
                print(" Oops! npack = %d  BUFLEN = %d\n", npack, self.bufferLen)
                raise IndexError
          
            
            # New data opcode
            struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
            struct.pack_into('B', self.buffer, self.bufferIndex+3, 3)
                
#             struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#             struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#             struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#             struct.pack_into('B', self.buffer, self.bufferIndex+3, 3)
            
            self.bufferIndex += 4
        
            #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
            struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
            #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+2, vflag)
            struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype)
            
            self.bufferIndex += 4
            
            self.__addStringToBuffer(gp.name, gp.nameLen)
            
            struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.stripes[i].nsVerts)
            self.bufferIndex += 4
            
            for j in range(3*gp.stripes[i].nsVerts):
                struct.pack_into('f', self.buffer, self.bufferIndex, gp.stripes[i].vertices[j])
                self.bufferIndex += 4
                
            
            if (gp.stripes[i].nsIndices != 0) and (gp.stripes[i].sIndice2  != NULL):
                
                struct.pack_into('i', self.buffer, self.bufferIndex, gp.stripes[i].nsIndices)
                self.bufferIndex += 4
                
                for j in range(gp.stripes[i].nsIndices):
                    struct.pack_into('H', self.buffer, self.bufferIndex, gp.stripes[i].sIndice2[j])
                    self.bufferIndex += 2
                    
                if (gp.stripes[i].nsIndices%2) != 0:
                  struct.pack_into('H', self.buffer, self.bufferIndex, 0)
                  self.bufferIndex += 2
     
            if gp.stripes[i].colors != NULL:
          
                struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.stripes[i].nsVerts)
                self.bufferIndex += 4
                
                for j in range(3*gp.stripes[i].nsVerts):
                    # This may need to be 'c' 
                    struct.pack_into('B', self.buffer, self.bufferIndex, gp.stripes[i].colors[j])
                    self.bufferIndex += 1
                    
                if ((3*gp.stripes[i].nsVerts)%4) != 0:
                    
                    for j in range(4 - (3*gp.stripes[i].nsVerts)%4):
                        struct.pack_into('B', self.buffer, self.bufferIndex, 0)
                        self.bufferIndex += 1
                    
           
            if gp.stripes[i].normals != NULL:
                
                struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.stripes[i].nsVerts)
                self.bufferIndex += 4
                
                for j in range(3*gp.stripes[i].nsVerts):
                    struct.pack_into('f', self.buffer, self.bufferIndex, gp.stripes[i].normals[j])
                    self.bufferIndex += 4
                
            # Line decorations -- no stripes 
            if (gp.gtype == cWV.WV_LINE) and (gp.normals != NULL) and (i == 0): 
                struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.nlIndex)
                self.bufferIndex += 4
                
                for j in range( 3*gp.nlIndex):
                    struct.pack_into('f', self.buffer, self.bufferIndex, gp.normals[j])
                    self.bufferIndex += 4
                    
            #print "Npack - ", npack, " self.bufferIndex - ", self.bufferIndex
      
            if (gp.stripes[i].npIndices != 0) and \
               (gp.stripes[i].pIndice2  != NULL):
            
                npack = 12+gp.nameLen + 2*gp.stripes[i].npIndices
            
                if (gp.stripes[i].npIndices%2) != 0:
                    npack += 2
            
                self.__sendBuffer(client, npack)
                
                # New data opcode
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 3)
                
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 3)
            
                self.bufferIndex += 4
        
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_INDICES)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 0) # local gtype 
            
                self.bufferIndex += 4
            
                self.__addStringToBuffer(gp.name, gp.nameLen)
            
                struct.pack_into('i', self.buffer, self.bufferIndex, gp.stripes[i].npIndices)
                self.bufferIndex += 4
                
                for j in range(gp.stripes[i].npIndices): 
                    struct.pack_into('H', self.buffer, self.bufferIndex, gp.stripes[i].pIndice2[j])
                    self.bufferIndex += 2
                    
       
                if (gp.stripes[i].npIndices%2) != 0:
                    struct.pack_into('H', self.buffer, self.bufferIndex, 0)
                    self.bufferIndex += 2
                    
            if (gp.stripes[i].nlIndices != 0) and \
               (gp.stripes[i].lIndice2  != NULL):
                
                npack = 12+gp.nameLen + 2*gp.stripes[i].nlIndices
            
                if (gp.stripes[i].nlIndices%2) != 0:
                    npack += 2
                    
                self.__sendBuffer(client, npack)
                
                # New data opcode
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 3) 
#                     
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 3)
            
                self.bufferIndex += 4
        
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_INDICES)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 1) # local gtype 
            
                self.bufferIndex += 4
            
                self.__addStringToBuffer(gp.name, gp.nameLen)
                
                struct.pack_into('i', self.buffer, self.bufferIndex, gp.stripes[i].nlIndices)
                self.bufferIndex += 4
           
                for j in range(gp.stripes[i].nlIndices):
                    struct.pack_into('H', self.buffer, self.bufferIndex, gp.stripes[i].lIndice2[j])
                    self.bufferIndex += 2
                    
       
                if (gp.stripes[i].nlIndices%2) != 0:
                    struct.pack_into('H', self.buffer, self.bufferIndex, 0)
                    self.bufferIndex += 2
                    
        return
    
    # Update gPrim
    cdef __updateGPrim(self, client, cWV.wvGPrim *gp):
        
       
        if (gp.updateFlg&cWV.WV_VERTICES) != 0:
            
            for i in range(gp.nStripe):
                
                if (gp.stripes[i].nsVerts  == 0) or \
                   (gp.stripes[i].vertices == NULL):
                    continue;
                
                npack = 12 + gp.nameLen + 3*4*gp.stripes[i].nsVerts
                self.sendBuffer(client, npack)
               
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                    
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                
                self.bufferIndex += 4
                
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_VERTICES)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype) 
                
                self.bufferIndex += 4
                
                self.__addStringToBuffer(client, gp.name, gp.nameLen)
                
                struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.stripes[i].nsVerts)
                self.bufferIndex += 4
           
                for j in range( 3*gp.stripes[i].nsVerts):
                    struct.pack_into('f', self.buffer, self.bufferIndex, gp.stripes[i].vertices[j])
                    self.bufferIndex += 4
          


        if (gp.updateFlg&cWV.WV_INDICES) != 0:
            for i in range(gp.nStripe):
                
                if (gp.stripes[i].nsIndices == 0) or \
                   (gp.stripes[i].sIndice2  == NULL):
                    continue
                
                npack = 12 + gp.nameLen + 2*gp.stripes[i].nsIndices
              
                if (gp.stripes[i].nsIndices%2) != 0:
                    npack += 2;
              
                self.sendBuffer(client, npack)
               
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                    
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode

                self.bufferIndex += 4
                
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_INDICES)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype) 
                
                self.bufferIndex += 4
                
                self.__addStringToBuffer(client, gp.name, gp.nameLen)
                
                struct.pack_into('i', self.buffer, self.bufferIndex, gp.stripes[i].nsIndices)
                self.bufferIndex += 4
                    
                for j in range(gp.stripes[i].nsIndices):
                    struct.pack_into('H', self.buffer, self.bufferIndex, gp.stripes[i].sIndice2[j])
                    self.bufferIndex += 2
                    
                if (gp.stripes[i].nsIndices%2) != 0:
                  struct.pack_into('H', self.buffer, self.bufferIndex, 0)
                  self.bufferIndex += 2


        if (gp.updateFlg&cWV.WV_COLORS) != 0:
            
            for i in range(gp.nStripe):
                
                if (gp.stripes[i].nsVerts == 0) or \
                   (gp.stripes[i].colors  == NULL): 
                    
                    continue
                
                npack = 12 + gp.nameLen + 3*gp.stripes[i].nsVerts
              
                if ((3*gp.stripes[i].nsVerts)%4) != 0:
                    
                    npack += 4 - (3*gp.stripes[i].nsVerts)%4
                           
                self.sendBuffer(client, npack)
                
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                    
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode

                self.bufferIndex += 4
                
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_COLORS)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype) 
                
                self.bufferIndex += 4
                
                self.__addStringToBuffer(client, gp.name, gp.nameLen)
                
                struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.stripes[i].nsVerts)
                self.bufferIndex += 4
      
                for j in range(3*gp.stripes[i].nsVerts):
                    # This may need to be 'c' 
                    struct.pack_into('B', self.buffer, self.bufferIndex, gp.stripes[i].colors[j])
                    self.bufferIndex += 1
         
                if ((3*gp.stripes[i].nsVerts)%4) != 0:
                    
                    for j in range(4 - (3*gp.stripes[i].nsVerts)%4):
                        struct.pack_into('B', self.buffer, self.bufferIndex, 0)
                        self.bufferIndex += 1
            
        if (gp.updateFlg&cWV.WV_NORMALS) != 0:
            
            if (gp.gtype == cWV.WV_TRIANGLE):
                
                for i in range(gp.nStripe):
                    
                    if (gp.stripes[i].nsVerts  == 0) or \
                       (gp.stripes[i].vertices == NULL):
                        continue
                    
                    npack = 12 + gp.nameLen + 3*4*gp.stripes[i].nsVerts
                    
                    self.sendBuffer(client, npack)
                    
                    struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                    struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                    
#                     struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                     struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                     struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                     struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode

                    self.bufferIndex += 4
                
                    #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                    struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                    #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                    struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_NORMALS)
                    struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype) 
                
                    self.bufferIndex += 4
                
                    self.__addStringToBuffer(client, gp.name, gp.nameLen)
                
                    struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.stripes[i].nsVerts)
                    self.bufferIndex += 4
           
                    for j in range( 3*gp.stripes[i].nsVerts):
                        struct.pack_into('f', self.buffer, self.bufferIndex, gp.stripes[i].normals[j])
                        self.bufferIndex += 4
                
            elif gp.gtype == cWV.WV_LINE: # line decorations 
                npack = 12 + gp.nameLen + 3*4*gp.nlIndex
                
                self.sendBuffer(client, npack)
                
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                    
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode

                self.bufferIndex += 4
                
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_NORMALS)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, gp.gtype) 
                
                self.bufferIndex += 4
                    
                self.__addStringToBuffer(client, gp.name, gp.nameLen)
       
                struct.pack_into('i', self.buffer, self.bufferIndex, 3*gp.nlIndex)
                self.bufferIndex += 4
           
                for j in range( 3*gp.nlIndex):
                    struct.pack_into('f', self.buffer, self.bufferIndex, gp.normals[j])
                    self.bufferIndex += 4
       
        if (gp.updateFlg&cWV.WV_PINDICES) != 0:
            
            for i in range(gp.nStripe):
                  
                if (gp.stripes[i].npIndices == 0) or \
                   (gp.stripes[i].pIndice2  == NULL):
                    continue
                
                npack = 12 + gp.nameLen + 2*gp.stripes[i].npIndices
                if (gp.stripes[i].npIndices%2) != 0:
                    npack += 2
                
                self.sendBuffer(client, npack)
                
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode

                self.bufferIndex += 4
                
                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_INDICES)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 0) 
                
                self.bufferIndex += 4
                    
                self.__addStringToBuffer(client, gp.name, gp.nameLen)

                struct.pack_into('i', self.buffer, self.bufferIndex, gp.stripes[i].npIndices)
                self.bufferIndex += 4
                
                for j in range(gp.stripes[i].npIndices): 
                    struct.pack_into('H', self.buffer, self.bufferIndex, gp.stripes[i].pIndice2[j])
                    self.bufferIndex += 2
                    
       
                if (gp.stripes[i].npIndices%2) != 0:
                    struct.pack_into('H', self.buffer, self.bufferIndex, 0)
                    self.bufferIndex += 2
       

        if (gp.updateFlg&cWV.WV_LINDICES) != 0:
            
            for i in range(gp.nStripe):
                
                if (gp.stripes[i].nlIndices == 0) or \
                   (gp.stripes[i].lIndice2  == NULL):
                    continue
                
                npack = 12 + gp.nameLen + 2*gp.stripes[i].nlIndices
                if (gp.stripes[i].nlIndices%2) != 0:
                    npack += 2
                
                self.sendBuffer(client, npack)
                
                struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode
                
#                 struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#                 struct.pack_into('B', self.buffer, self.bufferIndex+3, 4) # Edit opcode

                self.bufferIndex += 4

                #struct.pack_into('B', self.buffer, self.bufferIndex+0, gp.nameLen)
                struct.pack_into('i', self.buffer, self.bufferIndex+0, gp.nameLen)
                #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
                struct.pack_into('B', self.buffer, self.bufferIndex+2, cWV.WV_INDICES)
                struct.pack_into('B', self.buffer, self.bufferIndex+3, 1) 
                
                self.bufferIndex += 4

                self.__addStringToBuffer(client, gp.name, gp.nameLen)
                
                struct.pack_into('i', self.buffer, self.bufferIndex, gp.stripes[i].nlIndices)
                self.bufferIndex += 4

                for j in range(gp.stripes[i].nlIndices): 
                    struct.pack_into('H', self.buffer, self.bufferIndex, gp.stripes[i].lIndice2[j])
                    self.bufferIndex += 2
                
                if (gp.stripes[i].nlIndices%2) != 0:
                    struct.pack_into('H', self.buffer, self.bufferIndex, 0)
                    self.bufferIndex += 2

    def __addStringToBuffer(self, string, length):
        count= 0
        for i in string:
            #print "Character - ", i
            struct.pack_into('c', self.buffer, self.bufferIndex, i)
            self.bufferIndex +=1
            count += 1

        if count < length:
            for i in range(length-count):
                struct.pack_into('c', self.buffer, self.bufferIndex, " ")
                self.bufferIndex +=1
        
        elif count > length:
            print("Our string got longer!")
        
        return 
        
    # Send inital gPrim message (i.e. FoV, eye, center, etc.)
    def __sendInitMessage(self, client): 
        
        print("Initial message")
        
        struct.pack_into('B', self.buffer, 0, 0)
        struct.pack_into('B', self.buffer, 1, 0)
        struct.pack_into('B', self.buffer, 2, 0)
        struct.pack_into('B', self.buffer, 3, 8)
        
        struct.pack_into('f', self.buffer, 4, self.context.fov)
        struct.pack_into('f', self.buffer, 8, self.context.zNear)
        struct.pack_into('f', self.buffer, 12, self.context.zFar)
        struct.pack_into('3f', self.buffer, 16, self.context.eye[0], self.context.eye[1], self.context.eye[2])
        struct.pack_into('3f', self.buffer, 28, self.context.center[0], self.context.center[1], self.context.center[2])
        struct.pack_into('3f', self.buffer, 40, self.context.up[0], self.context.up[1], self.context.up[2])
        
        struct.pack_into('B', self.buffer, 52, 0)
        struct.pack_into('B', self.buffer, 53, 0)
        struct.pack_into('B', self.buffer, 54, 0)
        struct.pack_into('B', self.buffer, 55, 7)
        
        client.sendMessage(self.buffer[:56])
            
        print("Done - sending initial message")
    
    # Send thumbnail data 
    def __sendThumbnail(self, client):
            
        if self.context.thumbnail:
            
            print("Thumbnail data")
            #Init + width + height + 4*width*height + close

            self.bufferIndex = 0 
            
            struct.pack_into('B', self.buffer, 0, 0)
            struct.pack_into('B', self.buffer, 1, 0)
            struct.pack_into('B', self.buffer, 2, 0)
            struct.pack_into('B', self.buffer, 3, 10)
             
            struct.pack_into('i', self.buffer, 4, self.context.tnWidth)
            struct.pack_into('i', self.buffer, 8, self.context.tnHeight)
            
            self.bufferIndex = 12
            for i in range(self.context.tnWidth*self.context.tnHeight):
                
                struct.pack_into('4c', self.buffer, self.bufferIndex, self.context.thumbnail[4*i+0], 
                                                                      self.context.thumbnail[4*i+1], 
                                                                      self.context.thumbnail[4*i+2], 
                                                                      self.context.thumbnail[4*i+3])
                self.bufferIndex +=4

            struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+3, 7)
            
            self.bufferIndex += 4
            client.sendMessage(self.buffer[:self.bufferIndex])
            
            self.bufferIndex = 0
            
            print("Done - sending thumbnail data")
    
    # Send the key if there is one or if it has changed
    def __sendKey(self, client, flag):
        
        if (self.context.titleLen == -1) or \
           ((flag == -1) and (self.context.nColor != 0)) or \
           (self.context.nColor > 0): 

            print("Key data")
                     
            self.bufferIndex = 0
            i = self.context.titleLen
            
            if (i == -1):
                i = 0
            
            struct.pack_into('i', self.buffer, self.bufferIndex+0, i)
            #struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+3, 9) 
  
#             struct.pack_into('B', self.buffer, self.bufferIndex+0, i)
#             struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
#             struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
#             struct.pack_into('B', self.buffer, self.bufferIndex+3, 9)

            self.bufferIndex += 4
            
            if (self.context.titleLen == -1) or \
               (self.context.nColor == 0) or \
               (self.context.colors == NULL):
                
                struct.pack_into('i', self.buffer, self.bufferIndex, 0)
                self.bufferIndex += 4
            else:
                
                i = self.context.nColor
                if i < 0:
                    i = -i
                
                struct.pack_into('i', self.buffer, self.bufferIndex, i)
                self.bufferIndex += 4
                
                self.__addStringToBuffer(self.context.title, self.context.titleLen)
                
                struct.pack_into('f', self.buffer, self.bufferIndex, self.context.lims[0])
                self.bufferIndex += 4
               
                struct.pack_into('f', self.buffer, self.bufferIndex, self.context.lims[1])
                self.bufferIndex += 4
               
                for i in range(3*self.context.nColor):
                    struct.pack_into('f', self.buffer, self.bufferIndex, self.context.colors[i])
                    self.bufferIndex += 4
                   
       
            # End of frame marker
            struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
            struct.pack_into('B', self.buffer, self.bufferIndex+3, 7)
                
            self.bufferIndex+4
            
            client.sendMessage(self.buffer[:self.bufferIndex])
            
            self.bufferIndex = 0
            print("Done - sending key data")
            
    # Send delete all code
    def __sendDeleteAll(self):
        return
#         if self.context.cleanAll != 0:
#             npack = 8
#             wv_writeBuf(wsi, buf, npack, &iBuf)
#         buf[iBuf  ] = 0
#         buf[iBuf+1] = 0
#         buf[iBuf+2] = 0
#         buf[iBuf+3] = 2      #  /* delete opcode for all */
#         buf[iBuf+4] = 0
#         buf[iBuf+5] = 0
#         buf[iBuf+6] = 0
#         buf[iBuf+7] = 0
#         iBuf += npack
#          
#         self.context.cleanAll = 0
#  
    # Write buffer
    def __sendBuffer(self, client, nPack):
                
        if self.bufferIndex+nPack  <= self.bufferLen-4:
            return

        # /*
        #   printf("  iBuf = %d, BUFLEN = %d\n", *iBuf, BUFLEN)
        # */
        
        # Continue opcode
        struct.pack_into('B', self.buffer, self.bufferIndex+0, 0)
        struct.pack_into('B', self.buffer, self.bufferIndex+1, 0)
        struct.pack_into('B', self.buffer, self.bufferIndex+2, 0)
        struct.pack_into('B', self.buffer, self.bufferIndex+3, 0)
       
        self.bufferIndex += 4
        
#         print "Sending from __sendBuffer"
    
        client.sendMessage(self.buffer[:self.bufferIndex])
#         print "Done sending"
        
        self.bufferIndex = 0
        
    # Check and unify the limits for each color in the primitive 
    def checkUnifyLimits(self):
        
        maxVal = [] 
        minVal = [] 
        
        numColor = self.checkNumColor()
        print("Checking for a common max/min color limit all primitives...")
        
        for j in range(numColor):
            
            maxVal.append(-1E6)
            minVal.append(1E6)
     
            # Get the max value for all the current colors 
            for i in self.graphicPrim:
    
                item = i.colors[j]
                
                if item.minColor < minVal[j]:
                    minVal[j] = item.minColor
                    
                if item.maxColor > maxVal[j]:
                    maxVal[j] = item.maxColor
        
        for j in range(numColor):
            for i in self.graphicPrim:
                
                item = i.colors[j]
                
                item.updateColor(minVal[j], maxVal[j], item.numContour, item.reverseMap, item.colorMap)
    
    # Unify the bounding box for each primitive
    def unifyBoundingBox(self):
        
        print("Checking for a common bounding box for all primitives...")
    
        box = [1E6, -1E6, 1E6, -1E6, 1E6, -1E6]
        for i in self.graphicPrim:
            temp = i.vertex.boundingBox
            
            #print "Bounding Box = ", temp
            
            if temp[0] < box[0]:
                box[0] = temp[0]
            
            if temp[1] > box[1]:
                box[1] = temp[1]
            
            if temp[2] < box[2]:
                box[2] = temp[2]
            
            if temp[3] > box[3]:
                box[3] = temp[3]
                
            if temp[4] < box[4]:
                box[4] = temp[4]
            
            if temp[5] > box[5]:
                box[5] = temp[5]
        
        #print box
        
        focus = [0.5*(box[0]+box[1]), 0.5*(box[2]+box[3]), 0.5*(box[4]+box[5])]
        
        size = [box[1]-box[0], box[3]-box[2], box[5]-box[4]]
        focus.append(max(size))
        
        for i in self.graphicPrim:
            item = i.vertex
            item.setBoundingBox(box)
            item.adjustVertex(focus)
            
            i.updateVertex()
    
    # Check to make sure all the primitives have the same number colors 
    def checkNumColor(self):
        print("Checking to make sure all the primitives have the same number of colors...")
        
        numColor = None
        for i in self.graphicPrim:
    
            if not numColor:
                numColor = len(i.colors)
            else:
                
                if numColor != len(i.colors):
                    
                    print("The total number of colors do NOT match for each primitive")
                    print("Primitive", i.name , "has", len(i.colors), "colors, while", numColor, \
                           "colors were found for another primitive!")
                    
                    self.status = cCAPS.CAPS_MISMATCH
                    raise CAPSError(self.status, msg="while checking the number of colors for each primitive")
        
        
        return numColor
    
    def startServer(self, enableCheck = True, unifyBoundingBox = True):
        
        if enableCheck:
            # Check limits for primitives
            self.checkUnifyLimits()
        
        if unifyBoundingBox:
            # Check bouding box and adjust focus 
            self.unifyBoundingBox()
        
        # dissable the viewer (used for testing)
        if os.getenv("CAPS_BATCH"):
            cWV.wv_destroyContext(&self.context)
            self.context = NULL
            return
        
        # Start server   
        self.status = cWV.wv_startServer(<int> self.portNumber, NULL, NULL, NULL, 0, self.context)
        if self.status < 0 :
            raise CAPSError(self.status-1000, msg="while starting server, code error "+str(self.status))
        
        self.serverInstance = self.status
        
        CAPS_START = os.getenv("CAPS_START")
        
        initial = True
        while cWV.wv_statusServer(0):

            time.sleep(0.5)
            
            if initial:
                initial = False

                # Open browser
                # TODO: Ryan the controller is on capsViewer, not _wvContext...
                #if self.controller:
                #    self.controller.open(self.viewerHTML)
                if CAPS_START:
                    os.system(CAPS_START)

        #cWV.wv_removeAll(self.context)

        time.sleep(0.5)
        
        cWV.wv_cleanupServers()
    
    # Remove all graphic primitives from context 
    def removeAll(self):
        cWV.wv_removeAll(self.context)
            
        
        
## Defines a CAPS viewer object.
# A capsViewer object is a Pythonized version of Bob Haimes's "wv: A General Web-based 3D Viewer" API. wv's goal is "to
# generate a visual tool targeted for the 3D needs found within the MDAO process. A WebBrowser-based approach 
# is considered, in that it provides the broadest possible platform for deployment."
cdef class capsViewer:
    ## \example webviewer.py 
    # Basic example use case for creating visualizing data using CAPS's web viewer \ref pyCAPS.capsViewer
    
    cdef:
     
        int status 
        
        object oneBias
        
        readonly object context
       
        object controller
        
        object viewerHTML
        
        # List of items that have been added 
        object items # Must store all items
        object currentItems
        
    ## Initialize the viewer object. 
    # See \ref webviewer.py for an example use case.
    # \param browserName Name of browser to load (default - None). If left as None the system level default browser will 
    # be used. 
    #
    # \param portNumber Port number to start the server listening on (default - 7681).
    #
    # \param oneBias Flag to indicate the index biasing of the data. 
    # Options: 1 bias (default - True) or 0 bias (False) indexed.
    #
    # \param html Specify an alternative HTML file to launch when starting the server, 
    # see \ref startServer (default - None or env $CAPS_START)
    def __cinit__(self, browserName = None, portNumber=7681, oneBias=True, 
                  html=None):
                    
        #if html == "file:/$ESP_ROOT/lib/capsViewer.html":
        #    try:
        #        html = os.path.join("file:", os.environ["ESP_ROOT"])
        #        html = os.path.join(html, "lib", "capsViewer.html")
        #    except:
        #        
        #        print "Unable to determine $ESP_ROOT environment variable. Browser will NOT start automatically!"
        #        html = None
                
        self.viewerHTML = html
        
        if self.viewerHTML:
            if not os.path.isfile(self.viewerHTML.replace("file:", "")):
                print("Unable to locate file -", self.viewerHTML, ". Browser will NOT start automatically!")
                self.viewerHTML = None
        
        print("Starting context")
        self.context = _wvContext(portNumber=portNumber, oneBias = oneBias)
    
        self.items = []
        self.currentItems = []

        self.oneBias = oneBias 

        self.controller = None 
        
        if self.viewerHTML:
           
            try: 
                self.controller = webbrowser.get(browserName)
            
                if not self.controller:
                    self.status = cCAPS.CAPS_NOTFOUND
                    raise CAPSError(self.status, msg="while trying to find browser - " + browserName)
            except webbrowser.Error:
               self.viewerHTML = None
            
        
        # May need to remove this/ change in the future 
        clearGlobal()
        
    
    @property
    def __context(self):
        return <object> self.context
    
    def __dealloc__(self):
        
        self.clearItems()
        
        if self.items:
            self.items = []
            
        clearGlobal()

        self.context = None
        
        self.controller = None        
        
    def __del__(self):
        self.clearItems()
        
        if self.items:
            self.items = []
        
        self.context = None
        
        self.controller = None
            
    ## Start the server. The port number specified when initializing the viewer object will be 
    # used (see \ref capsViewer.__init__ ). 
    # Note that the server will continue to run as long as the browser is still connected, once the connection is broken 
    # the capViewer object should be deleted!
    # See \ref webviewer.py for an example use case.
    #
    # \param enableCheck Enable checks to ensure all primitives have the same number colors and that the limits
    # for each color are initially the same. 
    # \param unifyBoundingBox Unify the bounding boxes of all primitives (so all images are within the viewering area)
    def startServer(self, enableCheck = True, unifyBoundingBox = True):
        
        self.context.startServer(enableCheck=enableCheck, unifyBoundingBox=unifyBoundingBox)
    
    # Adjust the focus for all items.
    # \param unifyFocus Unify the focus for all items provided 
    # \param focus =  If a focal 
    #def adjustFocus(self, unifyFocus = True, focus=None):
    #    
    #    newFocus = [0, 0, 0, 0] 
    #    if unifyFocus:
# 
#             for i in self.items:
#                 
#                 if focus:
#                     newFocus = focus
#                 else:
#                     newFocus[0] = max(newFocus[0], i.focus[0])
#                     newFocus[1] = max(newFocus[1], i.focus[1])
#                     newFocus[2] = max(newFocus[2], i.focus[2])
#                     newFocus[3] = max(newFocus[3], i.focus[3])
#                
#                 i.setFocus(newFocus)
#                 i.adjustVertex()
#         else:
#             for i in self.items:
#                 
#                 if focus:
#                     i.setFocus(focus)
#              
#                 i.adjustVertex()
                
    
    def _determineBoundingBox(self, xyz):
        xyzMin = [ 99999999.0,  99999999.0,  99999999.0]
        xyzMax = [-99999999.0, -99999999.0, -99999999.0]
        if not isinstance(xyz, list):
            print("xyz should be a list")
            raise TypeError 
        
        if not isinstance(xyz[0], list):
            xyz = [xyz[:]]
            
        i = 0
        while i < len(xyz):
            temp = xyz[i]
            
            if temp[0] < xyzMin[0]:
                xyzMin[0] = temp[0]
            
            if temp[0] > xyzMax[0]:
                xyzMax[0] = temp[0]
            
            if temp[1] < xyzMin[1]:
                xyzMin[1] = temp[1]
                
            if temp[1] > xyzMax[1]:
                xyzMax[1] = temp[1]
                
            if temp[2] < xyzMin[2]:
                xyzMin[2] = temp[2]
            
            if temp[2] > xyzMax[2]:
                xyzMax[2] = temp[2]
            
            i += 1 
                    
        return [xyzMin[0], xyzMax[0], xyzMin[1], xyzMax[1], xyzMin[2], xyzMax[2]]
        
    ## Add a vertex set.
    # See \ref webviewer.py for an example use case.
    # \param xyz A list of lists of x,y, z values (e.g. [ [x1, y1, z1], [x2, y2, z2],[x3, y3, z3], etc. ] )
    # \param name Name of vertex set (default - None).
    # 
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored. 
    def addVertex(self, xyz, name=None):

        if not isinstance(xyz, list):
            print("xyz should be a list")
            raise TypeError 
        
         # Make sure a single vertex occurs twice due to a know bug in wvCleint (known to Bob anyways...)
        if not isinstance(xyz[0], list):
            xyz = [xyz[:], xyz[:]]
            
        if len(xyz) == 1:
            xyz = [xyz[0], xyz[0]]

        if name:
            self.items.append(wvData().createVertex(xyz, name = _byteify(str(name))))
        else:
            self.items.append(wvData().createVertex(xyz, name = NULL))
        self.currentItems.append(self.items[-1])
        
        self.items[-1].setBoundingBox(self._determineBoundingBox(xyz))
        
        return self.items[-1]
    
    ## Add a element connectivity set. 
    # See \ref webviewer.py for an example use case.
    # \param connectivity A list of lists of connectivity information 
    # (e.g. [ [node1, node2, node3], [node2, node3, node7, node8], etc. ]). Elements that >3 edges will internally
    # be decomposed into triangles. Note that whether or not the connectivity information is one biased was specified 
    # during capsViewer initialization (see \ref capsViewer.__init__ ), no further checks are currently made.
    # \param name Name of connectivity set (default - None). 
    #
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored. 
    def addIndex(self, connectivity, name=None):
        
        return self.__addIndex(None, connectivity, name)
    
    ## Add a element connectivity set with the intention of creating graphical lines on triangular primitives.
    # See \ref webviewer.py for an example use case.
    # \param connectivity A list of lists of connectivity information 
    # (e.g. [ [node1, node2, node3], [node2, node3, node7, node8], etc. ]). 
    # Note that whether or not the connectivity information is one biased was specified 
    # during capsViewer initialization (see \ref capsViewer.__init__ ), no further checks are currently made.
    # \param name Name of line set (default - None). 
    #
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored. 
    def addLineIndex(self, connectivity, name = None):
        
        return self.__addIndex("line", connectivity, name)

    ## Add a point connectivity set with the intention of creating graphical points on line or triangular primitives.
    # \param connectivity A list of connectivity information (e.g. [ node1, node2, node3, etc. ]). 
    # Note that whether or not the connectivity information is one biased was specified 
    # during capsViewer initialization (see \ref capsViewer.__init__ ), no further checks are currently made.
    # \param name Name of point set (default - None). 
    #
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored.     
    def addPointIndex(self, connectivity, name = None):
        
        return self.__addIndex("point", connectivity, name)
    
    # Sorting index function
    cdef __addIndex(self, type, connectivity, name):
    
        if type == "point":
            if not isinstance(connectivity, list):
                connectivity = [connectivity]
        else:
            if not isinstance(connectivity, list):
                print("connectivity should be a list")
                raise TypeError 
            
            if not isinstance(connectivity[0], list):
                connectivity = [connectivity[:]]
        
        if type == "point":
            if name:
                self.items.append(wvData().createPointIndex(connectivity, name = _byteify(str(name))))
            else:
                self.items.append(wvData().createPointIndex(connectivity, name=NULL))
        elif type == "line":
            if name:
                self.items.append(wvData().createLineIndex(connectivity, name = _byteify(str(name))))
            else:
                self.items.append(wvData().createLineIndex(connectivity, name=NULL))
        else:
            if name:
                self.items.append(wvData().createIndex(connectivity, name = _byteify(str(name))))
            else:
                self.items.append(wvData().createIndex(connectivity, name=NULL))

        self.currentItems.append(self.items[-1])
        
        return self.items[-1]

    ## Add color/scalar data set. 
    # See \ref webviewer.py for an example use case. 
    # \param colorData A list of scalar data values or Hexedecimal strings to be applied at each node
    # (e.g. [ node1_Color, node2_Color, node3_Color, etc.]). Alternatively,
    # a single Hexadecimal string (e.g. '#ff0000' == the color red) may be used to set all nodes to a constant color.
    #
    # \param minColor Minimum color value to use for scaling (default - None). If None, the minimum value is determined
    # automatically from the colorData provided. When plotting multiple color/data sets a minimum value should be provided
    # to ensure that all data is scaled the same. 
    #
    # \param maxColor Maximum color value to use for scaling (default - None). If None, the maximum value is determined
    # automatically from the colorData provided. When plotting multiple color/data sets a maximum value should be provided
    # to ensure that all data is scaled the same. 
    #
    # \param numContour Number of contour levels to use in the color map (default - 0 for a continuous color map).
    #
    # \param reverseMap Reverse or invert the color map (default - False).
    #
    # \param colorMap  Name of color map to use of visualization. 
    # Options: "blues" (default), "reds", "greys" (or "grays"), "blue_red" (or "bwr"), "jet", "rainbow", "hot", "cool".  
    #
    # \param name Name of color set (default - None). 
    #
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored. 
    def addColor(self, colorData, minColor=None, maxColor=None, numContour = 0, reverseMap = False, colorMap ="blues", name = None):
        
        return self.__addColor(None, colorData, minColor, maxColor, numContour, reverseMap, colorMap, name)   
    
    ## Add color/scalar data set for the graphical lines on triangular primitives. (\ref addLineIndex ). 
    # See \ref webviewer.py for an example use case. 
    # \param colorData A list of scalar data values or Hexedecimal strings to be applied at each node
    # (e.g. [ node1_Color, node2_Color, node3_Color, etc.]). Alternatively,
    # a single Hexadecimal string (e.g. '#ff0000' == the color red) may be used to set all nodes to a constant color.
    #
    # \param minColor Minimum color value to use for scaling (default - None). If None, the minimum value is determined
    # automatically from the colorData provided. When plotting multiple color/data sets a minimum value should be provided
    # to ensure that all data is scaled the same. 
    #
    # \param maxColor Maximum color value to use for scaling (default - None). If None, the maximum value is determined
    # automatically from the colorData provided. When plotting multiple color/data sets a maximum value should be provided
    # to ensure that all data is scaled the same. 
    #
    # \param numContour Number of contour levels to use in the color map (default - 0 for a continuous color map).
    #
    # \param reverseMap Reverse or invert the color map (default - False).
    # 
    # \param colorMap  Name of color map to use of visualization. 
    # Options: "blues" (default), "reds", "greys" (or "grays"), "blue_red" (or "bwr"), "jet", "rainbow", "hot", "cool".  
    #
    # \param name Name of color set (default - None). 
    #
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored. 
    def addLineColor(self, colorData, minColor=None, maxColor=None, numContour = 0, reverseMap = False, colorMap ="blues", name=None):
        
        return self.__addColor("line", colorData, minColor, maxColor, numContour, reverseMap, colorMap, name)   
    
    ## Add color/scalar data set for the graphical points on line or triangular primitives. (\ref addPointIndex ). 
    # See \ref webviewer.py for an example use case. 
    # \param colorData A list of scalar data values or Hexedecimal strings to be applied at each node
    # (e.g. [ node1_Color, node2_Color, node3_Color, etc.]). Alternatively,
    # a single Hexadecimal string (e.g. '#ff0000' == the color red) may be used to set all nodes to a constant color.
    #
    # \param minColor Minimum color value to use for scaling (default - None). If None, the minimum value is determined
    # automatically from the colorData provided. When plotting multiple color/data sets a minimum value should be provided
    # to ensure that all data is scaled the same. 
    #
    # \param maxColor Maximum color value to use for scaling (default - None). If None, the maximum value is determined
    # automatically from the colorData provided. When plotting multiple color/data sets a maximum value should be provided
    # to ensure that all data is scaled the same. 
    #
    # \param numContour Number of contour levels to use in the color map (default - 0 for a continuous color map).
    #
    # \param reverseMap Reverse or invert the color map (default - False).
    # 
    # \param colorMap  Name of color map to use of visualization. 
    # Options: "blues" (default), "reds", "greys" (or "grays"), "blue_red" (or "bwr"), "jet", "rainbow", "hot", "cool".  
    #
    # \param name Name of color set (default - None). 
    #
    # \return Optionally returns a reference to the "wvData item (object)" in which the data is stored. 
    def addPointColor(self, colorData, minColor=None, maxColor=None, numContour = 0, reverseMap = False, colorMap ="blues", name=None):
        
        return self.__addColor("point", colorData, minColor, maxColor, numContour, reverseMap, colorMap, name)
    
    # Sorting color function
    cdef __addColor(self, type, colorData, minColor, maxColor, numContour, reverseMap, colorMap, name): 
        
        if isinstance(colorData, str):
            # Check for valid hex string 
            match = regSearch(r'^#(?:[0-9a-fA-F]{3}){1,2}$', colorData) # This isn't fool proof

            if not match:
                raise ValueError("Invalid Hex string")
            
            minColor = 0
            maxColor = 0
                
        else:
            if not isinstance(colorData, list):
                colorData = [colorData]
            
            if isinstance(colorData[0], str): #Array of strings 
                minColor = 0
                maxColor = 0
            else:
                if not minColor:
                    minColor = min(colorData)
                if not maxColor:
                    maxColor = max(colorData)
        
        if type == "point":
            if name:
                self.items.append(wvData().createPointColor(colorData, float(minColor), float(maxColor), 
                                                            int(numContour), int(reverseMap),
                                                            colorMap = _byteify(colorMap.lower()), name = _byteify(str(name))))
            else:
                self.items.append(wvData().createPointColor(colorData, float(minColor), float(maxColor), 
                                                            int(numContour), int(reverseMap),
                                                            colorMap = _byteify(colorMap.lower()), name = NULL))
                
        elif type == "line":
            if name:
                self.items.append(wvData().createLineColor(colorData, float(minColor), float(maxColor), 
                                                        int(numContour), int(reverseMap),
                                                        colorMap = _byteify(colorMap.lower()), name = _byteify(str(name))))
            else:
                self.items.append(wvData().createLineColor(colorData, float(minColor), float(maxColor), 
                                                           int(numContour), int(reverseMap),
                                                           colorMap = _byteify(colorMap.lower()), name = NULL))
                
        else:
            if name:
                self.items.append(wvData().createColor(colorData, float(minColor), float(maxColor), 
                                                       int(numContour), int(reverseMap),
                                                       colorMap = _byteify(colorMap.lower()), name = _byteify(str(name))))
            else:
                self.items.append(wvData().createColor(colorData, float(minColor), float(maxColor), 
                                                       int(numContour), int(reverseMap),
                                                       colorMap = _byteify(colorMap.lower()), name = NULL))
        
        self.currentItems.append(self.items[-1])
        
        return self.items[-1]  
    
    ## Create a generic graphic primitive.
    # \param type Type of primitive to create. Options: "point", "line", or "triangle" (case insensitive)
    # \param items List of "wvData items (objects)" used to create the primitive (default - None). 
    # If none are provided all current items added will be used! 
    #
    # \param name Name of graphic primitive (default - GPrim_#, where # is 1 + number of primitives previously loaded). 
    # Important: If specified, the name must be unique.
    # \param lines Turn lines in the primitive on or off.
    #
    def createPrimitive(self, type, items=None, name=None, lines=True ):
        
        if items == None:
            items = self.currentItems #self.items
        
        gPrim = self.context.addPrimative(type, items, name = name, lines = lines)
        #gPrim = graphicPrimitive(self)
        
        if not isinstance(name, bytes) or isinstance(name, str):
            name = str(name)
    
        global globalInstance
        globalInstance.numPrim += 1
        
        globalInstance.gprim = <void **> realloc(globalInstance.gprim, globalInstance.numPrim*sizeof(void *))
        globalInstance.gprim[globalInstance.numPrim-1] = <void *> gPrim
        
#         return gPrim
        
    ## Create a point graphic primitive.
    # \param items List of "wvData items (objects)" used to create the point (default - None). 
    # If none are provided all current items added will be used! 
    #
    # \param name Name of graphic primitive (default - GPrim_#, where # is 1 + number of primitives previously loaded). 
    # Important: If specified, the name must be unique.
    #
    # \return Optionally returns a reference to the graphicPrimitive object. 
    def createPoint(self, items=None, name=None):
    # \return Index of graphic primitive
        
        return self.createPrimitive("point", items=items, name=name)
        
    ## Create a line graphic primitive.
    # \param items List of "wvData items (objects)" used to create the line (default - None). 
    # If none are provided all current items added will be used! 
    #
    # \param name Name of graphic primitive (default - GPrim_#, where # is 1 + number of primitives previously loaded). 
    # Important: If specified, the name must be unique.
    #
    # \return Optionally returns a reference to the graphicPrimitive object. 
    def createLine(self, items=None, name=None):
    # \return Index of graphic primitive
        
        return self.createPrimitive("line", items=items, name=name)
        
    ## Create a triangular graphic primitive.
    # See \ref webviewer.py for an example use case.
    # \param items List of "wvData items (objects)" used to create the triangle (default - None). 
    # If none are provided all current items added will be used! 
    #
    # \param name Name of graphic primitive (default - GPrim_#, where # is 1 + number of primitives previously loaded). 
    # Important: If specified, the name must be unique.
    # \param lines Turn lines in the primitive on or off.
    #
    # \return Optionally returns a reference to the graphicPrimitive object. 
    def createTriangle(self, items=None, name=None, lines=True):
   
        
        return self.createPrimitive("triangle", items=items, name=name, lines=lines)
        
    ## Clear the "current" list of "wvData items (objects)".
    def clearItems(self):
        
        if self.currentItems:
            self.currentItems = []

        #if self.items:
        #   self.items = []
    
    # Remove all graphic primitives
    def removeAll(self):
        self.context.removeAll()
         
    ## Add a _capsDataSet (see \ref _capsDataSet) object(s). A graphic primitive will be created automatically for
    # each data set. Note, however, if multiple data sets share the same vertex set (see \ref _capsVertexSet) 
    # the data in the repeated vertex sets will be appended into a single primitive.
    # 
    # \param dataSet A single or list of instances of _caspDataSet objects. 
    #
    def addDataSet(self, dataSet):
        
        if not isinstance(dataSet, list):
            dataSet = [dataSet]
        
        # Clear out list of current items
        self.clearItems()
        
        #xyzMin = [ 99.0,  99.0,  99.0]
        #xyzMax = [-99.0, -99.0, -99.0]
        
        boundVertexDict = {}
            
        colorRange = {}
        
        # Setup colorRange dictionary
        for i in dataSet:
            
            # Check to make sure data sets are actually data sets 
            if not isinstance(i, _capsDataSet): 
                print("Value is not a _capsDataSet")
                raise TypeError
            
            for k in range(i.dataRank):
                colorRange["Name="+i.dataSetName+"Rank="+str(k+1)] = {"Min" : 1E6, "Max" : -1E6}
            
        # Get max and min XYZ and max and min color
        for i in dataSet:
             
            if not isinstance(i, _capsDataSet): # This was already checked, supposedly 
                print("Value is not a _capsDataSet")
                raise TypeError
       
            #xyz = i.getDataXYZ()
            
            dataVal = i.getData()
            dataRank = i.dataRank
            
            #for j in range(len(xyz)):
            #    temp = xyz[j]
            #    if temp[0] < xyzMin[0]:
            #        xyzMin[0] = temp[0]
            #    
            #    if temp[0] > xyzMax[0]:
            #        xyzMax[0] = temp[0]
            #    
            #    if temp[1] < xyzMin[1]:
            #        xyzMin[1] = temp[1]
            #        
            #    if temp[1] > xyzMax[1]:
            #        xyzMax[1] = temp[1]
            #        
            #    if temp[2] < xyzMin[2]:
            #        xyzMin[2] = temp[2]
            #    
            #    if temp[2] > xyzMax[2]:
            #       xyzMax[2] = temp[2]
        
            if dataRank == 1:
                k = 0
                for val in dataVal:
                    if val < colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"]:
                        
                        colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"] = val 
                    
                    if val > colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"]:
                    
                        colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"] = val
                    
            else:
                
                for j in dataVal:
                    for k, val in enumerate(j):
                    
                        if val < colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"]:
                        
                            colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"] = val 
                    
                        if val > colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"]:
                    
                            colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"] = val
        
        #print xmin, ymin, zmin, xmax, ymax, zmax
        
        #focus = [0.5*(xyzMin[0]+xyzMax[0]), 0.5*(xyzMin[1]+xyzMax[1]), 0.5*(xyzMin[2]+xyzMax[2])]
        #focus.append(max(focus))
        
        # Add triangles 
        for i in dataSet:
            
            if not isinstance(i, _capsDataSet): # This was already checked, supposedly 
                print("Value is not a _capsDataSet")
                raise TypeError
       
            self.clearItems()
            
            boundVertex =  "Bound=" + i.capsBound.boundName + "Vertex=" + i.capsVertexSet.vertexSetName
            # Create the primitive, but only if it is a unique vertex set  
            if boundVertex not in boundVertexDict.keys():
            
                boundVertexDict[boundVertex] = None
                     
                self.addVertex(i.getDataXYZ())
            
                #item.adjustVertex(focus)
            
                self.addIndex(i.getDataConnect()[0])
                self.addLineIndex(i.getDataConnect()[0])
                self.addLineColor(1, colorMap = "grey")
            
                if i.dataRank == 1:
                    k = 0
                    self.addColor(i.getData(), 
                                  minColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"], 
                                  maxColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"], 
                                  name = _byteify(str(i.dataSetName)))

                else:
                    
                    dataVal = i.getData()
                
                    for k in range(i.dataRank):
                        temp = []
                        
                        # Lets create a color map for all ranks, createPrimative will only show the first one
                        for val in dataVal:  
                            temp.append(val[k])
                            
                        name = _byteify(i.dataSetName+"Rank="+str(k+1))
                        #print "Initial", name 
                        self.addColor(temp, 
                                      minColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"], 
                                      maxColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"], 
                                      name = name)
                
                k=0
                minColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"]
                maxColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"]

                name = _byteify(str(i.capsBound.boundName) + "|Range:%10.6g,%10.6g" % (minColor, maxColor))
                #print "Name = " , name
                boundVertexDict[boundVertex] = self.createTriangle(name=name)
                #self.context.setKey() gives a corrupt key message in the javascript...

            else: # Lets just create the color for the data set and append it to the primitive that already exists
                
                if i.dataRank == 1:
                    k = 0
                    item = self.addColor(i.getData(), 
                                         minColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"], 
                                         maxColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"], 
                                         name = _byteify(str(i.dataSetName)))
                    
                    boundVertexDict[boundVertex].appendColor(item)
    
                else:
                    
                    dataVal = i.getData()
                
                    for k in range(i.dataRank):
                        temp = []
                        
                        # Lets create a color map for all ranks, createPrimative will only show the first one
                        for val in dataVal:  
                            temp.append(val[k])
                            
                        name = _byteify(i.dataSetName+"Rank="+str(k+1))
                        #print "Append", name 
                        item = self.addColor(temp, 
                                             minColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Min"], 
                                             maxColor=colorRange["Name="+i.dataSetName+"Rank="+str(k+1)]["Max"], 
                                             name = name)
                            
                        boundVertexDict[boundVertex].appendColor(item)
                        
    ## Alias to \ref addDataBound
    def addBound(self, bound, dataSetType="both"):
        self.addDataBound( bound, dataSetType="both")

    ## Add a capsBound (see \ref capsBound) object(s). 
    #
    # \param bound A single or list of instances of capsBound objects. 
    # \param dataSetType Specifies which type of data sets in the bound should be added, source or destination.
    # Options: "source", "destination, or "both" (default) - values are case insensitive. 
    def addDataBound(self, bound, dataSetType="both"):  
        
        # Covert to list if not already one
        if not isinstance(bound, list):
            bound = [bound]
            
        dataSetType = dataSetType.lower()

        # Check data set type
        if dataSetType not in ["both", "source", "destination"]:
            
            print("Invalid option for dataSetType")
            raise ValueError
        
        # Initialize date set list
        dataSetList = []

        for i in bound:
            
            # Check to make sure the bound is actually a bound 
            if not isinstance(i, capsBound):
                print("Value is not a capsBound")
                raise TypeError
            
            if dataSetType == "both" or dataSetType == "source":
                
                dataSetList += i.dataSetSrc.values()
                
            if dataSetType == "both" or dataSetType == "destination":
                
                dataSetList += i.dataSetDest.values()
                
                    
        # Pass the list of data sets to addDataSet - it is supposed to do all the sorting 
        self.addDataSet(dataSetList)
    
    ## Add a Tecplot file.
    # 
    # \param filename Name of Tecplot file to load. 
    def addTecplot(self, filename):
        
        def determineXYZKey(variables, xyzName):
            
            for i in variables:
                
                if xyzName == i.lower():
                    return i
                elif xyzName + " " in i.lower():
                    return i 
                else:
                    pass
            
            print("Unable to determine which variable is x, y, or z ")
            raise ValueError
            
        tp = capsTecplot(filename)
        tp.loadZone()
        
        # Check for types 
        for i in tp.zone.values():
            
            if i.zoneType == "ORDERED":
                if i.k != 1:
                    print("3-dimensional ordered data isn't supported, zone - ", i.zoneName)
                    raise TypeError
                
            elif ( i.zoneType != "FELINESEG"  and 
                   i.zoneType != "FETRIANGLE" and 
                   i.zoneType != "FEQUADRILATERAL"):
                
                print("Unsupported zone type - ", i.zoneType)
                raise TypeError
        
        # Clear out list of current items
        self.clearItems()
        
        colorRange = {}
        
        # Setup colorRange dictionary
        for k in tp.variable:
            colorRange[k] = {"Min" : 1E6, "Max" : -1E6}
        
        xyzName = [determineXYZKey(tp.variable, "x"), 
                   determineXYZKey(tp.variable, "y"),
                   determineXYZKey(tp.variable, "z")]
        
        # Get max and min color
        for i in tp.zone.values():
            
            for j in tp.variable:
                
                # Skip x, y, z variables 
                if j in xyzName:
                    continue
                
                for val in i.data[j]:
                    if val < colorRange[j]["Min"]:
                        
                        colorRange[j]["Min"] = val 
                    
                    if val > colorRange[j]["Max"]:
                    
                        colorRange[j]["Max"] = val
            
        # Add triangles 
        for i in tp.zone.values():

            self.clearItems()

            name = _byteify(str(i.zoneName.rstrip().lstrip().replace("'", "").replace('"', "").replace(" ", "")))
            
            if i.zoneType == "ORDERED": # Ordered-zones
                
                # Copy over xyz to put in [[x1, y1, z1], [x2, y2, z2], etc... ] format
                xyz = []
                for j in range(i.i*i.j*i.k):
                    xyz.append([i.data[xyzName[0]][j], 
                                i.data[xyzName[1]][j], 
                                i.data[xyzName[2]][j]])
            
                self.addVertex(xyz)
                
                if i.i == 1 and  i.j == 1:
                    print("Single point in zone is not currently suppored!")
                    raise ValueError
                
                # Create connecitivy 
                connect = []
                
                if i.j == 1:
                    for row in range(i.i-1):
                        connect.append([row+1, row+2]);
                else:    
                    
                    node = range(1, i.i*i.j +1)      
                    for col in range(i.j-1):
                        for row in range(i.i-1):
                            connect.append( [node[ (col+0)*i.i + (row+0)],  
                                             node[ (col+0)*i.i + (row+1)], 
                                             node[ (col+1)*i.i + (row+1)],  
                                             node[ (col+1)*i.i + (row+0)]])
  
                self.addLineIndex(connect)
            
                self.addLineColor(1, colorMap = "grey")
                
                if i.j == 1:
                
                    print("Warning - color data for lines is not currently supported within the addTecplot function!")
                    self.createLine(name=name)
       
                else:
                    self.addIndex(connect)
                    
                    for j in tp.variable:
                        
                        if j in xyzName:
                            continue
                        
                        self.addColor(i.data[j], 
                                      minColor=colorRange[j]["Min"], 
                                      maxColor=colorRange[j]["Max"], 
                                      name = _byteify(str(j)))
                    
    #                 print "Name = " , name
                    self.createTriangle(name=name)
            
            else: # FE-zones
                # Copy over xyz to put in [[x1, y1, z1], [x2, y2, z2], etc... ] format
                xyz = []
                for j in range(i.numNode):
                    xyz.append([i.data[xyzName[0]][j], 
                                i.data[xyzName[1]][j], 
                                i.data[xyzName[2]][j]])
                
                self.addVertex(xyz)
            
                #item.adjustVertex(focus)
                
                self.addLineIndex(i.connectivity)
            
                self.addLineColor(1, colorMap = "grey")
                
                if i.zoneType == "FELINESEG":
                    print("Warning - color data for lines is not currently supported within the addTecplot function!")
                    self.createLine(name=name)
            
                else:
                    self.addIndex(i.connectivity)
                    
                    for j in tp.variable:
                        
                        if j in xyzName:
                            continue
                        
                        self.addColor(i.data[j], 
                                      minColor=colorRange[j]["Min"], 
                                      maxColor=colorRange[j]["Max"], 
                                      name = _byteify(str(j)))
                    
    #                 print "Name = " , name
                    self.createTriangle(name=name)
        
    ## Add a unstructured mesh file. Note that in most cases 
    # the mesh's name (see \ref capsUnstructMesh.$name) will be used for the name of the primitive.
    # 
    # \param meshFile Name of unstructured mesh file to load or an instance of a 
    # capsUnstructMesh(see \ref capsUnstructMesh) object.
    #  
    def addUnstructMesh(self, meshFile, **kwargs):
        
        if isinstance(meshFile, str):
            mesh = capsUnstructMesh(meshFile)
        elif isinstance(meshFile, capsUnstructMesh):
            mesh = meshFile 
        
        if mesh.numNode == 0 or mesh.numElement == 0:    
            mesh.loadFile(**kwargs)
        
        if mesh.oneBias != self.oneBias:
            print("Index biasing between mesh and server mismatch!")
            raise ValueError
        
        # Clear out list of current items
        self.clearItems()
        
        colorRange = {}
        
        # Setup colorRange dictionary
        for i in mesh.data.keys():
            colorRange[i] = {"Min" : 1E6, "Max" : -1E6}
        
        # Get max and min color
        for i in mesh.data.keys():
            
            for val in mesh.data[i]:
                if val < colorRange[i]["Min"]:
                        
                    colorRange[i]["Min"] = val 
                    
                if val > colorRange[i]["Max"]:
                    
                    colorRange[i]["Max"] = val
            
        # Origami mesh
        if mesh.meshType == "origami":
            
            # Add lines - create a new primitive for each line so we can individually select the color
            for index, i in enumerate(mesh.connectivity):
                
                self.clearItems()
                
                xyz = [mesh.xyz[i[0]-1], mesh.xyz[i[1]-1]]
                 
                self.addVertex(xyz)
                
                #item.adjustVertex(focus)
                
                self.addIndex([1, 2])
                
                for j in mesh.data.keys():
                    self.addLineColor(mesh.data[j][index], 
                                      minColor=colorRange[j]["Min"], 
                                      maxColor=colorRange[j]["Max"], 
                                      name = _byteify(str(j)),
                                      colorMap = "rainbow")
                    
                name = _byteify(str(mesh.name + str(index)))
                #print "Name = " , name
                self.createLine(name=name)
        
        elif mesh.meshType == "aflr3": # We can only handle surface meshes
            
            # Count how many triangle and quadrilaterals there are 
            
            count = len(mesh.data["BCID"])
            
            #print count, mesh.numElement
            
            if count != mesh.numElement:
                print("capsViewer can only visualize surface data - volume elements will be ignored!")
                
                # Determine what the last node in the connectivty should be 
                maxIndex = 0
                for i in range(count):
                    
                    temp = max(mesh.connectivity[i])
                    if  temp > maxIndex:
                        maxIndex = temp
                        
                #print "Maxindex", maxIndex
                
                # Assume first n-nodes just belong to the surface
                self.addVertex(mesh.xyz[0:maxIndex])
                
                # First e- elements are triangles and quadrilaterals
                self.addIndex(mesh.connectivity[0:count])
                self.addLineIndex(mesh.connectivity[0:count])
            
            else:
                
                self.addVertex(mesh.xyz)
             
                self.addIndex(mesh.connectivity)
                self.addLineIndex(mesh.connectivity)
            
            self.addLineColor(1, colorMap = "grey")
            self.addColor('#ffffa3') # Temporary - need to push boundary ids to nodes
            #for i in mesh.data.keys():
            #    self.addColor(mesh.data[i], 
            #                  minColor=colorRange[i]["Min"], 
            #                  maxColor=colorRange[i]["Max"], 
            #                  name = _byteify(str(i)))
                 
            name = _byteify(str(mesh.name.split(".")[0]))
            #print "Name = " , name
            self.createTriangle(name=name)
        
        elif mesh.meshType == "stl":
            
            self.addVertex(mesh.xyz)
            self.addIndex(mesh.connectivity)
            self.addLineIndex(mesh.connectivity)
            self.addLineColor(1, colorMap = "grey")
            self.addColor('#ffffa3')

            # Would need to convert element data to nodes 
#             for i in mesh.data.keys():
#                 self.addColor(mesh.data[i], 
#                               minColor=colorRange[i]["Min"], 
#                               maxColor=colorRange[i]["Max"], 
#                               name = _byteify(str(i)))
#                 
            name = _byteify(str(mesh.name.split(".")[0]))
            #print "Name = " , name
            self.createTriangle(name=name)

        # Default mesh 
        else:
            
            print("Mesh type", mesh.meshType, "isn't reconzied!")
            raise TypeError

#             self.addVertex(mesh.xyz)
#             self.addIndex(mesh.connectivity)
#             self.addLineIndex(mesh.connectivity)
#             self.addLineColor(1, colorMap = "grey")
#                 
#             for i in mesh.data.keys():
#                 self.addColor(mesh.data[i], 
#                               minColor=colorRange[i]["Min"], 
#                               maxColor=colorRange[i]["Max"], 
#                               name = _byteify(str(i)))
#                 
#             name = _byteify(str(mesh.name))
#             #print "Name = " , name
#             self.createTriangle(name=name)
    
    ## Add a geometry file. Files are loaded by extension: *.igs, *.iges, *.stp, 
    # *.step, *.brep (for native OpenCASCADE files), *.egads (for native files with persistent Attributes, splits ignored), 
    # or *.csm (ESP file)
    # 
    # \param geomFile Name of geometry file to load.
    # \param flags Options (default - None) : 
    #              1 - Don’t split closed and periodic entities
    #              2 - Split to maintain at least C1 in BSPLINEs
    #              4 - Don’t try maintaining Units on STEP read (always millimeters)
    # \param dictFile Dictionary file name when loading *.csm files, if applicable (default - None)
    def addGeometry(self, geomFile, flags=0, dictFile = None):
        
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
               
            raise CAPSError(self.status, msg=msg)
                        
        if not isinstance(geomFile, str):
            raise TypeError
        
        fileName = _byteify(geomFile)
        
        if ".csm" in geomFile: # If a ESP file - use OpenCSM loader
        
            # Set OCSMs output level
            cOCSM.ocsmSetOutLevel(<int> 1)
            self.status = cOCSM.ocsmLoad(fileName, &modl)
            if self.status < 0: 
                raiseError("while getting bodies from model (during a call to ocsmLoad)")
            
            # Cast void to model 
            MODL = <cOCSM.modl_T *> modl
            
            # If there was a dictionary, load the definitions from it now */
            if dictFile:
                dictName = _byteify(dictFile)
                self.status = cOCSM.ocsmLoadDict(modl, dictName)
                          
                if self.status < 0: 
                    raiseError("while loading dictionary (during a call to ocsmLoadDict)")
                
        
            self.status = cOCSM.ocsmBuild( modl,     # (in)  pointer to MODL
                                           <int> 0,  # (in)  last Branch to execute (or 0 for all, or -1 for no recycling)
                                           &builtTo, # (out) last Branch executed successfully
                                           &numBody, # (in)  number of entries allocated in body[], (out) number of Bodys on the stack
                                           NULL)     # (out) array  of Bodys on the stack (LIFO) (at least *nbody long)
            if self.status < 0: 
                raiseError("while building model (during a call to ocsmBuild)")
            
            numBody = 0
            for i  in range(MODL.nbody):
           
                if MODL.body[i+1].onstack != 1:
                    continue

                numBody += 1
                    
                body = <cEGADS.ego *> cEGADS.EG_reall(body, numBody*sizeof(cEGADS.ego))
                     
                if body == NULL:
                    self.status = cEGADS.EGADS_MALLOC
                    raiseError("while allocating bodies  (during a call to EG_reall)")
                     
                freeBody = True
                
                self.status = cEGADS.EG_copyObject(MODL.body[i+1].ebody, NULL, &body[numBody -1])
                if self.status != cEGADS.EGADS_SUCCESS:
                    raiseError("while getting bodies from model (during a call to EG_copyObject)") 
              
        else: # else use a EGADS model loader
        
            self.status = cEGADS.EG_open(&context)
            if self.status != cEGADS.EGADS_SUCCESS:
                raise CAPSError(self.status, msg="while opening context  (during a call to EG_open)")
            
            self.status = cEGADS.EG_loadModel(context, <int> flags, <char *>fileName, &model)
            if self.status != cEGADS.EGADS_SUCCESS: 
                raiseError("while loading model from geometry file (during a call to EG_loadModel)")
            
            self.status = cEGADS.EG_getTopology(model, &ref, &oclass, &mtype, data, &numBody, &body, &senses)
            if self.status != cEGADS.EGADS_SUCCESS or oclass != cEGADS.MODEL:
                raiseError("while getting bodies from model (during a call to EG_getTopology)")
        
        try:
            self.addEgadsBody(numBody, body)
        except CAPSError as e:
            self.status = e.errorCode
            raiseError(e.args[0])
        except:
            raise
            
        cleanup()
    
    ## Display an array of EGADS bodies 
    # \param numBody Number of bodies in array.
    # \param body Array of EGADS body objects.
    cdef addEgadsBody(self, int numBody, cEGADS.ego *body): 
        self.__loadEGADS(numBody, body, NULL) 
        
    ## Display an array of EGADS tessellations 
    # \param numTess Number of tessellations in array.
    # \param tess Array of EGADS tessellation objects.
    cdef addEgadsTess(self, int numTess, cEGADS.ego *tess): 
        self.__loadEGADS(numTess, NULL, tess, globalTess=True)
    
    # Tessellate an array of EGADS bodies and add them to the primitives or add an array
    # of tessellated bodies to primitives. One of 'bodies' or 'tessellations' should be null. 
    cdef __loadEGADS(self, int numBody, cEGADS.ego *bodies, cEGADS.ego *tessellations, globalTess=False):
        cdef:
            # EG_statusTessBody
            cEGADS.ego tempBody
            int dummyInt
            int numPoint = 0

            # EG_getTopology
            cEGADS.ego ref
            int oclass, mtype
            double[4] data
            int *senses
            int numChildren
            cEGADS.ego *children

            # EG_makeTessBody
            double box[6]
            double size = 0
            double params[3]

            cEGADS.ego tess
            cEGADS.ego body

        def raiseError(msg):
            raise CAPSError(self.status, msg=msg)
        
        nameList = []
        # Loop through bodies - tesselate and create gPrims
        for ibody in range(numBody):
            
            # Clear out list of current items
            self.clearItems()
            
            if not bodies:
                if tessellations:
                    if tessellations[ibody] == <cEGADS.ego> NULL: 
                        continue # skip bodies without tessellations
                    
                    self.status = cEGADS.EG_statusTessBody(tessellations[ibody], &body, &dummyInt, &numPoint)
                    if self.status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting the status of the tessellation (during a call to EG_statusTessBody)")
                else: 
                    self.status = cCAPS.CAPS_BADVALUE
                    raiseError("Neiteher bodies nor tessellations provided!");
            else:
                body = bodies[ibody]
                        
            if not tessellations:
                # Get tessellation
                self.status = cEGADS.EG_getBoundingBox(body, box)
                if self.status != cEGADS.EGADS_SUCCESS: 
                    raiseError("while getting bounding box of body (during a call to (during a call to EG_makeTessBody)")
                
                size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
                            (box[1]-box[4])*(box[1]-box[4]) +
                            (box[2]-box[5])*(box[2]-box[5]))
                
                params[0] = 0.0250*size
                params[1] = 0.0010*size
                params[2] = 15.0

                print("--> Tessellating Body %6d     (%12.5e %12.5e %7.3f)" % (ibody, params[0], params[1], params[2]))
                self.status = cEGADS.EG_makeTessBody(body, params, &tess)
                if self.status != cEGADS.EGADS_SUCCESS: 
                     raiseError("while making tesselation (during a call to EG_makeTessBody)")
            else: 
                tess = tessellations[ibody]

            self.status = cEGADS.EG_getTopology(body, &ref, &oclass, &mtype, data, &numChildren, &children, &senses)
            if self.status != cEGADS.EGADS_SUCCESS:
                raiseError("while getting body info (during a call to EG_getTopology)")

            if mtype == cEGADS.WIREBODY:
                nameList = self.__loadWireMesh(ibody, nameList, body, tess)
            else:
                nameList = self.__loadSurfaceMesh(ibody, nameList, body, tess, globalTess)

        return

    # Constructs surfaces meshs with quad/tri elements
    cdef __loadSurfaceMesh(self, object ibody, object nameList, cEGADS.ego body, cEGADS.ego tess, object globalTess):
        cdef:
            # EG_getBodyTopos
            int numFace
            cEGADS.ego *face = NULL

            # EG_statusTessBody
            cEGADS.ego tempBody
            int dummyInt
            int numPoint = 0
            
            # EG_getGlobal
            int pointType
            int pointIndex
            double[3] xyzs 
            
            # EG_localToGlobal
            int globalID = 0
            
            #EGADS tess return variables
            int numTri = 0
            const int *tris = NULL
            const int *triNeighbor = NULL
            const int *ptype = NULL
            const int *pindex = NULL
            const double *points = NULL
            const double *uv = NULL
            
            #Patch quadding variables
            int patch = 0, qlen = 0
            int numPatch = 0, n1, n2;
            const int *pvindex = NULL
            const int *pbounds = NULL
            
            # Attributes
            int atype
            int alen
            const int *ints
            const double *reals 
            const char *string
            
        def raiseError(msg):
            cEGADS.EG_free(face)
            raise CAPSError(self.status, msg=msg)

        xyz = []
        connectivity = []
        
        useColorArray = False
        colorArray = []
        colorDefault = '#ffffa3' # yellow - Default color for the body/face
        
        # How many faces do we have? 
        self.status = cEGADS.EG_getBodyTopos(body, 
                                             <cEGADS.ego> NULL, 
                                             cEGADS.FACE, 
                                             &numFace, 
                                             &face)
        if self.status != cEGADS.EGADS_SUCCESS: 
            raiseError("while getting number of faces (during a call to EG_getBodyTopos)")
        
        if globalTess:
            self.status = cEGADS.EG_statusTessBody(tess, &tempBody, &dummyInt, &numPoint)
            if self.status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting the status of the tessellation (during a call to EG_statusTessBody)")
            
            # Set color array
            colorArray = [colorDefault]*numPoint
            
            # Get the points 
            for j in range(int(numPoint)):
      
                self.status = cEGADS.EG_getGlobal(tess, j+1, &pointType, &pointIndex, xyzs)
                if self.status != cEGADS.EGADS_SUCCESS: 
                    raiseError("while getting coordinates for the tessellation (during a call to EG_getGlobal)")
                
                xyz.append([xyzs[0], xyzs[1], xyzs[2]])
               
        # check if the tessellation is globally quadded
        quadded = False
        attr = _byteify(".tessType")
        self.status = cEGADS.EG_attributeRet(tess, <char *>attr, &atype, &alen, &ints, &reals, &string)
        if self.status == cEGADS.EGADS_NOTFOUND:
            pass
        elif self.status != cEGADS.EGADS_SUCCESS:
            raiseError("while getting trying to retrieve .tessType (during a call to EG_attributeRet)")
        else: 
            if atype == cEGADS.ATTRSTRING and _strify(string) == "Quad":
                quadded = True;

        # Loop through face by face and get tessellation information 
        for faceIndex in range(numFace):
            
            attr = _byteify("capsIgnore")
            self.status = cEGADS.EG_attributeRet(face[faceIndex], <char *>attr, &atype, &alen, &ints, &reals, &string)
            if self.status == cEGADS.EGADS_NOTFOUND:
                pass
            elif self.status != cEGADS.EGADS_SUCCESS:
                raiseError("while getting trying to retrieve capsIgnore (during a call to EG_attributeRet)")
            else: 
                if atype == cEGADS.ATTRSTRING:
                    #print "capsIgnore attribute found for face - ", faceIndex+1, "!!"
                    continue # skip faces that are marked with capsIgnore

            # Get the color from the face
            try:
                temp = __getColor(face[faceIndex])
            except CAPSError as e: 
                self.status = e.errorCode
                raiseError("while getting color")
            except:
                raise
         
            if temp:
                color = temp
                useColorArray = True
            else:
                color = colorDefault
                
            # first check if the faces is a quad with classic patching quads
            self.status = cEGADS.EG_getQuads(tess, faceIndex+1, &qlen, &points, &uv, &ptype, &pindex, &numPatch);
            if self.status == cEGADS.EGADS_SUCCESS and numPatch != 0:
                if numPatch != 1:
                    self.status = cCAPS.CAPS_NOTIMPLEMENT
                    raiseError("EG_localToGlobal accidentally only works for a single quad patch!\n")

                patch = 1
                self.status = cEGADS.EG_getPatch(tess, faceIndex+1, patch, &n1, &n2, &pvindex, &pbounds)
                if self.status != cEGADS.EGADS_SUCCESS:
                    raiseError("while getting the quad patch (during a call to EG_getPatch)")

                if not globalTess:
                    for j in range(qlen):
                        xyz.append([points[3*j+0], points[3*j+1], points[3*j+2]])
                        colorArray.append(color)
                        
                    globalID += qlen
                    
                for j in range(1, n2):
                    for i in range(1, n1):

                        elem = []
                        
                        if globalTess:
                            self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, pvindex[(i-1)+n1*(j-1)], &globalID)
                            if self.status != cEGADS.EGADS_SUCCESS: 
                                raiseError("while getting local-to-global indexing [connectivity - 0] (during a call to EG_localToGlobal)")
                        
                            elem.append(globalID)
                            colorArray[globalID-1] = color
                         
                            self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, pvindex[(i  )+n1*(j-1)], &globalID)
                            if self.status != cEGADS.EGADS_SUCCESS: 
                                raiseError("while getting local-to-global indexing [connectivity - 1] (during a call to EG_localToGlobal)")
                        
                            elem.append(globalID)
                            colorArray[globalID-1] = color
                        
                            self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, pvindex[(i  )+n1*(j  )], &globalID)
                            if self.status != cEGADS.EGADS_SUCCESS: 
                                raiseError("while getting local-to-global indexing [connectivity - 2] (during a call to EG_localToGlobal)")                 
                            
                            elem.append(globalID)
                            colorArray[globalID-1] = color
                        
                            self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, pvindex[(i-1)+n1*(j  )], &globalID)
                            if self.status != cEGADS.EGADS_SUCCESS: 
                                raiseError("while getting local-to-global indexing [connectivity - 4] (during a call to EG_localToGlobal)")                 
                        
                            elem.append(globalID)
                            colorArray[globalID-1] = color
                        
                        else:
                            elem.append(pvindex[(i-1)+n1*(j-1)] + globalID - qlen)
                            elem.append(pvindex[(i  )+n1*(j-1)] + globalID - qlen)
                            elem.append(pvindex[(i  )+n1*(j  )] + globalID - qlen)
                            elem.append(pvindex[(i-1)+n1*(j  )] + globalID - qlen)
                          
                            colorArray[pvindex[(i-1)+n1*(j-1)] + globalID - qlen - 1] = color
                            colorArray[pvindex[(i  )+n1*(j-1)] + globalID - qlen - 1] = color
                            colorArray[pvindex[(i  )+n1*(j  )] + globalID - qlen - 1] = color
                            colorArray[pvindex[(i-1)+n1*(j  )] + globalID - qlen - 1] = color
                          
                        connectivity.append(elem)
                        
                continue # next face
            
            # Regular trinagles or quadded triangles
            self.status = cEGADS.EG_getTessFace(tess, faceIndex+1, &numPoint, &points, &uv, &ptype, 
                                                &pindex, &numTri, &tris, &triNeighbor)
            if self.status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting the face tessellation (during a call to EG_getTessFace)")

            if not globalTess:
                for j in range(numPoint):
                    xyz.append([points[3*j+0], points[3*j+1], points[3*j+2]])
                    colorArray.append(color)
                        
                globalID += numPoint
                
            # set the number of elements and stide based on complete quadding
            numElem = int(numTri/2) if quadded else numTri
            stride = 6 if quadded else 3

            for j in range(numElem):
                elem = []
                
                if globalTess:
                    self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[stride*j + 0], &globalID)
                    if self.status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting local-to-global indexing [connectivity - 0] (during a call to EG_localToGlobal)")
                
                    elem.append(globalID)
                    colorArray[globalID-1] = color
                
                    self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[stride*j + 1], &globalID)
                    if self.status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting local-to-global indexing [connectivity - 1] (during a call to EG_localToGlobal)")
                
                    elem.append(globalID)
                    colorArray[globalID-1] = color
                
                    self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[stride*j + 2], &globalID)
                    if self.status != cEGADS.EGADS_SUCCESS: 
                        raiseError("while getting local-to-global indexing [connectivity - 2] (during a call to EG_localToGlobal)")                 
                    
                    elem.append(globalID)
                    colorArray[globalID-1] = color
                
                    if quadded:
                        self.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[stride*j + 5], &globalID)
                        if self.status != cEGADS.EGADS_SUCCESS: 
                            raiseError("while getting local-to-global indexing [connectivity - 5] (during a call to EG_localToGlobal)")                 
                    
                        elem.append(globalID)
                        colorArray[globalID-1] = color
                else:
                    elem.append(tris[stride*j + 0] + globalID - numPoint)
                    elem.append(tris[stride*j + 1] + globalID - numPoint)
                    elem.append(tris[stride*j + 2] + globalID - numPoint)
                    if quadded:
                        elem.append(tris[stride*j + 5] + globalID - numPoint)
                          
                    colorArray[tris[stride*j + 0] + globalID - numPoint - 1] = color
                    colorArray[tris[stride*j + 1] + globalID - numPoint - 1] = color
                    colorArray[tris[stride*j + 2] + globalID - numPoint - 1] = color
                    
                    if quadded:
                        colorArray[tris[stride*j + 5] + globalID - numPoint - 1] = color
                    

                connectivity.append(elem)
                
        
        cEGADS.EG_free(face); face = NULL

        # Get color on body
        try: 
             color = __getColor(body) or color or '#ffffa3' # yellow - Default color for the body/face
        except CAPSError as e: 
            self.status = e.errorCode
            raiseError("while getting color")
        except:
            raise

        #print xyz
        #print connectivity
        
        if useColorArray:
            self.addColor(colorArray)
        else:
            self.addColor(color)
            
        self.addVertex(xyz)
        self.addIndex(connectivity)
        self.addLineIndex(connectivity)
        self.addLineColor(1, colorMap = "grey")
        
        name, nameList = self.__getBodyName(ibody, nameList, body)
        
#         cud = capsUnstructData()
#         cud.xyz = xyz
#         cud.connectivity = connectivity
#         if color: cud.metaData = {"color": color}
#         cud.writeFile(name)
        
        self.createTriangle(name=name,lines=globalTess)
    
        return nameList

    # Constructs a wire body mesh
    cdef __loadWireMesh(self, object ibody, object nameList, cEGADS.ego body, cEGADS.ego tess):
        cdef:
            # EG_getBodyTopos
            int numEdge
            cEGADS.ego *edge = NULL

            # EG_getTessEdge
            int plen = 0
            const double *points = NULL
            const double *t = NULL

            # EG_statusTessBody
            cEGADS.ego tempBody
            int dummyInt
            int numPoint = 0

            # EG_getGlobal
            int pointType
            int pointIndex
            double[3] xyzs 
            
            # EG_localToGlobal
            int globalID = 0

            # Attributes
            int atype
            int alen
            const int *ints
            const double *reals 
            const char *string

        def raiseError(msg):
            cEGADS.EG_free(edge)
            raise CAPSError(self.status, msg=msg)

        xyz = []
        connectivity = []

        # How many edges do we have? 
        self.status = cEGADS.EG_getBodyTopos(body, 
                                             <cEGADS.ego> NULL, 
                                             cEGADS.EDGE, 
                                             &numEdge, 
                                             &edge)
        if self.status != cEGADS.EGADS_SUCCESS: 
            raiseError("while getting number of edges (during a call to EG_getBodyTopos)")
        
        self.status = cEGADS.EG_statusTessBody(tess, &tempBody, &dummyInt, &numPoint)
        if self.status != cEGADS.EGADS_SUCCESS: 
            raiseError("while getting the status of the tessellation (during a call to EG_statusTessBody)")
        
        # Get the points 
        for j in range(numPoint):
        
            self.status = cEGADS.EG_getGlobal(tess, j+1, &pointType, &pointIndex, xyzs)
            if self.status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting coordinates for the tessellation (during a call to EG_getGlobal)")
            
            xyz.append([xyzs[0], xyzs[1], xyzs[2]])
        
        color = '#ffffa3' # yellow - Default color for the body/face
        
        # Loop through edge by edge and get tessellation information 
        for edgeIndex in range(numEdge):
            
            attr = _byteify("capsIgnore")
            self.status = cEGADS.EG_attributeRet(edge[edgeIndex], <char *>attr, &atype, &alen, &ints, &reals, &string)
            if self.status == cEGADS.EGADS_NOTFOUND:
                pass
            elif self.status != cEGADS.EGADS_SUCCESS:
                raiseError("while getting trying to retrieve capsIgnore (during a call to EG_attributeRet)")
            else: 
                if atype == cEGADS.ATTRSTRING:
                    #print "capsIgnore attribute found for edge - ", edgeIndex+1, "!!"
                    continue # skip edges that are marked with capsIgnore
            
            # get the edge tessellation
            self.status = cEGADS.EG_getTessEdge(tess, edgeIndex+1, &plen, &points, &t)
            if self.status == cEGADS.EGADS_DEGEN: continue
            if self.status != cEGADS.EGADS_SUCCESS: 
                raiseError("while getting the edge tessellation (during a call to EG_getTessEdge)")
        
            for j in range(plen-1):
                elem = []
                self.status = cEGADS.EG_localToGlobal(tess, -(edgeIndex+1), j+1, &globalID)
                if self.status != cEGADS.EGADS_SUCCESS: 
                    raiseError("while getting local-to-global edge indexing [connectivity - 0] (during a call to EG_localToGlobal)")
            
                elem.append(globalID)
            
                self.status = cEGADS.EG_localToGlobal(tess, -(edgeIndex+1), j+2, &globalID)
                if self.status != cEGADS.EGADS_SUCCESS: 
                    raiseError("while getting local-to-global edge indexing [connectivity - 1] (during a call to EG_localToGlobal)")
            
                elem.append(globalID)
            
                connectivity.append(elem)
            
            # Get the color from the edge
            try:
                temp = __getColor(edge[edgeIndex])
            except CAPSError as e: 
                self.status = e.errorCode
                raiseError("while getting color")
            except:
                raise
         
            if temp:
                color = temp
        
        cEGADS.EG_free(edge); edge = NULL
        
        # Get color
        try: 
             color = __getColor(body) or color or '#ffffa3' # yellow - Default color for the body/face
        except CAPSError as e: 
            self.status = e.errorCode
            raiseError("while getting color")
        except:
            raise
        
        name, nameList = self.__getBodyName(ibody, nameList, body)
        
        #print xyz
        #print connectivity
        
        self.addVertex(xyz)
        self.addColor(color)
        
        if numPoint == 1:
            self.addPointColor(1, colorMap = "black")
            self.createPoint(name=name)
            ####################################
            # Ryan: how do I set pSize = 6 ?!?!?
            ####################################
        else:
            self.addIndex(connectivity)
            self.addLineColor(1, colorMap = "grey")
            self.createLine(name=name)
        
        return nameList

    cdef __getBodyName(self, object ibody, object nameList, cEGADS.ego body):
        cdef:
            # Attributes
            int atype
            int alen
            const int *ints
            const double *reals 
            const char *string
            
        name = "Body_" + str(ibody+1)
        
        attr = _byteify("_name")
        # Get body _name - if available 
        self.status = cEGADS.EG_attributeRet(body, <char *>attr, &atype, &alen, &ints, &reals, &string)
        if self.status == cEGADS.EGADS_NOTFOUND:
            pass
        elif self.status != cEGADS.EGADS_SUCCESS:
            raise CAPSError(self.status, "while getting trying to retrieve body name (during a call to EG_attributeRet)")
        else:
            if atype == cEGADS.ATTRSTRING:
                name = _strify(string)
        if name:
            if name in nameList:
                name = "Body_" + str(ibody+1)
                print("Body", ibody, ": Name '", name, "' found more than once. Changing name to: '", name, "'")
            
            nameList.append(name)

        return name, nameList
