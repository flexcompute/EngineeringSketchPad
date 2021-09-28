#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.
from _ast import Or

cdef object sortAttrType(int atype, int alen, const int *ints, const double *reals, const char *string):
    cdef: 
        object value
        
    value = []    
    if atype == cEGADS.ATTRINT:
        if alen > 1:
            for i in range(alen):
                value.append(ints[i])
            return value
        else:
            return <object> ints[0]
    
    elif atype == cEGADS.ATTRREAL:
        if alen > 1:
            for i in range(alen):
                value.append(reals[i])
            return value
        else:
            return <object> reals[0]
        
    elif atype == cEGADS.ATTRCSYS:
        
        value = []
        for i in range(alen):
            value.append(reals[i])
        return value

    elif atype == cEGADS.ATTRSTRING:
        return <object> _strify(string)

    raise CAPSError(cEGADS.EGADS_ATTRERR, msg="Unknown attribte type : " + str(atype))

cdef object createAttributeDict(cEGADS.ego obj, object getInternalAttr):
    cdef: 
        int status
    
        int indexAttr = 0, numAttr = 0
        
        # EGADS attribute return variables
        int atype, alen 
        const int    *ints
        const double *reals
        const char *string
        const char *attrName
    
    attrDict = {}
    
    status = cEGADS.EG_attributeNum(obj, &numAttr)
    if status != cEGADS.EGADS_SUCCESS:
        return status, attrDict
    
    for indexAttr in range(numAttr):
        status = cEGADS.EG_attributeGet(obj, indexAttr+1, &attrName, &atype, &alen, &ints, &reals, &string)
        
        if (status != cCAPS.CAPS_SUCCESS and 
            status != cEGADS.EGADS_NOTFOUND):

            return status, attrDict
    
        elif status != cEGADS.EGADS_NOTFOUND:
        
            if getInternalAttr is False and attrName[0] == b'_':
                continue
            
            attrDict[_strify(attrName)] = sortAttrType(atype, alen, ints, reals, string)
    
    return cEGADS.EGADS_SUCCESS, attrDict

cdef object reverseAttributeMap( object attrMap):
    
    reverseMap = {}
    
    attribute = {}
    
    for body in attrMap.keys():
        reverseMap[body] = {}
        
        # Geometry level - Set up blank arrays - "Body", "Face", ...
        for geomLevel in attrMap[body].keys():
            
            if geomLevel is "Body":
                
                # Get attribute names 
                for attrName in attrMap[body][geomLevel].keys():
                    attribute[attrName] = [] 
                
            else:
                # Geometry index value 
                for geomIndex in attrMap[body][geomLevel].keys():
                
                    # Get attribute names 
                    for attrName in attrMap[body][geomLevel][geomIndex].keys():
                        attribute[attrName] = []
        
        # Geometry level - Append attribute values - "Body", "Face", ...
        for geomLevel in attrMap[body].keys():
             
            if geomLevel is "Body":
                # Set attribute names - create list of attribute values 
                for attrName in attrMap[body][geomLevel].keys(): 
                    attribute[attrName].append(attrMap[body][geomLevel][attrName])
        
            else:    
                # Geometry index value 
                for geomIndex in attrMap[body][geomLevel].keys():
                
                    # Set attribute names - create list of attribute values 
                    for attrName in attrMap[body][geomLevel][geomIndex].keys(): 
                        attribute[attrName].append(attrMap[body][geomLevel][geomIndex][attrName])
        
        # Set attribute names and values in a dictionary for the body  
        for attrName in attribute.keys():
            reverseMap[body][attrName] = {}
            
            # Set up geometry level for the attribute values 
            for attrValue in attribute[attrName]:
                
                if isinstance(attrValue, list): # Convert value to string if a list
                    attrValue = str(attrValue) 
                reverseMap[body][attrName][attrValue] = {"Body" : None, "Faces" : [], "Edges" : [], "Nodes" : []} 
       
        # Re-loop through the geometry and count the indexes
        # Geometry level - "Body", "Face", ...
        for geomLevel in attrMap[body].keys():
            
            if geomLevel is "Body":
                
                for attrName in attrMap[body][geomLevel].keys(): 
                    attrValue = attrMap[body][geomLevel][attrName]
                    
                    if isinstance(attrValue, list): # Convert value to string if a list
                        attrValue = str(attrValue) 
                        
                    reverseMap[body][attrName][attrValue][geomLevel] = True
    
            else:
                # Geometry index value 
                for geomIndex in attrMap[body][geomLevel].keys():
                
                    # Get attribute names and values 
                    for attrName in attrMap[body][geomLevel][geomIndex].keys():
                    
                        attrValue = attrMap[body][geomLevel][geomIndex][attrName]
                        
                        if isinstance(attrValue, list): # Convert value to string if a list
                            attrValue = str(attrValue)
                        
                        reverseMap[body][attrName][attrValue][geomLevel].append(geomIndex)
    
    return reverseMap
    
cdef object createAttributeMap(cEGADS.ego body, object getInternalAttr):
    
    cdef: 
        int status
        int i, numGeom
        
        cEGADS.ego *geom=NULL
        
        int bodyIndex
        
        # EGADS Body type return variables 
        cEGADS.ego eref
        cEGADS.ego *echilds=NULL
        int oclass, bodySubType, nchild
        int *senses=NULL
        double data[4]
 
    bodyDict = {}
    faceDict = {}
    edgeDict = {}
    nodeDict = {}
    
    # Fill the geometry dictionaray based on type 
    def fillGeomDict(geomType):
        
        geomDict = {}
        
        status = cEGADS.EG_getBodyTopos(body, 
                                        <cEGADS.ego> NULL, 
                                        geomType, 
                                        &numGeom, 
                                        &geom)
        if status != cEGADS.EGADS_SUCCESS:
            if geom != NULL:
                cEGADS.EG_free(geom)
                
            return (status, geomDict)
    
        # Loop through geometry 
        for i in range(numGeom):
        
            bodyIndex = cEGADS.EG_indexBodyTopo(body, geom[i])
        
            # Create attribution dictionary on face  
            status, geomDict[bodyIndex] = createAttributeDict(geom[i], getInternalAttr)
            
            if status != cEGADS.EGADS_SUCCESS:
            
                if geom != NULL:
                    cEGADS.EG_free(geom)
                   
                return (status, geomDict)
        
        if geom != NULL:
            cEGADS.EG_free(geom)
        
        return (cEGADS.EGADS_SUCCESS, geomDict)   
    
    def returnValue():
        return (status, {"Body" : bodyDict, 
                         "Faces" : faceDict, 
                         "Edges" : edgeDict, 
                         "Nodes" : nodeDict})
        
    # Get body type in case it is a node 
    status = cEGADS.EG_getTopology(body, &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses)
    if status < cEGADS.EGADS_SUCCESS:
        return returnValue()

    # Create attribution dictionary on body 
    status, bodyDict = createAttributeDict(body, getInternalAttr)
    if status != cEGADS.EGADS_SUCCESS:
        return returnValue()
    
    if oclass == cEGADS.NODE: # Just go down the body level 
        return returnValue()
        
    # Get Face dictionary 
    status, faceDict = fillGeomDict(cEGADS.FACE)
    if status != cEGADS.EGADS_SUCCESS:
        return returnValue()
    
    # Get Edge dictionary 
    status, edgeDict = fillGeomDict(cEGADS.EDGE)
    if status != cEGADS.EGADS_SUCCESS:
        return returnValue()
    
    # Get Node dictionary 
    status, nodeDict = fillGeomDict(cEGADS.NODE)
    if status != cEGADS.EGADS_SUCCESS:
        return returnValue()
    
    return returnValue()
        
cdef object createAttributeValList( cEGADS.ego body, int attrLevel, char *attributeKey):
    
    cdef: 
        int status
        
        int i, numGeom
        
        cEGADS.ego *geom=NULL
        
        # EGADS Body type return variables 
        cEGADS.ego eref
        cEGADS.ego *echilds=NULL
        int oclass, bodySubType, nchild
        int *senses=NULL
        double data[4]
        
        # EGADS attribute return variables
        int atype, alen 
        const int    *ints
        const double *reals
        const char *string
        
        object attrList
    
    # Initiate attribute list     
    attrList = []
    
    def getGeomAttr(geomType):
        geomAttr = []
        
        # Determine the number of geometry entities 
        status = cEGADS.EG_getBodyTopos(body, 
                                        <cEGADS.ego> NULL, 
                                        geomType, 
                                        &numGeom, 
                                        &geom)
        if status != cEGADS.EGADS_SUCCESS:

            if geom != NULL:
                cEGADS.EG_free(geom)
                
            return (status, geomAttr)

        # Loop through geometry entities
        for i in range(numGeom):
            
            # Get value following attribute key
            status = cEGADS.EG_attributeRet(geom[i], 
                                            <const char *> attributeKey, 
                                            &atype, 
                                            &alen, 
                                            &ints, 
                                            &reals, 
                                            &string)
            if (status != cEGADS.EGADS_SUCCESS and 
                status != cEGADS.EGADS_NOTFOUND):

                if geom != NULL:
                    cEGADS.EG_free(geom)
                    
                return (status, geomAttr)
            
            elif status != cEGADS.EGADS_NOTFOUND:
                
                geomAttr.append(sortAttrType(atype, alen, ints, reals, string))
        
        if geom != NULL:
            cEGADS.EG_free(geom)
        
        return (cEGADS.EGADS_SUCCESS, geomAttr)
    
    # If body is just a node body - change attrLevel to just the body 
    status = cEGADS.EG_getTopology(body, &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses)
    if status < cEGADS.EGADS_SUCCESS:
        return (status, attrList)
    
    if oclass == cEGADS.NODE:
        attrLevel = 0
    
    # Get value following attributeKey on the body 
    status = cEGADS.EG_attributeRet(body, 
                                    <const char *> attributeKey, 
                                    &atype, 
                                    &alen, 
                                    &ints, 
                                    &reals, 
                                    &string)
    if (status != cEGADS.EGADS_SUCCESS and 
        status != cEGADS.EGADS_NOTFOUND):

        return (status, attrList)
    
    elif status != cEGADS.EGADS_NOTFOUND:
        
        attrList.append(sortAttrType(atype, alen, ints, reals, string))
        
    # Search through faces
    if attrLevel > 0:
        status, tempList = getGeomAttr(cEGADS.FACE)
        
        attrList = attrList + tempList
        
        if status != cEGADS.EGADS_SUCCESS:
            return (status, attrList)
            
        
    # Search through edges
    if attrLevel > 1: 
        status, tempList = getGeomAttr(cEGADS.EDGE)
        
        attrList = attrList + tempList
        
        if status != cEGADS.EGADS_SUCCESS:
            return (status, attrList)
        
            
    # Search through nodes
    if attrLevel > 2:
        status, tempList = getGeomAttr(cEGADS.NODE)
        
        attrList = attrList + tempList
        
        if status != cEGADS.EGADS_SUCCESS:
            return (status, attrList)
        
    return (cEGADS.EGADS_SUCCESS, attrList)

cdef object setElementofCAPSTuple(tupleValue, tupleKey, elementValue, dictKey = None):
    
    if not isinstance(tupleValue, list):
        tupleValue = [tupleValue]
        
    output = []
    for tupleElement in tupleValue:
        
        key = tupleElement[0]
        value = tupleElement[1]
        
        if tupleKey == key:
            
            if isinstance(value, dict): # Modify element of a dictionary 
                value[dictKey] = elementValue
            else: # Else just modify the value 
                value = elementValue        
        
        output.append((key, value))
    
    return output