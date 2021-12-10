#
# Written by Dr. Ryan Durscher AFRL/RQVC and Dr. Marshall Galbraith MIT
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

from pyEGADS import egads

def _createAttributeValList( body, attrLevel, attributeKey ):

    # Initiate attribute list
    attrList = []
    
    def getTopoAttr(topoType):
        topoAttr = []
        
        # Determine the number of geometry entities 
        topologies = body.getBodyTopos(topoType)

        # Loop through geometry entities
        for topo in topologies:
            
            # Get value following attribute key
            attr = topo.attributeRet(attributeKey)
            if attr is not None:
                if isinstance(attr,list) and not isinstance(attr,egads.csystem) \
                   and len(attr) == 1: attr = attr[0]
                attrList.append(attr)
        
        return topoAttr
    
    # If body is just a node body - change attrLevel to just the body 
    (eref, oclass, bodySubType, lim, children, sens) = body.getTopology()
    
    if oclass == egads.NODE:
        attrLevel = 0
    
    # Get value following attributeKey on the body 
    attr = body.attributeRet(attributeKey)
    if attr is not None:
        if isinstance(attr,list) and len(attr) == 1: attr = attr[0]
        attrList.append(attr)
        
    # Search through faces
    if attrLevel > 0:
        attrList += getTopoAttr(egads.FACE)
             
    # Search through edges
    if attrLevel > 1: 
        attrList += getTopoAttr(egads.EDGE)

    # Search through nodes
    if attrLevel > 2:
        attrList += getTopoAttr(egads.NODE)
         
    return attrList


def _createAttributeDict(obj, getInternalAttr):
    
    attrDict = {}
    
    numAttr = obj.attributeNum()
    
    for indexAttr in range(numAttr):
        attrName, data = obj.attributeGet(indexAttr+1)
        
        if data is None: continue
        if getInternalAttr is False and attrName[0] == '_':
            continue
        
        attrDict[attrName] = data
    
    return attrDict

def _createAttributeMap(body, getInternalAttr):

    bodyDict = {}
    faceDict = {}
    edgeDict = {}
    nodeDict = {}
    
    # Fill the topology dictionaray based on type 
    def fillTopoDict(topoType):
        
        topoDict = {}
        
        topologies = body.getBodyTopos(topoType)
    
        # Loop through geometry 
        for topo in topologies:
        
            bodyIndex = body.indexBodyTopo(topo)
        
            # Create attribution dictionary on face  
            topoDict[bodyIndex] = _createAttributeDict(topo, getInternalAttr)
        
        return topoDict
    
    def returnValue():
        return {"Body"  : bodyDict, 
                "Faces" : faceDict, 
                "Edges" : edgeDict, 
                "Nodes" : nodeDict}
        
    # Get body type in case it is a node 
    (eref, oclass, bodySubType, lim, children, sens) = body.getTopology()

    # Create attribution dictionary on body 
    bodyDict = _createAttributeDict(body, getInternalAttr)
    
    if oclass == egads.NODE: # Just go down the body level 
        return returnValue()
        
    # Get Face dictionary 
    faceDict = fillTopoDict(egads.FACE)
    
    # Get Edge dictionary 
    edgeDict = fillTopoDict(egads.EDGE)
    
    # Get Node dictionary 
    nodeDict = fillTopoDict(egads.NODE)
    
    return returnValue()

def _reverseAttributeMap( attrMap ):
    
    reverseMap = {}
    
    attribute = {}
    
    for body in attrMap.keys():
        reverseMap[body] = {}
        
        # Geometry level - Set up blank arrays - "Body", "Face", ...
        for geomLevel in attrMap[body].keys():
            
            if geomLevel == "Body":
                
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
             
            if geomLevel == "Body":
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
            
            if geomLevel == "Body":
                
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
