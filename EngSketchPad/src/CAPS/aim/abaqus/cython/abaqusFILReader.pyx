
# Cython declarations
# cython: language_level=3

# Import Cython *.pxds 
cimport cEGADS

# Python version
from cpython.version cimport PY_MAJOR_VERSION

# Import Python modules
import sys
import os

# Handles all byte to uni-code and and byte conversion. 
def _str(data,ignore_dicts = False):
        
    if PY_MAJOR_VERSION >= 3:
        if isinstance(data, bytes):
            return data.decode("utf-8")
   
        if isinstance(data,list):
            return [_str(item,ignore_dicts=True) for item in data]
        if isinstance(data, dict) and not ignore_dicts:
            return {
                    _str(key,ignore_dicts=True): _str(value, ignore_dicts=True)
                    for key, value in data.iteritems()
                    }
    return data

# where variables are Grid Id, T1, T2, T3, R1, R2, R3
cdef public int abaqus_parseFILDisplacement(const char *filename, int *numGridOut, double ***dataMatrixOut) except -1:
    
    cdef:
        int numGrid 
        int numVariable = 7
        double **dataMatrix=NULL

    _, _, data = extractFilRecords(_str(filename))
    _, _, data = extractDispFil(data)
    
    numGrid = int(len(data.keys()))
    
    dataMatrix = <double **> cEGADS.EG_alloc(numGrid*sizeof(double *))
    if not dataMatrix: 
        raise MemoryError()

    i = 0
    while i < numGrid: 
        dataMatrix[i] = <double *> cEGADS.EG_alloc(numVariable*sizeof(double))       
        if not dataMatrix[i]:
            j = 0
            while j < i:
                if dataMatrix[j]:
                    cEGADS.EG_free(dataMatrix[j])
                j += 1
            cEGADS.EG_free(dataMatrix)
        
            raise MemoryError()
        
        i += 1
    
    j = 0  
    for i in data.keys(): 
        dataMatrix[j][0] = float(i)
        
        dataMatrix[j][1] = data[i][1]
        dataMatrix[j][2] = data[i][2]
        dataMatrix[j][3] = data[i][3]
        dataMatrix[j][4] = data[i][4]
        dataMatrix[j][5] = data[i][5]
        dataMatrix[j][6] = data[i][6]
        
        j += 1
    
    #print(data)
    
    numGridOut[0] = numGrid
    dataMatrixOut[0] = dataMatrix
    
    return 0


# where variables are Element Id, vonMises
cdef public int abaqus_parseFILvonMises(const char *filename, int *numGridOut, double ***dataMatrixOut) except -1:
    
    cdef:
        int numGrid 
        int numVariable = 2
        double **dataMatrix=NULL

    _, _, data = extractFilRecords(_str(filename))
    _, _, data = extractMisesFil(data)
    
    numGrid = int(len(data.keys()))
    
    dataMatrix = <double **> cEGADS.EG_alloc(numGrid*sizeof(double *))
    if not dataMatrix: 
        raise MemoryError()

    i = 0
    while i < numGrid: 
        dataMatrix[i] = <double *> cEGADS.EG_alloc(numVariable*sizeof(double))       
        if not dataMatrix[i]:
            j = 0
            while j < i:
                if dataMatrix[j]:
                    cEGADS.EG_free(dataMatrix[j])
                j += 1
            cEGADS.EG_free(dataMatrix)
        
            raise MemoryError()
        
        i += 1
    
    j = 0  
    for i in data.keys(): 
        dataMatrix[j][0] = float(i)
        
        dataMatrix[j][1] = data[i][1]

        j += 1
    
    #print(data)
    
    numGridOut[0] = numGrid
    dataMatrixOut[0] = dataMatrix
    
    return 0

#--------------------------------
# pull out displacement or von mises stress from file.

def extractDispFil(ParsedData, debug=False):
    """Pulls out the node displacement values from a node output file
    and returns a dictionary of the values.

    Returns [status, message, displacements]. Status of 0 is good.
        Displacements is a dictionary with node numbers and U and UR
        status 2 means there was no displacement in files
        Status 3 means the file doesn't exist
    """
    #Get data in a single string
    output=[]
    if ParsedData != []:
        for record in ParsedData:
            if record[0]==101:
                output.append(record[1:])
        if output==[]:
            return [2,"no data found",{}]
        else:

            data={}

            for record in output:
                if record==None:
                    pass
                else:
                    #NodeNum=record[0]
                    displacements={}
                    k=0
                    if debug:
                        print("record:", record)
                    while k < len(record)-1:
                        displacements[k+1]=record[k+1]
                        k+=1
                    data[record[0]]=displacements

            return [0,"found data",data]
    else:
        return [3,"file missing, abq probably didn't finish",{}]

def extractMisesFil(ParsedData, debug=False):
    """Pulls out the element mises values from a element output
    file and returns a dictionary of the values at the Section
    points.

    Returns [status, message, mises]. Status of 0 is good.
       mises is a dictionary with element numbers and Mises at the
       section points
        status 2 means there was no displacement in files
        Status 3 means the file doesn't exist
    """

    if ParsedData != []:
        data={}
        for record,data2 in enumerate((ParsedData)):
            if debug:
                print("data2: ",data2)
            if ParsedData[record][0]==1:
                elNum=ParsedData[record][1]
                secNum=ParsedData[record][3]
                if ParsedData[record+1][0]==12:
                    mises=ParsedData[record+1][1]
                    if elNum in data:
                        data[elNum].update({secNum:mises})
                    else:
                        data[elNum]={secNum:mises}

        if data==[]:
            return [2,"no data found",{}]
        else:
             return [0,"found data",data]
    else:
        return [3,"file missing, abq probably didn't finish",{}]


def extractFilRecords(abqFile,debug=False):
    """Parses the FIL file into a list of records.

    Returns [status, message, ParsedData].
        Status of 0 is good.
        ParsedData is a list of each record parsed into the correct
            format
        status 2 means there was no data in files
        Status 3 means the file doesn't exist
    """
    #Get data in a single string
    Rec=""
    if os.path.isfile(abqFile[0:-4]+".fil"):
        with open(abqFile[0:-4]+".fil") as f:
            for line in f:
                Rec+=line
        Rec=Rec.replace("\n","")
        Rec=Rec.replace("\r","")
        if debug:
            print('record',Rec) #debug

        ind = list(findAll(Rec, "*"))

        if debug:
            print('index:',ind)

        if ind==None:
            return 2,"empty file",[]
        RecSplit=Rec.split('*')
        del RecSplit[0] #first entry is blank

        if debug:
            print(RecSplit[0])

        ParsedData=[]

        for record in RecSplit:
            pos=0
            temp=[]

            if debug:
                print("record length: ",len(record))
                print("record: ",record)
            #numberRec=int(record[3:3+int(record[2])]) #pull the # of records
            if debug:
                print("record length: ",len(record))
            temp2=""
            while pos < len(record):
                if record[pos]=='I':
                    if pos==0: #at the beginning and need to get past intro
                        pos=pos+3+int(record[pos+1:pos+3])
                        if debug:
                            print("position",pos)
                    else:
                        temp.append(int(record[pos+3:pos+3+
                                    int(record[pos+1:pos+3])]))

                        pos=pos+3+int(record[pos+2])
                        if debug:
                            print("temp data: ", temp)
                            print("position", pos)
                elif record[pos]=='D':
                    temp.append(float(record[pos+1:pos+23].replace("D","E")))
                    pos+=22
                elif record[pos]=="A":
                    temp2+=record[pos+1:pos+9]
                    pos+=9
                    if pos < len(record):
                        if record[pos]!="A":
                            temp.append(temp2)
                            temp2==""
                    else:
                        temp.append(temp2)
                        temp2==""
                else:
                    pos+=1
            if debug:
                print("final temp data:", temp)
            ParsedData.append(temp)

        return 0, "found data", ParsedData 
    else:
        return 3, "file missing, abq probably didn't finish", []

def extractRecord(ParsedData,recnumber):
    """Pulls out the values from a parsed fil for just the records desired.

    Returns [status, message, data]. Status of 0 is good.
        data is records of interest
        status 2 means there no records of that type
        Status 3 means the data was empty
    """
    #Get data in a single string
    output=[]
    if ParsedData != []:
        for record in ParsedData:
            if record[0]==recnumber:
                output.append(record[1:])
        if output==[]:
            return [2,"no data found",[]]
        else:
            return [0,"found data",output]
    else:
        return [3,"data missing, abq probably didn't finish",[]]

def findAll(a_str, sub):
    """looks in a string for all occurances of a string
    """
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches
