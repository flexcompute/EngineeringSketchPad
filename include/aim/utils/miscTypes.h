// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Structures for miscellaneous utilities - Written by Dr. Ryan Durscher AFRL/RQVC

#ifndef _AIM_UTILS_MISCTYPES_H_
#define _AIM_UTILS_MISCTYPES_H_

// General container for to map attribute names to an assigned index
typedef struct {

    char *mapName;

    int numAttribute;
    char **attributeName;
    int  *attributeIndex;

} mapAttrToIndexStruct;

#endif // _AIM_UTILS_MISCTYPES_H_
