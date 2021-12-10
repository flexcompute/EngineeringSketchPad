// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef CARDTYPES_H
#define CARDTYPES_H

#include "feaTypes.h" // Bring in fea structures

#define CARD_SMALLWIDTH 8
#define CARD_LARGEWIDTH 16


typedef struct {

    // format
    feaFileTypeEnum formatType;

    int nameWidth;
    int fieldWidth;
    int contWidth;

    const char *delimiter;
    int delimWidth;
    int leftOrRight;


    // card name
    char *name;

    // fields
    size_t size;
    int capacity;
    char **fields;

    // number of fields to increase by when capacity is reached
    size_t allocStep;

} cardStruct;

#endif
