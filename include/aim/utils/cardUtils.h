// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_CARDUTILS_H_
#define _AIM_UTILS_CARDUTILS_H_

#include "cardTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// initialize card
int card_initiate(cardStruct *card, const char *name, feaFileTypeEnum formatType);

// destroy card
int card_destroy(cardStruct *card);

// resize card field array capacity
int card_resize(cardStruct *card, int size);

// add field to card
int card_addField(cardStruct *card, const char *fieldValue, int fieldSpan);

// add blank field
int card_addBlank(cardStruct *card);

// add blank fields
int card_addBlanks(cardStruct *card, int numBlanks);

// add blank fields until end of line (continuation)
int card_continue(cardStruct *card);

// add character string field
int card_addString(cardStruct *card, /*@null@*/ const char *fieldValue);

// add character string field, if null add blank
int card_addStringOrBlank(cardStruct *card, /*@null@*/ const char *fieldValue);

// add character string spanning `fieldSpan` fields
int card_addLongString(cardStruct *card, /*@null@*/ const char *fieldValue, int fieldSpan);

// add character string fields
int card_addStringArray(cardStruct *card, int numFieldValues, char *fieldValues[]);

// add integer field
int card_addInteger(cardStruct *card, int fieldValue);

// add integer fields
int card_addIntegerArray(cardStruct *card, int numFieldValues, const int fieldValues[]);

// add integer fields
int card_addIntegerOrBlankArray(cardStruct *card, int numFieldValues, const int *fieldValues[]);

// // add integer field from pointer, if null add blank
int card_addIntegerOrBlank(cardStruct *card, /*@null@*/ const int *fieldValue);

// add real field
int card_addDouble(cardStruct *card, double fieldValue);

// add real fields
int card_addDoubleArray(cardStruct *card, int numFieldValues, const double fieldValues[]);

// add real field from pointer, if null add blank
int card_addDoubleOrBlank(cardStruct *card, /*@null@*/ const double *fieldValue);

// does `field` represent a blank field ?
int card_isBlankField(const char *field);

// convert card name and fields to formatted string representing the card
/*@null@*/
char * card_toString(cardStruct *card);

// write the string representation of the card to file `fp`
void card_write(cardStruct *card, FILE *fp);

// write the string representation of the card to stdout
void card_print(cardStruct *card);


#ifdef __cplusplus
}
#endif

#endif // _AIM_UTILS_CARDUTILS_H_
