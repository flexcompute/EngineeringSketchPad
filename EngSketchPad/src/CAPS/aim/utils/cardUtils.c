// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>
// #include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include "miscUtils.h"
#include "cardUtils.h"


#ifdef WIN32
#define strcasecmp  stricmp
#define strtok_r   strtok_s
#endif



static int _initField(cardStruct *card, int fieldIndex) {

    if (card == NULL) return CAPS_NULLVALUE;

    if (fieldIndex < card->capacity) return CAPS_RANGEERR;

    card->fields[fieldIndex] = NULL;
    return CAPS_SUCCESS;
}

static int _allocFields(cardStruct *card, int capacity) {

    int i;

    if (card == NULL) return CAPS_NULLVALUE;

    // must have at least 1 field
    if (capacity < 1) {
        PRINT_ERROR("Card capacity must be greater than 0, "
                   "received capacity = %d", capacity);
        return CAPS_BADVALUE;
    }

    card->fields = EG_alloc(sizeof(char *) * capacity);
    if (card->fields == NULL) return EGADS_MALLOC;

    for (i = 0; i < capacity; i++) {
        _initField(card, i);
    }

    // update capacity
    card->capacity = capacity;

    return CAPS_SUCCESS;
}

static int _shrinkFields(cardStruct *card, int lessCapacity) {

    int i, capacity;
    char *field = NULL;

    if (card == NULL) return CAPS_NULLVALUE;

    for (i = 0; i < lessCapacity; i++) {
        field = card->fields[card->capacity - i - 1];
        if (field != NULL) {
            EG_free(field);
        }
    }

    capacity = card->capacity - lessCapacity;

    card->fields = EG_reall(card->fields,
                            sizeof(char *) * capacity);
    if (card->fields == NULL) return EGADS_MALLOC;

    // update capacity, and size if needed
    card->capacity = capacity;
    if (card->size > card->capacity) {
        card->size = card->capacity;
    }

    return CAPS_SUCCESS;
}

static int _growFields(cardStruct *card, int moreFields) {

    int i, capacity;

    if (card == NULL) return CAPS_NULLVALUE;

    capacity = card->capacity + moreFields;

    card->fields = EG_reall(card->fields,
                            sizeof(char *) * capacity);
    if (card->fields == NULL) return EGADS_MALLOC;

    for (i = moreFields; i > 0; i--) {
        _initField(card, capacity - i - 1);
    }

    // update capacity
    card->capacity = capacity;

    return CAPS_SUCCESS;
}

static const char *none = "";
static int _setFormatSmall(cardStruct *card) {

    card->nameWidth = 8;
    card->fieldWidth = 8;
    card->contWidth = 8;
    card->delimiter = none;
    card->delimWidth = 0;
    card->leftOrRight = 1;

    return CAPS_SUCCESS;
}

static const char *comma = ",";
static int _setFormatFree(cardStruct *card) {

    card->nameWidth = 8;
    card->fieldWidth = 7;
    card->contWidth = 8;
    card->delimiter = comma;
    card->delimWidth = 1;
    card->leftOrRight = 1;

    return CAPS_SUCCESS;
}

static const char *space = " ";
static int _setFormatLarge(cardStruct *card) {

    card->nameWidth = 8;
    card->fieldWidth = 15;
    card->contWidth = 8;
    card->delimiter = space;
    card->delimWidth = 1;
    card->leftOrRight = 1;

    return CAPS_SUCCESS;
}

static int _setFormat(cardStruct *card, const feaFileTypeEnum fileType) {

    if (card == NULL) return CAPS_NULLVALUE;

    card->formatType = fileType;

    if (fileType == SmallField) {
        return _setFormatSmall(card);
    }
    else if (fileType == FreeField) {
        return _setFormatFree(card);
    }
    else if (fileType == LargeField) {
        return _setFormatLarge(card);
    }
    else {
        PRINT_ERROR("Unrecognized format type: %d", card->formatType);
        return CAPS_BADVALUE;
    }
}

static int _setName(cardStruct *card, const char *name) {

    if (card == NULL) return CAPS_NULLVALUE;
    if (name == NULL) {
        PRINT_ERROR("Card name cannot be NULL");
        return CAPS_NULLVALUE;
    }

    card->name = EG_strdup(name);

    if (card->name == NULL) return EGADS_MALLOC;

    return CAPS_SUCCESS;
}

static int _setNull(cardStruct *card) {

    if (card == NULL) return CAPS_NULLVALUE;

    card->formatType = UnknownFileType;
    card->delimiter = NULL;
    card->fieldWidth = 0;
    card->delimWidth = 0;
    card->leftOrRight = 0;
    card->size = 0;
    card->name = NULL;
    card->capacity = 0;
    card->fields = NULL;
    card->allocStep = 0;

    return CAPS_SUCCESS;
}

/*@null@*/
static char * _formatField(const char *fieldValue, const char *pad, int fieldWidth) {

    int length, padlen;
    char *fieldFormatted = NULL;

    if (pad == NULL) {
        pad = "";
    }
    padlen = strlen(pad);

    // if has not been formatted
    if (strlen(fieldValue) != fieldWidth) {

        fieldFormatted = EG_alloc(sizeof(char) * (fieldWidth + 1));
        if (fieldFormatted == NULL) return NULL;

        // TODO: currently always pads right
        length = sprintf(fieldFormatted, "%*.*s%s", fieldWidth-padlen, fieldWidth-padlen, fieldValue, pad);
        assert(length == fieldWidth); // ensure `sprintf` did not overflow

        fieldFormatted[length] = 0;

    } else {
        fieldFormatted = EG_strdup(fieldValue);
    }

    return fieldFormatted;
}

static int _addFormattedField(cardStruct * card, const char *fieldFormatted,
                              int fieldWidth, int fieldSpan) {

    int i;
    char *field = NULL;

    if (fieldFormatted == NULL) return CAPS_NULLVALUE;

    if (fieldSpan > 1) {

        for ( i = 0; i < fieldSpan; i++) {

            // get i-th cell of the full formatted field

            field = EG_alloc(sizeof(char) * fieldWidth + 1);
            if (field == NULL) return EGADS_MALLOC;

            strncpy(field, fieldFormatted + i * fieldWidth, fieldWidth);
            card->fields[card->size++] = field;
        }
    }
    else {
        card->fields[card->size++] = EG_strdup(fieldFormatted);
    }

    return CAPS_SUCCESS;
}

static int _addFieldSmall(cardStruct *card, const char *fieldValue, int fieldSpan) {

    int status;
    char *fieldFormatted = NULL;

    if (card->size >= card->capacity) return CAPS_RANGEERR;
    if (fieldSpan <= 0) return CAPS_BADVALUE;

    fieldFormatted = _formatField(fieldValue, card->delimiter, CARD_SMALLWIDTH * fieldSpan);

    status = _addFormattedField(card, fieldFormatted, CARD_SMALLWIDTH, fieldSpan);

    if (fieldFormatted != NULL) {
        EG_free(fieldFormatted);
    }

    return status;
}

static int _addFieldFree(cardStruct *card, const char *fieldValue, int fieldSpan) {

    int status;
    char *fieldFormatted = NULL;

    if (card->size >= card->capacity) return CAPS_RANGEERR;
    if (fieldSpan <= 0) return CAPS_BADVALUE;

    fieldFormatted = _formatField(fieldValue, NULL, card->fieldWidth * fieldSpan);

    status = _addFormattedField(card, fieldFormatted, card->fieldWidth, fieldSpan);

    if (fieldFormatted != NULL) {
        EG_free(fieldFormatted);
    }

    return status;
}

static int _addFieldLarge(cardStruct *card, const char *fieldValue, int fieldSpan) {

    int status;
    char *fieldFormatted = NULL;

    if (card->size >= card->capacity) return CAPS_RANGEERR;
    if (fieldSpan <= 0) return CAPS_BADVALUE;

    fieldFormatted = _formatField(fieldValue, card->delimiter, CARD_LARGEWIDTH * fieldSpan);

    status = _addFormattedField(card, fieldFormatted, CARD_LARGEWIDTH, fieldSpan);

    if (fieldFormatted != NULL) {
        EG_free(fieldFormatted);
    }

    return status;
}
static int _countTrailingBlankFields(cardStruct *card) {

    int i, numTrailingBlanks = 0;

    for (i = card->size - 1; i >= 0; i--, numTrailingBlanks++) {
        if (!card_isBlankField(card->fields[i])) break;
    }

    return numTrailingBlanks;
}

static inline int _sprintEnd(char *buffer) {
    return sprintf(buffer, "\n");
}

static int _calcTotalChars(cardStruct *card, int numCells, int numCont) {

    int totalChars = card->nameWidth;

    if (card->formatType == SmallField) {
        totalChars += numCells * CARD_SMALLWIDTH; // chars in fields
        totalChars += numCont * (card->contWidth * 2 + 1); // chars in continuations
    }
    else if (card->formatType == FreeField) {
        totalChars += numCells * card->fieldWidth; // chars in fields
        totalChars += numCells * card->delimWidth;
        totalChars += numCont * (card->contWidth * 2 + 1); // chars in continuations
        totalChars += numCont * card->delimWidth;
    }
    else if (card->formatType == LargeField) {
        totalChars += numCells * CARD_LARGEWIDTH; // chars in fields
        totalChars += numCont * (card->contWidth * 2 + 1); // chars in continuations
    }

    totalChars += 1; // newline char

    return totalChars;
}

static void _removeTrailingDecimalZeros(char *doubleString) {

    int i, count;
    char *c, *decPoint = NULL, *exp = NULL;

    decPoint = strchr(doubleString, '.'); // find decimal point
    exp = strchr(doubleString, 'E');

    // if has decimal point 
    if (decPoint != NULL) {

        // if no exponent
        if (exp == NULL) {

            // overwrite zeros with null terminator
            for (i = strlen(doubleString)-1; i >= 0; i--) {

                if (doubleString[i] == '0') {
                    doubleString[i] = '\0';
                }
                else {
                    break;
                }
            }
        }
        else {

            // count trail zeros between decPoint and exp
            count = 0;
            for (c = exp-1; c != decPoint; c--) {
                
                if (*c == '0') {
                    count++;
                }
                else {
                    break;
                }
            }

            // shift exponent chars left by count, overwriting zeros
            memmove(exp - count, exp, strlen(exp));
            doubleString[strlen(doubleString) - count] = '\0';
        }

    }
}

static inline int _sprintName(char *buffer, char *name, int nameWidth) {
    int length = 0;
    length += sprintf(buffer, "%-*s", nameWidth, name);
    return length;
}

static inline int _sprintNameLarge(char *buffer, char *name, int nameWidth) {
    int length = 0;
    length += sprintf(buffer, "%-*s", nameWidth, name);
    assert(strlen(name) < nameWidth);
    buffer[strlen(name)] = '*';
    return length;
}

static inline int _sprintFieldFixed(char *buffer, const char *field) {
    return sprintf(buffer, "%s", field);
}

static inline int _sprintFieldFree(char *buffer, const char *field, const char *delimiter) {
    return sprintf(buffer, "%s%s", delimiter, field);
}

static inline int _sprintContSmall(char *buffer, int contCount,
                      /*@unused@*/ const char *delimiter,
                                   int contWidth) {
    int length = 0;
    length += sprintf(buffer + length, "+%-*d\n", contWidth-1, contCount);
    length += sprintf(buffer + length, "+%-*d", contWidth-1, contCount);
    return length;
}

static inline int _sprintContFree(char *buffer, int contCount,
                                  const char *delimiter, int contWidth) {
    int length = 0;
    length += sprintf(buffer, "%s", delimiter);
    length += sprintf(buffer + length, "+%-*d\n", contWidth-1, contCount);
    length += sprintf(buffer + length, "+%-*d", contWidth-1, contCount);
    return length;
}

static inline int _sprintContLarge(char *buffer, int contCount,
                      /*@unused@*/ const char *delimiter,
                                   int contWidth) {
    int length = 0;
    length += sprintf(buffer + length, "*%-*d\n", contWidth-1, contCount);
    length += sprintf(buffer + length, "*%-*d", contWidth-1, contCount);
    return length;
}

static inline int _requiresContSmall(int cellIndex) {
    return (cellIndex != 0) && (cellIndex % 8 == 0);
}

static inline int _requiresContFree(int cellIndex) {
    return (cellIndex != 0) && (cellIndex % 8 == 0);
}

static inline int _requiresContLarge(int cellIndex) {
    return (cellIndex != 0) && (cellIndex % 4 == 0);
}


static int _concatName(cardStruct *card, char *buffer) {

    if (card->formatType == SmallField) {
        return _sprintName(buffer, card->name, card->nameWidth);
    }
    else if (card->formatType == FreeField) {
        return _sprintName(buffer, card->name, card->nameWidth);
    }
    else if (card->formatType == LargeField) {
        return _sprintNameLarge(buffer, card->name, card->nameWidth);
    }
    else {
        return CAPS_BADVALUE;
    }
}

static int _concatField(cardStruct *card, const int i, char *buffer) {

    if (card->formatType == SmallField) {
        return _sprintFieldFixed(buffer, card->fields[i]);
    }
    else if (card->formatType == FreeField) {
        return _sprintFieldFree(buffer, card->fields[i], card->delimiter);
    }
    else if (card->formatType == LargeField) {
        return _sprintFieldFixed(buffer, card->fields[i]);
    }
    else {
        return CAPS_BADVALUE;
    }
}

static int _concatContinuation(cardStruct *card, int contIndex, char *buffer) {

    if (card->formatType == SmallField) {
        return _sprintContSmall(buffer, contIndex, card->delimiter, card->contWidth);
    }
    else if (card->formatType == FreeField) {
        return _sprintContFree(buffer, contIndex, card->delimiter, card->contWidth);
    }
    else if (card->formatType == LargeField) {
        return _sprintContLarge(buffer, contIndex, card->delimiter, card->contWidth);
    }
    else {
        return CAPS_BADVALUE;
    }
}

static int _requiresContinuation(cardStruct *card, int cellIndex) {

    if (card->formatType == SmallField) {
        return _requiresContSmall(cellIndex);
    }
    else if (card->formatType == FreeField) {
        return _requiresContFree(cellIndex);
    }
    else if (card->formatType == LargeField) {
        return _requiresContLarge(cellIndex);
    }
    else {
        return CAPS_BADVALUE;
    }
}

int card_initiate(cardStruct *card, const char *name, feaFileTypeEnum formatType) {

    int status = CAPS_SUCCESS;

    if (card == NULL) return CAPS_NULLVALUE;

    // initialize card variables to null
    _setNull(card);

    card->allocStep = 8;

    // set formatting options
    status = _setFormat(card, formatType);
    if (status != CAPS_SUCCESS) return status;

    // set card name
    status = _setName(card, name);
    if (status != CAPS_SUCCESS) return status;

    // allocate fields
    status = _allocFields(card, card->allocStep);
    if (status != CAPS_SUCCESS) return status;

    return CAPS_SUCCESS;
}

int card_destroy(cardStruct *card) {

    int i, status = CAPS_SUCCESS;

    if (card == NULL) return CAPS_NULLVALUE;

    if (card->name != NULL) EG_free(card->name);

    for (i = 0; i < card->size; i++) {
        if (card->fields[i] != NULL) EG_free(card->fields[i]);
    }
    if (card->fields != NULL) EG_free(card->fields);

    status = _setNull(card);
    if (status != CAPS_SUCCESS) return status;

    return CAPS_SUCCESS;
}

int card_resize(cardStruct *card, int capacity) {

    if (card == NULL) return CAPS_NULLVALUE;

    if (card->capacity == 0) {
        return _allocFields(card, capacity);
    }
    else if (card->capacity < capacity) {
        return _growFields(card, capacity - card->capacity);
    }
    else if (card->capacity > capacity) {
        return _shrinkFields(card, card->capacity - capacity);
    }
    else { // card->capacity == capacity
        return CAPS_SUCCESS;
    }
}

// int card_strip(cardStruct *card) { // TODO: this should only affect size, not capacity

//     int numTrailBlanks = _countTrailingBlankFields(card);

//     card->size -= numTrailBlanks;

//     // shrink fields if necessary, removing any trailing blank fields also
//     if (card->size > card->capacity) {
//         return CAPS_RANGEERR;
//     }
//     else if (card->size < card->capacity) {
//         return _shrinkFields(card, card->capacity - card->size);
//     }
//     return CAPS_SUCCESS;
// }

int card_addField(cardStruct *card, const char *fieldValue, int fieldSpan) {

    int status = CAPS_SUCCESS;

    if (fieldValue == NULL) return CAPS_NULLVALUE;

    if (card->size == card->capacity) {
        status = _growFields(card, card->allocStep);
        if (status != CAPS_SUCCESS) return status;
    }

    if (card->formatType == SmallField) {
        status = _addFieldSmall(card, fieldValue, fieldSpan);
    }
    else if (card->formatType == FreeField) {
        status = _addFieldFree(card, fieldValue, fieldSpan);
    }
    else if (card->formatType == LargeField) {
        status = _addFieldLarge(card, fieldValue, fieldSpan);
    }
    else {
        status = CAPS_BADVALUE;
    }

    return status;
}

int card_addBlank(cardStruct *card) {

    return card_addField(card, " ", 1);
}

int card_addBlanks(cardStruct *card, int numBlanks) {

    int i, status = CAPS_SUCCESS;

    for (i = 0; i < numBlanks; i++) {
        status = card_addBlank(card);
        if (status != CAPS_SUCCESS) return status;
    }

    return status;
}

int card_continue(cardStruct *card) {

    int fieldsPerLine;

    if (card->formatType == SmallField) {
        fieldsPerLine = 8;
    }
    else if (card->formatType == FreeField) {
        fieldsPerLine = 8;
    }
    else if (card->formatType == LargeField) {
        fieldsPerLine = 4;
    }
    else {
        return CAPS_BADVALUE;
    }

    if ((card->size % fieldsPerLine) != 0)
        return card_addBlanks(card, fieldsPerLine - (card->size % fieldsPerLine));

    return CAPS_SUCCESS;
}

int card_addString(cardStruct *card, const char *fieldValue) {
    return card_addLongString(card, fieldValue, 1);
}

int card_addLongString(cardStruct *card, const char *fieldValue, int fieldSpan) {

    if (fieldValue != NULL) {
        return card_addField(card, fieldValue, fieldSpan);
    }
    else {
        return card_addBlank(card);
    }
}

int card_addStringArray(cardStruct *card, int numFieldValues, char *fieldValues[]) {

    int i, status = CAPS_SUCCESS;

    for (i = 0; i < numFieldValues; i++) {
        status = card_addString(card, fieldValues[i]);
        if (status != CAPS_SUCCESS) return status;
    }
    return status;
}

// TODO: handle long int ?
int card_addInteger(cardStruct *card, int fieldValue) {

    int status = CAPS_SUCCESS;
    char *integerString = NULL;

    integerString = convert_integerToString(fieldValue, card->fieldWidth, card->leftOrRight);
    status = card_addField(card, integerString, 1);
    EG_free(integerString);

    return status;
}

int card_addIntegerArray(cardStruct *card, int numFieldValues, int fieldValues[]) {

    int i, status = CAPS_SUCCESS;

    for (i = 0; i < numFieldValues; i++) {
        status = card_addInteger(card, fieldValues[i]);
        if (status != CAPS_SUCCESS) return status;
    }
    return status;
}

int card_addIntegerOrBlank(cardStruct *card, int *fieldValue) {
    
    if (fieldValue == NULL) {
        return card_addBlank(card);
    }
    else {
        return card_addInteger(card, *fieldValue);
    }
}

int card_addDouble(cardStruct *card, double fieldValue) {
    int status = CAPS_SUCCESS;
    char *doubleString = NULL;

    doubleString = convert_doubleToString(fieldValue, card->fieldWidth, card->leftOrRight);
    _removeTrailingDecimalZeros(doubleString);
    status = card_addField(card, doubleString, 1);
    EG_free(doubleString);

    return status;
}

int card_addDoubleArray(cardStruct *card, int numFieldValues, double fieldValues[]) {

    int i, status = CAPS_SUCCESS;

    for (i = 0; i < numFieldValues; i++) {
        status = card_addDouble(card, fieldValues[i]);
        if (status != CAPS_SUCCESS) return status;
    }
    return status;
}

int card_addDoubleOrBlank(cardStruct *card, double *fieldValue) {
    
    if (fieldValue == NULL) {
        return card_addBlank(card);
    }
    else {
        return card_addDouble(card, *fieldValue);
    }
}

int card_isBlankField(const char *field) {
    int i;
    for (i = 0; i < strlen(field); i++) {
        if (!isspace(field[i])) return false;
    }
    return true;
}

char * card_toString(cardStruct *card) {

    int i;

    int fieldCount, contCount, numFields;
    int totalChars, cardLength;
    char *cardString = NULL;

    numFields = card->size - _countTrailingBlankFields(card);

    fieldCount = 0;
    contCount = 0;
    for (i = 0; i < numFields; i++) {
        if (_requiresContinuation(card, fieldCount)) {
            contCount++;
        }
        fieldCount++;
    }

    totalChars = _calcTotalChars(card, fieldCount, contCount);

    cardString = EG_alloc(sizeof(char) * (totalChars + 1));
    if (cardString == NULL) return NULL;

    cardLength = 0;
    cardLength += _concatName(card, cardString);

    fieldCount = 0;
    contCount = 0;
    for (i = 0; i < numFields; i++) {

        if (_requiresContinuation(card, fieldCount)) {
            // write continuation
            cardLength += _concatContinuation(card, contCount, cardString + cardLength);
            contCount++;
        }

        // write field
        cardLength += _concatField(card, i, cardString + cardLength);
        fieldCount++;
    }
    cardLength += _sprintEnd(cardString + cardLength);
    // printf("totalChars: %d, cardLength: %d\n", totalChars, cardLength);
    assert(cardLength == totalChars);

    cardString[cardLength] = 0;

    return cardString;
}

void card_write(cardStruct *card, FILE *fp) {

    char *cardString = NULL;

    cardString = card_toString(card);
    fprintf(fp, "%s", cardString);

    if (cardString != NULL) {
        EG_free(cardString);
    }
}

void card_print(cardStruct *card) {
    card_write(card, stdout);
}
