// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <stdio.h>
#include <string.h>

#include "aimUtil.h"
#include "jsonUtils.h"
#include "miscUtils.h"


// JSON

// Return whether the string value represents a JSON string
// Currently only checks whether the first character is '{'
int json_isDict(const char *string)
{
    return strncmp(string, "{", 1) == 0;
}

// Get raw string `value` with given `key` in `jsonDict` string
// Simple wrapper over search_jsonDictionary function
int json_get(const char *jsonDict, const char *key, char **value)
{

    if (!json_isDict(jsonDict)) {
        return CAPS_BADVALUE;
    }

    return search_jsonDictionary(jsonDict, key, value);
}

// Get string `value` with given `key` in `jsonDict` string
// NOTE: overwrites value pointer with newly allocated string
int json_getString(const char *jsonDict, const char *key, char **value)
{
    int status;
    char *valueStr = NULL;

    if (*value != NULL) {
        AIM_FREE(*value);
    }

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        *value = string_removeQuotation(valueStr);
    }

    if (valueStr != NULL) AIM_FREE(valueStr);

    return status;
}

// Get array of strings `value` with given `key` in `jsonDict` string
int json_getStringArray(const char *jsonDict, const char *key, int size, char *value[])
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toStringArray(valueStr, size, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get dynamic array of strings `value` with given `key` in `jsonDict` string
int json_getStringDynamicArray(const char *jsonDict, const char *key,
                               int *size, char **value[])
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toStringDynamicArray(valueStr, size, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get integer `value` with given `key` in `jsonDict` string
int json_getInteger(const char *jsonDict, const char *key, int *value)
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toInteger(valueStr, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get array of integers `value` with given `key` in `jsonDict` string
int json_getIntegerArray(const char *jsonDict, const char *key, int size, int value[])
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toIntegerArray(valueStr, size, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get dynamic array of integers `value` with given `key` in `jsonDict` string
int json_getIntegerDynamicArray(const char *jsonDict, const char *key,
                                int *size, int *value[])
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toIntegerDynamicArray(valueStr, size, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get double `value` with given `key` in `jsonDict` string
int json_getDouble(const char *jsonDict, const char *key, double *value)
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toDouble(valueStr, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get array of doubles `value` with given `key` in `jsonDict` string
int json_getDoubleArray(const char *jsonDict, const char *key, int size, double value[])
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toDoubleArray(valueStr, size, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}

// Get dynamic array of doubles `value` with given `key` in `jsonDict` string
int json_getDoubleDynamicArray(const char *jsonDict, const char *key,
                               int *size, double *value[])
{
    int status;
    char *valueStr = NULL;

    status = json_get(jsonDict, key, &valueStr);
    if (status == CAPS_SUCCESS) {
        status = string_toDoubleDynamicArray(valueStr, size, value);
    }

    if (valueStr != NULL) EG_free(valueStr);

    return status;
}
