// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef __JSON_UTILS_H__
#define __JSON_UTILS_H__


// JSON

// Whether the string value represents a JSON string
int json_isDict(char *string);

// Get raw string `value` with given `key` in `jsonDict` string
int json_get(char *jsonDict, char *key, char **value);

// Get string `value` with given `key` in `jsonDict` string
int json_getString(char *jsonDict, char *key, char **value);

// Get array of strings `value` with given `key` in `jsonDict` string
int json_getStringArray(char *jsonDict, char *key, int size, char *value[]);

// Get dynamic array of strings `value` with given `key` in `jsonDict` string
int json_getStringDynamicArray(char *jsonDict, char *key, int *size,
                               char **value[]);

// Get integer `value` with given `key` in `jsonDict` string
int json_getInteger(char *jsonDict, char *key, int *value);

// Get array of integers `value` with given `key` in `jsonDict` string
int json_getIntegerArray(char *jsonDict, char *key, int size, int value[]);

// Get dynamic array of integers `value` with given `key` in `jsonDict` string
int json_getIntegerDynamicArray(char *jsonDict, char *key, int *size,
                                int *value[]);

// Get double `value` with given `key` in `jsonDict` string
int json_getDouble(char *jsonDict, char *key, double *value);

// Get array of doubles `value` with given `key` in `jsonDict` string
int json_getDoubleArray(char *jsonDict, char *key, int size, double value[]);

// Get dynamic array of doubles `value` with given `key` in `jsonDict` string
int json_getDoubleDynamicArray(char *jsonDict, char *key, int *size,
                               double *value[]);


#endif // __JSON_UTILS_H__
