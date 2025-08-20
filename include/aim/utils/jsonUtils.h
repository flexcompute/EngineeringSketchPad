// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#ifndef _AIM_UTILS_JSON_UTILS_H_
#define _AIM_UTILS_JSON_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

// JSON

// Whether the string value represents a JSON string
int json_isDict(const char *string);

// Get raw string `value` with given `key` in `jsonDict` string
int json_get(const char *jsonDict, const char *key, char **value);

// Get string `value` with given `key` in `jsonDict` string
int json_getString(const char *jsonDict, const char *key, char **value);

// Get array of strings `value` with given `key` in `jsonDict` string
int json_getStringArray(const char *jsonDict, const char *key, int size, char *value[]);

// Get dynamic array of strings `value` with given `key` in `jsonDict` string
int json_getStringDynamicArray(const char *jsonDict, const char *key,
                               int *size, char **value[]);

// Get integer `value` with given `key` in `jsonDict` string
int json_getInteger(const char *jsonDict, const char *key, int *value);

// Get array of integers `value` with given `key` in `jsonDict` string
int json_getIntegerArray(const char *jsonDict, const char *key, int size, int value[]);

// Get dynamic array of integers `value` with given `key` in `jsonDict` string
int json_getIntegerDynamicArray(const char *jsonDict, const char *key,
                                int *size, int *value[]);

// Get double `value` with given `key` in `jsonDict` string
int json_getDouble(const char *jsonDict, const char *key, double *value);

// Get array of doubles `value` with given `key` in `jsonDict` string
int json_getDoubleArray(const char *jsonDict, const char *key, int size, double value[]);

// Get dynamic array of doubles `value` with given `key` in `jsonDict` string
int json_getDoubleDynamicArray(const char *jsonDict, const char *key,
                               int *size, double *value[]);

#ifdef __cplusplus
}
#endif

#endif // _AIM_UTILS_JSON_UTILS_H_
