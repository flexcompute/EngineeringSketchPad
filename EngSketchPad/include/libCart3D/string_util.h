/*
 * general string utilities.  these are overflow-safe.
 */

#ifndef __STRING_UTIL_H_
#define __STRING_UTIL_H_

/*
 * $Id: string_util.h,v 1.3 2004/11/30 00:05:28 smurman Exp $
 */

#include <string.h>

/* get STRING_LEN from c3d_global */
#include "c3d_global.h"
#include "basicTypes.h"

#define EQUAL_STRINGS(a, b) (0 == strcmp((a), (b)))
#define NEQUAL_STRINGS(a, b) (0 != strcmp((a), (b)))

/* 
 * replace all occurences of substr in string with the new_substr
 * return the result in new_string.  
 * new_string contains the original string if substr not found.
 * buflen is the length of the new_string buffer 
 */
void string_substitute(char *string, const char *substr,
		       char *new_string, const char *new_substr, const int buflen);

/*
 * append the addstr to string, and return result in new_string
 * buflen is the length of the new_string buffer 
 */ 
void string_append(char *string, const char *addstr,
		   char *new_string, const int buflen);

/*
 * remove any leading and trailing whitespace from the string and
 * return the result in new_string
 */
void string_remove_ws_ends(const char *string, char *new_string, const int buflen);

/*
 * conversion routines
 * convert the string to a number.  returns FALSE if conversion fails.
 */
bool string_convert2double(const char *str, double *val);

bool string_convert2int(const char *str, int *val);

#endif /* __STRING_UTIL_H_ */
