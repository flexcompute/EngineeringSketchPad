/*
 *	The Web Viewer
 *
 *		Default CallBack
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <string.h>

#include "wsserver.h"


void browserMessage(/*@unused@*/ void *userPtr, void *wsi, char *text,
                    /*@unused@*/ int len)
{
  int stat;
  
  if (wsi == NULL) {
    printf(" BuiltIn browserMessage (from server): %s\n", text);
  } else {
    printf(" BuiltIn browserMessage: %s\n", text);
    if (text != NULL)
      if (strcmp(text, "bounce") == 0) {
        stat = wv_makeMessage(wsi, text);
        if (stat != 0) printf(" ERROR: wv_makeMessage = %d\n", stat);
      }
  }
}
