/*
 * $Id: myLList.h,v 1.3 2004/04/07 15:18:14 aftosmis Exp $
 *
 */
/* ---------
  myllist.h
  ----------
                            ...a high performance array based linked
			                         list structure for general use
*/
#ifndef __MY_LLIST_H__
#define __MY_LLIST_H__

#define TAIL                  -999      /* flag for Tail of linked list      */

typedef  struct      LinkedListStructure  tsLList;
typedef  tsLList*    p_tsLList;
struct LinkedListStructure{              /* ...in loops, this is faster than */
  int index,next;                        /*       pointer based Linked Lists */
};
                                         /*--  and Linked List Info struct - */
typedef  struct      LinkedListInfoStruct tsLLinfo;
typedef  tsLLinfo*   p_tsLLinfo;
struct LinkedListInfoStruct{
  int head,nObjs;                     /*  ...index of the thread head and #  */
};                                    /*  ...of objects on this thread       */

/* ------- */
#endif
