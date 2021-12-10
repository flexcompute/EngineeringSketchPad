/*
 ************************************************************************
 *                                                                      *
 * RedBlackTree.h - header for RedBlackTree.c                           *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2021  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

#ifndef _REDBLACKTREE_H_
#define _REDBLACKTREE_H_

/*
 ************************************************************************
 *                                                                      *
 * Structures                                                           *
 *                                                                      *
 ************************************************************************
 */

/*
 * Tree    - a list of Nodes arranged into a balanced binary Tree
 *              that is collored Red and Black
 *
 * subTree - a Tree that starts other than at the Root
 *
 * Root    - the top of the Tree
 *
 * Node    - the elements on the Tree
 */

#define LONG  long long

typedef struct {
    int    nnode;                       /* current number of Nodes */
    int    mnode;                       /* maximum number of Nodes */
    int    root;                        /* index of root Node */
    int    chunk;                       /* chunk size */
    LONG*  key1;                        /* array of primary   keys */
    LONG*  key2;                        /* array of secondary keys */
    LONG*  key3;                        /* array of tertiary  keys */
    int*   left;                        /* array of left children */
    int*   rite;                        /* array of rite children */
    int*   prnt;                        /* array of parents */
    int*   colr;                        /* array of colors */
} rbt_T;

/*
 ************************************************************************
 *                                                                      *
 * Defined constants                                                    *
 *                                                                      *
 ************************************************************************
 */

#define RBT_BLACK 0
#define RBT_RED   1

/*
 ************************************************************************
 *                                                                      *
 * Callable routines                                                    *
 *                                                                      *
 ************************************************************************
 */

/* create an enpty red-black Tree */
int rbtCreate(int    chunk,            /* (in)  chunk size (for allocations ) */
              rbt_T  **tree);          /* (out) pointer to new RBT */

/* delete a Tree */
int rbtDelete(rbt_T  *tree);           /* (in)  pointer to RBT */

/* insert a Node into the Tree (and return its index) */
int rbtInsert(rbt_T  *tree,            /* (in)  pointer to RBT */
              LONG   key1,             /* (in)  first  key */
              LONG   key2,             /* (in)  second key */
              LONG   key3);            /* (in)  third  key */

/* find the ritemost Node in a subTree */
int rbtMaximum(rbt_T  *tree,           /* (in)  pointer to RBT */
               int    istart);         /* (in)  Node at which to start search */

/* find the leftmost Node in a subTree */
int rbtMinimum(rbt_T  *tree,           /* (in)  pointer to RBT */
               int    istart);         /* (in)  Node at which to start search */

/* find the next (on the right) Node in a Tree */
int rbtNext(rbt_T  *tree,              /* (in)  pointer to BRT */
            int    istart);            /* (in)  Node at which to start */

/* find the previous (on the left) Node in a Tree */
int rbtPrev(rbt_T  *tree,              /* (in)  pointer to RBT */
            int    istart);            /* (in)  Node at which to start */

/* find a Node in a Tree (or -1) */
int rbtSearch(rbt_T  *tree,            /* (in)  pointer to RBT */
              LONG   key1,             /* (in)  first  key */
              LONG   key2,             /* (in)  second key */
              LONG   key3);            /* (in)  third  key */

#endif  /* _REDBLACKTREE_H_ */
