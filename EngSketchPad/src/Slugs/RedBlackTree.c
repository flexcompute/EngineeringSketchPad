/*
 ******************************************************************************
 *                                                                            *
 * create and manage Red-Black Trees                                          *
 *                                                                            *
 * Reference: "Introduction to Algorithms" by Thomas Cormen, Charles Leiserson*
 *            and Ronald Rivest, McGraw-Hill, 1991, pp 244-280.               *
 *                                                                            *
 *            Written by John Dannenhoffer @ Syracuse University              *
 *                                                                            *
 ******************************************************************************
 */

/*
 * Copyright (C) 2013/2022  John F. Dannenhoffer, III (Syracuse University)
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "common.h"
#include "RedBlackTree.h"

/* forward declarations of internal routines defined below */
static int  compareKeys(LONG ikey1, LONG jkey1, LONG ikey2, LONG jkey2, LONG ikey3, LONG jkey3);
static void rotateLeft( rbt_T *tree, int inode);
static void rotateRite( rbt_T *tree, int inode);


/*
 ******************************************************************************
 *                                                                            *
 * rbtCreate - create a new Tree                                              *
 *                                                                            *
 ******************************************************************************
 */
int
rbtCreate(int    chunk,
          rbt_T  **tree)
{

    /* --------------------------------------------------------------- */

    /* default return value */
    *tree = NULL;

    /* get a new Tree and a default number of Nodes */
    *tree = malloc(sizeof(rbt_T));
    if(*tree != NULL) {
        (*tree)->mnode = chunk;
        (*tree)->chunk = chunk;

        (*tree)->key1 = (LONG *) malloc(((*tree)->mnode) * sizeof(LONG));
        (*tree)->key2 = (LONG *) malloc(((*tree)->mnode) * sizeof(LONG));
        (*tree)->key3 = (LONG *) malloc(((*tree)->mnode) * sizeof(LONG));
        (*tree)->left = (int  *) malloc(((*tree)->mnode) * sizeof(int ));
        (*tree)->rite = (int  *) malloc(((*tree)->mnode) * sizeof(int ));
        (*tree)->prnt = (int  *) malloc(((*tree)->mnode) * sizeof(int ));
        (*tree)->colr = (int  *) malloc(((*tree)->mnode) * sizeof(int ));

        /* initialize the Tree */
        (*tree)->nnode =  0;
        (*tree)->root  = -1;
    }

    return 0;
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtDelete - delete an entire Tree                                          *
 *                                                                            *
 ******************************************************************************
 */
int
rbtDelete(rbt_T   *tree)
{

    /* --------------------------------------------------------------- */

    /* check for valid tree */
    if (tree == NULL) {
        printf("rbtDelete called without valid Tree\n");
        exit(0);
    }

    /* free all elements in the Tree and then the Tree itself */
    free(tree->key1);
    free(tree->key2);
    free(tree->key3);
    free(tree->left);
    free(tree->rite);
    free(tree->prnt);
    free(tree->colr);

    return 0;
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtInsert - insert a Node into a Tree                                      *
 *                                                                            *
 ******************************************************************************
 */
int
rbtInsert(rbt_T   *tree,
          LONG    key1,
          LONG    key2,
          LONG    key3)
{
    int ix, iy, iz, inode, ipz, ippz;

    /* --------------------------------------------------------------- */

    /* check for valid Tree */
    if (tree == NULL) {
        printf("rbtInsert called without valid Tree\n");
        exit(0);
    }

    /* expand size of Tree if we have used up all available Nodes */
    if (tree->nnode >= tree->mnode) {
        tree->mnode += tree->chunk;
        tree->key1 = realloc(tree->key1, (tree->mnode) * sizeof(LONG));
        tree->key2 = realloc(tree->key2, (tree->mnode) * sizeof(LONG));
        tree->key3 = realloc(tree->key3, (tree->mnode) * sizeof(LONG));
        tree->left = realloc(tree->left, (tree->mnode) * sizeof(int ));
        tree->rite = realloc(tree->rite, (tree->mnode) * sizeof(int ));
        tree->prnt = realloc(tree->prnt, (tree->mnode) * sizeof(int ));
        tree->colr = realloc(tree->colr, (tree->mnode) * sizeof(int ));
    }

    /* put the new Node at the bottom of the Tree */
    tree->key1[tree->nnode] = key1;
    tree->key2[tree->nnode] = key2;
    tree->key3[tree->nnode] = key3;
    tree->left[tree->nnode] = -1;
    tree->rite[tree->nnode] = -1;
    tree->prnt[tree->nnode] = -1;       /* over-written below */
    tree->colr[tree->nnode] = -1;       /* over-written below */
    tree->nnode ++;

    /* find out where the new Node will fall into the current Tree */
    iz = tree->nnode - 1;
    iy = -1;
    ix = tree->root;

    while (ix >= 0) {
        iy = ix;
        if (compareKeys(key1, tree->key1[ix],
                        key2, tree->key2[ix],
                        key3, tree->key3[ix]) < 0) {
            ix = tree->left[ix];
        } else {
            ix = tree->rite[ix];
        }
    }

    /* link the new Node to its parent and vice versa */

    tree->prnt[iz] = iy;

    if (iy == -1) {
        tree->root = iz;
    } else if (compareKeys(key1, tree->key1[iy],
                           key2, tree->key2[iy],
                           key3, tree->key3[iy]) < 0) {
        tree->left[iy] = iz;
    } else {
        tree->rite[iy] = iz;
    }

    /* make inode (the returned value) point to the new Node */
    inode = iz;

    /* now start reordering the Tree following the red-black algorithm
          so that the Tree is a balanced as possible */
    tree->colr[iz] = RBT_RED;

    /* determine what violations of the red-black properties were
          introduced above */
    while ((iz != tree->root) && (tree->colr[tree->prnt[iz]] == RBT_RED)) {

        /* move a violation of the red-child-black-parent violation up
              the Tree while maintaining that every simple path from
              a Nde to a descendent leaf contains the same number
              of black Nodes */

        ipz  = tree->prnt[iz ];
        ippz = tree->prnt[ipz];

        if (ipz == tree->left[ippz]) {
            iy = tree->rite[ippz];

            if ((iy >= 0) && (tree->colr[iy] == RBT_RED)) {
                tree->colr[ipz ] = RBT_BLACK;
                tree->colr[iy  ] = RBT_BLACK;
                tree->colr[ippz] = RBT_RED;
                iz               = ippz;
            } else {
                if (iz == tree->rite[ipz]) {
                    iz = ipz;
                    rotateLeft(tree, iz);
                }

                ipz  = tree->prnt[iz ];
                ippz = tree->prnt[ipz];

                tree->colr[ipz ] = RBT_BLACK;
                tree->colr[ippz] = RBT_RED;
                rotateRite(tree, ippz);
            }
        } else {
            iy = tree->left[ippz];

            if ((iy >= 0) && (tree->colr[iy] == RBT_RED)) {
                tree->colr[ipz ] = RBT_BLACK;
                tree->colr[iy  ] = RBT_BLACK;
                tree->colr[ippz] = RBT_RED;
                iz               = ippz;
            } else {
                if (iz == tree->left[ipz]) {
                    iz = ipz;
                    rotateRite(tree, iz);
                }

                ipz  = tree->prnt[iz ];
                ippz = tree->prnt[ipz];

                tree->colr[ipz ] = RBT_BLACK;
                tree->colr[ippz] = RBT_RED;
                rotateLeft(tree, ippz);
            }
        }
    }

    /* finally color the Root of the tree black */
    tree->colr[tree->root] = RBT_BLACK;

    return inode;
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtMaximum - find rite-most Node in subTree                                *
 *                                                                            *
 ******************************************************************************
 */
int
rbtMaximum(rbt_T   *tree,
           int     istart)
{
    int ix;

    /* --------------------------------------------------------------- */

    /* check for valid Tree */
    if (tree == NULL) {
        printf("rbtMaximum called without valid Tree\n");
        exit(0);
    }

    /* start at given istart or root if istart < 0 */
    if (istart < 0) {
        ix = tree->root;
    } else {
        ix = istart;
    }

    /* find the ritemost Node by following the rite children */
    while (tree->rite[ix] >= 0) {
        ix = tree->rite[ix];
    }

    return ix;
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtMinimum - find left-most Node in subTree                                *
 *                                                                            *
 ******************************************************************************
 */
int
rbtMinimum(rbt_T   *tree,
           int     istart)
{
    int ix;

    /* --------------------------------------------------------------- */

    /* check for valid Tree */
    if (tree == NULL) {
        printf("rbtMinimum called without valid Tree\n");
        exit(0);
    }

    /* start at given istart or root if istart < 0 */
    if (istart < 0) {
        ix = tree->root;
    } else {
        ix = istart;
    }

    /* find the leftmost Node by following the left children */
    while (tree->left[ix] >= 0) {
        ix = tree->left[ix];
    }

    return ix;
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtNext - find Node immediately to the rite                                *
 *                                                                            *
 ******************************************************************************
 */
int
rbtNext(rbt_T   *tree,
        int     istart)
{
    int ix, iy;

    /* --------------------------------------------------------------- */

    /* check for valid Tree */
    if (tree == NULL) {
        printf("rbtNext called without valid Tree\n");
        exit(0);
    }

    /* start at given istart or root if istart < 0 */
    if (istart < 0) {
        ix = tree->root;
    } else {
        ix = istart;
    }

    /* if the rite child is not empty, then find the minimum of the
          Tree starting at the rite child */
    if (tree->rite[ix] >= 0) {
        return rbtMinimum(tree, tree->rite[ix]);

    /* otherwise, iy is the lowest ancestor of ix whose left child
          is also an ancestor of ix */
    } else {
        iy = tree->prnt[ix];

        while ((iy >= 0) && (ix == tree->rite[iy])) {
            ix = iy;
            iy = tree->prnt[iy];
        }

        return iy;
    }
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtPrev - find Node immediately to the left                                *
 *                                                                            *
 ******************************************************************************
 */
int
rbtPrev(rbt_T   *tree,
        int     istart)
{
    int ix, iy;

    /* --------------------------------------------------------------- */

    /* check for valid Tree */
    if (tree == NULL) {
        printf("rbtPrev called without valid Tree\n");
        exit(0);
    }

    /* start at given istart or root if istart < 0 */
    if (istart < 0) {
        ix = tree->root;
    } else {
        ix = istart;
    }

    /* if the left child is not empty, then find the maximum of the
          Tree starting at the left child */
    if (tree->left[ix] >= 0) {
        return rbtMaximum(tree, tree->left[ix]);

    /* otherwise, iy is the lowest ancestor of ix whose rite child
          is also an ancestor of ix */
    } else {
        iy = tree->prnt[ix];

        while ((iy >= 0) && (ix == tree->left[iy])) {
            ix = iy;
            iy = tree->prnt[iy];
        }

        return iy;
    }
}


/*
 ******************************************************************************
 *                                                                            *
 * rbtSearch - search for Node in Tree                                        *
 *                                                                            *
 ******************************************************************************
 */
int
rbtSearch(rbt_T   *tree,
          LONG    key1,
          LONG    key2,
          LONG    key3)
{
    int ix, ians;

    /* --------------------------------------------------------------- */

    /* check for valid Tree */
    if (tree == NULL) {
        printf("rbtSearch called without valid Tree\n");
        exit(0);
    }

    /* start at the Root of the Tree */
    ix = tree->root;

    /* iteratively descend the Tree, moving left or rite depending on
          the relative location of the key with the keys in the Tree */
    while(ix >= 0) {
        ians = compareKeys(key1, tree->key1[ix],
                           key2, tree->key2[ix],
                           key3, tree->key3[ix]);
        if        (ians < 0) {
            ix = tree->left[ix];
        } else if (ians > 0) {
            ix = tree->rite[ix];
        } else {
            break;
        }
    }

    return ix;
}


/*
 ******************************************************************************
 *                                                                            *
 * compareKeys - compare (ikey1,ikey2,ikey3) with (jkey1,jkey2,jkey3) and return -1/0/+1 *
 *                                                                            *
 ******************************************************************************
 */
static int
compareKeys(LONG   ikey1,
            LONG   jkey1,
            LONG   ikey2,
            LONG   jkey2,
            LONG   ikey3,
            LONG   jkey3)
{

    /* --------------------------------------------------------------- */

    /* compare primary keys */
    if        (ikey1 < jkey1) {
        return -1;
    } else if (ikey1 > jkey1) {
        return +1;

    /* compare secondary keys */
    } else if (ikey2 < jkey2) {
        return -1;
    } else if (ikey2 > jkey2) {
        return +1;

    /* compare tertiary keys */
    } else if (ikey3 < jkey3) {
        return -1;
    } else if (ikey3 > jkey3) {
        return +1;

    /* all keys are the same */
    } else {
        return 0;
    }
}


/*
 ******************************************************************************
 *                                                                            *
 * rotateLeft - rotate Nodes to the left                                      *
 *                                                                            *
 ******************************************************************************
 */
static void
rotateLeft(rbt_T   *tree,
           int     inode)
{
    int ix, iy;

    /* --------------------------------------------------------------- */

    ix = inode;
    iy = tree->rite[ix];

    /* turn iy's left subTree into ix's rite subTree */
    tree->rite[ix] = tree->left[iy];

    if (tree->left[iy] >= 0) {
        tree->prnt[tree->left[iy]] = ix;
    }

    /* link ix's parent to iy */
    tree->prnt[iy] = tree->prnt[ix];

    if (tree->prnt[ix] == -1) {
        tree->root = iy;
    } else if (ix == tree->left[tree->prnt[ix]]) {
        tree->left[tree->prnt[ix]] = iy;
    } else {
        tree->rite[tree->prnt[ix]] = iy;
    }

    /* put ix on iy's left */
    tree->left[iy] = ix;
    tree->prnt[ix] = iy;
}


/*
 ******************************************************************************
 *                                                                            *
 * rotateRite - rotate Nodes to the rite                                      *
 *                                                                            *
 ******************************************************************************
 */
static void
rotateRite(rbt_T   *tree,
           int     inode)
{
    int ix, iy;

    /* --------------------------------------------------------------- */

    ix = inode;
    iy = tree->left[ix];

    /* turn iy's rite subTree into ix's left subTree */
    tree->left[ix] = tree->rite[iy];

    if (tree->rite[iy] >= 0) {
        tree->prnt[tree->rite[iy]] = ix;
    }

    /* link ix's parent to iy */
    tree->prnt[iy] = tree->prnt[ix];

    if (tree->prnt[ix] == -1) {
        tree->root = iy;
    } else if (ix == tree->rite[tree->prnt[ix]]) {
        tree->rite[tree->prnt[ix]] = iy;
    } else {
        tree->left[tree->prnt[ix]] = iy;
    }

    /* put ix on iy's rite */
    tree->rite[iy] = ix;
    tree->prnt[ix] = iy;
}
