
extern int  VH_init(int bias, int nVert, int startIndex, void **vhash);
/*
    bias       - 0 for C counting, 1 for FORTRAN indexing
    nVert      - the number of vertices to hash
    startIndex - the first index to return when adding something new
    vhash      - returned blind pointer to be used for the following
    icode      - 0 OK, 2 malloc error, 4 bad input(s)
 */

extern int  VH_index(void *vhash, int index1, int index2, int *newIndex);
/*
    vhash      - blind pointer from VH_init
    index1     - the first  index in the pair - not larger than nVert-(1-bias)
    index2     - the second index in the pair - not larger than nVert-(1-bias)
    newIndex   - the returned assigned index (starting with startIndex)
    icode      - 0 OK, 1 just added, 2 malloc error, 3 bad blind pointer,
                 4 bad input(s)
 */

extern int  VH_populate(void *vhash, int index1, int index2, int index3);
/*
    vhash      - blind pointer from VH_init
    index1     - the first  index in the pair - not larger than nVert-(1-bias)
    index2     - the second index in the pair - not larger than nVert-(1-bias)
    index3     - the index stored for pair (index1,index2)
    icode      - 0 OK, 1 just added, 2 malloc error, 3 bad blind pointer,
                 4 bad input(s)
 */

extern void VH_free(void **vhash);
/*
    vhash      - blind pointer from VH_init to free up
 */
