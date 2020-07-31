#ifndef VECTORS_H
#define VECTORS_H

// int sptNewIndexVector(sptIndexVector *vec, sptNnzIndex len, sptNnzIndex cap);
// int sptNewValueVector(sptValueVector *vec, sptNnzIndex len, sptNnzIndex cap);
#include "errors.h"
#include "macros.h"
/*
 * Initialize a new sptIndex vector
 *
 * @param vec a valid pointer to an uninitialized sptIndex variable,
 * @param len number of values to create
 * @param cap total number of values to reserve
 *
 * Vector is a type of one-dimentional array with dynamic length
 */
int sptNewIndexVector(sptIndexVector *vec, sptNnzIndex len, sptNnzIndex cap) {
  if(cap < len) {
    cap = len;
  }
  if(cap < 2) {
    cap = 2;
  }
  vec->len = len;
  vec->cap = cap;
  vec->data = (sptIndex *) malloc(cap * sizeof *vec->data);
  spt_CheckOSError(!vec->data, "IdxVec New");
  memset(vec->data, 0, cap * sizeof *vec->data);
  return 0;
}

/**
 * Initialize a new value vector
 *
 * @param vec a valid pointer to an uninitialized sptValueVector variable,
 * @param len number of values to create
 * @param cap total number of values to reserve
 *
 * Vector is a type of one-dimentional array with dynamic length
 */
int sptNewValueVector(sptValueVector *vec, sptNnzIndex len, sptNnzIndex cap) {
  if(cap < len) {
    cap = len;
  }
  if(cap < 2) {
    cap = 2;
  }
  vec->len = len;
  vec->cap = cap;
  vec->data = (sptValue *) malloc(cap * sizeof *vec->data);
  spt_CheckOSError(!vec->data, "ValVec New");
  memset(vec->data, 0, cap * sizeof *vec->data);
  return 0;
}

/**
 * Fill an existed dense value vector with a randomized values
 *
 * @param vec   a valid pointer to an existed sptVector variable,
 * @param val   a given value constant
 *
 * Vector is a type of one-dimentional array with dynamic length
 */
int sptRandomValueVector(sptValueVector * const vec) {
    // srand(time(NULL));
    for(sptNnzIndex i=0; i<vec->len; ++i)
        vec->data[i] = rand() % 10 + 1;
    return 0;
}

/*
 * Initialize a new sptNnzIndexVector vector
 *
 * @param vec a valid pointer to an uninitialized sptNnzIndex variable,
 * @param len number of values to create
 * @param cap total number of values to reserve
 *
 * Vector is a type of one-dimentional array with dynamic length
 */

int sptNewNnzIndexVector(sptNnzIndexVector *vec, sptNnzIndex len, sptNnzIndex cap) {
    if(cap < len) {
        cap = len;
    }
    if(cap < 2) {
        cap = 2;
    }
    vec->len = len;
    vec->cap = cap;
    vec->data = (sptNnzIndex *) malloc(cap * sizeof *vec->data);
    spt_CheckOSError(!vec->data, "NnzIdxVec New");
    memset(vec->data, 0, cap * sizeof *vec->data);
    return 0;
}

/*
 * Initialize a new sptBlockIndexVector vector
 *
 * @param vec a valid pointer to an uninitialized sptBlockIndex variable,
 * @param len number of values to create
 * @param cap total number of values to reserve
 *
 * Vector is a type of one-dimentional array with dynamic length
 */

int sptNewBlockIndexVector(sptBlockIndexVector *vec, sptNnzIndex len, sptNnzIndex cap) {
    if(cap < len) {
        cap = len;
    }
    if(cap < 2) {
        cap = 2;
    }
    vec->len = len;
    vec->cap = cap;
    vec->data = (sptBlockIndex *) malloc(cap * sizeof *vec->data);
    spt_CheckOSError(!vec->data, "BlkIdxVec New");
    memset(vec->data, 0, cap * sizeof *vec->data);
    return 0;
}

/*
 * Initialize a new sptElementIndexVector vector
 *
 * @param vec a valid pointer to an uninitialized sptElementIndex variable,
 * @param len number of values to create
 * @param cap total number of values to reserve
 *
 * Vector is a type of one-dimentional array with dynamic length
 */

int sptNewElementIndexVector(sptElementIndexVector *vec, sptNnzIndex len, sptNnzIndex cap) {
    if(cap < len) {
        cap = len;
    }
    if(cap < 2) {
        cap = 2;
    }
    vec->len = len;
    vec->cap = cap;
    vec->data = (sptElementIndex *) malloc(cap * sizeof *vec->data);
    spt_CheckOSError(!vec->data, "EleIdxVec New");
    memset(vec->data, 0, cap * sizeof *vec->data);
    return 0;
}

/**
 * Release the memory buffer a sptNnzIndexVector is holding
 *
 * @param vec a pointer to a valid long nnz vector
 *
 */
void sptFreeNnzIndexVector(sptNnzIndexVector *vec) {
    free(vec->data);
    vec->len = 0;
    vec->cap = 0;
}

/**
 * Release the memory buffer a sptBlockIndexVector is holding
 *
 * @param vec a pointer to a valid size vector
 *
 */
void sptFreeBlockIndexVector(sptBlockIndexVector *vec) {
    free(vec->data);
    vec->len = 0;
    vec->cap = 0;
}

/**
 * Release the memory buffer a value vector is holding
 *
 * @param vec a pointer to a valid value vector
 *
 */
void sptFreeValueVector(sptValueVector *vec) {
    vec->len = 0;
    vec->cap = 0;
    free(vec->data);
}

/**
 * Release the memory buffer a sptIndexVector is holding
 *
 * @param vec a pointer to a valid size vector
 *
 */
void sptFreeIndexVector(sptIndexVector *vec) {
    free(vec->data);
    vec->len = 0;
    vec->cap = 0;
}

/**
 * Release the memory buffer a sptElementIndexVector is holding
 *
 * @param vec a pointer to a valid size vector
 *
 */
void sptFreeElementIndexVector(sptElementIndexVector *vec) {
    free(vec->data);
    vec->len = 0;
    vec->cap = 0;
}


/**
 * Add a value to the end of a sptNnzIndexVector
 *
 * @param vec   a pointer to a valid long nnz index vector
 * @param value the value to be appended
 *
 * The length of the long nnz index vector will be changed to contain the new value.
 */
int sptAppendNnzIndexVector(sptNnzIndexVector *vec, sptNnzIndex const value) {
    if(vec->cap <= vec->len) {
#ifndef MEMCHECK_MODE
        sptNnzIndex newcap = vec->cap + vec->cap/2;
#else
        sptNnzIndex newcap = vec->len+1;
#endif
        sptNnzIndex *newdata = (sptNnzIndex *) realloc(vec->data, newcap * sizeof *vec->data);
        spt_CheckOSError(!newdata, "NnzIdxVec Append");
        vec->cap = newcap;
        vec->data = newdata;
    }
    vec->data[vec->len] = value;
    ++vec->len;
    return 0;
}

/**
 * Add a value to the end of a sptBlockIndexVector
 *
 * @param vec   a pointer to a valid block index vector
 * @param value the value to be appended
 *
 * The length of the block index vector will be changed to contain the new value.
 */
int sptAppendBlockIndexVector(sptBlockIndexVector *vec, sptBlockIndex const value) {
    if(vec->cap <= vec->len) {
#ifndef MEMCHECK_MODE
        sptNnzIndex newcap = vec->cap + vec->cap/2;
#else
        sptNnzIndex newcap = vec->len+1;
#endif
        sptBlockIndex *newdata = (sptBlockIndex *) realloc(vec->data, newcap * sizeof *vec->data);
        spt_CheckOSError(!newdata, "BlkIdxVec Append");
        vec->cap = newcap;
        vec->data = newdata;
    }
    vec->data[vec->len] = value;
    ++vec->len;
    return 0;
}

/**
 * Add a value to the end of a sptElementIndexVector
 *
 * @param vec   a pointer to a valid element index vector
 * @param value the value to be appended
 *
 * The length of the element index vector will be changed to contain the new value.
 */
int sptAppendElementIndexVector(sptElementIndexVector *vec, sptElementIndex const value) {
    if(vec->cap <= vec->len) {
#ifndef MEMCHECK_MODE
        sptNnzIndex newcap = vec->cap + vec->cap/2;
#else
        sptNnzIndex newcap = vec->len+1;
#endif
        sptElementIndex *newdata = (sptElementIndex *) realloc(vec->data, newcap * sizeof *vec->data);
        spt_CheckOSError(!newdata, "EleIdxVec Append");
        vec->cap = newcap;
        vec->data = newdata;
    }
    vec->data[vec->len] = value;
    ++vec->len;
    return 0;
}

/**
 * Add a value to the end of a value vector
 *
 * @param vec   a pointer to a valid value vector
 * @param value the value to be appended
 *
 * The length of the value vector will be changed to contain the new value.
 */
int sptAppendValueVector(sptValueVector *vec, sptValue const value) {
    if(vec->cap <= vec->len) {
#ifndef MEMCHECK_MODE
        sptNnzIndex newcap = vec->cap + vec->cap/2;
#else
        sptNnzIndex newcap = vec->len+1;
#endif
        sptValue *newdata = (sptValue *) realloc(vec->data, newcap * sizeof *vec->data);
        spt_CheckOSError(!newdata, "ValVec Append");
        vec->cap = newcap;
        vec->data = newdata;
    }
    vec->data[vec->len] = value;
    ++vec->len;
    return 0;
}

void sptFreeSparseMatrix(sptSparseMatrix *mtx) {
    sptFreeIndexVector(&mtx->rowind);
    sptFreeIndexVector(&mtx->colind);
    sptFreeValueVector(&mtx->values);
    mtx->nrows = 0;
    mtx->ncols = 0;
    mtx->nnz = 0;
}

int sptNewSparseMatrixHiCOO(
    sptSparseMatrixHiCOO *himtx, 
    const sptIndex nrows, 
    const sptIndex ncols,
    const sptNnzIndex nnz,
    const sptElementIndex sb_bits,
    const sptElementIndex sk_bits)
{
    int result;

    himtx->nrows = nrows;
    himtx->ncols = ncols;
    himtx->nnz = nnz;

    /* Parameters */
    himtx->sb_bits = sb_bits; // block size by nnz
    himtx->sk_bits = sk_bits; // superblock size by nnz

    result = sptNewNnzIndexVector(&himtx->bptr, 0, 0);
    spt_CheckError(result, "HiSpMtx New", NULL);
    result = sptNewBlockIndexVector(&himtx->bindI, 0, 0);
    spt_CheckError(result, "HiSpMtx New", NULL);
    result = sptNewBlockIndexVector(&himtx->bindJ, 0, 0);
    spt_CheckError(result, "HiSpMtx New", NULL);

    result = sptNewElementIndexVector(&himtx->eindI, 0, 0);
    spt_CheckError(result, "HiSpMtx New", NULL);
    result = sptNewElementIndexVector(&himtx->eindJ, 0, 0);
    spt_CheckError(result, "HiSpMtx New", NULL);

    result = sptNewValueVector(&himtx->values, 0, 0);
    spt_CheckError(result, "HiSpMtx New", NULL);

    /* Allocate superblock scheduler */
    sptIndex sk = (sptIndex)pow(2, sk_bits);
    sptIndex kernel_ndim = (nrows + sk - 1)/sk;
    himtx->kschr = (sptIndexVector*)malloc(kernel_ndim * sizeof(*(himtx->kschr)));
    spt_CheckOSError(!himtx->kschr, "HiSpTns New");
    for(sptIndex i = 0; i < kernel_ndim; ++i) {
        result = sptNewIndexVector(&(himtx->kschr[i]), 0, 0);
        spt_CheckError(result, "HiSpTns New", NULL);
    }
    himtx->nkiters = 0;

    result = sptNewNnzIndexVector(&himtx->kptr, 0, 0);
    spt_CheckError(result, "HiSpTns New", NULL);

    return 0;
}


void sptFreeSparseMatrixHiCOO(sptSparseMatrixHiCOO *himtx)
{
    sptFreeNnzIndexVector(&himtx->bptr);
    sptFreeBlockIndexVector(&himtx->bindI);
    sptFreeBlockIndexVector(&himtx->bindJ);
    sptFreeElementIndexVector(&himtx->eindI);
    sptFreeElementIndexVector(&himtx->eindJ);
    sptFreeValueVector(&himtx->values);

    sptFreeNnzIndexVector(&himtx->kptr);
    sptIndex sk = (sptIndex)pow(2, himtx->sk_bits);
    sptIndex kernel_ndim = (himtx->nrows + sk - 1)/sk;
    for(sptIndex i = 0; i < kernel_ndim; ++i) {
        sptFreeIndexVector(&(himtx->kschr[i]));
    }
    free(himtx->kschr);

    himtx->nnz = 0;
    himtx->nrows = 0;
    himtx->ncols = 0;
    himtx->sb_bits = 0;
    himtx->sk_bits = 0;
    himtx->nkiters = 0;
}

static int sptEqualWithTwoCoordinates(
    const sptIndex * item1,
    const sptIndex * item2,
    const sptIndex nmodes)
{
    sptIndex i1, i2;
    for(sptIndex m=0; m<nmodes; ++m) {
        i1 = item1[m];
        i2 = item2[m];
        if(i1 != i2) {
            return 0;
            break;
        }
    }
    return 1;
}

static int sptBlockEnd(
    sptIndex * out_item,
    sptIndex nmodes,
    sptIndex nrows,
    sptIndex ncols,
    const sptIndex * in_item,
    const sptElementIndex sb)
{
    sptAssert(in_item[0] < nrows);
    out_item[0] = in_item[0]+sb < nrows ? in_item[0]+sb : nrows;    // exclusive
    sptAssert(in_item[1] < ncols);
    out_item[1] = in_item[1]+sb < ncols ? in_item[1]+sb : ncols;    // exclusive

    return 0;
}

static int sptLocateBeginCoord(
    sptIndex * out_item,
    sptIndex nmodes,
    const sptIndex * in_item,
    const sptElementIndex bits)
{
    for(sptIndex m=0; m<nmodes; ++m) {
        out_item[m] = in_item[m] >> bits;
    }

    return 0;
}

int sptSparseMatrixPartition(
    sptSparseMatrixHiCOO *himtx,
    sptNnzIndex *max_nnzb,
    sptSparseMatrix *mtx, 
    const sptElementIndex sb_bits)
{
    int result;
    sptNnzIndex nnz = mtx->nnz;
    sptIndex nrows = mtx->nrows;
    sptIndex ncols = mtx->ncols;
    sptIndex nmodes = 2;
    sptElementIndex sb = pow(2, sb_bits);

    /* Temporary storage */
    sptIndex * block_begin = (sptIndex *)malloc(nmodes * sizeof(*block_begin));
    sptIndex * block_end = (sptIndex *)malloc(nmodes * sizeof(*block_end));
    sptIndex * block_begin_prior = (sptIndex *)malloc(nmodes * sizeof(*block_begin_prior));
    sptIndex * block_coord = (sptIndex *)malloc(nmodes * sizeof(*block_coord));

    sptNnzIndex k_begin, k_end; // #Nonzeros locations
    sptNnzIndex nk = 0; // #Kernels  
    sptNnzIndex nb = 1; // #Blocks  // counting from the first nnz
    sptNnzIndex nb_tmp = 0;
    sptNnzIndex ne = 0; // #Nonzeros per block
    sptIndex eindex = 0;

    /* different appending methods:
     * elements: append every nonzero entry
     * blocks: append when seeing a new block.
     * chunks: appending when seeting a new chunk. Notice the boundary of kernels and the last chunk of the whole tensor may be larger than the sc.
     * kernels: append when seeing a new kernel. Not appending a vector, just write data into an allocated array.
     */
    /* Process first nnz */
    block_coord[0] = mtx->rowind.data[0];    // first nonzero indices
    block_coord[1] = mtx->colind.data[0];    // first nonzero indices
    result = sptLocateBeginCoord(block_begin_prior, nmodes, block_coord, sb_bits);
    spt_CheckError(result, "HiSpMtx Convert", NULL);
    sptAppendBlockIndexVector(&himtx->bindI, (sptBlockIndex)block_begin_prior[0]);
    sptAppendBlockIndexVector(&himtx->bindJ, (sptBlockIndex)block_begin_prior[1]);
    sptAppendNnzIndexVector(&himtx->bptr, 0);

    /* Loop for all kernels, 0 - himtx->kptr.len - 1 */
    for(sptNnzIndex k=0; k<himtx->kptr.len - 1; ++k) {
        k_begin = himtx->kptr.data[k];
        k_end = himtx->kptr.data[k+1]; // exclusive
        nb_tmp = k == 0 ? 0: nb;
        /* Modify kptr pointing to block locations */
        himtx->kptr.data[k] = nb_tmp;
        ++ nk;

        /* Loop nonzeros in each kernel */
        for(sptNnzIndex z = k_begin; z < k_end; ++z) {
            // printf("z: %"PARTI_PRI_NNZ_INDEX "\n", z);

            block_coord[0] = mtx->rowind.data[z];    // first nonzero indices
            block_coord[1] = mtx->colind.data[z];    // first nonzero indices
            // printf("block_coord:\n");
            // sptAssert(sptDumpIndexArray(block_coord, nmodes, stdout) == 0);

            result = sptLocateBeginCoord(block_begin, nmodes, block_coord, sb_bits);
            spt_CheckError(result, "HiSpMtx Convert", NULL);
            // printf("block_begin_prior:\n");
            // sptAssert(sptDumpIndexArray(block_begin_prior, nmodes, stdout) == 0);
            // printf("block_begin:\n");
            // sptAssert(sptDumpIndexArray(block_begin, nmodes, stdout) == 0);

            result = sptBlockEnd(block_end, nmodes, nrows, ncols, block_begin, sb);  // exclusive
            spt_CheckError(result, "HiSpMtx Convert", NULL);

            /* Append einds and values */
            eindex = mtx->rowind.data[z] < (block_begin[0] << sb_bits) ? mtx->rowind.data[z] : mtx->rowind.data[z] - (block_begin[0] << sb_bits);
            sptAssert(eindex < sb);
            sptAppendElementIndexVector(&himtx->eindI, (sptElementIndex)eindex);
            eindex = mtx->colind.data[z] < (block_begin[1] << sb_bits) ? mtx->colind.data[z] : mtx->colind.data[z] - (block_begin[1] << sb_bits);
            sptAssert(eindex < sb);
            sptAppendElementIndexVector(&himtx->eindJ, (sptElementIndex)eindex);
            sptAppendValueVector(&himtx->values, mtx->values.data[z]);


            /* z in the same block with last z */
            if (sptEqualWithTwoCoordinates(block_begin, block_begin_prior, nmodes) == 1)
            {
                /* ne: #Elements in current block */
                ++ ne;
            } else { /* New block */
                /* ne: #Elements in the last block */
                /* Append block bptr and bidx */
                sptAppendNnzIndexVector(&himtx->bptr, (sptBlockIndex)z);
                sptAppendBlockIndexVector(&himtx->bindI, (sptBlockIndex)block_begin[0]);
                sptAppendBlockIndexVector(&himtx->bindJ, (sptBlockIndex)block_begin[1]);
                for(sptIndex m=0; m<nmodes; ++m)
                    block_begin_prior[m] = block_begin[m];

                ++ nb;
                ne = 1;              
            } // End new block
            // printf("nb: %"PARTI_PRI_NNZ_INDEX ", ne: %"PARTI_PRI_NNZ_INDEX "\n\n", nb, ne);

        }   // End z loop
    }   // End k loop
    
    sptAssert(nb <= nnz);
    sptAssert(nb == himtx->bindI.len); 
    sptAssert(nk == himtx->kptr.len - 1);

    /* Last element for kptr, bptr */
    himtx->kptr.data[himtx->kptr.len - 1] = himtx->bptr.len;    // change kptr pointing to blocks
    sptAppendNnzIndexVector(&himtx->bptr, nnz);

    *max_nnzb = himtx->bptr.data[1] - himtx->bptr.data[0];
    sptNnzIndex sum_nnzb = 0;
    for(sptIndex i=0; i < himtx->bptr.len - 1; ++i) {
        sptNnzIndex nnzb = himtx->bptr.data[i+1] - himtx->bptr.data[i];
        sum_nnzb += nnzb;
        if(*max_nnzb < nnzb) {
          *max_nnzb = nnzb;
        }
    }
    sptAssert(sum_nnzb == himtx->nnz);

    free(block_begin);
    free(block_end);
    free(block_begin_prior);
    free(block_coord);
    return 0;
}


static const uint32_t MASKS[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
static const uint32_t SHIFTS[] = {1, 2, 4, 8};

void spt_SwapValues(sptSparseMatrix *mtx, sptNnzIndex ind1, sptNnzIndex ind2) {

    sptIndex eleind1;
    eleind1 = mtx->rowind.data[ind1];
    mtx->rowind.data[ind1] = mtx->rowind.data[ind2];
    mtx->rowind.data[ind2] = eleind1;
    eleind1 = mtx->colind.data[ind1];
    mtx->colind.data[ind1] = mtx->colind.data[ind2];
    mtx->colind.data[ind2] = eleind1;

    sptValue val1 = mtx->values.data[ind1];
    mtx->values.data[ind1] = mtx->values.data[ind2];
    mtx->values.data[ind2] = val1;
}

/* Compare functions */
int spt_SparseMatrixCompareIndicesMorton2D(
    sptSparseMatrix * const mtx1, 
    uint64_t loc1, 
    sptSparseMatrix * const mtx2, 
    uint64_t loc2,
    sptElementIndex sb_bits)
{
    uint64_t mkey1 = 0, mkey2 = 0;
    
    /* Only support 3-D tensors, with 32-bit indices. */
    uint32_t x1 = mtx1->rowind.data[loc1];
    uint32_t y1 = mtx1->colind.data[loc1];
    uint32_t x2 = mtx2->rowind.data[loc2];
    uint32_t y2 = mtx2->colind.data[loc2];

    /* Compare block indices */
    sptIndex blk_x1 = x1 >> sb_bits;
    sptIndex blk_y1 = y1 >> sb_bits;
    sptIndex blk_x2 = x2 >> sb_bits;
    sptIndex blk_y2 = y2 >> sb_bits;

    if(blk_x1 < blk_x2) {
        return -1;
    } else if(blk_x1 > blk_x2) {
        return 1;
    } else if(blk_y1 < blk_y2) {  // if blk_x1 == blk_x2
        return -1;
    } else if(blk_y1 > blk_y2) {  // if blk_x1 == blk_x2
        return 1;
    }

    /* blk_x1 == blk_x2, blk_y1 == blk_y2, sort inside a block in Z-Morton order */
    uint64_t x = x1 - (blk_x1 << sb_bits);
    uint64_t y = y1 - (blk_y1 << sb_bits);
    x = (x | (x << SHIFTS[3])) & MASKS[3];
    x = (x | (x << SHIFTS[2])) & MASKS[2];
    x = (x | (x << SHIFTS[1])) & MASKS[1];
    x = (x | (x << SHIFTS[0])) & MASKS[0];
    y = (y | (y << SHIFTS[3])) & MASKS[3];
    y = (y | (y << SHIFTS[2])) & MASKS[2];
    y = (y | (y << SHIFTS[1])) & MASKS[1];
    y = (y | (y << SHIFTS[0])) & MASKS[0];
    mkey1 = y | (x << 1);

    x = x2 - (blk_x2 << sb_bits);
    y = y2 - (blk_y2 << sb_bits);
    x = (x | (x << SHIFTS[3])) & MASKS[3];
    x = (x | (x << SHIFTS[2])) & MASKS[2];
    x = (x | (x << SHIFTS[1])) & MASKS[1];
    x = (x | (x << SHIFTS[0])) & MASKS[0];
    y = (y | (y << SHIFTS[3])) & MASKS[3];
    y = (y | (y << SHIFTS[2])) & MASKS[2];
    y = (y | (y << SHIFTS[1])) & MASKS[1];
    y = (y | (y << SHIFTS[0])) & MASKS[0];
    mkey2 = y | (x << 1);

    if(mkey1 < mkey2) {
        return -1;
    } else if(mkey1 > mkey2) {
        return 1;
    } else {
        return 0;
    }
    
}


int spt_SparseTensorCompareIndicesSingleMode(sptSparseMatrix * const mtx1, sptNnzIndex loc1, sptSparseMatrix * const mtx2, sptNnzIndex loc2, sptIndex const mode)
{
    sptIndex eleind1, eleind2;
    if (mode == 0) {
        eleind1 = mtx1->rowind.data[loc1];
        eleind2 = mtx2->rowind.data[loc2];
    } else if (mode == 1) {
        eleind1 = mtx1->colind.data[loc1];
        eleind2 = mtx2->colind.data[loc2];
    }
    if(eleind1 < eleind2) {
        return -1;
    } else if(eleind1 > eleind2) {
        return 1;
    }

    return 0;
}


int spt_SparseMatrixCompareIndicesRowBlock(
    sptSparseMatrix * const mtx1, 
    sptNnzIndex loc1, 
    sptSparseMatrix * const mtx2, 
    sptNnzIndex loc2,
    sptElementIndex sk_bits)
{
    sptIndex eleind1 = mtx1->rowind.data[loc1];
    sptIndex eleind2 = mtx2->rowind.data[loc2];
    sptIndex blkind1 = eleind1 >> sk_bits;
    sptIndex blkind2 = eleind2 >> sk_bits;
    // printf("blkind1: %lu, blkind2: %lu\n", blkind1, blkind2);

    if(blkind1 < blkind2) {
        return -1;
    } else if(blkind1 > blkind2) {
        return 1;
    } 

    eleind1 = mtx1->colind.data[loc1];
    eleind2 = mtx2->colind.data[loc2];
    blkind1 = eleind1 >> sk_bits;
    blkind2 = eleind2 >> sk_bits;

    if(blkind1 < blkind2) {
        return -1;
    } else if(blkind1 > blkind2) {
        return 1;
    } 

    return 0;
}


/* Quick sort functions */
static void spt_QuickSortIndexMorton2D(sptSparseMatrix *mtx, sptNnzIndex l, sptNnzIndex r, sptElementIndex sb_bits)
{

    uint64_t i, j, p;
    if(r-l < 2) {
        return;
    }
    p = (l+r) / 2;
    for(i = l, j = r-1; ; ++i, --j) {
        while(spt_SparseMatrixCompareIndicesMorton2D(mtx, i, mtx, p, sb_bits) < 0) {
            // printf("(%lu, %lu) result: %d\n", i, p, spt_SparseMatrixCompareIndicesMorton2D(mtx, i, mtx, p));
            ++i;
        }
        while(spt_SparseMatrixCompareIndicesMorton2D(mtx, p, mtx, j, sb_bits) < 0) {
            // printf("(%lu, %lu) result: %d\n", p, j,spt_SparseMatrixCompareIndicesMorton2D(mtx, p, mtx, j));
            --j;
        }
        if(i >= j) {
            break;
        }
        spt_SwapValues(mtx, i, j);
        if(i == p) {
            p = j;
        } else if(j == p) {
            p = i;
        }
    }
    #pragma omp task firstprivate(l,i) shared(mtx)
    {
        spt_QuickSortIndexMorton2D(mtx, l, i, sb_bits);
    }
    spt_QuickSortIndexMorton2D(mtx, i, r, sb_bits);
    #pragma omp taskwait 
    
}


static void spt_QuickSortIndexSingleMode(sptSparseMatrix *mtx, sptNnzIndex l, sptNnzIndex r, sptIndex mode) 
{
    sptNnzIndex i, j, p;
    if(r-l < 2) {
        return;
    }
    p = (l+r) / 2;
    for(i = l, j = r-1; ; ++i, --j) {
        while(spt_SparseTensorCompareIndicesSingleMode(mtx, i, mtx, p, mode) < 0) {
            ++i;
        }
        while(spt_SparseTensorCompareIndicesSingleMode(mtx, p, mtx, j, mode) < 0) {
            --j;
        }
        if(i >= j) {
            break;
        }
        spt_SwapValues(mtx, i, j);
        if(i == p) {
            p = j;
        } else if(j == p) {
            p = i;
        }
    }
    #pragma omp task firstprivate(l,i) shared(mtx, mode)
    {
        spt_QuickSortIndexSingleMode(mtx, l, i, mode);
    }
    spt_QuickSortIndexSingleMode(mtx, i, r, mode);
    #pragma omp taskwait
}


static void spt_QuickSortIndexRowBlock(sptSparseMatrix *mtx, sptNnzIndex l, sptNnzIndex r,  sptElementIndex sk_bits) 
{
    sptNnzIndex i, j, p;
    if(r-l < 2) {
        return;
    }
    p = (l+r) / 2;
    for(i = l, j = r-1; ; ++i, --j) {
        while(spt_SparseMatrixCompareIndicesRowBlock(mtx, i, mtx, p, sk_bits) < 0) {
            ++i;
        }
        while(spt_SparseMatrixCompareIndicesRowBlock(mtx, p, mtx, j, sk_bits) < 0) {
            --j;
        }
        if(i >= j) {
            break;
        }
        spt_SwapValues(mtx, i, j);
        if(i == p) {
            p = j;
        } else if(j == p) {
            p = i;
        }
    }
    #pragma omp task firstprivate(l,i) shared(mtx, sk_bits)
    {
        spt_QuickSortIndexRowBlock(mtx, l, i, sk_bits);
    }
    spt_QuickSortIndexRowBlock(mtx, i, r, sk_bits);
    #pragma omp taskwait
}


/****************************
 * Sorting functions
 ****************************/
void sptSparseMatrixSortIndexMorton(
    sptSparseMatrix *mtx, 
    int force,
    sptNnzIndex begin,
    sptNnzIndex end,
    sptElementIndex sb_bits) 
{
    if(force) {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                spt_QuickSortIndexMorton2D(mtx, begin, end, sb_bits);
            }
        }
    }
}

void sptSparseMatrixSortIndexSingleMode(sptSparseMatrix *mtx, int force, sptIndex mode, int tk) 
{
    if(force) {
        #pragma omp parallel num_threads(tk) 
        {
            #pragma omp single nowait 
            {
                spt_QuickSortIndexSingleMode(mtx, 0, mtx->nnz, mode);
            }
        }
    }
}

/**
 * Reorder the elements in a COO sparse matrix lexicographically, sorting by row major order.
 * @param mtx  the sparse matrix to operate on
 */
void sptSparseMatrixSortIndexRowBlock(
    sptSparseMatrix *mtx, 
    int force,
    sptNnzIndex begin,
    sptNnzIndex end,
    sptElementIndex sk_bits) 
{
    if(force) {
        #pragma omp parallel 
        {
            #pragma omp single nowait
            {
                spt_QuickSortIndexRowBlock(mtx, begin, end, sk_bits);
            }
        }
    }
}

/**
 * Randomly shuffle all indices.
 *
 * @param[in] mtx matrix to be shuffled
 * @param[out] map_inds records the randomly generated mapping
 *
 */
void sptGetRandomShuffledIndices(sptSparseMatrix *mtx, sptIndex ** map_inds)
{
    /* Get randomly renumbering indices */
    for(sptIndex m = 0; m < 2; ++m) {
        sptIndex dim_len;
        if (m == 0) dim_len = mtx->nrows;
        else if (m == 1) dim_len = mtx->ncols;
        
        for(long int i = dim_len - 1; i > 0; --i) {
            srand(m+i+1+time(NULL));
            sptIndex new_loc = (sptIndex) (rand() % (i+1));            
            /* Swap i <-> new_loc */
            sptIndex tmp = map_inds[m][i];
            map_inds[m][i] = map_inds[m][new_loc];
            map_inds[m][new_loc] = tmp;
        }
    }
}

#endif
