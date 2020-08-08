/*
    This file is part of HiParTI!.

    HiParTI! is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    HiParTI! is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with HiParTI!.
    If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARTI_STRUCTS_H
#define PARTI_STRUCTS_H

#include <stdio.h>
#include <cinttypes>
#include "types.h"
#include "macros.h"
#include "errors.h"

/**
 * Dense dynamic array of specified type of scalars
 */
typedef struct {
    sptNnzIndex    len;   /// length
    sptNnzIndex    cap;   /// capacity
    sptValue    *data; /// data
} sptValueVector;

/**
 * Dense dynamic array of different types of integers
 */
typedef struct {
    sptNnzIndex len;   /// length
    sptNnzIndex cap;   /// capacity
    sptIndex *data; /// data
} sptIndexVector;

typedef struct {
    sptNnzIndex len;   /// length
    sptNnzIndex cap;   /// capacity
    sptElementIndex *data; /// data
} sptElementIndexVector;

typedef struct {
    sptNnzIndex len;   /// length
    sptNnzIndex cap;   /// capacity
    sptBlockIndex *data; /// data
} sptBlockIndexVector;

typedef struct {
    sptNnzIndex len;   /// length
    sptNnzIndex cap;   /// capacity
    sptNnzIndex *data; /// data
} sptNnzIndexVector;


/**
 * Dense matrix type
 */
typedef struct {
    sptIndex nrows;   /// # rows
    sptIndex ncols;   /// # columns
    sptIndex cap;     /// # of allocated rows
    sptIndex stride;  /// ncols rounded up to 8
    sptValue *values; /// values, length cap*stride
} sptMatrix;

/**
 * Sparse matrix type, COO format
 */
typedef struct {
    sptIndex nrows;  /// # rows
    sptIndex ncols;  /// # colums
    sptNnzIndex nnz;    /// # non-zeros
    sptIndexVector rowind; /// row indices, length nnz
    sptIndexVector colind; /// column indices, length nnz
    sptValueVector values; /// non-zero values, length nnz
} sptSparseMatrix;


/**
 * Sparse matrix type, CSR format
 */
typedef struct {
    sptIndex nrows;  /// # rows
    sptIndex ncols;  /// # colums
    sptNnzIndex nnz;    /// # non-zeros
    sptNnzIndexVector rowptr; /// row indices, length nnz
    sptIndexVector colind; /// column indices, length nnz
    sptValueVector values; /// non-zero values, length nnz
} sptSparseMatrixCSR;


/**
 * Sparse tensor type, Hierarchical COO format (HiCOO)
 */
typedef struct {
    /* Basic information */
    sptIndex            nrows;  /// # rows
    sptIndex            ncols;  /// # columns
    sptNnzIndex         nnz;         /// # non-zeros

    /* Parameters */
    sptElementIndex       sb_bits;         /// block size by nnz
    sptElementIndex       sk_bits;         /// superblock size by nnz

    /* Index data arrays */
    sptNnzIndexVector         bptr;      /// Block pointers to all nonzeros, nb = bptr.length - 1
    sptBlockIndexVector       bindI;    /// Block indices for rows, length nb
    sptBlockIndexVector       bindJ;    /// Block indices for columns, length nb
    sptElementIndexVector     eindI;    /// Element indices within each block for rows, length nnz
    sptElementIndexVector     eindJ;    /// Element indices within each block for columns, length nnz
    sptValueVector            values;      /// non-zero values, length nnz

    /* Scheduling information */    /// TODO: move scheduler out of HiCOO format
    sptNnzIndexVector         kptr;      /// Nonzero kernel pointers in 1-D array, indexing blocks. sptIndexVector may be enough
    sptIndexVector            *kschr;    /// Kernel scheduler
    sptIndex                  nkiters;     /// max-length of iterations


} sptSparseMatrixHiCOO;


/**
 * Key-value pair structure
 */
typedef struct 
{
  sptIndex key;
  sptIndex value;
} sptKeyValuePair;

char * sptBytesString(uint64_t const bytes)
{
  double size = (double)bytes;
  int suff = 0;
  const char *suffix[5] = {"B", "KiB", "MiB", "GiB", "TiB"};
  while(size > 1024 && suff < 4) {
    size /= 1024.;
    ++suff;
  }
  char * ret = NULL;
  if(asprintf(&ret, "%0.2f %s", size, suffix[suff]) == -1) {
    fprintf(stderr, "SPT: asprintf failed with %" PRIu64 " bytes.\n", bytes);
    ret = NULL;
  }
  return ret;
}

static double sptSparseMatrixDensity(sptSparseMatrix const * const mtx)
{
  double density = (double)mtx->nnz / ((double)mtx->nrows * mtx->ncols);
  return density;
}

void sptSparseMatrixStatus(sptSparseMatrix *mtx, FILE *fp)
{
  fprintf(fp, "COO Sparse Matrix information (use sptIndex, sptValue))---------\n");
  fprintf(fp, " DIMS=%"PARTI_PRI_INDEX "x%"PARTI_PRI_INDEX "\n", mtx->nrows, mtx->ncols);
  fprintf(fp, " NNZ=%"PARTI_PRI_NNZ_INDEX, mtx->nnz);
  fprintf(fp, " DENSITY=%e\n" , sptSparseMatrixDensity(mtx));

  fprintf(fp, " Average row length: %.2lf\n", (double)mtx->nnz / mtx->nrows);
  fprintf(fp, " Average column length: %.2lf\n", (double)mtx->nnz / mtx->ncols);

  char * bytestr = sptBytesString(mtx->nnz * (sizeof(sptIndex) * 2 + sizeof(sptValue)));
  fprintf(fp, " COO-STORAGE=%s\n", bytestr);
  fprintf(fp, "\n");
  free(bytestr);
}

int sptDumpSparseMatrix(const sptSparseMatrix *mtx, sptIndex start_index, FILE *fp) 
{
    int iores;
    sptNnzIndex i;
    iores = fprintf(fp, "%"PARTI_PRI_INDEX, mtx->nrows);
    spt_CheckOSError(iores < 0, "SpMtx Dump");
    iores = fprintf(fp, "x%"PARTI_PRI_INDEX ", ", mtx->ncols);
    spt_CheckOSError(iores < 0, "SpMtx Dump");
    iores = fprintf(fp, "%"PARTI_PRI_NNZ_INDEX "\n", mtx->nnz);
    spt_CheckOSError(iores < 0, "SpMtx Dump");

    for(i = 0; i < mtx->nnz; ++i) {
        iores = fprintf(fp, "%"PARTI_PRI_INDEX "\t", mtx->rowind.data[i] + start_index);
        spt_CheckOSError(iores < 0, "SpMtx Dump");
        iores = fprintf(fp, "%"PARTI_PRI_INDEX "\t", mtx->colind.data[i] + start_index);
        spt_CheckOSError(iores < 0, "SpMtx Dump");
        iores = fprintf(fp, "%"PARTI_PRI_VALUE "\n", (double) mtx->values.data[i]);
        spt_CheckOSError(iores < 0, "SpMtx Dump");
    }
    return 0;
}

/**
 * Dump a dense sptNnzIndex array to file
 *
 * @param array a pointer to a valid sptNnzIndex array
 * @param size of the array
 * @param fp a file pointer
 *
 */
int sptDumpNnzIndexArray(sptNnzIndex const *array, sptNnzIndex const n, FILE *fp) {
    int iores;
    iores = fprintf(fp, "sptNnzIndex array length: %"PARTI_PRI_NNZ_INDEX "\n", n);
    spt_CheckOSError(iores < 0, "NnzIdxArray Dump");
    for(sptNnzIndex i=0; i < n; ++i) {
        iores = fprintf(fp, "%"PARTI_PRI_NNZ_INDEX "\t", array[i]);
        spt_CheckOSError(iores < 0, "NnzIdxArray Dump");
    }
    iores = fprintf(fp, "\n");

    return 0;
}

#endif
