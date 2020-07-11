#ifndef _MMIO_HIGHLEVEL_
#define _MMIO_HIGHLEVEL_

#include "mmio.h"

#define STATS_ON
#define PRINT_MAP

#ifdef STATS_ON	
#include <stdint.h>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#endif

// read matrix infomation from mtx file
int mmio_info(int *m, int *n, int *nnz, int *isSymmetric, const char *filename,
	      uint64_t blkSize = 16)
{
    int m_tmp, n_tmp, nnz_tmp;

    int ret_code;
    MM_typecode matcode;
    FILE *f;

    int nnz_mtx_report;
    int isInteger = 0, isReal = 0, isPattern = 0, isSymmetric_tmp = 0, isComplex = 0;

    // load matrix
    if ((f = fopen(filename, "r")) == NULL)
        return -1;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return -2;
    }

    if ( mm_is_pattern( matcode ) )  { isPattern = 1; /*printf("type = Pattern\n");*/ }
    if ( mm_is_real ( matcode) )     { isReal = 1; /*printf("type = real\n");*/ }
    if ( mm_is_complex( matcode ) )  { isComplex = 1; /*printf("type = real\n");*/ }
    if ( mm_is_integer ( matcode ) ) { isInteger = 1; /*printf("type = integer\n");*/ }

    /* find out size of sparse matrix .... */
    ret_code = mm_read_mtx_crd_size(f, &m_tmp, &n_tmp, &nnz_mtx_report);
    if (ret_code != 0)
        return -4;

    if ( mm_is_symmetric( matcode ) || mm_is_hermitian( matcode ) )
    {
        isSymmetric_tmp = 1;
        //printf("input matrix is symmetric = true\n");
    }
    else
    {
        //printf("input matrix is symmetric = false\n");
    }

    int *csrRowPtr_counter = (int *)malloc((m_tmp+1) * sizeof(int));
    memset(csrRowPtr_counter, 0, (m_tmp+1) * sizeof(int));

    int *csrRowIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));
    int *csrColIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));
    double *csrVal_tmp    = (double *)malloc(nnz_mtx_report * sizeof(double));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    std::map<std::pair<uint64_t, uint64_t>, uint64_t > nonzerosPerBlockMap;
    std::map<std::pair<uint64_t, uint64_t>, double > PerBlockDensityMap;
    std::map<std::pair<uint64_t, uint64_t>, std::pair<int, int> >
      nnzRowColCountPerBlockMap;
    std::map<std::pair<uint64_t, uint64_t>, std::set<int> > nnzRowsPerBlockMap;
    std::map<std::pair<uint64_t, uint64_t>, std::set<int> > nnzColsPerBlockMap;
    auto blockSize = blkSize;
    
    for (int i = 0; i < nnz_mtx_report; i++)
    {
        int idxi, idxj;
        double fval, fval_im;
        int ival;
        int returnvalue;

        if (isReal)
        {
            returnvalue = fscanf(f, "%d %d %lg\n", &idxi, &idxj, &fval);
        }
        else if (isComplex)
        {
            returnvalue = fscanf(f, "%d %d %lg %lg\n", &idxi, &idxj, &fval, &fval_im);
        }
        else if (isInteger)
        {
            returnvalue = fscanf(f, "%d %d %d\n", &idxi, &idxj, &ival);
            fval = ival;
        }
        else if (isPattern)
        {
            returnvalue = fscanf(f, "%d %d\n", &idxi, &idxj);
            fval = 1.0;
        }

        // adjust from 1-based to 0-based
        idxi--;
        idxj--;

        csrRowPtr_counter[idxi]++;
        csrRowIdx_tmp[i] = idxi;
        csrColIdx_tmp[i] = idxj;
        csrVal_tmp[i] = fval;

#ifdef STATS_ON	
	auto blockIdxRowNum = floor (idxi / blockSize);
	auto blockIdxColNum = floor (idxj / blockSize);

	auto search = nonzerosPerBlockMap.find(std::make_pair(blockIdxRowNum,
							      blockIdxColNum));
	if (search != nonzerosPerBlockMap.end()) {
	  auto currentBlockNNZ = search->second;
	  search->second = currentBlockNNZ + 1;
	  // Needs c++17
	  // nonzerosPerBlockMap.insert_or_assign(std::make_pair(blockIdxRowNum,
	  //				       blockIdxColNum),	currentBlockNNZ + 1);
	} else {
	  nonzerosPerBlockMap.emplace(std::make_pair(blockIdxRowNum, blockIdxColNum),
				     1);
	}

	//record nnzRowIdx and nnzColIdx
	auto foundRow = nnzRowsPerBlockMap.find(std::make_pair(blockIdxRowNum,
							    blockIdxColNum));
	if (foundRow != nnzRowsPerBlockMap.end()) {
	  auto rowIdxSet = foundRow->second;
	  rowIdxSet.insert(blockIdxRowNum);
	} else {
	  std::set<int> tempSet;
	  tempSet.insert(blockIdxRowNum);
	  nnzRowsPerBlockMap.emplace(std::make_pair(blockIdxRowNum, blockIdxColNum), tempSet);
	}
	
	auto foundCol = nnzColsPerBlockMap.find(std::make_pair(blockIdxRowNum,
							    blockIdxColNum));
	if (foundCol != nnzColsPerBlockMap.end()) {
	  auto colIdxSet = foundCol->second;
	  colIdxSet.insert(blockIdxColNum);
	} else {
	  std::set<int> tempSet;
	  tempSet.insert(blockIdxColNum);
	  nnzColsPerBlockMap.emplace(std::make_pair(blockIdxRowNum, blockIdxColNum), tempSet);
	}
    }

    for (auto &p : nonzerosPerBlockMap) {
      auto r  = p.first.first;
      auto c = p.first.second;
      // Get nnzRowCount for this block
      auto nnzRowCount = nnzRowsPerBlockMap.find(std::make_pair(r, c))->second.size();
      // Get nnzColCount for this block
      auto nnzColCount = nnzColsPerBlockMap.find(std::make_pair(r, c))->second.size();
      nnzRowColCountPerBlockMap.emplace(std::make_pair(r,c ),
    					std::make_pair(nnzRowCount, nnzColCount));
    }

#ifdef PRINT_MAP
    std::cout << "Printing total number of non-zero rows and cols per block" << std::endl;
    for (const auto &p : nnzRowColCountPerBlockMap) {
      std::cout << "(" << p.first.first << "," << p.first.second << ")" << "-->(" << p.second.first << "," << p.second.second << ")" << std::endl;
    }
#endif

    

#ifdef PRINT_MAP
    std::cout << "Printing total number of non-zeros per block" << std::endl;
    for (const auto &p : nonzerosPerBlockMap) {
      std::cout << "(" << p.first.first << "," << p.first.second << ")" << "-->" << p.second << std::endl;
    }
#endif

    // compute per-block density
    for (const auto &p : nonzerosPerBlockMap) {
      PerBlockDensityMap.emplace(std::make_pair(p.first.first, p.first.second),
				 (double) p.second / (double) (blockSize * blockSize));
    }

#ifdef PRINT_MAP
    std::cout << "Printing per block density" << std::endl;
    for (const auto &p : PerBlockDensityMap) {
      std::cout << "(" << p.first.first << "," << p.first.second << ")" << "-->" << p.second << std::endl;
    }
#endif    
    auto maxBlkNum = ceil (m_tmp / blockSize);
    auto totalNumZeroBlocks = 0;
    std::vector<std::pair<uint64_t, uint64_t> > zeroBlkIndices;
    for (auto i = 0; i < maxBlkNum; i ++) {
      for (auto j = 0; j < maxBlkNum; j++) {
	auto search = nonzerosPerBlockMap.find(std::make_pair(i, j));
	if (search == nonzerosPerBlockMap.end()) {
	  zeroBlkIndices.push_back(std::make_pair(i,j));
	  totalNumZeroBlocks++;
	}
      }
    }

    std::cout << "Total number of non-zero blocks: " << nonzerosPerBlockMap.size() << " and total number of zero-blocks: " << totalNumZeroBlocks << std::endl;

    std::vector<uint64_t> histogram((blockSize + 1) * (blockSize + 1), 0);
    for (const auto &blockStat : nonzerosPerBlockMap) {
      histogram.at(blockStat.second) = histogram.at(blockStat.second) + 1;
    }

    for (auto i = 0; i < histogram.size(); i++) {
      if (histogram.at(i) != 0) {
	std::cout << "No of blocks containing elements " << i << " is "  << histogram[i] << std::endl;
      }
    }
#endif
    
    if (f != stdin)
        fclose(f);

    if (isSymmetric_tmp)
    {
        for (int i = 0; i < nnz_mtx_report; i++)
        {
            if (csrRowIdx_tmp[i] != csrColIdx_tmp[i])
                csrRowPtr_counter[csrColIdx_tmp[i]]++;
        }
    }

    // exclusive scan for csrRowPtr_counter
    int old_val, new_val;

    old_val = csrRowPtr_counter[0];
    csrRowPtr_counter[0] = 0;
    for (int i = 1; i <= m_tmp; i++)
    {
        new_val = csrRowPtr_counter[i];
        csrRowPtr_counter[i] = old_val + csrRowPtr_counter[i-1];
        old_val = new_val;
    }

    nnz_tmp = csrRowPtr_counter[m_tmp];

    *m = m_tmp;
    *n = n_tmp;
    *nnz = nnz_tmp;
    *isSymmetric = isSymmetric_tmp;

    // free tmp space
    free(csrColIdx_tmp);
    free(csrVal_tmp);
    free(csrRowIdx_tmp);
    free(csrRowPtr_counter);

    return 0;
}

// read matrix infomation from mtx file
int mmio_data(int *csrRowPtr, int *csrColIdx, double *csrVal, const char *filename)
{
    int m_tmp, n_tmp, nnz_tmp;

    int ret_code;
    MM_typecode matcode;
    FILE *f;

    int nnz_mtx_report;
    int isInteger = 0, isReal = 0, isPattern = 0, isSymmetric_tmp = 0, isComplex = 0;

    // load matrix
    if ((f = fopen(filename, "r")) == NULL)
    {
        printf("Error loading matrix file.\n");
        return -1;
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return -2;
    }

    if ( mm_is_pattern( matcode ) )  { isPattern = 1; /*printf("type = Pattern\n");*/ }
    if ( mm_is_real ( matcode) )     { isReal = 1; /*printf("type = real\n");*/ }
    if ( mm_is_complex( matcode ) )  { isComplex = 1; /*printf("type = real\n");*/ }
    if ( mm_is_integer ( matcode ) ) { isInteger = 1; /*printf("type = integer\n");*/ }

    /* find out size of sparse matrix .... */
    ret_code = mm_read_mtx_crd_size(f, &m_tmp, &n_tmp, &nnz_mtx_report);
    if (ret_code != 0)
        return -4;


    if ( mm_is_symmetric( matcode ) || mm_is_hermitian( matcode ) )
    {
        isSymmetric_tmp = 1;
        //printf("input matrix is symmetric = true\n");
    }
    else
    {
        //printf("input matrix is symmetric = false\n");
    }

    
    int *csrRowPtr_counter = (int *)malloc((m_tmp+1) * sizeof(int));
    memset(csrRowPtr_counter, 0, (m_tmp+1) * sizeof(int));

    int *csrRowIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));
    int *csrColIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));
    double *csrVal_tmp    = (double*)malloc(nnz_mtx_report * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (int i = 0; i < nnz_mtx_report; i++)
    {
        int idxi, idxj;
        double fval, fval_im;
        int ival;
        int returnvalue;

        if (isReal)
        {
            returnvalue = fscanf(f, "%d %d %lg\n", &idxi, &idxj, &fval);
        }
        else if (isComplex)
        {
            returnvalue = fscanf(f, "%d %d %lg %lg\n", &idxi, &idxj, &fval, &fval_im);
        }
        else if (isInteger)
        {
            returnvalue = fscanf(f, "%d %d %d\n", &idxi, &idxj, &ival);
            fval = ival;
        }
        else if (isPattern)
        {
            returnvalue = fscanf(f, "%d %d\n", &idxi, &idxj);
            fval = 1.0;
        }

        // adjust from 1-based to 0-based
        idxi--;
        idxj--;

        csrRowPtr_counter[idxi]++;
        csrRowIdx_tmp[i] = idxi;
        csrColIdx_tmp[i] = idxj;
        csrVal_tmp[i] = fval;
    }


    if (f != stdin)
        fclose(f);

    if (isSymmetric_tmp)
    {
        for (int i = 0; i < nnz_mtx_report; i++)
        {
            if (csrRowIdx_tmp[i] != csrColIdx_tmp[i])
                csrRowPtr_counter[csrColIdx_tmp[i]]++;
        }
    }


    // exclusive scan for csrRowPtr_counter
    int old_val, new_val;

    old_val = csrRowPtr_counter[0];
    csrRowPtr_counter[0] = 0;
    for (int i = 1; i <= m_tmp; i++)
    {
        new_val = csrRowPtr_counter[i];
        csrRowPtr_counter[i] = old_val + csrRowPtr_counter[i-1];
        old_val = new_val;
    }


    nnz_tmp = csrRowPtr_counter[m_tmp];
    memcpy(csrRowPtr, csrRowPtr_counter, (m_tmp+1) * sizeof(int));
    memset(csrRowPtr_counter, 0, (m_tmp+1) * sizeof(int));

    if (isSymmetric_tmp)
    {
        for (int i = 0; i < nnz_mtx_report; i++)
        {
            if (csrRowIdx_tmp[i] != csrColIdx_tmp[i])
            {
                int offset = csrRowPtr[csrRowIdx_tmp[i]] + csrRowPtr_counter[csrRowIdx_tmp[i]];
                csrColIdx[offset] = csrColIdx_tmp[i];
                csrVal[offset] = csrVal_tmp[i];
                csrRowPtr_counter[csrRowIdx_tmp[i]]++;

                offset = csrRowPtr[csrColIdx_tmp[i]] + csrRowPtr_counter[csrColIdx_tmp[i]];
                csrColIdx[offset] = csrRowIdx_tmp[i];
                csrVal[offset] = csrVal_tmp[i];
                csrRowPtr_counter[csrColIdx_tmp[i]]++;
            }
            else
            {
                int offset = csrRowPtr[csrRowIdx_tmp[i]] + csrRowPtr_counter[csrRowIdx_tmp[i]];
                csrColIdx[offset] = csrColIdx_tmp[i];
                csrVal[offset] = csrVal_tmp[i];
                csrRowPtr_counter[csrRowIdx_tmp[i]]++;
            }
        }
    }
    else
    {
        for (int i = 0; i < nnz_mtx_report; i++)
        {
            int offset = csrRowPtr[csrRowIdx_tmp[i]] + csrRowPtr_counter[csrRowIdx_tmp[i]];
            csrColIdx[offset] = csrColIdx_tmp[i];
            csrVal[offset] = csrVal_tmp[i];
            csrRowPtr_counter[csrRowIdx_tmp[i]]++;
        }
    }


    // free tmp space
    free(csrColIdx_tmp);
    free(csrVal_tmp);
    free(csrRowIdx_tmp);
    free(csrRowPtr_counter);

    return 0;
}

#endif
