#include "types.h"
#include "structs.h"
#include "matrix.h"
#include "sblas.h"

#include <map>
#include <vector>
#include <set>

int main(int argc, char** argv) {
  string matrixFile = "";
  sptElementIndex sb_bits = 6;    // block size
  uint64_t blockSize = 2 << sb_bits;
  if( argc == 3 ) {
    matrixFile = argv[1];
    blockSize = (2 << std::stoull(argv[2]));
  } else   if( argc == 2 ) {
    matrixFile = argv[1];
  } else {
    std::cout << "Usage: ./cppfile matrixFile blocksize\n";
    return 1;
  }

  // CsrSparseMatrix<unsigned,double> csrMtx2("./ash85.mtx");
  auto usehicoo = 1;

  sptSparseMatrix mtx;
  FILE *fi;
  if ((fi = fopen(matrixFile.c_str(), "r")) == NULL)
    return;
  sptAssert(sptLoadSparseMatrix(&mtx, 1, fi) == 0);
  fclose(fi);
  sptRandomValueVector(&(mtx.values));    // to better compare results
  sptSparseMatrixStatus(&mtx, stdout);
  // printf("Input COO Matrix\n");
  // sptAssert(sptDumpSparseMatrix(&mtx, 0, stdout) == 0);

  sptNnzIndex max_nnzb = 0;
  sptSparseMatrixHiCOO himtx;
  sptNnzIndex * bptr; // non-zero block pointers
  // sptStartTimer(timer);
  sptNewSparseMatrixHiCOO(&himtx, mtx.nrows, mtx.ncols, mtx.nnz, sb_bits, sb_bits);


  // Natural block sorting
  sptSparseMatrixSortIndexRowBlock(&mtx, 1, 0, mtx.nnz, sb_bits);  // OMP-Parallelized inside
  // Morton-order block sorting
  // sptSparseMatrixSortIndexMorton(&mtx, 1, 0, mtx.nnz, sb_bits);  // OMP-Parallelized inside
  sptAppendNnzIndexVector(&himtx.kptr, 0);
  sptAppendNnzIndexVector(&himtx.kptr, mtx.nnz);
  sptSparseMatrixPartition(&himtx, &max_nnzb, &mtx, sb_bits);   // Create a HiCOO matrix underneath
  bptr = himtx.bptr.data; // Extract block pointers from HiCOO matrix struct
  printf("Block pointers:\n");
  sptAssert(sptDumpNnzIndexArray(bptr, himtx.bptr.len, stdout) == 0);


  std::map<std::pair<uint64_t, uint64_t>, uint64_t > totalNnzsPerBlockMap;
  std::map<std::pair<uint64_t, uint64_t>, double > perBlockDensityMap;
  std::map<std::pair<uint64_t, uint64_t>, std::set<int> > nnzRowsPerBlockMap;
  std::map<std::pair<uint64_t, uint64_t>, std::set<int> > nnzColsPerBlockMap;
  std::map<std::pair<uint64_t, uint64_t>, std::pair<int, int> >
    nnzRowColCountPerBlockMap;

  size_t maxRowIdx = 0;
  size_t maxColIdx = 0;
  sptNnzIndex nb = himtx.bptr.len - 1;
  // sptElementIndex sb_bits = himtx.sb_bits;
  for(sptNnzIndex b = 0; b < nb; ++b) {   // Loop blocks
    sptNnzIndex bptr_begin = himtx.bptr.data[b];
    sptNnzIndex bptr_end = himtx.bptr.data[b+1];

    // store total num of nnzs in the current block
    auto blkIdxI = himtx.bindI.data[b];
    auto blkIdxJ = himtx.bindJ.data[b];

    if (blkIdxI > maxRowIdx) maxRowIdx = blkIdxI; 
    if (blkIdxJ > maxColIdx) maxColIdx = blkIdxJ;
    
    auto totalNnzCurBlk = bptr_end - bptr_begin;
    totalNnzsPerBlockMap.emplace(std::make_pair(blkIdxI, blkIdxJ), totalNnzCurBlk);
    
    // sptValue * restrict blocked_yvals = y->data + (himtx->bindI.data[b] << sb_bits);
    // sptValue * restrict blocked_xvals = x->data + (himtx->bindJ.data[b] << sb_bits);

    auto blkIdxICOO = himtx.bindI.data[b] << sb_bits;
    auto blkIdxJCOO = himtx.bindJ.data[b] << sb_bits;
    
    for(sptNnzIndex z = bptr_begin; z < bptr_end; ++z) {   // Loop entries
      sptElementIndex ei = himtx.eindI.data[z];
      sptElementIndex ej = himtx.eindJ.data[z];
      auto ICOO = blkIdxICOO + ei;
      auto JCOO = blkIdxJCOO + ej;
      //   // blocked_yvals[ei] += himtx->values.data[z] * blocked_xvals[ej];
      std::cout << "i: " << ICOO <<  " j: " << JCOO  << std::endl;

      //record nnzRowIdx and nnzColIdx
      auto foundBlk = nnzRowsPerBlockMap.find(std::make_pair(blkIdxI,
							     blkIdxJ));
      if (foundBlk != nnzRowsPerBlockMap.end()) {
	auto rowNum = ICOO % blockSize;
	foundBlk->second.insert(rowNum);
	/* for (auto &p: rowIdxSet) {std::cout << "set r: " << p << " ";} */
	/* std::cout << std::endl; */
      } else {
	std::set<int> tempSet;
	auto rowNum = ICOO % blockSize;
	tempSet.insert(rowNum);
	nnzRowsPerBlockMap.emplace(std::make_pair(blkIdxI, blkIdxJ), tempSet);
      }
	
      foundBlk = nnzColsPerBlockMap.find(std::make_pair(blkIdxI,
							blkIdxJ));
      if (foundBlk != nnzColsPerBlockMap.end()) {
	auto colNum = JCOO % blockSize;
	foundBlk->second.insert(colNum);
	/* for (auto &p: colIdxSet) {std::cout << "set c: " << p << " ";} */
	/* std::cout << std::endl; */
      } else {
	std::set<int> tempSet;
	auto colNum = JCOO % blockSize;
	tempSet.insert(colNum);
	nnzColsPerBlockMap.emplace(std::make_pair(blkIdxI, blkIdxJ), tempSet);
      }
    }
  }

  std::cout << "Printing total number of non-zeros per block" << std::endl;
  for (const auto &p : totalNnzsPerBlockMap) {
    std::cout << "(" << p.first.first << "," << p.first.second << ")" << "-->"
	      << p.second << std::endl;
  }

  // compute per-block density
  for (const auto &p : totalNnzsPerBlockMap) {
    perBlockDensityMap.emplace(std::make_pair(p.first.first, p.first.second),
			       (double) p.second / (double) (blockSize * blockSize));
  }

  std::cout << "Printing per block density" << std::endl;
  for (const auto &p : perBlockDensityMap) {
    std::cout << "(" << p.first.first << "," << p.first.second << ")" << "-->" << p.second << std::endl;
  }

  for (auto &p : totalNnzsPerBlockMap) {
    auto r  = p.first.first;
    auto c = p.first.second;
    // Get nnzRowCount for this block
    auto nnzRowCount = nnzRowsPerBlockMap.find(std::make_pair(r, c))->second.size();
    /* for (auto &q : nnzRowsPerBlockMap.find(std::make_pair(r, c))->second) { */
    /* 	//for (auto &r : q) */
    /* 	  std::cout << q << " "; */
    /* 	std::cout << std::endl; */
    /* } */
    // Get nnzColCount for this block
    auto nnzColCount = nnzColsPerBlockMap.find(std::make_pair(r, c))->second.size();
    // std::cout << nnzRowCount << " " << nnzColCount << std::endl;
    nnzRowColCountPerBlockMap.emplace(std::make_pair(r,c ),
				      std::make_pair(nnzRowCount, nnzColCount));
  }

  std::cout << "Printing total number of non-zero rows and cols per block" << std::endl;
  for (const auto &p : nnzRowColCountPerBlockMap) {
    std::cout << "(" << p.first.first << "," << p.first.second << ")" << "-->(" << p.second.first << "," << p.second.second << ")" << std::endl;
  }


  // Determine whether the blk should be represented as dense or sparse.
  int *block_sd = (int *)malloc((nb + 1) * sizeof(int));
  memset(block_sd, 0, (nb + 1) * sizeof(int));
  int cnt = 0;
  int nnzblk = 0;
  for (const auto &p : perBlockDensityMap) {
    if (p.second > 0) {
      block_sd[cnt] = 1;
      nnzblk++;
    } else {
      block_sd[cnt] = 0;
    }
    cnt++;
  }

  std::cout << "cnt " << cnt << std::endl;
  
  // Host part of the dense blocks
  int *denseMatrix = (int*) malloc (nnzblk * blockSize * blockSize * sizeof(int));
  auto cnt_ = 0;
  //densify some blocks
  for(sptNnzIndex b = 0; b < nb; ++b) {   // Loop blocks
    if (block_sd[b] == 1) {
      sptNnzIndex bptr_begin = himtx.bptr.data[b];
      sptNnzIndex bptr_end = himtx.bptr.data[b+1];
      auto current_start = cnt_ * blockSize * blockSize;
      for(sptNnzIndex z = bptr_begin; z < bptr_end; ++z) {   // Loop entries
	sptElementIndex ei = himtx.eindI.data[z];
	sptElementIndex ej = himtx.eindJ.data[z];
	*(denseMatrix + current_start + (ei * blockSize + ej)) = 1;
      }
      cnt_++;
    }
  }

  // print the dense matrix for verification.
  std::cout << "printing dense matrix" << std::endl;
  for (auto i = 0; i < nnzblk; i++) {
    auto current_start = denseMatrix + (i * blockSize * blockSize );
    for (auto j = 0; j < blockSize; j++) {
       for (auto k = 0; k < blockSize; k++) {
	 std::cout <<  *(current_start + (j * blockSize + k)) << " ";
       }
       std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
  }

  
  // copy this information back to GPU
  sptNnzIndexVector         *block_sd_gpu;
  cudaMalloc((void **) &block_sd_gpu, (nb + 1) * sizeof(int)); // TODO: convert to bool
  cudaMemcpy(block_sd_gpu, block_sd, (nb + 1) * sizeof(int),
	     cudaMemcpyHostToDevice);

  // copy the dense blocks to gpu
  int *denseMatrixGpu;
  cudaMalloc((void **) &denseMatrixGpu,  nnzblk * blockSize * blockSize * sizeof(int));
  cudaMemcpy(denseMatrixGpu, denseMatrix, nnzblk * blockSize * blockSize * sizeof(int),
	     cudaMemcpyHostToDevice);
  
  /*GPU-related data structures*/
  sptNnzIndexVector         *bptr_gpu;      /// Block pointers to all nonzeros on the GPU, nb = bptr.length - 1
  sptBlockIndexVector       *bindI_gpu;    /// Block indices for rows on the GPU, length nb
  sptBlockIndexVector       *bindJ_gpu;    /// Block indices for columns on the GPU, length nb
  sptElementIndexVector     *eindI_gpu;    /// Element indices within each block for rows on the GPU, length nnz
  sptElementIndexVector     *eindJ_gpu;    /// Element indices within each block for columns on the GPU, length nnz
  sptValueVector            values_gpu;      /// non-zero values on the GPU, length nnz

  sptNnzIndex *bptr_len_gpu;
  sptNnzIndex *bindI_len_gpu;
  sptNnzIndex *bindJ_len_gpu;
  sptNnzIndex *eindI_len_gpu;
  sptNnzIndex *eindJ_len_gpu;
  // sptNnzIndex val_len_gpu;
  
  cudaMalloc((void **) &bptr_gpu, himtx.bptr.len * sizeof(sptNnzIndex));
  cudaMemcpy(bptr_gpu, himtx.bptr.data, himtx.bptr.len * sizeof(sptNnzIndex),
	     cudaMemcpyHostToDevice);
  cudaMalloc((void **) &bptr_len_gpu, sizeof(sptNnzIndex));
  cudaMemcpy(bptr_len_gpu, &himtx.bptr.len, sizeof(sptNnzIndex),
	     cudaMemcpyHostToDevice);

  cudaMalloc((void **) &bindI_gpu, himtx.bindI.len * sizeof(sptBlockIndex));
  cudaMemcpy(bindI_gpu, himtx.bindI.data, himtx.bindI.len * sizeof(sptBlockIndex),
	     cudaMemcpyHostToDevice);
  cudaMalloc((void **) &bindI_len_gpu, sizeof(sptNnzIndex));
  cudaMemcpy(bindI_len_gpu, &himtx.bindI.len, sizeof(sptNnzIndex),
	     cudaMemcpyHostToDevice);

  cudaMalloc((void **) &bindJ_gpu, himtx.bindJ.len * sizeof(sptBlockIndex));
  cudaMemcpy(bindI_gpu, himtx.bindI.data, himtx.bindI.len * sizeof(sptBlockIndex),
	     cudaMemcpyHostToDevice);
  cudaMalloc((void **) &bindJ_len_gpu, sizeof(sptNnzIndex));
  cudaMemcpy(bindJ_len_gpu, &himtx.bindJ.len, sizeof(sptNnzIndex),
	     cudaMemcpyHostToDevice);

  cudaMalloc((void **) &eindI_gpu, himtx.eindI.len * sizeof(sptElementIndex));
  cudaMemcpy(eindI_gpu, himtx.eindI.data, himtx.eindI.len * sizeof(sptElementIndex),
	     cudaMemcpyHostToDevice);
  cudaMalloc((void **) &eindI_len_gpu, sizeof(sptNnzIndex));
  cudaMemcpy(eindI_len_gpu, &himtx.eindI.len, sizeof(sptNnzIndex),
	     cudaMemcpyHostToDevice);

  cudaMalloc((void **) &eindJ_gpu, himtx.eindJ.len * sizeof(sptElementIndex));
  cudaMemcpy(eindJ_gpu, himtx.eindJ.data, himtx.eindJ.len * sizeof(sptElementIndex),
	     cudaMemcpyHostToDevice);
  cudaMalloc((void **) &eindJ_len_gpu, sizeof(sptNnzIndex));
  cudaMemcpy(eindJ_len_gpu, &himtx.eindJ.len, sizeof(sptNnzIndex),
	     cudaMemcpyHostToDevice);


  //current frontier blocks
  int* currentFrontiersBlocks;
  cudaMalloc((void **) &currentFrontiersBlocks, sizeof(int) * maxRowIdx * blockSize);
  cudaMemset(currentFrontiersBlocks, 0, sizeof(int) * maxRowIdx * blockSize);

  //next frontier blocks
  int* nextFrontiersBlocks;
  cudaMalloc((void **) &nextFrontiersBlocks, sizeof(int) * nb * blockSize);
  cudaMemset(currentFrontiersBlocks, 0, sizeof(int) * nb * blockSize);

  
  // set source for multiple frontiers
  auto source = 1; // TODO: change to take it from the cmd prompt.
  auto sourceFrontierBlk = (int) (source / blockSize);
  auto sourceFrontierIdx = source % blockSize;
  int tmp = 1;
  cudaMemcpy(&currentFrontiersBlocks[(sourceFrontierBlk * blockSize) + sourceFrontierIdx],
	     &tmp, sizeof(int), cudaMemcpyHostToDevice);
  
  // launch kernel

  //
  
  sptFreeSparseMatrixHiCOO(&himtx);
  sptFreeSparseMatrix(&mtx);


  // CsrSparseMatrix<unsigned,double> csrMtx(matrixFile.c_str(), blockSize);
  // CsrSparseMatrix<unsigned,double> csrMtx(matrixFile.c_str(), blockSize, usehicoo);
  // CsrSparseMatrix<unsigned,double> csrMtx("./data/citHepPh/citHepPh.mtx", blockSize);
  return 0;
}
