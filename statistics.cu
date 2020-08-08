#include "types.h"
#include "structs.h"
#include "matrix.h"
#include "sblas.h"



int main(int argc, char** argv) {
  string matrixFile = "";
  uint64_t blockSize = 16;
  if( argc == 3 ) {
    matrixFile = argv[1];
    blockSize = std::stoull(argv[2]);
  } else   if( argc == 2 ) {
    matrixFile = argv[1];
  } else {
    std::cout << "Usage: ./cppfile matrixFile \n";
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
  printf("Input COO Matrix\n");
  sptAssert(sptDumpSparseMatrix(&mtx, 0, stdout) == 0);

  sptNnzIndex max_nnzb = 0;
  sptSparseMatrixHiCOO himtx;
  sptNnzIndex * bptr; // non-zero block pointers
  sptElementIndex sb_bits = 6;    // block size
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
  
  sptNnzIndex nb = himtx.bptr.len - 1;
  // sptElementIndex sb_bits = himtx.sb_bits;
  for(sptNnzIndex b = 0; b < nb; ++b) {   // Loop blocks
    sptNnzIndex bptr_begin = himtx.bptr.data[b];
    sptNnzIndex bptr_end = himtx.bptr.data[b+1];
    // sptValue * restrict blocked_yvals = y->data + (himtx->bindI.data[b] << sb_bits);
    // sptValue * restrict blocked_xvals = x->data + (himtx->bindJ.data[b] << sb_bits);

    for(sptNnzIndex z = bptr_begin; z < bptr_end; ++z) {   // Loop entries
      sptElementIndex ei = himtx.eindI.data[z];
      sptElementIndex ej = himtx.eindJ.data[z];
      // blocked_yvals[ei] += himtx->values.data[z] * blocked_xvals[ej];
      std::cout << "i: " << b*64+ei <<  "j: " << b*64+ej  << std::endl;
    }
  }

  
  sptFreeSparseMatrixHiCOO(&himtx);
  sptFreeSparseMatrix(&mtx);


  // CsrSparseMatrix<unsigned,double> csrMtx(matrixFile.c_str(), blockSize);
  // CsrSparseMatrix<unsigned,double> csrMtx(matrixFile.c_str(), blockSize, usehicoo);
  // CsrSparseMatrix<unsigned,double> csrMtx("./data/citHepPh/citHepPh.mtx", blockSize);
  return 0;
}
