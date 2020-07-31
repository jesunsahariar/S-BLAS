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
  // CsrSparseMatrix<unsigned,double> csrMtx(matrixFile.c_str(), blockSize);
  CsrSparseMatrix<unsigned,double> csrMtx(matrixFile.c_str(), blockSize, usehicoo);
  // CsrSparseMatrix<unsigned,double> csrMtx("./data/citHepPh/citHepPh.mtx", blockSize);
  return 0;
}
