#include "BlockMatrix.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

// Tests the construction and debug printing of the BlockMatrix class for
// dense matrices
int test_1()
{
   unsigned nRowBlocks(4);
   unsigned nColBlocks(3);

   BlockMatrix bm(nRowBlocks,nColBlocks);

   // Used to fill the matrix entries with index information
   DebugFunctor debugFunctor;
  
   // Initialize the structure of our BlockMatrix.  Note this does not
   // allocate any memory.
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           // These are the dimensions of the current tile
           unsigned nRows(2*(row+1));
           unsigned nCols(2*(col+1)+1);
           bm(row,col).init(nRows,nCols,VMatrix::Dense,&debugFunctor);
       }
   }

   // Bind our data, this allocates the memory and sets the entries
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           bm(row,col).bind();
       }
   }

   bm.info();
   bm.print();

   return 0;
}


// Tests the different storage formats for the tiles.
int test_2()
{
   unsigned nRowBlocks(4);
   unsigned nColBlocks(3);

   BlockMatrix bm(nRowBlocks,nColBlocks);

   // Used to fill the matrix entries with index information
   DebugFunctor    debugFunctor;
   ZeroFunctor     zeroFunctor;
   DiagonalFunctor diagonalFunctor;
  
   // Initialize the structure of our BlockMatrix.  Note in this case
   // we bind the data at the same time.
   bm(0,0).init(2,3,VMatrix::Diagonal,   &diagonalFunctor).bind();
   bm(0,1).init(2,7,VMatrix::Zero,       &zeroFunctor).bind();
   bm(0,2).init(2,5,VMatrix::Dense,      &debugFunctor).bind();

   bm(1,0).init(6,3,VMatrix::Diagonal,   &diagonalFunctor).bind();
   bm(1,1).init(6,7,VMatrix::Tridiagonal,&debugFunctor).bind();
   bm(1,2).init(6,5,VMatrix::Diagonal,   &diagonalFunctor).bind();

   bm(2,0).init(2,3,VMatrix::Diagonal,   &diagonalFunctor).bind();
   bm(2,1).init(2,7,VMatrix::Zero,       &zeroFunctor).bind();
   bm(2,2).init(2,5,VMatrix::Diagonal,   &diagonalFunctor).bind();

   bm(3,0).init(5,3,VMatrix::Diagonal,   &debugFunctor).bind();
   bm(3,1).init(5,7,VMatrix::Dense,       &debugFunctor).bind();
   bm(3,2).init(5,5,VMatrix::Diagonal,   &diagonalFunctor).bind();

   bm.info();
   bm.print();

   return 0;
}


// Tests matrix multiplication for VMatrix class
int test_3()
{
   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init( 9,10, VMatrix::Dense, &testFunctor).bind();
   b.init(10, 8, VMatrix::Dense, &testFunctor).bind();
   c.init( 9, 8, VMatrix::Dense, &zeroFunctor).bind();
   VMatrix::matrix_product(c, a, b); 

   c.print();

   std::string fname("test_3.dat");

   std::ifstream ifs(fname.c_str(), std::ios::in);

   if (ifs.is_open()) {
      std::string line;
      double x, res(0);
      double* data(c.data());
      unsigned k(0);

      while (ifs >> x) {
         res += std::abs(x-*(data+k));
         ++k;
      }
      std::cout << "Matrix residue: " << res << std::endl;

      ifs.close();
      if (res > 1e-12) return 1;
   }else {
      std::cerr << "Failed to open flie " << fname << std::endl;
      return 1;
   }

  return 0;
}


// Test the block matrix multiplication C = A.A
int test_4()
{
   std::cout << "Entering test_4()" << std::endl;
   BlockMatrix A(3,3);
   BlockMatrix C(3,3);
   ZeroFunctor zeroFunctor;
   
   TestFunctor fun1(0,0);
   A(0,0).init(4,4, VMatrix::Dense, &fun1).bind();
   C(0,0).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun2(4,0);
   A(1,0).init(4,4, VMatrix::Dense, &fun2).bind();
   C(1,0).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun3(8,0);
   A(2,0).init(4,4, VMatrix::Dense, &fun3).bind();
   C(2,0).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun4(0,4);
   A(0,1).init(4,4, VMatrix::Dense, &fun4).bind();
   C(0,1).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun5(4,4);
   A(1,1).init(4,4, VMatrix::Dense, &fun5).bind();
   C(1,1).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun6(8,4);
   A(2,1).init(4,4, VMatrix::Dense, &fun6).bind();
   C(2,1).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun7(0,8);
   A(0,2).init(4,4, VMatrix::Dense, &fun7).bind();
   C(0,2).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun8(4,8);
   A(1,2).init(4,4, VMatrix::Dense, &fun8).bind();
   C(1,2).init(4,4, VMatrix::Dense, &zeroFunctor).bind();

   TestFunctor fun9(8,8);
   A(2,2).init(4,4, VMatrix::Dense, &fun9).bind();
   C(2,2).init(4,4, VMatrix::Dense, &zeroFunctor).bind();


   BlockMatrix::matrix_product(C, A, A);

   VMatrix a, b, c;
   a.init(12,12, VMatrix::Dense, &fun1).bind();
   c.init(12,12, VMatrix::Dense, &zeroFunctor).bind();

   VMatrix::matrix_product(c,a,a);

   C.print();
   C.toDense(&b);

   c.print(); 
   b.print();
   
   double* data_c(c.data());
   double* data_b(c.data());
   double  res(0);
   for (unsigned i = 0; i < 144; ++i) {
       res += std::abs(*(data_c+i) - *(data_b+i));
   }

   std::cout << "Matrix residue: " << res << std::endl;
   return  (res > 1e-12) ? 1 : 0;
}


int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok = 
      test_1() + 
      test_2() + 
      test_3() +
      test_4();

   std::cout << "Test run finished: ";

   if (ok) {
      std::cout << "FAIL" << std::endl;
   }else {
      std::cout << "PASS" << std::endl;
   }

   return ok;
}
