#include "BlockMatrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sys/time.h>

// Tests the construction and debug printing of the BlockMatrix class for
// dense matrices
int test_1()
{
   std::cout << "=================================" << std::endl;
   std::cout << " test_1: Ctor and debug printing" << std::endl;
   std::cout << "=================================" << std::endl;

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
   std::cout << "=======================" << std::endl;
   std::cout << " test_2: Block formats " << std::endl;
   std::cout << "=======================" << std::endl;

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

   std::vector<int> stripes{-4,-2,-1,1};

   bm(3,0).init(9,3,VMatrix::Diagonal,   &debugFunctor).bind();
   bm(3,1).init(9,7,VMatrix::Striped,    &debugFunctor).bindStripes(stripes);
   bm(3,2).init(9,5,VMatrix::Diagonal,   &diagonalFunctor).bind();

   bm.info();
   bm.print();

   return 0;
}


int test_3()
{
   std::cout << "================================" << std::endl;
   std::cout << " test_3: dense <- dense x dense " << std::endl;
   std::cout << "================================" << std::endl;

   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init( 9,10, VMatrix::Dense, &testFunctor).bind();
   b.init(10, 8, VMatrix::Dense, &testFunctor).bind();
   c.init( 9, 8, VMatrix::Dense, &zeroFunctor).bind();
   VMatrix::matrix_product(c, a, b); 

   std::string fname("test_3.dat");
   std::ifstream ifs(fname.c_str(), std::ios::in);

   if (ifs.is_open()) {
      std::string line;
      double x, res(0);
      double* data(c.data());
      unsigned k(0);

      while (ifs >> x) {
         res += std::abs(x-data[k]);
         ++k;
      }
      std::cout << "Matrix residue: " << res << std::endl;

      ifs.close();
      if (res > 1e-12) return 1;
   }else {
      std::cerr << "Failed to open flie " << fname << std::endl;
      return 1;
   }

   std::cout << "PASS" << std::endl << std::endl;
  return 0;
}


// Test the block matrix multiplication C = A.A
int test_4()
{
   std::cout << "=====================================" << std::endl;
   std::cout << " test_4: Block matrix multiplication" << std::endl;
   std::cout << "=====================================" << std::endl;

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

   C.toDense(&b);

   double* data_c(c.data());
   double* data_b(b.data());
   double  res(0);
   unsigned n(c.nCols()*c.nRows());

   for (unsigned i = 0; i < n; ++i) {
       res += std::abs(data_c[i] - data_b[i]);
   }


   std::cout << "Matrix residue:   " << res << std::endl;
   if (res > n*1e-12) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


//Timing test
int test_5()
{
   std::cout << "=====================" << std::endl;
   std::cout << " test_5: Timing test " << std::endl;
   std::cout << "=====================" << std::endl;

   const unsigned dim(1024);
   const unsigned blocks(2);

   BlockMatrix A(blocks,blocks);
   BlockMatrix C(blocks,blocks);

   ZeroFunctor zeroFunctor;
   DebugFunctor fun1;

   for (unsigned bi = 0; bi < blocks; ++bi) {
       for (unsigned bj = 0; bj < blocks; ++bj) {
           //if (((bi + bj) % 2) == 0) {
           if (true ) {
           //if (bi ==  bj) {
              A(bi,bj).init(dim,dim, VMatrix::Dense   , &fun1).bind();
              C(bi,bj).init(dim,dim, VMatrix::Dense   , &zeroFunctor).bind();
           }else {
              A(bi,bj).init(dim,dim, VMatrix::Zero    , &fun1).bind();
              C(bi,bj).init(dim,dim, VMatrix::Dense   , &zeroFunctor).bind();
           }
       }
   }

//   A.print("A Matrix");
   
   VMatrix a, c;
   C.toDense(&c);
   A.toDense(&a);

   struct timeval tv1, tv2;
   struct timezone tz;
   double elapsed;
   
   gettimeofday(&tv1, &tz);
   BlockMatrix::matrix_product(C,A,A);
   gettimeofday(&tv2, &tz);
   elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
   std::cout << "BlockMatrix time: " << elapsed << std::endl;;

 //  C.print("Product from BlockMatrix multiplication");

   gettimeofday(&tv1, &tz);
   VMatrix::matrix_product(c,a,a);
   gettimeofday(&tv2, &tz);
   elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
   std::cout << "VMatrix time:     " << elapsed << std::endl;;


   std::cout << "-----------------" << std::endl;

   VMatrix b;
   C.toDense(&b);


   double* data_c(c.data());
   double* data_b(b.data());
   double  res(0);
   unsigned n(c.nCols()*c.nRows());

   for (unsigned i = 0; i < n; ++i) {
       res += std::abs( data_c[i] - data_b[i]);
   }

   std::cout << "Matrix residue:   " << res << std::endl;
   if (res > n*1e-12) {
//      b.print("b VMatrix obtained from densifying C");
//      c.print("c VMatrix obtained from VMatrix multiplication");
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}



int test_6()
{
   std::cout << "===============================" << std::endl;
   std::cout << " test_6: dense <- diag x dense" << std::endl;
   std::cout << "===============================" << std::endl;

   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Diagonal, &testFunctor).bind();
   b.init(10,10, VMatrix::Dense,    &testFunctor).bind();
   c.init(10,10, VMatrix::Dense,    &zeroFunctor).bind();

   VMatrix::matrix_product(c, a, b); 

   std::string fname("test_6.dat");
   std::ifstream ifs(fname.c_str(), std::ios::in);

   if (ifs.is_open()) {
      std::string line;
      double x, res(0);
      double* data(c.data());
      unsigned k(0);

      while (ifs >> x) {
         res += std::abs(x-data[k]);
         ++k;
      }
      std::cout << "Matrix residue: " << res << std::endl;

      ifs.close();
      if (res > 1e-12) return 1;

   }else {
      std::cerr << "Failed to open flie " << fname << std::endl;
      return 1;
   }

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


int test_7()
{
   std::cout << "===============================" << std::endl;
   std::cout << " test_7: dense <- dense x diag " << std::endl;
   std::cout << "===============================" << std::endl;

   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Dense,    &testFunctor).bind();
   b.init(10,10, VMatrix::Diagonal, &testFunctor).bind();
   c.init(10,10, VMatrix::Dense,    &zeroFunctor).bind();

   VMatrix::matrix_product(c, a, b); 

   std::string fname("test_7.dat");
   std::ifstream ifs(fname.c_str(), std::ios::in);

   if (ifs.is_open()) {
      std::string line;
      double x, res(0);
      double* data(c.data());
      unsigned k(0);

      while (ifs >> x) {
         res += std::abs(x-data[k]);
         ++k;
      }
      std::cout << "Matrix residue: " << res << std::endl;

      ifs.close();
      if (res > 1e-12) {
         std::cerr << "FAIL: residue too large" << std::endl<< std::endl;
         return 1;
      }

   }else {
      std::cerr << "FAIL: unable to open flie " << fname << std::endl;
      return 1;
   }

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


int test_8()
{
   std::cout << "==============================" << std::endl;
   std::cout << " test_8: dense <- diag x diag " << std::endl;
   std::cout << "==============================" << std::endl;

   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Diagonal, &testFunctor).bind();
   b.init(10,10, VMatrix::Diagonal, &testFunctor).bind();
   c.init(10,10, VMatrix::Dense,    &zeroFunctor).bind();

   VMatrix::matrix_product(c, a, b); 

   std::string fname("test_8.dat");
   std::ifstream ifs(fname.c_str(), std::ios::in);

   c.print();

   if (ifs.is_open()) {
      std::string line;
      double x, res(0);
      double* data(c.data());
      unsigned k(0);

      while (ifs >> x) {
         res += std::abs(x-data[k]);
         ++k;
      }
      std::cout << "Matrix residue: " << res << std::endl;

      ifs.close();
      if (res > 1e-12) {
         std::cerr << "FAIL: residue too large" << std::endl<< std::endl;
         return 1;
      }

   }else {
      std::cerr << "FAIL: unable to open flie " << fname << std::endl;
      return 1;
   }

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


int test_9()
{
   std::cout << "==================================" << std::endl;
   std::cout << " test_9: dense <- striped x dense " << std::endl;
   std::cout << "==================================" << std::endl;

   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;
   DebugFunctor debugFunctor;

   VMatrix a, b, c, d;

   //std::vector<int> stripes = {-3,-1,1,3};
   std::vector<int> stripes = {1,3};

   a.init(5,5, VMatrix::Striped, &debugFunctor).bindStripes(stripes);
   b.init(5,5, VMatrix::Dense,   &testFunctor).bind();
   c.init(5,5, VMatrix::Dense,   &zeroFunctor).bind();
   d.init(5,5, VMatrix::Dense,   &zeroFunctor).bind();

   a.print();
   b.print();


   VMatrix::matrix_product(c, a, b); 

   std::string fname("test_9.dat");
   std::ifstream ifs(fname.c_str(), std::ios::in);

   c.print();

   // Do a dense version for checking
   a.toDense();
   VMatrix::matrix_product(d, a, b); 

   double* data_c(c.data());
   double* data_d(d.data());
   double  res(0);
   unsigned n(c.nCols()*c.nRows());

   for (unsigned i = 0; i < n; ++i) {
       res += std::abs(data_c[i] - data_d[i]);
   }


   std::cout << "Matrix residue:   " << res << std::endl;
   if (res > n*1e-12) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;

}

int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok = 
        test_1()
/*
      + test_2()
      + test_3()
      + test_4()
      + test_5()
      + test_6()
      + test_7()
      + test_8()
*/
      + test_9()
      ;

   std::cout << std::endl;
   std::cout << "Test run finished: ";

   if (ok) {
      std::cout << "FAIL" << std::endl;
   }else {
      std::cout << "PASS" << std::endl;
   }

   std::cout << "With " << ok << " failure(s)" << std::endl;

   return ok;
}
