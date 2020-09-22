#include "BlockMatrix.h"
#include "MatMult.h"
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
   std::cout << " test_1: Ctor and debug printing " << std::endl;
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
           bm(row,col).init(nRows,nCols,VMatrix::Dense);
       }
   }

   // Bind our data, this allocates the memory and sets the entries
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           bm(row,col).bind(debugFunctor);
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
   bm(0,0).init(2,3,VMatrix::Diagonal).bind(diagonalFunctor);
   bm(0,1).init(2,7,VMatrix::Zero    ).bind(zeroFunctor);
   bm(0,2).init(2,5,VMatrix::Dense   ).bind(debugFunctor);

   bm(1,0).init(6,3,VMatrix::Diagonal).bind(diagonalFunctor);
   bm(1,1).init(6,7,VMatrix::Diagonal).bind(debugFunctor);
   bm(1,2).init(6,5,VMatrix::Diagonal).bind(diagonalFunctor);

   bm(2,0).init(2,3,VMatrix::Diagonal).bind(diagonalFunctor);
   bm(2,1).init(2,7,VMatrix::Zero    ).bind(zeroFunctor);
   bm(2,2).init(2,5,VMatrix::Diagonal).bind(diagonalFunctor);

   std::vector<int> stripes{-4,-2,-1,1};
   std::vector<int> stripes2{-1,1};

   bm(3,0).init(9,3,VMatrix::Zero    ).bind(debugFunctor);
   bm(3,1).init(9,7,stripes          ).bind(debugFunctor);
   bm(3,2).init(9,5,stripes2         ).bind(debugFunctor);

   bm.info();
   bm.print();

   return 0;
}


int test_3()
{
   std::cout << "================================" << std::endl;
   std::cout << " test_3: dense <- dense x dense " << std::endl;
   std::cout << "================================" << std::endl;

   TestFunctor testFunctor;
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init( 9,10).bind(testFunctor);
   b.init(10, 8).bind(testFunctor);
   c.init( 9, 8).bind(zeroFunctor);
   matrix_product(c, a, b); 

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
   std::cout << " test_4: Block matrix multiplication " << std::endl;
   std::cout << "=====================================" << std::endl;

   BlockMatrix A(3,3);
   BlockMatrix C(3,3);
   ZeroFunctor zeroFunctor;

   for (unsigned i = 0; i < 3; ++i) {
       for (unsigned j = 0; j < 3; ++j) {
           A(i,j).init(4,4);
           C(i,j).init(4,4).bind(zeroFunctor);
       }
   }

   // Need to init before binding A to use the row/colOffset functions
   for (unsigned i = 0; i < 3; ++i) {
       for (unsigned j = 0; j < 3; ++j) {
           A(i,j).bind(TestFunctor(A.rowOffset(i),A.rowOffset(j)));
       }
   }

   // Compute the matrix product via blocks
   matrix_product(C, A, A);

   VMatrix a, b, c;
   a.init(12,12, VMatrix::Dense).bind(TestFunctor());
   c.init(12,12, VMatrix::Dense).bind(zeroFunctor);

   // Compute the matrix product all as one
   matrix_product(c,a,a);

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

   const unsigned dim(124);
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
              A(bi,bj).init(dim,dim, VMatrix::Dense).bind(fun1);
              C(bi,bj).init(dim,dim, VMatrix::Dense).bind(zeroFunctor);
           }else {
              A(bi,bj).init(dim,dim, VMatrix::Zero ).bind(fun1);
              C(bi,bj).init(dim,dim, VMatrix::Dense).bind(zeroFunctor);
           }
       }
   }

   A.info("A Matrix:");
// A.print("A Matrix");
   
   VMatrix a, c;
   C.toDense(&c);
   A.toDense(&a);

   struct timeval tv1, tv2;
   struct timezone tz;
   double elapsed;
   
   gettimeofday(&tv1, &tz);
   matrix_product(C,A,A);
   gettimeofday(&tv2, &tz);
   elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
   std::cout << "BlockMatrix time: " << elapsed << std::endl;;

 //  C.print("Product from BlockMatrix multiplication");

   gettimeofday(&tv1, &tz);
   matrix_product(c,a,a);
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
       res = std::max(res, std::abs( data_c[i] - data_b[i]));
   }

   std::cout << "Max deviation:    " << res << std::endl;
   if (res > 1e-8) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}



int test_6()
{
   std::cout << "===============================" << std::endl;
   std::cout << " test_6: dense <- diag x dense " << std::endl;
   std::cout << "===============================" << std::endl;

   TestFunctor testFunctor(0,0);
   ZeroFunctor zeroFunctor;

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   b.init(10,10, VMatrix::Dense   ).bind(testFunctor);
   c.init(10,10, VMatrix::Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

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

   a.init(10,10, VMatrix::Dense   ).bind(testFunctor);
   b.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   c.init(10,10, VMatrix::Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

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

   a.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   b.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   c.init(10,10, VMatrix::Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

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

   TestFunctor  testFunctor;
   ZeroFunctor  zeroFunctor;
   DebugFunctor debugFunctor;

   VMatrix a, b, c, d;

   std::vector<int> stripes = {-8,-2,2,4};

   a.init(15,16, stripes).bind(debugFunctor);
   b.init(16,10, VMatrix::Dense).bind(debugFunctor);
   c.init(15,10, VMatrix::Dense).bind(zeroFunctor);
   d.init(15,10, VMatrix::Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 

   //a.print("A with stripes");
   //b.print("B is dense");
   //c.print("Product of striped x dense");

   // Do a dense version for checking
   a.toDense();
   matrix_product(d, a, b); 

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

int test_10()
{
   std::cout << "===================================" << std::endl;
   std::cout << " test_10: dense <- dense x striped " << std::endl;
   std::cout << "===================================" << std::endl;

   ZeroFunctor  zeroFunctor;
   DebugFunctor debugFunctor;

   VMatrix a, b, c, d;
   std::vector<int> stripes = {-8,-2,2,4};

   a.init(15,16, VMatrix::Dense).bind(debugFunctor);
   b.init(16,10, stripes).bind(debugFunctor);
   c.init(15,10, VMatrix::Dense).bind(zeroFunctor);
   d.init(15,10, VMatrix::Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 
   //a.print();
   //b.print();
   //c.print();

   // Do a dense version for checking
   b.toDense();
   matrix_product(d, a, b); 

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


int test_11()
{
   std::cout << "=====================================" << std::endl;
   std::cout << " test_11 dense <- striped x striped " << std::endl;
   std::cout << "=====================================" << std::endl;

   ZeroFunctor  zeroFunctor;
   DebugFunctor debugFunctor;

   VMatrix a, b, c, d;
   std::vector<int> stripesA = {-2};
   std::vector<int> stripesB = {-2};

   a.init(6,7, stripesA).bind(debugFunctor);
   b.init(7,8, stripesB).bind(debugFunctor);
   c.init(6,8, VMatrix::Dense).bind(zeroFunctor);
   d.init(6,8, VMatrix::Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 
   a.print();
   b.print();
   c.print();

   // Do a dense version for checking
   a.toDense();
   b.toDense();
   matrix_product(d, a, b); 

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


int test_12()
{
   std::cout << "=========================" << std::endl;
   std::cout << " test_12 Banded matrices " << std::endl;
   std::cout << "=========================" << std::endl;

   DebugFunctor debugFunctor;
   VMatrix a;

   a.init(10,10,4,4).bind(debugFunctor);

   double t(0);
   for (int i = 0; i < 5; ++i) {
       for (int j = 0; j < 5; ++j) {
           t+=a(i,j);
       }
   }

   a.print("banded matrix print:");

   a.toDense();

   a.print("banded to dense");

   return 0;
}



int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok = 
/*
        test_1()
      + test_2()
      + test_3()
      + test_4()
      + test_5()
      + test_6()
      + test_7()
      + test_8()
      + test_9()
      + test_10()
      + test_11()
*/
      + test_12()
      ;

   std::cout << std::endl;
   std::cout << "Test run finished: ";

   if (ok) {
      std::cout << "FAIL" << std::endl;
      std::cout << "With " << ok << " failure(s)" << std::endl;
   }else {
      std::cout << "ALL PASS" << std::endl;
   }


   return ok;
}
