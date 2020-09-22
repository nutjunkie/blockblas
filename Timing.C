#include "BlockMatrix.h"
#include "MatMult.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <stdio.h>
#include <ctime>


struct timeval  time_start;
struct timezone time_zone;

void timerStart() 
{
    gettimeofday(&time_start, &time_zone);
}

double timerStop() 
{
   struct timeval time_stop;
   gettimeofday(&time_stop, &time_zone);
   return (double) (time_stop.tv_sec -time_start.tv_sec) + 
          (double) (time_stop.tv_usec-time_start.tv_usec) * 1.e-6;
}


void makeDense(BlockMatrix& bm, unsigned dim, Functor const& functor)
{ 
   unsigned nRows(bm.nRowBlocks());
   unsigned nCols(bm.nColBlocks());

   ZeroFunctor  zeroFunctor;
   DebugFunctor debugFunctor;

   for (unsigned bi = 0; bi < nRows; ++bi) {
       for (unsigned bj = 0; bj < nCols; ++bj) {
           bm(bi,bj).init(dim,dim, VMatrix::Dense).bind(functor);
       }
   }
}


void makeDiagonal(BlockMatrix& bm, unsigned dim, Functor const& functor)
{ 
   unsigned nRows(bm.nRowBlocks());
   unsigned nCols(bm.nColBlocks());

   ZeroFunctor  zeroFunctor;
   DebugFunctor debugFunctor;

   for (unsigned bi = 0; bi < nRows; ++bi) {
       for (unsigned bj = 0; bj < nCols; ++bj) {
           if (bi ==  bj) {
              bm(bi,bj).init(dim,dim, VMatrix::Dense).bind(functor);
           }else {
           }
       }
   }
}


int test_1(unsigned blocks, unsigned dim)
{
   std::cout << "==============================" << std::endl;
   std::cout << " test_1: Dense block multiply " << std::endl;
   std::cout << "==============================" << std::endl;

   BlockMatrix A(blocks,blocks);
   BlockMatrix C(blocks,blocks);

   makeDense(A, dim, DebugFunctor());
   makeDense(C, dim, ZeroFunctor());

   VMatrix a, c;
   C.toDense(&c);
   A.toDense(&a);

   //A.info("A Matrix:");
   std::cout << "Blocks: " << blocks << " x " << blocks << std::endl;
   std::cout << "Tiles:  " << dim << " x " << dim << std::endl;
   std::cout << "Size:   " << A.nRows() << " x " << A.nCols() << std::endl;
   
   
   timerStart();
   matrix_product(C,A,A);
   std::cout << "BlockMatrix time: " << timerStop() << std::endl;;

   //  C.print("Product from BlockMatrix multiplication");
   timerStart();
   matrix_product(c,a,a);
   std::cout << "VMatrix time:     " << timerStop() << std::endl;

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

   std::cout << "-----------------" << std::endl;
   if (res > 1e-8) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


int test_2(unsigned dim)
{

   std::cout << "==================================" << std::endl;
   std::cout << " test_2: Striped - Dense multiply " << std::endl;
   std::cout << "==================================" << std::endl;

   VMatrix a, b, c, d;
   std::vector<int> stripes = {-3,-4, 0 ,4 ,3 };
   
   a.init(dim, dim, stripes).bind(TestFunctor());
   b.init(dim, dim, VMatrix::Dense).bind(TestFunctor());
   c.init(dim, dim, VMatrix::Dense).bind(ZeroFunctor());
   d.init(dim, dim, VMatrix::Dense).bind(ZeroFunctor());

   std::cout << "Size:   " << a.nRows() << " x " << a.nCols() << std::endl;
   std::cout << "Press Enter to begin" << std::endl;
   int ch = getchar();
   

   struct timeval  tv1, tv2;
   struct timezone tz;
   std::clock_t clock_start = std::clock();

   gettimeofday(&tv1, &tz);
      matrix_product(c,a,b);
   std::clock_t clock_end = std::clock();

   gettimeofday(&tv2, &tz);
   double elapsed = (double) (tv2.tv_sec -tv1.tv_sec) + 
                    (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
   std::cout << "Striped Matrix time: " << elapsed << std::endl;
   std::cout << "Clock time:          " << 1000.0 * (clock_end-clock_start) / CLOCKS_PER_SEC 
             << " ms\n";

   return 0;

   a.toDense();
   timerStart();
   //matrix_product(d,a,b);
   std::cout << "Dense Matrix time:   " << timerStop() << std::endl;

   double* data_c(c.data());
   double* data_d(d.data());
   double  res(0);
   unsigned n(c.nCols()*c.nRows());

   for (unsigned i = 0; i < n; ++i) {
       res = std::max(res, std::abs( data_c[i] - data_d[i]));
   }

   std::cout << "Max deviation:    " << res << std::endl;

   std::cout << "-----------------" << std::endl;
   if (res > 1e-8) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


int main()
{
   std::cout << "Running tests:" << std::endl;

   unsigned total(std::pow(2,10));
   int ok(0);

   for (unsigned p = 14; p < 15; ++p) {
       ok += test_2(std::pow(2,p));
   }

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
