#include "BlockMatrix.h"
#include "MatMult.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sys/time.h>


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


//Timing test
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
   
   struct timeval tv1, tv2;
   struct timezone tz;
   double elapsed;
   
   timerStart();
   matrix_product(C,A,A);
   elapsed = timerStop();
   std::cout << "BlockMatrix time: " << elapsed << std::endl;;

   //  C.print("Product from BlockMatrix multiplication");
   timerStart();
   matrix_product(c,a,a);
   elapsed = timerStop();
   std::cout << "VMatrix time:     " << elapsed << std::endl;;

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


int main()
{
   std::cout << "Running tests:" << std::endl;

   unsigned total(std::pow(2,10));
   int ok(0);

   for (unsigned p = 1; p < 8; ++p) {
       unsigned blocks(std::pow(2,p));
       ok += test_1(blocks, total/blocks);
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
