#include <iostream>
#include <veclib/veclib.h>

#include "VMatrix.h"
#include "Timer.h"
#include "util.h"


int test_dgemm(unsigned n, unsigned dim)
{
   print_header(1, "DGEMM");

   VMatrix A, B, C;
   A.init(dim, dim, VMatrix::Dense).bind(DebugFunctor());
   B.init(dim, dim, VMatrix::Dense).bind(DebugFunctor());
   C.init(dim, dim, VMatrix::Dense).bind(ZeroFunctor());

   std::cout << "Size:   " << A.nRows() << " x " << A.nCols() << std::endl;
   std::cout << "Performing " << n << " iterations for time-averaging" << std::endl;

   Timer timer;
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nCols(),
            B.data(), B.nRows(), 1.0, C.data(), C.nRows());
   }   
   double elapsed(timer.stop());
   unsigned flops = 2*A.nRows()*A.nCols()*B.nCols() - A.nRows()*B.nCols();
   std::cout << " Average time:  " << timer.format(elapsed/n) << std::endl;
   std::cout << " Average FLOPS: " << n*flops/elapsed << std::endl;
   return 0;
}


int test_dgemv(unsigned n, unsigned dim)
{
   print_header(1, "DGEMV");

   VMatrix A, B, C;
   A.init(dim, dim, VMatrix::Dense).bind(DebugFunctor());
   B.init(dim, 1,   VMatrix::Dense).bind(DebugFunctor());
   C.init(dim, 1,   VMatrix::Dense).bind(ZeroFunctor());

   std::cout << "Size:   " << A.nRows() << " x " << A.nCols() << std::endl;
   std::cout << "Performing " << n << " iterations for time-averaging" << std::endl;

   Timer timer;
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       cblas_dgemv(CblasRowMajor, CblasNoTrans,
             A.nRows(), A.nCols(), 1.0, A.data(), A.nCols(),
             B.data(), 1, 1.0, C.data(), 1);
   }   
   double elapsed(timer.stop());
   unsigned flops = 2*A.nRows()*A.nCols()*B.nCols() - A.nRows()*B.nCols();
   std::cout << " Average time:  " << timer.format(elapsed/n) << std::endl;
   std::cout << " Average FLOPS: " << n*flops/elapsed << std::endl;
   return 0;
}


int test_daxpy(unsigned n, unsigned dim)
{
   print_header(1, "DAXPY");

   VMatrix A, C;
   A.init(dim, 1, VMatrix::Dense).bind(DebugFunctor());
   C.init(dim, 1,   VMatrix::Dense).bind(ZeroFunctor());

   std::cout << "Size:   " << A.nRows() << " x " << A.nCols() << std::endl;
   std::cout << "Performing " << n << " iterations for time-averaging" << std::endl;

   Timer timer;
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       cblas_daxpy(dim, 1.1, A.data(), 1, C.data(), 1);
   }   
   double elapsed(timer.stop());
   unsigned flops = 2*A.nRows() - 1;
   std::cout << " Average time:  " << timer.format(elapsed/n) << std::endl;
   std::cout << " Average FLOPS: " << n*flops/elapsed << std::endl;
   return 0;
}


int main(int argc, char **argv)
{
   unsigned dim(1024);

   if (argc >= 2) {
      std::istringstream iss(argv[1]);
      if (iss >> dim) {
         // Conversion successful
         if (dim < 100) {
            dim = std::pow(2,dim);
          }
          //test_dgemv(10, dim);
          test_daxpy(10, dim);
        }else {
           std::cerr << "Invalid argument" << std::endl;
        }
   }else {
      std::cerr << "Feed me moar!!" << std::endl;
   }

   return 0;
}
