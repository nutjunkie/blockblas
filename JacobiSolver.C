#include "JacobiSolver.h"
#include "BlockMatrix.h"
#include "MatMult.h"
#include <iostream>
#include <veclib/veclib.h>


#define MAX_ITER  100

// Solves
//    A.x = b

void jacobi_solver(BlockMatrix<double>& x,  BlockMatrix<double> const& A, BlockMatrix<double> const& b)
{
   ZeroFunctor zeroFunctor;

   //A.print("A matrix");
   //b.print("b vector");

   x(0,0).set(0,0,1.0);

   // form the diagonal inverse matrices
   unsigned nBlocks(A.nRowBlocks());
   BlockMatrix<double> Aii(nBlocks, 1);

   for (int i = 0; i < nBlocks; ++i) {
       Aii(i) = A(i,i);
       Aii(i).invert();
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   BlockMatrix<double> work(x);
   BlockMatrix<double> lastx(x);

   for (unsigned iter = 0; iter < MAX_ITER; ++iter) {
       work.bind(zeroFunctor);
       matrix_product_sans_diagonal(work,A,x);
       -work;
       work += b;

	   // This needs to be swapped out for a linear solve using LU
	   // decomposition.
       x.bind(zeroFunctor);
       for (int i = 0; i < nBlocks; ++i) {
           matrix_product(x(i),Aii(i),work(i));
       }

       work = lastx;
       work -= x;

       double res(0.0);
       for (int i = 0; i < nBlocks; ++i) {
           res += work(i).norm2();
           lastx(i) = x(i);
       }

       std::cout << std::scientific;
       std::cout << "Iter: " << iter << " Vector residual = " << std::sqrt(res) << std::endl;
       if (std::sqrt(res) < 1e-8) {
          std::cout << "CONVERGED, iterations "<< iter << std::endl;
          break;
       }
   }
}
