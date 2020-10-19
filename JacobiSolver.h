#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H

#include <iostream>
#include <veclib/veclib.h>
#include "BlockMatrix.h"
#include "MatMult.h"

#define MAX_ITER      100

// Solves  A.x = b

template <class T>
void jacobi_solver(BlockMatrix<T>& x,  BlockMatrix<T> const& A, BlockMatrix<T> const& b)
{
   ZeroFunctor<T> zeroFunctor;

   //A.print("A matrix");
   //b.print("b vector");

   x(0,0).set(0,0,1.0);

   // form the diagonal inverse matrices
   unsigned nBlocks(A.nRowBlocks());
   BlockMatrix<T> Aii(nBlocks, 1);

   for (int i = 0; i < nBlocks; ++i) {
       Aii(i) = A(i,i);
       Aii(i).invert();
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   BlockMatrix<T> work(x);
   BlockMatrix<T> lastx(x);

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
#endif
