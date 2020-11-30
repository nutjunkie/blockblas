#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H

#include <iostream>

#include "BlockMatrix.h"
#include "MatMult.h"

#define MAX_ITER      100

// Solves  A.x = b

template <class T,LayoutT L>
void jacobi_solver(BlockMatrix<T,L>& x,  BlockMatrix<T,L> const& A, BlockMatrix<T,L> const& b)
{
   ZeroFunctor<T> zeroFunctor;

   x(0,0).set(0,0,T(1.0));

   //A.print("A matrix ---");
   //b.print("b vector ---");
   //x.print("x vector ---");
   //
   // form the diagonal inverse matrices
   unsigned nBlocks(A.nRowBlocks());
   BlockMatrix<T,L> Aii(nBlocks, 1);

   for (int i = 0; i < nBlocks; ++i) {
       Aii(i) = A(i,i);
       Aii(i).invert();
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   BlockMatrix<T,L> work(x);
   BlockMatrix<T,L> lastx(x);

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
          return;
          BlockMatrix<T,L> result(x);
          result.bind(zeroFunctor);
          matrix_product(result, A, x);
          result.print("Prodct A.x");
          b.print("RHS b");
          
       }
   }
   std::cout << "FAILED TO CONVERGED, iterations "<< iter << std::endl;
}

#endif
