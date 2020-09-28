#include "JacobiSolver.h"
#include "BlockMatrix.h"
#include "MatMult.h"
#include <iostream>
#include <veclib/veclib.h>


// Solves
//    A.x = b

void jacobi_solver(BlockMatrix& x,  BlockMatrix const& A, BlockMatrix const& b)
{

   ZeroFunctor zeroFunctor;

   A.print("A matrix");
   b.print("b vector");

   x(0,0).set(0,0,1.0);

   // form the diagonal inverse matrices
   unsigned nBlocks(A.nRowBlocks());
   BlockMatrix Aii(nBlocks, 1);

   for (int i = 0; i < nBlocks; ++i) {
       Aii(i) = A(i,i);
       Aii(i).invert();
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   BlockMatrix work(x);
   BlockMatrix lastx(x);
   lastx.bind(zeroFunctor);

   for (unsigned iter = 0; iter < 10; ++iter) {
       work.bind(zeroFunctor);
       matrix_product_sans_diagonal(work,A,x);
       -work;
       work += b;

       x.bind(zeroFunctor);
       for (int i = 0; i < nBlocks; ++i) {
           matrix_product(x(i),Aii(i),work(i));
       }

       x.print("updated x vector");
       lastx.print("lastx vector");

       work = lastx;
       work -= x;
       work.print("residual  vector");

       double res(0.0);
       for (int i = 0; i < nBlocks; ++i) {
           res += work(i).norm2();
           lastx(i) = x(i);
       }
       lastx.print("lastx now updated ");

       std::cout << "Vector residual = " << std::sqrt(res) << std::endl;

   }
}
