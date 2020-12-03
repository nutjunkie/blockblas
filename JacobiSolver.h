#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H

#include <iostream>

#include "BlockMatrix.h"
#include "MatMult.h"

#include "CMTile.h"
#include "TileArray.h"
#include "Log.h"

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

   unsigned iter(0);
   for (iter = 0; iter < MAX_ITER; ++iter) {
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




// Solves A.x = b using an interative Jacobi algorithm
template <class T>
int jacobi_solver(TileArray<T> const& A, TileArray<T>& x, TileArray<T> const& b)
{
   //A.print("A matrix ---");
   //b.print("b vector ---");
   //x.print("x vector ---");

   // form the diagonal inverse matrices
   unsigned nTiles(A.nRowTiles());
   TileArray<T> Aii(nTiles, 1);

   for (int i = 0; i < nTiles; ++i) {
       CMTile<T> const& t = dynamic_cast<CMTile<T> const&>(A(i,i));
       CMTile<T>* tile = new CMTile<T>(t);
       tile->invert();
       Aii.set(i,0,tile);
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   TileArray<T> work(x);
   TileArray<T> lastx(x);

   for (unsigned iter = 0; iter < MAX_ITER; ++iter) {
       work.fill();
       product_sans_diagonal(A, x, work);

       work.scale(-1.0);
       work += b;

	   // This needs to be swapped out for a linear solve using LU
	   // decomposition.
       for (int i = 0; i < nTiles; ++i) {
           tile_product(Aii(i), work(i), 0.0, x(i));
       }

       work  = lastx;
       work -= x;

       double res(0.0);
       for (int i = 0; i < nTiles; ++i) {
           res += work(i).norm2();
           lastx(i) = x(i);
       }

       std::cout << std::scientific;
       std::cout << "Iter: " << iter << " Vector residual = " << std::sqrt(res) << std::endl;
       if (std::sqrt(res) < 1e-8) {
          std::cout << "CONVERGED, iterations "<< iter << std::endl;
          return iter;
       }
   }

   return -MAX_ITER;
}


template <class T>
void lu_solve(Tile<T> const& A, Tile<T>& b, int* ipiv);


// Solves A.x = b using an interative Jacobi algorithm
template <class T>
int jacobi_solverLU(TileArray<T> const& A, TileArray<T>& x, TileArray<T> const& b)
{
   int rc = -MAX_ITER;
   //A.print("A matrix ---");
   //b.print("b vector ---");
   //x.print("x vector ---");

   // form the diagonal inverse matrices
   unsigned nTiles(A.nRowTiles());
   TileArray<T> Aii(nTiles, 1);

   // We need to save the ipiv values for lu_solve
   int** ipiv = new int*[nTiles]; 

   for (int i = 0; i < nTiles; ++i) {
       CMTile<T> const& t = dynamic_cast<CMTile<T> const&>(A(i,i));
       CMTile<T>* tile = new CMTile<T>(t);
       tile->print("before factorization");
       ipiv[i] = new int[tile->nRows()];
       tile->factorLU(ipiv[i]);
       tile->print("after factorization");
       Aii.set(i,0,tile);
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   TileArray<T> work(x);
   TileArray<T> lastx(x);

   for (unsigned iter = 0; iter < MAX_ITER; ++iter) {
       work.fill();
       // A.x = work
       product_sans_diagonal(A, x, work);

       work.scale(-1.0);
       work += b;

       for (int i = 0; i < nTiles; ++i) {
           lu_solve(Aii(i), work(i), ipiv[i]);
       }

       x    = work;
       work = lastx;
       work -= x;

       double res(0.0);
       for (int i = 0; i < nTiles; ++i) {
           res += work(i).norm2();
           lastx(i) = x(i);
       }

       std::cout << std::scientific;
       std::cout << "Iter: " << iter << " Vector residual = " << std::sqrt(res) << std::endl;
       if (std::sqrt(res) < 1e-8) {
          std::cout << "CONVERGED, iterations "<< iter << std::endl;
          rc = iter;
          goto cleanup;
       }
   }

   cleanup:
      for (int i = 0; i < nTiles; ++i) {
          delete [] ipiv[i];
      }
      delete [] ipiv;

   return rc;
}

#endif
