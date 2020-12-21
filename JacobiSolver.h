#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H

#include <iostream>
#include "CMTile.h"
#include "TileArray.h"
#include "TileProduct.h"


template <class T>
int diagonalDominance(TileArray<T> const& A)
{
   //A.info("A matrix", std::cout);
   //A.print("A matrix", std::cout);
   unsigned n(A.nRows());

   CMTile<T> a(A);

   for (unsigned j = 0; j < n; ++j) {
       double diag(std::abs(a(j,j)));
       double colsum(0.0); 

       for (unsigned i = 0; i < n; ++i) {
           colsum += std::abs(a(i,j));
       }

       colsum -= diag;
      
       if (diag < colsum) {
          std::cout << "Colsum " << j << "  " << diag << "  " << colsum;
          std::cout << "     !!!!" << std::endl;
       }
   }
}



// Solves A.x = b using an interative Jacobi algorithm
template <class T>
int jacobi_solver(TileArray<T> const& A, TileArray<T>& x, TileArray<T> const& b)
{
   //diagonalDominance(A);

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
   double res(0.0);
   unsigned iter(0);

   for (iter = 0; iter < MAX_ITER; ++iter) {
       work.fill();
       product_sans_diagonal(A, x, work);

       work.scale(-1.0);
       work += b;

	   // Swap this out for a linear solve using LU decomposition.
       for (int i = 0; i < nTiles; ++i) {
           tile_product(Aii(i), work(i), T(0.0), x(i));
       }

       work  = lastx;
       work -= x;

       res = 0.0;
       for (int i = 0; i < nTiles; ++i) {
           res += work(i).norm2();
           lastx(i) = x(i);
       }
       res = std::sqrt(res);

       //std::cout << std::scientific;
       //std::cout << "Iter: " << iter << " Vector residual = " << res << std::endl;
       if (res < 1e-8) break;
   }

   if (iter < MAX_ITER) {
      std::cout << "CONVERGED, iterations "<< iter << std::endl;
   }else {
      std::cout << "FAILED to converge: " << res << std::endl;
      iter *= -1;
   }

   if (false) {
      work.fill();
      product(A, x, work);
      work -= b;
      std::cout <<"product norm: " << work.norm2() << std::endl;
   }


   return iter;
}


template <class T>
void lu_solve(Tile<T> const& A, Tile<T>& b, int* ipiv);


// Solves A.x = b using an interative Jacobi algorithm
template <class T>
int jacobi_solverLU(TileArray<T> const& A, TileArray<T>& x, TileArray<T> const& b)
{
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
       //tile->print("before factorization");
       ipiv[i] = new int[tile->nRows()];
       tile->factorLU(ipiv[i]);
       //tile->print("after factorization");
       Aii.set(i,0,tile);
   }

   //Aii.info("diags Matrix");
   //Aii.print("diags Matrix");

   TileArray<T> work(x);
   TileArray<T> lastx(x);

   double res(0.0);
   unsigned iter(0);

   for (iter = 0; iter < MAX_ITER; ++iter) {
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

       res = 0.0;
       for (int i = 0; i < nTiles; ++i) {
           res += work(i).norm2();
           lastx(i) = x(i);
       }

       //std::cout << std::scientific;
       //std::cout << "Iter: " << iter << " Vector residual = " << std::sqrt(res) << std::endl;
       if (std::sqrt(res) < 1e-8) break;
   }

   cleanup:
      for (int i = 0; i < nTiles; ++i) {
          delete [] ipiv[i];
      }
      delete [] ipiv;

   if (iter < MAX_ITER) {
      std::cout << "CONVERGED, iterations "<< iter << std::endl;
   }else {
      std::cout << "FAILED to converge: " << res << std::endl;
      iter *= -1;
   }

   return iter;
}

#endif
