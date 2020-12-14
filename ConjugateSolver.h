#ifndef CONJUGATESOLVER_H
#define CONJUGATESOLVER_H

#include <iostream>
#include "CMTile.h"
#include "TileArray.h"
#include "TileProduct.h"

#define MAX_ITER      20


// Solves A.x = b using the conjugate gradient method.
template <class T>
int conjugate_gradient(TileArray<T> const& A, TileArray<T>& X, TileArray<T> const& B)
{
   //A.print("A matrix ---");
   //B.print("b vector ---");
   //X.print("x vector ---");
   
   unsigned N(B.nRows());
   unsigned nRHS(B.nCols());
   unsigned nBlocks(A.nRowTiles());

   T zero(0.0);
   T one(1.0);

   // S template
   TileArray<T> S(nBlocks,1);
   for (unsigned i = 0; i < nBlocks; ++i) {
       S.set(i,0, new CMTile<T>(B(i,0).nRows(),B(i,0).nCols()));
   }

   // P Array
   TileArray<T> P(S);
   CMTile<T>    p(N,nRHS);
   p.alloc();
   P.bind(p.data());

   // Q Array (AP)
   TileArray<T> Q(S);
   CMTile<T>    q(N,nRHS);
   q.alloc();
   Q.bind(q.data());

   // R Array (residual)
   TileArray<T> R(S);
   CMTile<T>    r(B);
   R.bind(r.data());

   // Lambda
   TileArray<T> Lambda(1,1);
   Lambda.set(0,0, new CMTile<T>(nRHS,nRHS));
   CMTile<T> lambda(nRHS,nRHS);
   lambda.alloc();
   Lambda.bind(lambda.data());
   
   // small arrays
   CMTile<T> rtr(nRHS,nRHS);
   CMTile<T> psi(nRHS,nRHS);

   rtr.alloc();
   psi.alloc();

   R.scale(-one);
   product(A,X,R);

   //std::cout << "Pointers: " << r.data() << " <-> " << R(0,0).data() << std::endl;
   //r.print("r");
   tile_product(r, r, zero, rtr, CblasTrans); // rtr <- r^t.r
   //rtr.print("RtR");

   p = r;
   p.scale(-one);

   double res(0.0);
   unsigned iter(0);

   for (iter = 0; iter < MAX_ITER; ++iter) {
       // Line 1
       Q.fill();                                    // Q <- 0
       product(A,P,Q);                              // Q <- A.P
       //Q.print("Q mat");
       //
       tile_product(p, q, zero, psi, CblasTrans); // psi <- p^t.q
       psi.invert();
       tile_product(psi, rtr, zero, lambda);
       //Lambda.print("lambda mat");

       // Line 2
       product(P,Lambda,X);  // X <- X + P.Lambda
       //X.print("X mat");

       // Line 3
       product(Q,Lambda,R);  // R <- R + Q.Lambda
       //R.print("R mat");
       //std::cout << "Pointers: " << r.data() << " <-> " << R(0,0).data() << std::endl;

       // Line 4
       lambda = rtr; 
       lambda.invert();
       tile_product(r, r, zero, rtr, CblasTrans); // rtr <- r^t.r
       //rtr.print("RtR mat");

       res = std::sqrt(r.norm2());
       //std::cout << iter << " Residue = " << res << std::endl;
       if (std::sqrt(res) < 1e-8) break;

       tile_product(lambda, rtr, zero, psi);
       //psi.print("psi mat");

       // Line 5
       tile_product(p, psi, zero, q);
       //q.print("q mat");

       p = q;
       p -= r;
       //p.print("p mat");
   }

   if (iter < MAX_ITER) {
      std::cout << "CONVERGED, iterations "<< iter << "  residue:" << res;
   }else {
      std::cout << "FAILED to converge: " << res;
      iter *= -1;
   }

   if (false) {
      R.fill();
      product(A, X, R);
      R -= B;
      std::cout <<"  product norm: " << R.norm2();
   }

   std::cout << std::endl;



   return iter;
}

#endif
