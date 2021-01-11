#ifndef CONJUGATESOLVER_H
#define CONJUGATESOLVER_H

#include <iostream>
#include "CMTile.h"
#include "TileArray.h"
#include "TileProduct.h"



// Solves (A+root*I).x = b using the conjugate gradient method.
// Allows for different types of A and X
template <class T, class U>
int conjugate_gradient(TileArray<U> const& A, TileArray<T>& X, TileArray<T> const& B, T const root)
{
   //A.print("A matrix ---");
   //B.print("b vector ---");
   //X.print("x vector ---");
   //
   
   unsigned N(B.nRows());
   unsigned nRHS(B.nCols());
   unsigned nBlocks(A.nRowTiles());

   T zero(0.0);
   T one(1.0);

   // S template
   TileArray<T> S(nBlocks,1);

   //CMTile<T>** M = new CMTile<T>*[nBlocks];
   DiagonalTile<T>** M = new DiagonalTile<T>*[nBlocks];

   for (unsigned i = 0; i < nBlocks; ++i) {
       //M[i] = new CMTile<T>(A(i,i));
       M[i] = new DiagonalTile<T>(A(i,i));
       M[i]->addToDiag(root);
       M[i]->invert();
       S.set(i,0, new CMTile<T>(B(i,0).nRows(),B(i,0).nCols()));
   }

   // P Array
   TileArray<T> P(S);
   CMTile<T>    p(N,nRHS);
   p.alloc();
   P.bind(p.data());

   // Q Array (AP)
   TileArray<T> Q(S);
   CMTile<T>    q(X);
   Q.bind(q.data());

   // R Array (residual)
   TileArray<T> R(S);
   CMTile<T>    r(B);
   R.bind(r.data());

   // Z Array (residual)
   TileArray<T> Z(S);
   CMTile<T>    z(N,nRHS);
   z.alloc();
   Z.bind(z.data());

   // Lambda
   TileArray<T> Lambda(1,1);
   Lambda.set(0,0, new CMTile<T>(nRHS,nRHS));
   CMTile<T> lambda(nRHS,nRHS);
   lambda.alloc();
   Lambda.bind(lambda.data());
   
   // small arrays
   CMTile<T> rtz(nRHS,nRHS);
   CMTile<T> psi(nRHS,nRHS);

   rtz.alloc();
   psi.alloc();

   Q.scale(root);
   R.scale(-one);
   r += q;
   product(A,X,R);
   //R.print("R mat");

   for (unsigned i = 0; i < nBlocks; ++i) {
       tile_product(*M[i], R(i,0), zero, Z(i,0));
       //M[i]->print("M inverse");
   }
   //Z.print("Z mat");

   //std::cout << "Pointers: " << r.data() << " <-> " << R(0,0).data() << std::endl;
   //r.print("r");
   tile_product(r, z, zero, rtz, CblasTrans); // rtz <- r^t.z
   //rtz.print("RtZ");

   p = z;
   p.scale(-one);
   //P.print("P mat");

   double res(0.0);
   unsigned iter(0);

   //A.addToDiag(root);
   for (iter = 0; iter < MAX_ITER; ++iter) {
       // Line 1:
       q = p;       
       q.scale(root);
       product(A,P,Q);                           // Q <- A.P
       //Q.print("Q mat");

       // Line 2:
       tile_product(p, q, zero, psi, CblasTrans); // psi <- p^t.q
       psi.invert();
       tile_product(psi, rtz, zero, lambda);
       //Lambda.print("lambda mat");

       // Line 3:
       product(P,Lambda,X);  // X <- X + P.Lambda
       //X.print("X mat");

       // Line 4:
       product(Q,Lambda,R);  // R <- R + Q.Lambda
       //R.print("R mat");

       res = std::sqrt(r.norm2());
       //std::cout << iter << " Residue = " << res << std::endl;
       if (res < 1e-10) break;


       // Line 5:
#pragma omp parallel for 
       for (unsigned i = 0; i < nBlocks; ++i) {
           tile_product(*M[i], R(i,0), zero, Z(i,0));
       }

       // Line 6:
       lambda = rtz; 
       tile_product(r, z, zero, rtz, CblasTrans); // rtz <- r^t.z
       //rtz.print("RtZ mat");

       // Line 7:
       lambda.invert();
       tile_product(lambda, rtz, zero, psi);
       //psi.print("psi mat");

       // Line 8:
       tile_product(p, psi, zero, q);
       p = q;
       p -= z;
       //q.print("q mat");
       //p.print("p mat");
   }

   if (iter < MAX_ITER) {
      std::cout << std::fixed << std::showpoint << std::setprecision(10);
      std::cout << "CONVERGED, iterations "<< iter << "  residue: " << res;
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

   cleanup: 
      for (unsigned i = 0; i < nBlocks; ++i) {
          delete M[i];
      }
      delete []  M;

   return iter;
}


#endif
