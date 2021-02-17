#ifndef CONJUGATESOLVER_H
#define CONJUGATESOLVER_H

#include <iostream>
#include "CMTile.h"
#include "TileArray.h"
#include "TileProduct.h"



// Solves (A+root*I).x = b using the conjugate gradient method.
// Allows for different types of A and X
template <template<class> class TT, class T, class U>
int conjugate_gradient(TT<U> const& A, TileArray<T>& X, TileArray<T> const& B, T const root)
{
   unsigned N(B.nRows());
   unsigned nRHS(B.nCols());
   unsigned nBlocks(A.nRowTiles());

   T zero(0.0);
   T one(1.0);

   // S template
   TileArray<T> S(nBlocks,1);

   DiagonalTile<T>** M = new DiagonalTile<T>*[nBlocks];

   #pragma omp parallel for 
   for (unsigned i = 0; i < nBlocks; ++i) {
       size_t n(A(i,i).nRows());
       M[i] = new DiagonalTile<T>(A(i,i));
       M[i]->addToDiag(root);
       M[i]->invert();
       S.set(i,0, new CMTile<T>(B(i,0).nRows(),B(i,0).nCols()));
   }


   // The product() function assumes the data in these arrays is contiguous
   // The CMTiles are just a convenient way of ensuring this.

   TileArray<T> P(S);
   CMTile<T>    p(N,nRHS);
   p.alloc();
   P.bind(p.data());

   TileArray<T> Q(S);
   CMTile<T>    q(X);
   Q.bind(q.data());

   TileArray<T> R(S);
   CMTile<T>    r(B);
   R.bind(r.data());

   TileArray<T> Z(S);
   CMTile<T>    z(N,nRHS);
   z.alloc();
   Z.bind(z.data());


   // Tiny tiles
   TileArray<T> Lambda(1,1);
   Lambda.set(0,0, new CMTile<T>(nRHS,nRHS));
   Lambda(0,0).alloc();

   TileArray<T> Psi(1,1);
   Psi.set(0,0, new CMTile<T>(nRHS,nRHS));
   Psi(0,0).alloc();
   
   TileArray<T> RtZ(1,1);
   RtZ.set(0,0, new CMTile<T>(nRHS,nRHS));
   RtZ(0,0).fill();


   Q.scale(root);
   R.scale(-one);
   R += Q;
   product(A,X,R);

   #pragma omp parallel for 
   for (unsigned i = 0; i < nBlocks; ++i) {
       tile_product(*M[i], R(i,0), zero, P(i,0));
   }

   product(R, P, RtZ, CblasTrans); // rtz <- r^t.z                        //zgemm

   P.scale(-one);

   double res(0.0);
   unsigned iter(0);

   for (iter = 0; iter < MAX_ITER; ++iter) {
       // Line 1:
       memcpy(q.data(), p.data(), N*nRHS*sizeof(T));
       Q.scale(root);
       product(A,P,Q);                        // Q <- (A-root I).P

       // Line 2:
       Psi(0,0).fill();
       product(P, Q, Psi, CblasTrans);        // psi <- p^t.q         //zgemm
       
       Psi(0,0).invert();
       Lambda(0,0).fill();
       product(Psi, RtZ, Lambda);             // Lambda = Psi.RTZ     //zgemm 

       // Line 3:
       product(P,Lambda,X);                   // X <- X + P.Lambda    //zgemm

       // Line 4:
       product(Q,Lambda,R);                   // R <- R + Q.Lambda    //zgemm

       res = std::sqrt(R.norm2());
       if (res < 1e-10) break;

       // Line 5:
       #pragma omp parallel for 
       for (unsigned i = 0; i < nBlocks; ++i) {
           tile_product(*M[i], R(i,0), zero, Z(i,0));
       }

       // Line 6:
       Lambda = RtZ; 
       RtZ(0,0).fill();
       product(R, Z, RtZ, CblasTrans);  // rtz <- r^t.z

       // Line 7:
       Lambda(0,0).invert();
       Psi(0,0).fill();
       product(Lambda, RtZ, Psi);

       // Line 8:
       Q.fill();
       product(P, Psi, Q);    // Q = P.Psi

       memcpy(p.data(), q.data(), N*nRHS*sizeof(T));
       P -= Z;
   }

   if (iter < MAX_ITER) {
      std::cout << std::fixed << std::showpoint << std::setprecision(10);
      std::cout << "CONVERGED, iterations "<< iter << "  residue: " << res << std::endl;
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

   cleanup: 
      for (unsigned i = 0; i < nBlocks; ++i) {
          delete M[i];
      }
      delete []  M;

   return iter;
}


#endif
