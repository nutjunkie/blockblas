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

   for (unsigned i = 0; i < nBlocks; ++i) {
       size_t n(A(i,i).nRows());
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


   // Lambda
   TileArray<T> Lambda(1,1);
   Lambda.set(0,0, new CMTile<T>(nRHS,nRHS));
   CMTile<T> lambda(nRHS,nRHS);
   lambda.alloc();
   Lambda.bind(lambda.data());

   // Psi
   TileArray<T> Psi(1,1);
   Psi.set(0,0, new CMTile<T>(nRHS,nRHS));
   CMTile<T> psi(nRHS,nRHS);
   psi.alloc();
   Psi.bind(psi.data());

   
   // small arrays
   CMTile<T> rtz(nRHS,nRHS);

   rtz.alloc();

   Q.scale(root);
   R.scale(-one);
   R += Q;
   product(A,X,R);

   //Z Array (residual)
   TileArray<T> Z(S);
   CMTile<T>    z(N,nRHS);
   z.alloc();
   Z.bind(z.data());

   for (unsigned i = 0; i < nBlocks; ++i) {
       tile_product(*M[i], R(i,0), zero, Z(i,0));
   }

   //std::cout << "Pointers: " << r.data() << " <-> " << R(0,0).data() << std::endl;
   //r.print("r");
   tile_product(r, z, zero, rtz, CblasTrans); // rtz <- r^t.z                        //zgemm
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

       // Line 2:
       //Psi(0,0).fill();
       //product(P, Q, Psi, CblasTrans); // psi <- p^t.q                       //zgemm
       //Psi.print("First Psi");
       
       Psi(0,0).fill();
       tile_product(p, q, zero, psi, CblasTrans); // psi <- p^t.q                       //zgemm
       //Psi.print("Second psi");

       Psi(0,0).invert();
       tile_product(psi, rtz, zero, lambda);                                            //zgemm

       // Line 3:
       product(P,Lambda,X);  // X <- X + P.Lambda                                       //zgemm

       // Line 4:
       product(Q,Lambda,R);  // R <- R + Q.Lambda

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

       // Line 7:
       Lambda(0,0).invert();
       tile_product(lambda, rtz, zero, psi);

       // Line 8:
       tile_product(p, psi, zero, q);
       p = q;
       p -= z;
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
