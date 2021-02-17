#include "Davidson.h"
#include "TileProduct.h"

#define MAX_ITERATIONS 50

typedef TileArray<double> TA;
typedef DiagonalTile<double> DT;



void axpy(double* a, double const beta, double* const b, size_t const n)
{
   for (size_t i = 0; i < n; ++i) {
       a[i] += beta*b[i];
   }
}


double dot(double const* a, double const* b, size_t const n)
{
   double d(0.0);
   for (size_t i = 0; i < n; ++i) {
       d += a[i]*b[i];
   }
   return d;
}



int DavidsonMethod(TA const& A, TA const& vec0, double lam0)
{
   int rv(0);

   unsigned const maxIter(20);
   unsigned N(A.nRows());
   unsigned nBlocks(A.nRowTiles());
   unsigned ku(0), kw(0);

   double* U = new double[maxIter*N];
   double* W = new double[maxIter*N];
   double* y = new double[N];
   double* r = new double[N];
   double* work = new double[maxIter];
   double  lam1;
   double  zero(0.0);


   CMTile<double> B(maxIter,maxIter);
   B.fill();

   // Digaonal Tiles:
   DT** D = new DT*[nBlocks];
   DT** M = new DT*[nBlocks];

   #pragma omp parallel for 
   for (unsigned i = 0; i < nBlocks; ++i) {
       D[i] = new DT(A(i,i));
   }



   // Set up the block structure for u and w (lazy)
   TA uj(vec0);
   TA wj(vec0);

   memcpy(U, vec0(0).data(), N*sizeof(double));

   for (int j = 0; j < maxIter-1; ++j) {
       // w_j = A u_j
       memset(W+j*N, 0, N*sizeof(double));
       uj.bind(U+j*N);
       wj.bind(W+j*N);
       product(A,uj,wj);

       //uj.print("Print u vec");
       //wj.print("Print w vec");

       // augment B(i,j)
       for (int k = 0; k < j; ++k) {
           B.set(k,j, dot(U+k*N, W+j*N, N));
           // std::cout << "setting B(" << k << "," << j << ") = " << B(k,j) <<  " = " << dot(U+k*N, W+j*N, N) << std::endl;
           B.set(j,k, dot(U+j*N, W+k*N, N));
       }
       B.set(j,j, dot(U+j*N, W+j*N, N));

       // copy to b(i,j), this will get overwritten with the eigenvectors
       CMTile<double> b(j+1,j+1);
       b.alloc();

       for (unsigned i = 0; i <= j; ++i) {
           for (unsigned k = 0; k <= j; ++k) {
               //std::cout << "setting b(" << i << "," << k << ") = " << B(i,k) << std::endl;
               b.set(i,k, B(i,k));
           }
       }

       b.print("b for eigen:");

       // Eigenvalues of B
       int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', b.nRows(), b.data(), b.nRows(), work);
       lam1 = work[0];
       for (unsigned i = 0; i <= j; ++i) {
           std::cout << "work array: " << work[i] << std::endl;
       }

       std::cout << "First eigenvalue: " << lam1 << std::endl;

       if (std::abs(lam1-lam0) < 1e-10) {
          std::cout << "Converged: " << lam1 << std::endl;
          goto cleanup;
       }
       lam0 = lam1;

       // form y
       memset(y, 0, N*sizeof(double));
       for (unsigned i = 0; i <= j; ++i) {
           axpy(y, b(i,0), U+i*N, N);
       }

       // form r in wj
       memcpy(r, y, N*sizeof(double));
       uj.bind(y);
       wj.bind(r);
       wj.scale(-lam1);
       product(A,uj,wj);
       //wj.print("r vector :");

       // form t in uj
       #pragma omp parallel for 
       for (unsigned i = 0; i < nBlocks; ++i) {
          D[i]->addToDiag(-lam1);
          D[i]->invert();
          tile_product(*D[i], wj(i), zero, uj(i));
          D[i]->invert();
          D[i]->addToDiag(lam1);
       }

       CMTile<double> Uj(N,j+1);
       Uj.bind(U);
       CMTile<double> t(N,1);
       t.bind(y);
       CMTile<double> tmp(j+1,1);
       tmp.fill();
       tile_product(Uj, t, zero, tmp, CblasTrans);
       tile_product(Uj, tmp, -1.0, t, CblasNoTrans);

       double n2(t.norm2());
       t.scale(-1.0/std::sqrt(n2));

       memcpy(U+(j+1)*N, y, N*sizeof(double));
   }

   std::cout << "Davidson barffed: " << std::endl;

   cleanup:
      delete [] U;
      delete [] W;
      delete [] y;
      delete [] r;
      delete [] work;

   return rv;
}



int DavidsonIteration(TA const& A, TA const& vec0, double lam0)
{
   int rv(0);
   unsigned N(A.nRows());
   unsigned nBlocks(A.nRowTiles());

   double const thresh(1e-10);

   DT** D = new DT*[nBlocks];
   DT** M = new DT*[nBlocks];

   #pragma omp parallel for 
   for (unsigned i = 0; i < nBlocks; ++i) {
       D[i] = new DT(A(i,i));
       M[i] = new DT(A(i,i));
       D[i]->scale(-1.0);
   }

   TA v0(vec0);
   TA v1(vec0);
   TA Av(vec0);

   // normalize v0;
   double norm(v0.norm2());
   norm = std::sqrt(norm);
   v0.scale(1.0/norm);

   TA Lam1(1,1);
   Lam1.set(0,0, new CMTile<double>(1,1));
   Lam1(0,0).set(0,0,0.0);
   
   double zero(0.0);
   double one(1.0);

   for (unsigned it = 0; it < MAX_ITERATIONS; ++it) {

       v1.fill();
       product(A,v0,v1);

       #pragma omp parallel for 
       for (unsigned i = 0; i < nBlocks; ++i) {
           tile_product(*D[i], v0(i), one, v1(i));
           *M[i] = *D[i];
           M[i]->addToDiag(lam0);
           M[i]->invert();
           tile_product(*M[i], v1(i), zero, v0(i));
       }

       // normalize v0;
       norm = v0.norm2();
       norm = std::sqrt(norm);
       v0.scale(1.0/norm);

       v1.fill();
       Lam1(0,0).set(0,0,0.0);
       product(A, v0, v1);

       //v0.print("v0");
       //v1.print("v1");
       //Lam1.print("lam");
       //double lam1(Lam1(0,0).get(0,0));
       double lam1(0.0);

       for (unsigned i = 0; i < nBlocks; ++i) {
           for (unsigned j = 0; j < v0(i).nRows(); ++j) {
               lam1 += v0(i).get(j)*v1(i).get(j);
           }
       } 

       std::cout << "Eigenvalue: " << lam1 << std::endl;
      
       if (std::abs(lam1-lam0) < thresh) {
          std::cout << "Eigenvalue converged in " << it 
                    << " iterations: " << lam1 << std::endl;
          goto cleanup;
       }

       lam0 = lam1;
   }

   std::cout << "Error: Max iterations reached in DavidsonIteration " << lam0 << std::endl;
   rv = 1;

   cleanup:
      for (unsigned i = 0; i < nBlocks; ++i) {
          delete D[i];
          delete M[i];
      }
      delete [] D;
      delete [] M;

  return rv;
}
