#include "Davidson.h"
#include "TileProduct.h"

#define DEBUG false

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


void normalize(double* vec, size_t const n)
{
   double d(dot(vec,vec,n));
   d = 1.0/std::sqrt(d);
   for (size_t i = 0; i < n; ++i) {
       vec[i] *= d;
   }
}

template<>
int DavidsonMethod(TileArray<double> const& A, TA const& vec0, double lam0)
{
    std::cout << "Need to convert to SymmetricTileArray" << std::endl;
    return 1;
}


template<>
int DavidsonMethod(SymmetricTileArray<double> const& A, TA const& vec0, double lam0)
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

   #pragma omp parallel for 
   for (unsigned i = 0; i < nBlocks; ++i) {
       D[i] = new DT(A(i,i));
   }

   // Set up the block structure for u and w (lazy)
   TA uj(vec0);
   TA wj(vec0);

   memcpy(U, vec0(0).data(), N*sizeof(double));

   normalize(U,N);

   for (int j = 0; j < maxIter-1; ++j) {
       // w_j = A u_j
       memset(W+j*N, 0, N*sizeof(double));
       uj.bind(U+j*N);
       wj.bind(W+j*N);
       product(A,uj,wj);

       // augment B(i,j)
       for (int k = 0; k < j; ++k) {
           B.set(k,j, dot(U+k*N, W+j*N, N));
           B.set(j,k, dot(U+j*N, W+k*N, N));
       }
       B.set(j,j, dot(U+j*N, W+j*N, N));

       // copy to b(i,j), this will get overwritten with the eigenvectors
       CMTile<double> b(j+1,j+1);
       b.alloc();

       for (unsigned i = 0; i <= j; ++i) {
           for (unsigned k = 0; k <= j; ++k) {
               b.set(i,k, B(i,k));
           }
       }

       if (DEBUG) b.print("b for eigen:");

       // Eigenvalues of B
       int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', b.nRows(), b.data(), b.nRows(), work);
       lam1 = work[0];

       std::cout << std::fixed << std::showpoint << std::setprecision(10);
       std::cout << "Iter: " << j+1 << " ev = " << lam1 << std::endl;

       if (std::abs(lam1-lam0) < 1e-10) {
          std::cout << std::fixed << std::showpoint << std::setprecision(10);
          std::cout << "Converged in " << j+1 << " iterations: " << lam1 << std::endl;
          goto cleanup;
       }else if(j == 0) {
          lam1 *= 1.00001;
       }
       lam0 = lam1;

       if (DEBUG) b.print("b from eigen:");


       // form y
       memset(y, 0, N*sizeof(double));
       for (unsigned i = 0; i <= j; ++i) {
           axpy(y, b(i), U+i*N, N);
       }

       // form r in wj
       memcpy(r, y, N*sizeof(double));
       uj.bind(y);
       wj.bind(r);
       wj.scale(-lam1);
       product(A,uj,wj);

       // form t in uj (12)
       #pragma omp parallel for 
       for (unsigned i = 0; i < nBlocks; ++i) {
          D[i]->addToDiag(-lam1);
          D[i]->invert();
          tile_product(*D[i], wj(i), zero, uj(i));
          D[i]->invert();
          D[i]->addToDiag(lam1);
       }

       // project out U from t (13)
       CMTile<double> Uj(N,j+1);
       Uj.bind(U);
       CMTile<double> T(N,1);
       T.bind(y);
       CMTile<double> Tmp(j+1,1);
       Tmp.alloc();

       tile_product(Uj, T, zero, Tmp, CblasTrans);
       if (DEBUG) Tmp.print("Projections");

       tile_product(Uj, Tmp, -1.0, T, CblasNoTrans);

       normalize(y,N);
       memcpy(U+(j+1)*N, y, N*sizeof(double));

       if (DEBUG) {
          Uj.resize(N,j+2);
          Uj.bind(U);
          Tmp.resize(j+2,j+2);
          Tmp.alloc();
          tile_product(Uj, Uj, zero, Tmp, CblasTrans);
          Tmp.print("Ortho mat");
       }
   }

   std::cout << "Davidson barffed: " << std::endl;

   cleanup:
      delete [] U;
      delete [] W;
      delete [] y;
      delete [] r;
      delete [] work;
      for (unsigned i = 0; i < nBlocks; ++i) {
          delete D[i];
      }

      delete [] D;


   return rv;
}
