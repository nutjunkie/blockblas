#ifndef TILEPRODUCT_H
#define TILEPRODUCT_H
/******************************************************************************
 * 
 *  Non-memeber function declarations for multiplying various tile types.
 * 
 *****************************************************************************/

#include <omp.h>
#include "CMTile.h"
#include "TileArray.h"
#include "StripedTile.h"
#include "SymmetricTileArray.h"


//Computes A.B = C
template <class T, class U>
void tile_product(CMTile<U> const& A, CMTile<T> const& B, 
   T const c, CMTile<T>& C, CBLAS_TRANSPOSE const Atrans = CblasNoTrans);


template <class T, class U>
void tile_product(DiagonalTile<U> const& A, CMTile<T> const& B, T const alpha, CMTile<T>& C, 
   CBLAS_TRANSPOSE const Atrans = CblasNoTrans) 
{
   unsigned nc(B.nCols());
   unsigned nr(A.nRows());

   unsigned ldb(B.leadingDim());
   unsigned ldc(C.leadingDim());

   U const* a(A.data());
   T const* b(B.data());
   T*       c(C.data());

   C.scale(alpha);

   for (unsigned j = 0; j < nc; ++j) {
       for (unsigned i = 0; i < nr; ++i) {
           c[i+ldc*j] += a[i] * b[i+j*ldb];
       }
   }
}


template <class T, class U>
void tile_product(StripedTile<U> const& A, CMTile<T> const& B, T const alpha, CMTile<T>& C,
   CBLAS_TRANSPOSE const Atrans = CblasNoTrans)
{
   unsigned rowsA(A.nRows());
   unsigned colsA(A.nCols());
   unsigned colsC(C.nCols());

   std::vector<int> stripes(A.stripes());

   if (Atrans == CblasTrans) {
      std::swap(rowsA,colsA);
      
      std::vector<int>::iterator iter;
      for (iter = stripes.begin(); iter != stripes.end(); ++iter) {
          *iter = -(*iter);
      }
   }

   unsigned ldb(B.leadingDim());
   unsigned ldc(C.leadingDim());

   unsigned m(std::min(rowsA,colsA));

   U const* a(A.data());
   T const* b(B.data());
   T*       c(C.data());

   C.scale(alpha);

   for (unsigned s = 0; s < stripes.size(); ++s) {
       int offset(stripes[s]);
       int offC(std::max(0,-offset));
       int offB(std::max(0, offset));
       // Contraction length
       int len = (offset < 0) ? std::min(rowsA+offset, colsA)
                              : std::min(rowsA, colsA-offset);
       U const* a0(&a[s*m]);
       for (unsigned j = 0; j < colsC; ++j) {
           T const* b0(&b[offB+j*ldb]);
           T*       c0(&c[offC+j*ldc]);
           for (unsigned i = 0; i < len; ++i) {
               c0[i] += a0[i]*b0[i];
           }   
       }   
   }   
}


template <class T, class U>
void tile_product(Tile<U> const& A, Tile<T> const& B, T const gamma, Tile<T>& C,
   CBLAS_TRANSPOSE const Atrans = CblasNoTrans)
{
   CMTile<T> const& b = dynamic_cast<CMTile<T> const&>(B);
   CMTile<T>& c = dynamic_cast<CMTile<T>&>(C);

   switch (A.storage()) {
      case Zero: {
         C.scale(gamma);
         return;
      } break;

      case Diagonal: {
         DiagonalTile<U> const& a = dynamic_cast<DiagonalTile<U> const&>(A);
         tile_product(a, b, gamma, c, Atrans);
      } break;

      case Striped: {
         StripedTile<U> const& a = dynamic_cast<StripedTile<U> const&>(A);
         tile_product(a, b, gamma, c, Atrans);
      } break;

      case CMDense: {
         CMTile<U> const& a = dynamic_cast<CMTile<U> const&>(A);
         tile_product(a, b, gamma, c, Atrans);
      } break;

      default:
        std::cerr << "ERROR: Unimplemented tile_product" << std::endl;
        break;
   }
}


//
// Note these product functions *accumulate* into C.
// C must be initialized if this is not what you want.
template <class T>
void product(Tile<T> const& A, Tile<T> const& B, Tile<T>& C)
{
   tile_product(A, B, T(1.0), C);
}


// Note these product functions *accumulate* into C.
// C must be initialized if this is not what you want.
template <class T, class U>
void product(TileArray<U> const& A, TileArray<T> const& B, TileArray<T>& C,
   CBLAS_TRANSPOSE const Atrans = CblasNoTrans)
{
   size_t nRowTiles(0);
   size_t nColTiles(0);

   if (Atrans == CblasNoTrans) {
      assert(A.nRowTiles() == C.nRowTiles());
      assert(A.nColTiles() == B.nRowTiles());
      assert(B.nColTiles() == C.nColTiles());
      nRowTiles = A.nRowTiles();
      nColTiles = B.nColTiles();
   }else {
      assert(A.nColTiles() == C.nRowTiles());
      assert(A.nRowTiles() == B.nRowTiles());
      assert(B.nColTiles() == C.nColTiles());
      nRowTiles = A.nColTiles();
      nColTiles = B.nColTiles();
   }

   unsigned bi, bj;
#pragma omp parallel for private(bj) collapse(2)
   for (bi = 0; bi < nRowTiles; ++bi) {
       for (bj = 0; bj < nColTiles; ++bj) {
           for (unsigned k = 0; k < A.nColTiles(); ++k) {
               // !!! Accumulate into C !!!
               tile_product(A(bi,k), B(k,bj), T(1.0), C(bi,bj), Atrans);
           }
       }
   }
}



template <class T>
void product_sans_diagonal(TileArray<T> const& A, TileArray<T> const& B, TileArray<T>& C)
{
   assert(A.nRowTiles() == C.nRowTiles());
   assert(A.nColTiles() == B.nRowTiles());
   assert(B.nColTiles() == C.nColTiles());
#pragma omp parallel for
   for (unsigned bi = 0; bi < A.nRowTiles(); ++bi) {
       for (unsigned bj = 0; bj < B.nColTiles(); ++bj) {
           for (unsigned k = 0; k < A.nColTiles(); ++k) {
               if (bi != k ) {
                  // !!! Accumulate into C !!!
                  tile_product(A(bi,k), B(k,bj), T(1.0), C(bi,bj));
               }
           }
       }
  }
}



template <class T, class U>
void productBlocking(SymmetricTileArray<U> const& A, TileArray<T> const& B, TileArray<T>& C)
{
   assert(A.nRowTiles() == C.nRowTiles());
   assert(A.nColTiles() == B.nRowTiles());
   assert(B.nColTiles() == C.nColTiles());

   size_t const nRowTiles(A.nRowTiles());
   size_t const nColTiles(B.nColTiles());
   unsigned i, j, k;

   std::vector<TileIndex> indices(A.sort());
   std::vector<TileIndex>::iterator iter;

   std::vector<TileIndex>::iterator const start(indices.begin());
   std::vector<TileIndex>::iterator const stop(indices.end());



   for (j = 0; j < nColTiles; ++j) {

/*
       for (i = 0; i < nRowTiles; ++i) {
           for (unsigned k = 0; k < A.nColTiles(); ++k) {
               // !!! Accumulate into C !!!
               if (i > k) {
                  tile_product(A(k,i), B(k,j), T(1.0), C(i,j), CblasTrans);
               }else {
                  tile_product(A(i,k), B(k,j), T(1.0), C(i,j), CblasNoTrans);
               }
           }
       }
*/
#pragma omp parallel for
       for (iter = start; iter != stop; ++iter) {
           i = iter->first;
           k = iter->second;
           if (i > k) {
              tile_product(A(k,i), B(k,j), T(1.0), C(i,j), CblasTrans);
           }else {
              tile_product(A(i,k), B(k,j), T(1.0), C(i,j), CblasNoTrans);
           }
       }
       

   }
}




template <class T, class U>
void product(SymmetricTileArray<U> const& A, TileArray<T> const& B, TileArray<T>& C)
{
   assert(A.nRowTiles() == C.nRowTiles());
   assert(A.nColTiles() == B.nRowTiles());
   assert(B.nColTiles() == C.nColTiles());

   size_t const nRowTiles(A.nRowTiles());
   size_t const nColTiles(B.nColTiles());

   std::vector<TileIndex> indices(A.sort());
   std::vector<TileIndex>::iterator iter;
   std::vector<TileIndex>::iterator const start(indices.begin());
   std::vector<TileIndex>::iterator const stop(indices.end());

   for (size_t j = 0; j < nColTiles; ++j) {

       size_t const nRows(C.nRows());
       size_t const nCols(C(0,j).nCols());
       size_t const nData(nRows*nCols);

       TileArray<T>** tileArrays;
       T* ptr;

       #pragma omp parallel
       {
          int const nThreads(omp_get_num_threads());

          #pragma omp single
          {
              // In general this reallocation is inefficient, but 
              // the loop over j is of length 1 for our case.
              tileArrays = new TileArray<T>*[nThreads];
              ptr = new T[nData*nThreads];
          }

          // Initialization
          #pragma omp for
          for (int n = 0; n < nThreads; ++n) {
              int    iThread(omp_get_thread_num());
              size_t offset(iThread*nData);

              TileArray<T>* pTA = new TileArray<T>(nRowTiles,1);
              tileArrays[iThread] = pTA;

              for (size_t i = 0; i < nRowTiles; ++i) {
                  CMTile<T>* pT = new CMTile<T>(C(i,j).nRows(), C(i,j).nCols());
                  pTA->set(i,0, pT);
              }

              memset(ptr+offset,0,nData*sizeof(T));
              pTA->bind(ptr+offset);
          }
 
          // Loop over tile products
          #pragma omp for schedule(dynamic,1) 
          for (iter = start; iter != stop; ++iter) {
              int iThread(omp_get_thread_num());
              size_t i(iter->first);
              size_t k(iter->second);

              if (i > k) {
                 tile_product(A(k,i), B(k,j), T(1.0), tileArrays[iThread]->tile(i,0), CblasTrans);
              }else {
                 tile_product(A(i,k), B(k,j), T(1.0), tileArrays[iThread]->tile(i,0), CblasNoTrans);
              }
          }

          // Reduce
/*
          #pragma omp single
          for (int iThread = 0; iThread < nThreads; ++iThread) {
              for (size_t i = 0; i < nRowTiles; ++i) {
                  //std::cout << "Reduction for thread/i = " << iThread << " / " << i << std::endl;
                  CMTile<T> const& t = dynamic_cast<CMTile<T> const&>(tileArrays[iThread]->tile(i,0));
                  C(i,j) += t;
              }
          }
*/
          T* cData(C(0,j).data());
          #pragma omp for 
          for (size_t i = 0; i < nData; ++i) {
              for (int n = 0; n < nThreads; ++n) {
                  cData[i] += ptr[i+n*nData];
              }
          }

          // Cleanup
          #pragma omp for
          for (int iThread = 0; iThread < nThreads; ++iThread) {
              delete tileArrays[omp_get_thread_num()];
          }
       }


       delete[] tileArrays;
       delete[] ptr;
   }
}


#endif
