#ifndef TILEPRODUCT_H
#define TILEPRODUCT_H
/******************************************************************************
 * 
 *  Non-memeber function declarations for multiplying various tile types.
 * 
 *****************************************************************************/

#include "CMTile.h"
#include "TileArray.h"
#include "StripedTile.h"


//Computes A.B = C
template <class T, class U>
void tile_product(CMTile<U> const& A, CMTile<T> const& B, 
   T const c, CMTile<T>& C, CBLAS_TRANSPOSE const Atrans = CblasNoTrans);


template <class T, class U>
void tile_product(DiagonalTile<U> const& A, CMTile<T> const& B, T const alpha, CMTile<T>& C)
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
void tile_product(StripedTile<U> const& A, CMTile<T> const& B, T const alpha, CMTile<T>& C)
{
   unsigned rowsA(A.nRows());
   unsigned colsA(A.nCols());
   unsigned colsC(C.nCols());

   unsigned ldb(B.leadingDim());
   unsigned ldc(C.leadingDim());

   unsigned m(std::min(rowsA,colsA));

   U const* a(A.data());
   T const* b(B.data());
   T*       c(C.data());

   std::vector<int> const& stripes(A.stripes());

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
void tile_product(Tile<U> const& A, Tile<T> const& B, T const gamma, Tile<T>& C)
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
         tile_product(a, b, gamma, c);
      } break;

      case Striped: {
         StripedTile<U> const& a = dynamic_cast<StripedTile<U> const&>(A);
         tile_product(a, b, gamma, c);
      } break;

      case CMDense: {
         CMTile<U> const& a = dynamic_cast<CMTile<U> const&>(A);
         tile_product(a, b, gamma, c);
      } break;

      default:
        std::cerr << "ERROR: Unimplemented tile_product" << std::endl;
        break;
   }
}


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
void product(TileArray<U> const& A, TileArray<T> const& B, TileArray<T>& C)
{
   assert(A.nRowTiles() == C.nRowTiles());
   assert(A.nColTiles() == B.nRowTiles());
   assert(B.nColTiles() == C.nColTiles());
#pragma omp parallel for 
   for (unsigned bi = 0; bi < A.nRowTiles(); ++bi) {
       for (unsigned bj = 0; bj < B.nColTiles(); ++bj) {
           for (unsigned k = 0; k < A.nColTiles(); ++k) {
               // !!! Accumulate into C !!!
               tile_product(A(bi,k), B(k,bj), T(1.0), C(bi,bj));
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


#endif
