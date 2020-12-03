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


template <class T>
void tile_product(CMTile<T> const& A, CMTile<T> const& B, T const c, CMTile<T>& C);


template <class T>
void tile_product(StripedTile<T> const& A, CMTile<T> const& B, T const c, CMTile<T>& C);


template <class T>
void tile_product(Tile<T> const& A, Tile<T> const& B, T const alpha, Tile<T>& C)
{
   CMTile<T> const& b = dynamic_cast<CMTile<T> const&>(B);
   CMTile<T>& c = dynamic_cast<CMTile<T>&>(C);

   switch (A.storage()) {
      case Zero: {
         return;
      } break;

      case Striped: {
         StripedTile<T> const& a = dynamic_cast<StripedTile<T> const&>(A);
         tile_product(a, b, alpha, c);
      } break;

      case CMDense: {
         CMTile<T> const& a = dynamic_cast<CMTile<T> const&>(A);
         tile_product(a, b, alpha, c);
      } break;

      default:
        Log::error("Unimplemented tile_product");
        break;
   }
}



template <class T>
void product(TileArray<T> const& A, TileArray<T> const& B, TileArray<T>& C)
{
   assert(A.nRowTiles() == C.nRowTiles());
   assert(A.nColTiles() == B.nRowTiles());
   assert(B.nColTiles() == C.nColTiles());
//#pragma omp parallel for
   for (unsigned bi = 0; bi < A.nRowTiles(); ++bi) {
       for (unsigned bj = 0; bj < B.nColTiles(); ++bj) {
           for (unsigned k = 0; k < A.nColTiles(); ++k) {
//             std::cout << "Multiplying block: C(" << bi << "," << bj << ") <- A(" 
//                       << bi << "," << k << ") x B(" << k << "," << bj << ")" << std::endl;
               tile_product(A(bi,k), B(k,bj), 1.0, C(bi,bj));
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
//#pragma omp parallel for
   for (unsigned bi = 0; bi < A.nRowTiles(); ++bi) {
       for (unsigned bj = 0; bj < B.nColTiles(); ++bj) {
           for (unsigned k = 0; k < A.nColTiles(); ++k) {
               if (bi != k ) {
                  tile_product(A(bi,k), B(k,bj), 1.0, C(bi,bj));
               }
           }
       }
  }
}

#endif
