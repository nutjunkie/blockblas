#include "CMTile.h"
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


