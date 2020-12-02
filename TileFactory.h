#ifndef TILEFACTORY_H
#define TILEFACTORY_H
/******************************************************************************
 * 
 *  Factory function for generating arbitrary Tile types
 *
 *****************************************************************************/

#include "ZeroTile.h"
#include "DiagonalTile.h"
#include "StripedTile.h"
#include "CMTile.h"


template <class T>
Tile<T>* TileFactory(Tile<T> const& that)
{
   Tile<T>* tile(0);

   switch (that.storage()) {
      case Zero: {
         ZeroTile<T> const& t = dynamic_cast<ZeroTile<T> const&>(that);
         tile = new ZeroTile<T>(t);
      } break;

      case Diagonal: {
         DiagonalTile<T>const& t = dynamic_cast<DiagonalTile<T> const&>(that);
         tile = new DiagonalTile<T>(t);
      } break;

      case Striped: {
         StripedTile<T>const& t = dynamic_cast<StripedTile<T> const&>(that);
         tile = new StripedTile<T>(t);
      } break;

      case CMDense: {
         CMTile<T>const& t = dynamic_cast<CMTile<T> const&>(that);
         tile = new CMTile<T>(t);
      } break;

      default: {
         Log::error("Unimplemented Storage in TileFactory");
      } break;
   }

   return tile;
};


#endif
