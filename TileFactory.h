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
Tile<T>* TileFactory(StorageT const storage)
{
    Tile<T>* tile;

    switch (storage) {
       case Zero:
          tile = new ZeroTile<T>();
          break;
       case Diagonal:
          tile = new DiagonalTile<T>();
          break;
       case Striped:
          tile = new StripedTile<T>();
          break;
       case Dense:
          tile = new CMTile<T>();
          break;
    }

    return tile;
};


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

      case Dense: {
         CMTile<T>const& t = dynamic_cast<CMTile<T> const&>(that);
         tile = new CMTile<T>(t);
      } break;
   }

   return tile;
};


#endif
