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
         std::cerr << "Unimplemented Storage in TileFactory" << std::endl;
      } break;
   }

   return tile;
};


template <class T, class U>
Tile<T>* TileFactory2(Tile<U> const& that)
{
   Tile<T>* tile(0);

   switch (that.storage()) {
      case Zero: {
         ZeroTile<U> const& t = dynamic_cast<ZeroTile<U> const&>(that);
         tile = new ZeroTile<T>();
         tile->from(t);
      } break;

      case Diagonal: {
         DiagonalTile<U>const& t = dynamic_cast<DiagonalTile<U> const&>(that);
         tile = new DiagonalTile<T>();
         tile->from(t);
      } break;

      case Striped: {
         StripedTile<U>const& t = dynamic_cast<StripedTile<U> const&>(that);
         tile = new StripedTile<T>();
         tile->from(t);
      } break;

      case CMDense: {
         CMTile<U>const& t = dynamic_cast<CMTile<U> const&>(that);
         tile = new CMTile<T>();
         tile->from(t);
      } break;

      default: {
         std::cerr << "Unimplemented Storage in TileFactory2" << std::endl;
      } break;
   }

   return tile;
};




#endif
