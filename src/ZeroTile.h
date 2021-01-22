#ifndef ZEROTILE_H
#define ZEROTILE_H
/******************************************************************************
 * 
 *  Tile that contains only zeros.
 * 
 *****************************************************************************/

#include "Tile.h"


template <class T>
class CMTile;

template <class T>
class ZeroTile : public Tile<T>
{
   public:
      ZeroTile(size_t const nRows = 0, size_t const nCols = 0) : Tile<T>(nRows, nCols)
      { }


      ZeroTile(ZeroTile<T> const& that)
      {
         copy(that);
      }


      ZeroTile(CMTile<T> const& that)
      {
         copy(that);
      }


      StorageT storage() const
      {
         return Zero;
      }


      size_t numData() const
      {
          return 1;  // allocate at least one to avoid malloc errors
      }


      size_t indexOf(size_t const i, size_t const j) const
      {
         return 0;
      }


      double norm2() const
      {
         return 0.0;
      }


      std::string id() const 
      {
          return ".";
      }


      Tile<T>& operator+=(Tile<T> const& that)
      {
         std::cerr << "operator+= NYI for ZeroTile" << std::endl;;
         return *this;
      }


      Tile<T>& operator-=(Tile<T> const& that)
      {
        std::cerr << "operator-= NYI for ZeroTile" << std::endl;;
         return *this;
      }


      Tile<T>& scale(T const t)
      {
         return *this;
      }


   protected:
      void copy(Tile<T> const& that)
      {
         resize(that.nRows(), that.nCols());

         if (that.isBound()) {
            this->alloc();
            this->m_data[0] = T();
         }
      }
};

#endif
