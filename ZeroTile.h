#ifndef ZEROTILE_H
#define ZEROTILE_H
/******************************************************************************
 * 
 *  Tile that contains only zeros.
 * 
 *****************************************************************************/

#include "Tile.h"


template <class T>
class ZeroTile : public Tile<T>
{
   public:
      ZeroTile(size_t const nRows = 0, size_t const nCols = 0) : Tile<T>(nRows, nCols)
      { }


      StorageT storage() const
      {
         return Zero;
      }


      size_t numData() const
      {
          return 1;  // allocate at least one to avoid malloc errors
      }


      size_t indexOf(unsigned const i, unsigned const j) const
      {
         return 0;
      }


      double norm2() const
      {
         return 0.0;
      }


   protected:
      void copy(Tile<T> const& that)
      {
         size_t nRows(that.nRows());
         size_t nCols(that.nCols());

         this->dealloc();
         this->m_nRows = nRows;
         this->m_nCols = nCols;

         if (that.isBound()) {
            this->alloc();
            this->m_data[0] = T();
         }
      }
};

#endif
