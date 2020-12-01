#ifndef DIGAONALTILE_H
#define DIGAONALTILE_H
/******************************************************************************
 * 
 *  Tile that contains only elements along the diagonal
 * 
 *****************************************************************************/

#include "Tile.h"


template <class T>
class DiagonalTile : public Tile<T>
{
    public:
       DiagonalTile(size_t const nRows = 0, size_t const nCols = 0) : Tile<T>(nRows, nCols)
       { }


       StorageT storage() const
       {
          return Diagonal;
       }


       size_t numData() const
       {
           return std::min(this->m_nRows, this->m_nCols);
       }


       size_t indexOf(unsigned const i, unsigned const j) const
       {
           return (i==j) ? i : this->m_nData; 
       }


      void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }
      void fill() { Tile<T>::fill(); }
};

#endif
