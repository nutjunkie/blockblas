#ifndef STRIPEDTILE_H
#define STRIPEDTILE_H
/******************************************************************************
 * 
 *  Tile that contains only elements along specified stripes
 * 
 *****************************************************************************/

#include <vector>
#include <algorithm>
#include "Tile.h"


template <class T>
class StripedTile : public Tile<T>
{
    public:
       StripedTile(size_t const nRows = 0, size_t const nCols = 0, 
         std::vector<int> const& stripes = std::vector<int>()) 
         : Tile<T>(nRows, nCols), m_stripes(stripes)
       { }


       StorageT storage() const
       {
          return Striped;
       }


       std::vector<int> const& stripes() const 
       {
          return m_stripes;
       }


       size_t numData() const
       {
          return std::min(this->m_nRows,this->m_nCols) * m_stripes.size();
       }


       size_t indexOf(unsigned const i, unsigned const j) const
       {
          size_t idx(this->m_nData);
          int stripe((int)j-(int)i);
          std::vector<int>::const_iterator it;
          it = std::find(m_stripes.begin(), m_stripes.end(), stripe);
          
          if (it != m_stripes.end()) {
             // We have hit a non-zero element
             unsigned m(std::min(this->m_nRows,this->m_nCols));
             unsigned index = std::distance(m_stripes.begin(), it);
             int ij = (stripe < 0) ? j : i;
             idx = ij + index*m;
          }

          return idx;
       }

       // Not sure why these pass-throughs are necessary
       void bind(T* data) { Tile<T>::bind(data); }
       void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }


   protected:
      void copy(StripedTile<T> const& that)
      {
          m_stripes = that.stripes;
          Tile<T>::copy(that);
      }


   private:
      std::vector<int> m_stripes;
};

#endif
