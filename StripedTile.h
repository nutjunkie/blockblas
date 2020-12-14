#ifndef STRIPEDTILE_H
#define STRIPEDTILE_H
/******************************************************************************
 * 
 *  Tile that contains only elements along specified diagonal stripes.  The
 *  stripes are determined by an integer  vector containing offsets from the
 *  main diagonal.  0 => main diagonal, -n => lower off-diagonal, +n => upper
 *  off-diagonal.
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


      StripedTile(StripedTile<T> const& that)
      {
         copy(that);
      }

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


      size_t indexOf(size_t const i, size_t const j) const
      {
         size_t idx(this->m_nData);
         int stripe((int)j-(int)i);
         std::vector<int>::const_iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);
         
         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            size_t m(std::min(this->m_nRows,this->m_nCols));
            size_t index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            idx = ij + index*m;
         }

         return idx;
      }


      std::string id() const 
      {
         return std::to_string(m_stripes.size());
      }


      Tile<T>& operator+=(Tile<T> const& that)
      {   
         std::cerr << "operator+= NYI for StripedTile" << std::endl;;
         return *this;
      }   


      Tile<T>& operator-=(Tile<T> const& that)
      {   
          std::cerr << "operator-= NYI for StripedTile" << std::endl;;
          return *this;
      }  
 
 
      Tile<T>& scale(T const t)
      {
         for (size_t i = 0; i < this->m_nData; ++i) {
             this->m_data[i] *= t;
         } 
      }


   protected:
      void copy(Tile<T> const& that)
      {
          switch (that.storage()) {
             case Striped: {
                StripedTile<T> const& t = dynamic_cast<StripedTile<T> const&>(that);
                copy(t);
             }
          }
      }


      void copy(StripedTile<T> const& that)
      {
         m_stripes = that.m_stripes;
         resize(that.nRows(), that.nCols());

         if (that.isBound()) {
            this->alloc();
            memcpy(this->m_data, that.m_data, this->m_nData*sizeof(T));
         }
      }


   private:
      std::vector<int> m_stripes;
};

#endif
