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
class CMTile;

template <class T>
class StripedTile : public Tile<T>
{
   public:
      StripedTile(size_t const nRows = 0, size_t const nCols = 0, 
        std::vector<int> const& stripes = std::vector<int>()) 
        : Tile<T>(nRows, nCols), m_stripes(stripes)
      { }


      StripedTile(size_t const nRows = 0, size_t const nCols = 0, size_t const nStripes = 0)
        : Tile<T>(nRows, nCols)
      { 
         // yuk
         for (size_t i = 0; i < nStripes; ++i)  m_stripes.push_back(0);
      }


      StripedTile(StripedTile<T> const& that)
      {
         copy(that);
      }


      StripedTile(CMTile<T> const& that, std::vector<int> const& stripes)
      {
         m_stripes = stripes;
         size_t nr(that.nRows());
         size_t nc(that.nCols());
         this->resize(nr, nc);

         if (that.isBound()) {
            this->alloc();
            size_t m(std::min(nr,nc));

            for (unsigned s = 0; s < stripes.size(); ++s) {
                int offset(stripes[s]);
                int len = (offset < 0) ? std::min(nr+offset, nc)
                                       : std::min(nr, nc-offset);
                int ioff = std::max(0,-offset);
                int joff = std::max(0, offset);
                
                T* d0 = &(this->m_data[s*m]);
                for (unsigned ij = 0; ij < len; ++ij) {
                    d0[ij] = that(ij+ioff,ij+joff);
                }
            }
         }
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


      void setStripes(std::vector<int> const& stripes) 
      {
         m_stripes = stripes;
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
         this->resize(that.nRows(), that.nCols());

         if (that.isBound()) {
            this->alloc();
            memcpy(this->m_data, that.m_data, this->m_nData*sizeof(T));
         }
      }


      void copy(CMTile<T> const& that)
      {
         m_stripes = that.m_stripes;
         this->resize(that.nRows(), that.nCols());

         if (that.isBound()) {
            this->alloc();
            memcpy(this->m_data, that.m_data, this->m_nData*sizeof(T));
         }
      }



   private:
      std::vector<int> m_stripes;
};

#endif
