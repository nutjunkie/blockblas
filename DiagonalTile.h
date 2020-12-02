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


      DiagonalTile(DiagonalTile<T> const& that)
      {
         copy(that);
      }


      StorageT storage() const
      {
         return Diagonal;
      }


      size_t numData() const
      {
          return std::min(this->m_nRows, this->m_nCols);
      }


      size_t indexOf(size_t const i, size_t const j) const
      {
          return (i==j) ? i : this->m_nData; 
      }


      std::string id() const 
      {
          return "\\";
      }


      Tile<T>& operator+=(Tile<T> const& that)
      {   
         Log::error("operator+= NYI for DiagonalTile");
         return *this;
      }   


      Tile<T>& operator-=(Tile<T> const& that)
      {   
         Log::error("operator-= NYI for DiagonalTile");
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
         resize(that.nRows(), that.nCols());

         if (that.isBound()) {
            this->alloc();
            for (size_t i = 0; i < this->m_nData; ++i) {
                this->m_data[i] = that(i,i);
            } 
         }
      }
      
     


      //void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }
      //void fill() { Tile<T>::fill(); }
};

#endif
