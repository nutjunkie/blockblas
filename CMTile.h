#ifndef CMTILE_H
#define CMTILE_H
/******************************************************************************
 * 
 *  Dense Tile with storage in column major order
 * 
 *****************************************************************************/

#include "Tile.h"
#include "StripedTile.h"


template <class T>
class CMTile : public Tile<T>
{
    public:
       CMTile(size_t const nRows = 0, size_t const nCols = 0) : Tile<T>(nRows, nCols)
       { 
          m_leadingDim = nRows;
       }


       StorageT storage() const
       {
          return Dense;
       }


       size_t numData () const
       {
          return this->m_leadingDim * this->m_nCols;
       }


       size_t leadingDim() const
       {
          return m_leadingDim;
       }


       void alloc()
       {
           if (m_leadingDim == 0) m_leadingDim = this->m_nRows;
           Tile<T>::alloc();
       }

       size_t indexOf(unsigned const i, unsigned const j) const
       {
          //std::cout << "(" << i << "," << j << ") -> " << i + j*this->m_leadingDimd << std::endl;
          return (i + j*this->m_leadingDim);
       }


       void bind(T* data, size_t leadingDim = 0)
       {
          m_leadingDim = (leadingDim == 0) ? this->m_nRows : leadingDim;
          Tile<T>::bind(data);
       }


       // Not sure why these pass-throughs are necessary
       void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }


       void info(const char* msg = 0, std::ostream& os = std::cout) const
       {
          Tile<T>::info(msg,os);
          std::cout << "Leading Dim: " << m_leadingDim <<  std::endl;
       }



    private:
       size_t m_leadingDim;
};


template <class T>
void tile_product(CMTile<T> const& A, CMTile<T> const& B, T const c, CMTile<T>& C);

template <class T>
void tile_product(StripedTile<T> const& A, CMTile<T> const& B, T const c, CMTile<T>& C);


#endif
