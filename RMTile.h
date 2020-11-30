#ifndef RMTILE_H
#define RMTILE_H
/******************************************************************************
 * 
 *  Dense Tile with storage in row-major order
 * 
 *****************************************************************************/

#include "Tile.h"


template <class T>
class RMTile : public Tile<T>
{
    public:
       RMTile(size_t const nRows = 0, size_t const nCols = 0) : Tile<T>(nRows, nCols)
       { 
          m_leadingDim = nRows;
       }


       StorageT storage() const
       {
          return Dense;
       }


       size_t numData () const
       {
          return this->m_leadingDim * this->m_nRows;
       }


       size_t leadingDim() const
       {
          return m_leadingDim;
       }


       size_t indexOf(unsigned const i, unsigned const j) const
       {
          //std::cout << "(" << i << "," << j << ") -> " << i + j*this->m_leadingDim << std::endl;
          return (j + i*this->m_leadingDim);
       }


       void bind(T* data, size_t leadingDim = 0)
       {
          m_leadingDim = (leadingDim == 0) ? this->m_nCols : leadingDim;
          Tile<T>::bind(data);
       }


       // Not sure why these pass-throughs are necessary
       void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }
          

    private:
       size_t m_leadingDim;
};


void tile_product(RMTile<double> const& A, RMTile<double> const& B, double const c, RMTile<double>& C);

#endif
