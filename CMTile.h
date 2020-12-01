#ifndef CMTILE_H
#define CMTILE_H
/******************************************************************************
 * 
 *  Dense Tile with storage in column major order
 * 
 *****************************************************************************/

#include "Tile.h"
#include "StripedTile.h"
#include "DiagonalTile.h"


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
         return m_leadingDim * this->m_nCols;
      }


      size_t leadingDim() const
      {
         return m_leadingDim;
      }


      size_t indexOf(unsigned const i, unsigned const j) const
      {
         return (i + j*leadingDim());
      }


      void bind(T* data, size_t leadingDim = 0)
      {
         m_leadingDim = (leadingDim == 0) ? this->m_nRows : leadingDim;
         Tile<T>::bind(data);
      }


      CMTile& operator+=(DiagonalTile<T> const& that)
      {
         size_t nRows(this->m_nRows);
         size_t nCols(this->m_nCols);
#ifdef DEBUG          
         assert(nCols == that.nCols());
         assert(nRows == that.nRows());
#endif
         unsigned m(std::min(nRows,nCols));

         T* a(this->data());
         T const* b(that.data());

         for (unsigned i = 0; i < m; ++i) {
             *a += b[i];
              a += m_leadingDim + 1;
         }
      }


      CMTile& operator+=(CMTile<T> const& that)
      {
         size_t nRows(this->m_nRows);
         size_t nCols(this->m_nCols);
         size_t lda(this->m_leadingDim);
         size_t ldb(that.leadingDim());
#ifdef DEBUG          
         assert(nCols == that.nCols());
         assert(nRows == that.nRows());
#endif

         T* a0(this->data());
         T const* b0(that.data());
#pragma omp parallel for
         for (unsigned j = 0; j < nCols; ++j) {
             T* a(a0+j*lda);            
             T const* b(b0+j*ldb);            

             for (unsigned i = 0; i < nRows; ++i) {
                 a[i] += b[i];
             }
         }
      }


      // Not sure why these pass-throughs are necessary
      void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }


      void info(const char* msg = 0, std::ostream& os = std::cout) const
      {
         Tile<T>::info(msg,os);
         std::cout << "Leading Dim: " << m_leadingDim <<  std::endl;
      }


   protected:
      void setRC(size_t nRows, size_t nCols)
      {
          this->m_nRows = nRows;
          this->m_nCols = nCols;
          m_leadingDim  = nRows;
      }


      void copy(CMTile<T> const& that)
      {
         size_t nRows(that.m_nRows());
         size_t nCols(that.m_nCols());

         this->dealloc();
         this->m_nRows = nRows;
         this->m_nCols = nCols;

         if (that.isBound()) {
            this->alloc();
            T* a0(this->data());
            T const* b0(that.data());
            size_t lda(m_leadingDim);
            size_t ldb(that.leadingDim());

#pragma omp parallel for
            for (unsigned j = 0; j < nCols; ++j) {
                T*       a(a0+j*lda);            
                T const* b(b0+j*ldb);            
                memcpy(a,b,nRows*sizeof(T));
            }
         }
      }


    private:
       size_t m_leadingDim;
};


template <class T>
void tile_product(CMTile<T> const& A, CMTile<T> const& B, T const c, CMTile<T>& C);

template <class T>
void tile_product(StripedTile<T> const& A, CMTile<T> const& B, T const c, CMTile<T>& C);


#endif
