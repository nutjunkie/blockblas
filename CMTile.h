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

      CMTile(CMTile<T> const& that)
      {
         copy(that);
      }


      StorageT storage() const
      {
         return Dense;
      }


      size_t numData() const
      {
         return m_leadingDim * this->m_nCols;
      }


      size_t leadingDim() const
      {
         return m_leadingDim;
      }


      size_t indexOf(unsigned const i, unsigned const j) const
      {
         return (i + j*m_leadingDim);
      }


      void invert();


      double norm2() const;


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

         return *this;
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

         return *this;
      }


      CMTile& operator-=(DiagonalTile<T> const& that)
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
             *a -= b[i];
              a += m_leadingDim + 1;
         }

         return *this;
      }


      CMTile& operator-=(CMTile<T> const& that)
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
                 a[i] -= b[i];
             }
         }

         return *this;
      }


      CMTile& scale(T const alpha)
      {
         T* a0(this->data());
         size_t lda(this->m_leadingDim);

#pragma omp parallel for
         for (unsigned j = 0; j < this->m_nCols; ++j) {
             T* a(a0+j*lda);            
             for (unsigned i = 0; i < this->m_nRows; ++i) {
                 a[i] *= alpha;
             }
         }

         return *this;
      }


      // Not sure why these pass-throughs are necessary
      void fill(Functor<T> const& functor) { Tile<T>::fill(functor); }
      void fill() { Tile<T>::fill(); }


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
         std::cout << "Calling CMTile::copy" << std::endl;
         size_t nRows(that.nRows());
         size_t nCols(that.nCols());

         this->dealloc();
         setRC(nRows,nCols);

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

template <class T>
void tile_product(Tile<T> const& A, Tile<T> const& B, T const alpha, Tile<T>& C)
{
   CMTile<T> const& b = dynamic_cast<CMTile<T> const&>(B);
   CMTile<T>& c = dynamic_cast<CMTile<T>&>(C);

   switch (A.storage()) {
      case Zero: {
         return;
      } break;

      case Striped: {
         StripedTile<T> const& a = dynamic_cast<StripedTile<T> const&>(A);
         tile_product(a, b, alpha, c);
      } break;

      case Dense: {
         CMTile<T> const& a = dynamic_cast<CMTile<T> const&>(A);
         tile_product(a, b, alpha, c);
      } break;

      default:
        std::cerr << "ERROR: unimplemented tile_product" << std::endl;
        break;
   }
}


#endif
