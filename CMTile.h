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

template <class U>
class TileArray;

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


      CMTile(TileArray<T> const& TA);


      StorageT storage() const
      {
         return CMDense;
      }


      size_t numData() const
      {
         return m_leadingDim * this->m_nCols;
      }


      size_t indexOf(size_t const i, size_t const j) const
      {
         return (i + j*m_leadingDim);
      }


      std::string id() const
      {
         return "X";
      }


      size_t leadingDim() const
      {
         return m_leadingDim;
      }


      void invert();


      void factorLU(int*);


      double norm2() const;


      void bind(T* data, size_t leadingDim = 0)
      {
         m_leadingDim = (leadingDim == 0) ? this->m_nRows : leadingDim;
         Tile<T>::bind(data);
      }


      Tile<T>& operator+=(Tile<T> const& that)
      {   
         switch (that.storage()) {
            case Diagonal: {
               DiagonalTile<T> const& t = dynamic_cast<DiagonalTile<T> const&>(that);
               *this += t;
            } break;
               
            case CMDense: {
               CMTile<T> const& t = dynamic_cast<CMTile<T> const&>(that);
               *this += t;
            } break;

            default: {
               size_t lda(this->m_leadingDim);
               size_t nr(this->m_nRows);
               size_t nc(this->m_nCols);
 
               for (unsigned j = 0; j < nc; ++j) {
                   for (unsigned i = 0; i < nr; ++i) {
                       this->m_data[i+lda*j] += that(i,j); 
                   }
               }
            } break;
         }
         return *this;
      }   


      Tile<T>& operator-=(Tile<T> const& that)
      {   
         switch (that.storage()) {
            case Diagonal: {
               DiagonalTile<T> const& t = dynamic_cast<DiagonalTile<T> const&>(that);
               *this -= t;
            } break;
               
            case CMDense: {
               CMTile<T> const& t = dynamic_cast<CMTile<T> const&>(that);
               *this -= t;
            } break;

            default: {
               size_t lda(this->m_leadingDim);
               size_t nr(this->m_nRows);
               size_t nc(this->m_nCols);
               for (unsigned j = 0; j < nc; ++j) {
                   for (unsigned i = 0; i < nr; ++i) {
                       this->m_data[i+lda*j] -= that(i,j); 
                   }
               }
            } break;
         }
         return *this;
      }  


      CMTile& operator+=(DiagonalTile<T> const& that)
      {
         size_t nRows(this->m_nRows);
         size_t nCols(this->m_nCols);
         assert(nCols == that.nCols());
         assert(nRows == that.nRows());
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
         assert(nCols == that.nCols());
         assert(nRows == that.nRows());

         T* a0(this->data());
         T const* b0(that.data());
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
         size_t m(std::min(nRows,nCols));

         assert(nCols == that.nCols());
         assert(nRows == that.nRows());

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

         assert(nCols == that.nCols());
         assert(nRows == that.nRows());

         T* a0(this->data());
         T const* b0(that.data());
         for (unsigned j = 0; j < nCols; ++j) {
             T* a(a0+j*lda);            
             T const* b(b0+j*ldb);            

             for (unsigned i = 0; i < nRows; ++i) {
                 a[i] -= b[i];
             }
         }

         return *this;
      }


      Tile<T>& scale(T const alpha)
      {
         T* a0(this->data());
         size_t lda(this->m_leadingDim);

         for (unsigned j = 0; j < this->m_nCols; ++j) {
             T* a(a0+j*lda);            
             for (unsigned i = 0; i < this->m_nRows; ++i) {
                 a[i] *= alpha;
             }
         }

         return *this;
      }


      void fill0()
      {
         if (!this->isBound()) this->alloc();

         T* a0(this->data());
         size_t lda(this->m_leadingDim);

         for (unsigned j = 0; j < this->m_nCols; ++j) {
             T* a(a0+j*lda);            
             memset(a,0,this->m_nRows*sizeof(T));
         }
      }


      void info(const char* msg = 0, std::ostream& os = std::cout) const
      {
         Tile<T>::info(msg,os);
         std::cout << "Leading Dim: " << m_leadingDim <<  std::endl;
      }


      Tile<T>* reduce()
      {
          unsigned const nr(this->m_nRows);
          unsigned const nc(this->m_nCols);

          int const nDiags(nr+nc-1);
          int const offset(nr-1);

          T* sums = new T[nDiags];

          memset(sums, 0, nDiags*sizeof(T));

          for (int j = 0; j < nc; ++j) {
              for (int i = 0; i < nr; ++i) {
                  sums[j-i+offset] += std::abs(this->m_data[i+j*m_leadingDim]);
              }
          }

          double thresh(1e-10);
          std::vector<int> stripes;
          for (int i = 0; i < nDiags; ++i) {
              if (sums[i] > thresh) stripes.push_back(i-offset);
          }

          Tile<T>* tile;
          size_t k(stripes.size());
          if (k == 0) {
             tile = new ZeroTile<T>(*this); 
          }else if (k == 1 && sums[offset] > thresh && nDiags != 1) {
             tile = new DiagonalTile<T>(*this); 
          }else if (k < nDiags/2) {
             tile = new StripedTile<T>(*this, stripes); 
          }else {
             tile = new CMTile<T>(*this);
          }

          //std::cout << "Storage reduced to " << k << " " << toString(tile->storage()) << std::endl;
          //this->print();
          delete [] sums;

          return tile;
      }



   protected:
      void resize(size_t nRows, size_t nCols)
      {
          this->dealloc();
          this->m_nRows = nRows;
          this->m_nCols = nCols;
          m_leadingDim  = nRows;
      }


      void copy(Tile<T> const& that)
      {
         switch (that.storage()) {
            case CMDense: {
               CMTile<T> const& t = dynamic_cast<CMTile<T> const&>(that);
               copy(t);
            }
               break;

            default:
               std::cerr <<"ERROR: CMTile::copy called with unimplemnted type: " 
                         << that.storage() << std::endl;
               break;
         }
      }


      void copy(CMTile<T> const& that)
      {
         size_t nRows(that.nRows());
         size_t nCols(that.nCols());

         this->dealloc();
         resize(nRows,nCols);

         if (that.isBound()) {
            this->alloc();
            T* a0(this->data());
            T const* b0(that.data());
            size_t lda(m_leadingDim);
            size_t ldb(that.leadingDim());

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

#endif
