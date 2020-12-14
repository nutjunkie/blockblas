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

            default:
               std::cerr << "ERROR: operator+= NYI for CMTile" << std::endl;
               break;
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

            default:
               std::cerr << "ERROR: operator-= NYI for CMTile" << std::endl;
               break;
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


      Tile<T>& scale(T const alpha)
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


      void fill0()
      {
         if (!this->isBound()) this->alloc();

         T* a0(this->data());
         size_t lda(this->m_leadingDim);

//#pragma omp parallel for
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

/*
      template <class T>
      void toDense(CMTile<T>& vm) const
      {
         vm.resize(nRows(), nCols());
         vm.alloc();

         for (unsigned bi = 0; bi < nRowTiles(); ++bi) {
             unsigned iOff(rowOffset(bi));
             for (unsigned bj = 0; bj < nColTiles(); ++bj) {
                 unsigned jOff(colOffset(bj));
//               std::cout << "Tile("<< bi << "," << bj << ") with offset (" 
//                         << iOff << "," << jOff << ")" << std::endl;
                 U const& m(tile(bi,bj));
                 for (unsigned i = 0; i < m.nRows(); ++i) {
                     for (unsigned j = 0; j < m.nCols(); ++j) {
                         vm.set(iOff+i, jOff+j, m(i,j));
                     }
                 }
             }
         }
      }
*/




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

#endif
