#ifndef TILEARRAY_H
#define TILEARRAY_H
/******************************************************************************
 * 
 *  Class for managing arrays of Tiles (i.e. a block matrix).  The interface 
 *  is desiged to homogenize the handling of both dense, banded and zero Tiles.
 * 
 *****************************************************************************/

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include "TileFactory.h"
#include "Functor.h"


template <class T>
class TileArray
{
   public:

      TileArray(size_t const nRowTiles = 0, size_t const nColTiles = 0) :
         m_nRowTiles(nRowTiles), m_nColTiles(nColTiles), m_tiles(0)
      {
         size_t nTiles(m_nRowTiles*m_nColTiles);
         if (nTiles > 0) {
            m_tiles = new Tile<T>*[nTiles];
            memset(m_tiles, 0, nTiles*sizeof(Tile<T>*));
         }
      }


      TileArray(TileArray const& that) : m_tiles(0)
      {
         copy(that);
      }


      TileArray(Tile<T> const& that) : m_tiles(0)
      {
         copy(that);
      }


      TileArray& operator=(TileArray const& that)
      {
         if (this != &that) copy(that);
         return *this;
      }


      ~TileArray()
      {
         destroy();
      }


      void set(unsigned const row, unsigned const col, Tile<T>* tile) 
      {
         unsigned idx(row + col*m_nRowTiles);

         if (idx < m_nRowTiles*m_nColTiles) {
            if (m_tiles[idx]) delete m_tiles[idx];
            m_tiles[idx] = tile;
         }
      }


      template <class U>
      void from(TileArray<U> const& that)
      {
         resize(that.nRowTiles(), that.nColTiles());

         for (unsigned col = 0; col < m_nColTiles; ++col) {
             for (unsigned row = 0; row < m_nRowTiles; ++row) {
                 m_tiles[row + col*m_nRowTiles] = TileFactory2<T,U>(that(row,col));
            }
         }
      }


      TileArray& operator+=(TileArray const& that)
      {
         assert(that.nRowTiles() == m_nRowTiles);
         assert(that.nColTiles() == m_nColTiles);

#pragma omp parallel for
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 tile(bi,bj) += that(bi,bj);
             }
         }
         return *this;
      }


      TileArray& operator-=(TileArray const& that)
      {
         assert(that.nRowTiles() == m_nRowTiles);
         assert(that.nColTiles() == m_nColTiles);

#pragma omp parallel for
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 tile(bi,bj) -= that(bi,bj);
             }
         }
         return *this;
      }


      void addToDiag(T const t)
      {
         unsigned m(std::min(m_nColTiles, m_nRowTiles));
         for (unsigned bi = 0; bi < m; ++bi) {
             tile(bi,bi).addToDiag(t);
         }
      }


      double norm2() const
      {
         double norm(0.0);
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 norm += tile(bi,bj).norm2();
             }
         }
         return norm;
      }
      

      // This implies the TileArrays are always stored in column-major form
      Tile<T>& tile(unsigned row, unsigned const col)
      {
          return *m_tiles[row + col*m_nRowTiles];
      }


      Tile<T> const& tile(unsigned row, unsigned const col) const
      {
          return *m_tiles[row + col*m_nRowTiles];
      }


      Tile<T>& operator()(unsigned const row, unsigned const col = 0)
      {
          return tile(row,col);
      }


      Tile<T> const& operator()(unsigned const row, unsigned const col = 0) const
      {
          return tile(row,col);
      }


      template <class T>
      void fill(Functor<T>& functor) 
      {
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 tile(bi,bj).fill(functor);
             }
         }
      }

 
      void fill() 
      {
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 tile(bi,bj).fill();
             }
         }
      }

      // This will blow up if the Tiles are not all dense
      void bind(T* data)
      {
          unsigned lda(nRows());

          for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
              T* d0(data+colOffset(bj)*lda);

              for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
                  CMTile<T>& cmt = dynamic_cast<CMTile<T>&>(tile(bi,bj));
                  cmt.bind(d0,lda);
                  d0 += cmt.nRows();
             }
          }
      }


      template <class T>
      void scale(T const t) 
      {
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 tile(bi,bj).scale(t);
             }
         }
      }


      // This will blow up if the Tiles are not all dense
      void reduce() 
      {
         for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
              for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
                 CMTile<T>& cmt = dynamic_cast<CMTile<T>&>(tile(bi,bj));
                 set(bi,bj,cmt.reduce());
             }
         }
      }


      size_t nRowTiles() const 
      { 
         return m_nRowTiles; 
      }


      size_t nColTiles() const 
      { 
         return m_nColTiles; 
      }


      unsigned nRows() const
      {
         if (m_nRowTiles == 0 || m_nColTiles == 0) return 0;

         unsigned n(0);
         for (unsigned i = 0; i < m_nRowTiles; ++i) {
             n += tile(i,0).nRows();
         }
         return n;
      }


      unsigned nCols() const
      {
         if (m_nRowTiles == 0 || m_nColTiles == 0) return 0;

         unsigned n(0);
         for (unsigned i = 0; i < m_nColTiles; ++i) {
             n += tile(0,i).nCols();
         }
         return n;
      }


      bool consistent() const
      {
         if (m_nRowTiles == 0 || m_nColTiles == 0) return true;

         unsigned nTiles(m_nRowTiles*m_nColTiles);

         for (unsigned i = 0; i < nTiles; ++i) {
             if (m_tiles[i] == 0) return false;
         }

         // Check row heights
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             unsigned n = tile(bi,0).nRows();

             for (unsigned bj = 1; bj < m_nColTiles; ++bj) {
                 if (n != tile(bi,bj).nRows()) return false;
             }
         }

         // Check column widths
         for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
             unsigned n = tile(0,bj).nCols();

             for (unsigned bi = 1; bi < m_nRowTiles; ++bi) {
                 if (n != tile(bi,bj).nCols()) return false;
             }
         }

         return true;
      }


      // Inefficient implmentation targeted for debugging.
      unsigned rowOffset(unsigned bi) const
      {
         unsigned offset(0);
         for (unsigned i = 0; i < bi; ++i) {
             offset += tile(i,0).nRows();
         }
         return offset;
      }


      // Inefficient implmentation targeted for debugging.
      unsigned colOffset(unsigned bj) const
      {
         unsigned offset(0);
         for (unsigned j = 0; j < bj; ++j) {
             offset += tile(0,j).nCols();
         }
         return offset;
      }


      void info(const char* msg = 0, std::ostream& os = std::cout) const
      {
         if (msg)  os << msg << std::endl;

         os << "Number of row blocks: :    " << nRowTiles() << std::endl;
         os << "Number of column blocks:   " << nColTiles() << std::endl;
         os << "Total number of rows:      " << nRows() << std::endl;
         os << "Total number of columns:   " << nCols() << std::endl;
         os << "TileArray is consistent: " << (consistent() ? "true" : "false") << std::endl;
         os << std::endl;

         os << "Tile sizes:" << std::endl;
         for (unsigned row = 0; row < m_nRowTiles; ++row) {
             for (unsigned col = 0; col < m_nColTiles; ++col) {
                 Tile<T> const& m(tile(row, col));
                 //os << "(" << m.nRows() <<","<< m.nCols() << ")  -> " << m.numData();
                 os << "(" << m.nRows() <<","<< m.nCols() << ")  ";
             }   
             os << std::endl;
         } 

         os << std::endl;
         os << "Data arrangement:" << std::endl;

         for (unsigned row = 0; row < m_nRowTiles; ++row) {
             for (unsigned col = 0; col < m_nColTiles; ++col) {
                 Tile<T> const& m(tile(row, col));
                 os << m.id() << "   ";
             }   
             os << std::endl;
         }
      }

      void print(const char* msg = 0, std::ostream& os = std::cout) const
      {
         if (msg) os << msg << std::endl;

         std::cout << std::fixed << std::showpoint << std::setprecision(2);
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             unsigned nr = tile(bi,0).nRows();
             for (unsigned i = 0; i < nr; ++i) {
                 for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                     Tile<T> const& m(tile(bi, bj));
                     unsigned nc = m.nCols();
                     for (unsigned j = 0; j < nc; ++j) {
                            os << std::setw(5) << m(i,j) << " ";
/*
                         if (m(i,j) == 0) {
                            os << "    . ";
                         }else {
                            os << std::setw(5) << m(i,j) << " ";
                         }
*/
                     }
                     os << " ";
                }
                os << std::endl;
             }
            os << std::endl;
         }
      }


      void resize(size_t nRowTiles, size_t nColTiles)
      {
         destroy();
         m_nRowTiles = nRowTiles;
         m_nColTiles = nColTiles;
         m_tiles = new Tile<T>*[m_nRowTiles*m_nColTiles];
         memset(m_tiles, 0, m_nRowTiles*m_nColTiles*sizeof(Tile<T>*));
         
      }


   private:
      size_t m_nRowTiles;
      size_t m_nColTiles;
      Tile<T>** m_tiles;


      void destroy()
      {
         if (m_tiles) {
            size_t nTiles(m_nRowTiles*m_nColTiles);
            for (size_t i = 0; i < nTiles; ++i) {
                delete m_tiles[i];
            }
            delete [] m_tiles;
         }

         m_tiles = 0;
      }


      void copy(TileArray<T> const& that) 
      {
         size_t const nRowTiles(that.m_nRowTiles);
         size_t const nColTiles(that.m_nColTiles);

         if (nRowTiles != m_nRowTiles ||
             nColTiles != m_nColTiles) resize(nRowTiles,nColTiles);

         for (unsigned col = 0; col < m_nColTiles; ++col) {
             for (unsigned row = 0; row < m_nRowTiles; ++row) {
                 Tile<T> const& TB(that(row,col));
/*
                 Tile<T> const& TA(tile(row,col));
                 if (TA.storage() == TB.storage() &&
                     TA.nRows()   == TB.nRows()   &&
                     TA.nCols()   == TB.nCols()) {
                 }else {
                 }
*/
                m_tiles[row + col*m_nRowTiles] = TileFactory(TB);
             }
         }
      }


      void copy(Tile<T> const& that) 
      {
         destroy();
         m_nRowTiles = 1;
         m_nColTiles = 1;

         m_tiles = new Tile<T>*[m_nRowTiles*m_nColTiles];
         m_tiles[0] = TileFactory(that);
      }
};
#endif
