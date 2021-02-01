#ifndef SYMMETRICTILEARRAY_H
#define SYMMETRICTILEARRAY_H
/******************************************************************************
 * 
 *  Class for managing arrays of Tiles (i.e. a block matrix).  The interface 
 *  is desiged to homogenize the handling of both dense, banded and zero Tiles.
 *  Only the upper-triangle of tiles is stored and no check is performed to
 *  ensure the diagonal Tiles are themselves symmetric.
 * 
 *****************************************************************************/

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include "TileFactory.h"
#include "Functor.h"


template <class T>
class SymmetricTileArray 
{
   public:

      SymmetricTileArray(size_t const nTiles = 0) :
         m_nColTiles(nTiles), m_tiles(0)
      {
         resize(nTiles);
      }


      SymmetricTileArray(SymmetricTileArray const& that) : m_tiles(0)
      {
         copy(that);
      }


      SymmetricTileArray(TileArray<T> const& that) : m_tiles(0)
      {
         copy(that);
      }


      SymmetricTileArray& operator=(SymmetricTileArray const& that)
      {
         if (this != &that) copy(that);
         return *this;
      }


      ~SymmetricTileArray()
      {
         destroy();
      }


      void set(size_t const row, size_t const col, Tile<T>* tile) 
      {
         size_t idx(index(row,col));
         if (idx < m_nTotTiles) {
            if (m_tiles[idx]) delete m_tiles[idx];
            m_tiles[idx] = tile;
         }
      }


      template <class U>
      void from(SymmetricTileArray<U> const& that)
      {
         resize(that.nRowTiles());

         size_t k(0);
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row, ++k) {
                 m_tiles[k] = TileFactory2<T,U>(that(row,col));
             }
         }
      }


      SymmetricTileArray& operator+=(SymmetricTileArray const& that)
      {
         assert(that.nColTiles() == m_nColTiles);

#pragma omp parallel for
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 tile(row,col) += that(row,col);
             }
         }

         return *this;
      }


      SymmetricTileArray& operator-=(SymmetricTileArray const& that)
      {
         assert(that.nColTiles() == m_nColTiles);

#pragma omp parallel for
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 tile(row,col) -= that(row,col);
             }
         }

         return *this;
      }


      void addToDiag(T const t)
      {
         unsigned m(std::min(m_nColTiles, m_nColTiles));
         for (unsigned bi = 0; bi < m; ++bi) {
             tile(bi,bi).addToDiag(t);
         }
      }


      double norm2() const
      {
         double norm(0.0);
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 double n2(tile(bi,bj).norm2());
                 norm += (row == col) ? n2 : 2*n2;
             }
         }
         return norm;
      }
      

      // This implies the TileArrays are always stored in column-major form
      Tile<T>& tile(unsigned row, unsigned const col)
      {
          return *m_tiles[index(row,col)];
      }


      Tile<T> const& tile(unsigned row, unsigned const col) const
      {
          return *m_tiles[index(row,col)];
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
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 tile(row,col).fill(functor);
             }
         }
      }

 
      void fill() 
      {
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 tile(row,col).fill0();
             }
         }
      }

/*
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
*/

      template <class T>
      void scale(T const t) 
      {
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 tile(bi,bj).scale(t);
             }
         }
      }


      // This will blow up if the Tiles are not all dense
      void reduce() 
      {
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row) {
                 CMTile<T>& cmt = dynamic_cast<CMTile<T>&>(tile(row,col));
                 set(row,col,cmt.reduce());
             }
         }
      }


      size_t nRowTiles() const 
      { 
         return m_nColTiles; 
      }


      size_t nColTiles() const 
      { 
         return m_nColTiles; 
      }


      unsigned nRows() const
      {
         if (m_nColTiles == 0) return 0;

         unsigned n(0);
         for (unsigned i = 0; i < m_nColTiles; ++i) {
             n += tile(i,m_nColTiles-1).nRows();
         }
         return n;
      }


      unsigned nCols() const
      {
         if (m_nColTiles == 0) return 0;

         unsigned n(0);
         for (unsigned i = 0; i < m_nColTiles; ++i) {
             n += tile(0,i).nCols();
         }
         return n;
      }


      bool consistent() const
      {
         if (m_nColTiles == 0) return true;

         for (unsigned i = 0; i < m_nTotTiles; ++i) {
             if (m_tiles[i] == 0) return false;
         }

         // Check row heights
         for (unsigned bi = 0; bi < m_nColTiles; ++bi) {
             unsigned n = tile(bi,m_nColTiles-1).nRows();

             for (unsigned bj = bi; bj < m_nColTiles; ++bj) {
                 if (n != tile(bi,bj).nRows()) return false;
             }
         }

         // Check column widths
         for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
             unsigned n = tile(0,bj).nCols();

             for (unsigned bi = 1; bi <= bj; ++bi) {
                 if (n != tile(bi,bj).nCols()) return false;
             }
         }

         // Check diagonals are all square
         for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
             if (tile(bj,bj).nRows() != tile(bj,bj).nCols()) return false;
         }

         return true;
      }


      // Inefficient implmentation targeted for debugging.
      unsigned rowOffset(unsigned bi) const
      {
         unsigned offset(0);
         for (unsigned i = 0; i < bi; ++i) {
             offset += tile(i,m_nColTiles-1).nRows();
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

         os << "Number of row/col blocks:  " << nRowTiles() << std::endl;
         os << "Total number of rows/cols: " << nRows() << std::endl;
         os << "TileArray is consistent: " << (consistent() ? "true" : "false") << std::endl;
         os << std::endl;

         os << "Tile sizes:" << std::endl;
         for (unsigned row = 0; row < m_nColTiles; ++row) {
             for (unsigned col = 0; col < m_nColTiles; ++col) {
                 if (col < row) {
                    os << "       ";
                 }else {
                    Tile<T> const& m(tile(row, col));
                    //os << "(" << m.nRows() <<","<< m.nCols() << ")  -> " << m.numData();
                    os << "(" << m.nRows() <<","<< m.nCols() << ")  ";
                 }
             }   
             os << std::endl;
         } 

         os << std::endl;
         os << "Data arrangement:" << std::endl;

         for (unsigned row = 0; row < m_nColTiles; ++row) {
             for (unsigned col = 0; col < m_nColTiles; ++col) {
                 if (col < row) {
                    os << "    ";
                 }else {
                    Tile<T> const& m(tile(row, col));
                    os << m.id() << "   ";
                 }
             }   
             os << std::endl;
         }
      }


      void print(const char* msg = 0, std::ostream& os = std::cout) const
      {
         if (msg) os << msg << std::endl;

         std::cout << std::fixed << std::showpoint << std::setprecision(2);
         for (unsigned bi = 0; bi < m_nColTiles; ++bi) {
             //std::cout << " bi/bj " << bi << " " << m_nColTiles-1 << std::endl;
             unsigned nr = tile(bi,m_nColTiles-1).nRows();

             for (unsigned i = 0; i < nr; ++i) {
                 for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                     if (bi > bj) {
                        Tile<T> const& m(tile(0, bj));
                        unsigned nc = m.nCols();
                        for (unsigned j = 0; j < nc; ++j) {
                               os << "      ";
                        }
                        os << " ";
                     }else {
                        Tile<T> const& m(tile(bi, bj));
                        unsigned nc = m.nCols();
                        for (unsigned j = 0; j < nc; ++j) {
                               os << std::setw(5) << m(i,j) << " ";
                        }
                        os << " ";
                     }
                }
                os << std::endl;
             }
            os << std::endl;
         }
      }


      void resize(size_t nTiles)
      {
         destroy();
         m_nColTiles = nTiles;
         m_nTotTiles = ((m_nColTiles+1)*m_nColTiles) / 2;

         if (m_nTotTiles > 0) {
            m_tiles = new Tile<T>*[m_nTotTiles];
            memset(m_tiles, 0, m_nTotTiles*sizeof(Tile<T>*));
         }
      }


   private:
      size_t m_nColTiles;
      size_t m_nTotTiles;
      Tile<T>** m_tiles;


      size_t index(size_t const row, size_t const col) const
      {
         size_t idx = (row > col) ? m_nTotTiles : ((col+1)*col)/2 + row;
         //std::cout << "Tile mapping  (" << row << "," << col << ") -> " << idx << std::endl;
         return idx;
      }


      void destroy()
      {
         if (m_tiles) {
            for (size_t i = 0; i < m_nTotTiles; ++i) {
                delete m_tiles[i];
            }
            delete [] m_tiles;
         }

         m_tiles = 0;
      }


      void copy(TileArray<T> const& that) 
      {
         assert(that.nRowTiles() == that.nColTiles());

         resize(that.nRowTiles());

         size_t k(0);
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row, ++k) {
                 m_tiles[k] = TileFactory(that(row,col));
             }
         }
      }


      void copy(SymmetricTileArray<T> const& that) 
      {
         assert(that.nRowTiles() == that.nColTiles());

         resize(that.nRowTiles());

         size_t k(0);
         for (size_t col = 0; col < m_nColTiles; ++col) {
             for (size_t row = 0; row <= col; ++row, ++k) {
                 m_tiles[k] = TileFactory(that(row,col));
             }
         }
     }


     void copy(Tile<T> const& that) 
     {
         resize(1);
         m_tiles[0] = TileFactory(that);
     }
};
#endif
