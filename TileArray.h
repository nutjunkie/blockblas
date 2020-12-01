#ifndef TILEARRAY_H
#define TILEARRAY_H
/******************************************************************************
 * 
 *  Class declarations for managing arrays of Tiles.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero Tiles.
 * 
 *****************************************************************************/

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include "Tile.h"
#include "TileFactory.h"
#include "Functor.h"


template <class T>
class TileArray
{
   public:
      TileArray(size_t const nRowTiles, size_t const nColTiles) :
         m_nRowTiles(nRowTiles), m_nColTiles(nColTiles)
      {
         size_t nTiles(m_nRowTiles*m_nColTiles);
         m_tiles = new Tile<T>*[nTiles];

         for (unsigned i = 0; i < nTiles; ++i) {
             m_tiles[i] = new ZeroTile<T>();
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


      template <class U>
      void from(TileArray<U> const& that)
      {
         std::cerr << "TileArray::from NYI" << std::endl;
      } 

      
      void set(unsigned const row, unsigned const col, Tile<T>* tile) 
      {
         unsigned idx(row + col*m_nRowTiles);

         if (idx < m_nRowTiles*m_nColTiles) {
            delete m_tiles[idx];
            m_tiles[idx] = tile;
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
         return *this;
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


      void bind(Functor<T>& functor) 
      {
         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                 tile(bi,bj).bind(functor);
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


      unsigned nRowTiles() const 
      { 
         return m_nRowTiles; 
      }


      unsigned nColTiles() const 
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
                 Tile<T> const& m(tile(bi,bj));
                 for (unsigned i = 0; i < m.nRows(); ++i) {
                     for (unsigned j = 0; j < m.nCols(); ++j) {
                         vm.set(iOff+i, jOff+j, m(i,j));
                     }
                 }
             }
         }
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


      void info(const char* msg = 0) const
      {
         if (msg) {
            std::cout << msg << std::endl;
         }
         std::cout << "Number of row blocks: :    " << nRowTiles() << std::endl;
         std::cout << "Number of column blocks:   " << nColTiles() << std::endl;
         std::cout << "Total number of rows:      " << nRows() << std::endl;
         std::cout << "Total number of columns:   " << nCols() << std::endl;
         std::cout << "TileArray is consistent: " << (consistent() ? "true" : "false") << std::endl;
         std::cout << std::endl;

         std::cout << "Tile sizes:" << std::endl;
         for (unsigned row = 0; row < m_nRowTiles; ++row) {
             for (unsigned col = 0; col < m_nColTiles; ++col) {
                 Tile<T> const& m(tile(row, col));
                 std::cout << "(" << m.nRows() <<","<< m.nCols() << ")  ";
             }   
             std::cout << std::endl;
         } 

         std::cout << std::endl;
         std::cout << "Data arrangement:" << std::endl;

         for (unsigned row = 0; row < m_nRowTiles; ++row) {
             for (unsigned col = 0; col < m_nColTiles; ++col) {
                 StorageT s(tile(row, col).storage());
                 switch (s) {
                    case Zero:     std::cout << '.';  break;
                    case Diagonal: std::cout << '\\'; break;
                    case Banded:   std::cout << 'b';  break;
                    case Dense:    std::cout << 'X';  break;

                    case Striped: {
                       StripedTile<T> const& st = dynamic_cast<StripedTile<T> const&>(tile(row, col));
                       unsigned n(st.stripes().size());
                       std::cout << n;
                    } break;
                 }
                 std::cout << "   ";
             }   
             std::cout << std::endl;
         }
      }

      void print(const char* msg = 0) const
      {
         if (msg) std::cout << msg << std::endl;

         for (unsigned bi = 0; bi < m_nRowTiles; ++bi) {
             unsigned nr = tile(bi,0).nRows();
             for (unsigned i = 0; i < nr; ++i) {
                 for (unsigned bj = 0; bj < m_nColTiles; ++bj) {
                     Tile<T> const& m(tile(bi, bj));
                     unsigned nc = m.nCols();
                     for (unsigned j = 0; j < nc; ++j) {
                         if (m(i,j) == 0) {
                            std::cout << "    . ";
                         }else {
                            std::cout << std::setw(5) << m(i,j) << " ";
                         }
                     }
                     std::cout << " ";
                }
                std::cout << std::endl;
             }
            std::cout << std::endl;
         }
      }


   protected:
      size_t m_nRowTiles;
      size_t m_nColTiles;

      Tile<T>** m_tiles;
       

   private:
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
         destroy();
         m_nRowTiles = that.m_nRowTiles;
         m_nColTiles = that.m_nColTiles;

         m_tiles = new Tile<T>*[m_nRowTiles*m_nColTiles];

        for (unsigned col = 0; col < m_nColTiles; ++col) {
            for (unsigned row = 0; row < m_nRowTiles; ++row) {
                 m_tiles[row + col*m_nRowTiles] = TileFactory(that(row,col));
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
