#ifndef TILE_H
#define TILE_H
/******************************************************************************
 * 
 *  Abstract base class for an individual tile of a BlockMatrix.
 * 
 *****************************************************************************/

#include <assert.h>

#include <cstddef>
#include <cstring>   // memset
#include <iostream>
#include <iomanip>
#include "Functor.h"



template <class T>
class Tile
{
   public:
      template <class T>
      friend Tile<T>* TileFactory(StorageT const storage);
       
      template <class T>
      friend  Tile<T>* TileFactory(Tile<T> const& that);

   public:
      // On construction, the Tile does not have any data associated with it.
      // Actual data can be either bound externally (using bind(T*)) or allocated
      // using alloc().
      Tile(size_t const nRows = 0, size_t const nCols = 0) : m_data(0), m_nData(0), m_ownData(false)
      { 
         resize(nRows,nCols);
      }


      Tile(Tile<T> const& that) : m_data(0)
      {
         copy(that);
      }


      ~Tile() 
      { 
         dealloc(); 
      }


      void resize(size_t nRows, size_t nCols) 
      {
         dealloc();
         m_nRows = nRows;
         m_nCols = nCols;
      }


      Tile<T>& operator=(Tile<T> const& that)
      {
         if (this != &that) copy(that);
         return *this;
      }


      virtual StorageT storage() const = 0;


      size_t nRows() const 
      { 
         return m_nRows; 
      }


      size_t nCols() const 
      { 
         return m_nCols; 
      }


      // Allow for raw access to the memory
      T const* data() const 
      { 
         return m_data; 
      }


      T* data() 
      {
         return m_data; 
      }


      // Allocates sufficient space for the Tile based on its storage type.
      virtual void alloc(size_t nData = 0) 
      {
         if (m_data) dealloc();
         m_nData   = (nData == 0) ? numData() : nData;
         m_data    = new T[m_nData];
         m_ownData = true;
      }


      void dealloc() 
      {
         if (m_ownData && m_data) delete [] m_data;
         m_data    = 0;
         m_nData   = 0;
         m_ownData = false;
      }


      bool isBound() const 
      { 
         return m_data != 0; 
      }


      // This 'borrows' the memory (shared pointer?)
      Tile& bind(T* data)
      {
         if (isBound()) dealloc();
         m_ownData = false;
         m_nData = numData();
         m_data  = data;
         return *this;
      }


      void fill()
      {
         if (!isBound()) alloc();
         memset(m_data, 0, m_nData*sizeof(T));
      }


      // The following functions are inefficient and are for convenience only.
      void fill(Functor<T> const& functor)
      {
         if (!isBound()) alloc();

         for (unsigned j = 0; j < m_nCols; ++j) {
             for (unsigned i = 0; i < m_nRows; ++i) {
                 set(i,j,functor(i,j));
             }
         }
      }


      T operator()(unsigned const i, unsigned const j) const
      {
          size_t idx(indexOf(i,j));
          return idx < m_nData ? m_data[idx] : T();
      }


      void set(unsigned const i, unsigned const j, T const value)
      {
          size_t idx(indexOf(i,j));
          if (idx < m_nData) m_data[idx] = value;
      }


      // This is general, but very inefficient.  The copy() functions
      // should be better when Tile classes and data types match.
      template <class U>
      void from(Tile<U> const& that)
      {
         dealloc();
         setRC(that.nRows(),that.nCols());

         if (that.isBound()) {
            alloc();

            for (unsigned j = 0; j < m_nCols; ++j) {
                for (unsigned i = 0; i < m_nRows; ++i) {
                    set(i,j,that(i,j));
                }
            }
         }
      }

       
      virtual Tile<T>& operator+=(Tile<T> const& that) 
      {
         std::cerr << "WARNING: unoptimized operator += called" << std::endl;
         assert(that.m_nCols == m_nCols);
         assert(that.m_nRows == m_nRows);

         T val;
         for (unsigned j = 0; j < m_nCols; ++j) {
             for (unsigned i = 0; i < m_nRows; ++i) {
                 val = (*this)(i,j);
                 set(i,j,val + that(i,j));
             }
         }

         return *this;
      }


      virtual Tile<T>& operator-=(Tile<T> const& that) 
      {
         std::cerr << "WARNING: unoptimized operator += called" << std::endl;
         assert(that.m_nCols == m_nCols);
         assert(that.m_nRows == m_nRows);

         T val;
         for (unsigned j = 0; j < m_nCols; ++j) {
             for (unsigned i = 0; i < m_nRows; ++i) {
                 val = (*this)(i,j);
                 set(i,j,val - that(i,j));
             }
         }

         return *this;
      }


      virtual Tile<T>& scale(T const alpha)
      {
         for (unsigned i = 0; i < m_nData; ++i) {
             m_data[i] = alpha*m_data[i];
         }

         return *this;
      }

 
      // Adds the value t to the diagonals
      Tile<T>& addToDiag(T const t)
      {
         unsigned k(std::min(m_nCols, m_nRows));
         T val;
         for (unsigned i = 0; i < k; ++i) {
             val = (*this)(i,i);
             set(i,i,val+t);
         }
         return *this;
      }


      virtual void info(const char* msg = 0, std::ostream& os = std::cout) const
      {
         if (msg) std::cout << msg << std::endl;
         std::cout << "Storage:     " << toString(storage()) << std::endl;
         std::cout << "Num data:    " << m_nData << std::endl;
         std::cout << "Type size:   " << sizeof(T) << std::endl;
         std::cout << "Dimensions:  " << m_nRows << "x" << m_nCols << std::endl;
      }


      void print(const char* msg = 0, std::ostream& os = std::cout) const
      {
         if (isBound()) {
            if (msg)  os << msg << std::endl;

            std::cout << std::fixed << std::showpoint << std::setprecision(2);
            for (unsigned i = 0; i < m_nRows; ++i) {
                for (unsigned j = 0; j < m_nCols; ++j) {
                    os << std::setw(5) << (*this)(i,j) << " ";
                }   
                os << std::endl;
            }   
            os << std::endl;
         }else {
            os << "WARNING: print called on unbound Tile" << std::endl;
         }
      }


      virtual size_t numData() const = 0;

      virtual size_t indexOf(unsigned const i, unsigned const j) const = 0;


   protected:
      size_t m_nRows;       
      size_t m_nCols;       
      size_t m_nData;       
      bool   m_ownData;
      T*     m_data;


      virtual void setRC(size_t nRows, size_t nCols)
      {
         m_nRows = nRows;
         m_nCols = nCols;
      }


      virtual void copy(Tile<T> const& that) 
      {
         dealloc();
         setRC(that.nRows(),that.nCols());

         if (that.isBound()) {
            alloc(that.m_nData);
            memcpy(m_data, that.m_data, m_nData*sizeof(T));
         }
      }
};

#endif
