#ifndef TILE_H
#define TILE_H
/******************************************************************************
 * 
 *  Abstract base class for an individual tile of a TileArray.
 * 
 *****************************************************************************/

#include <assert.h>

#include <cstddef>
#include <cstring>   // memset
#include <iostream>
#include <iomanip>
#include "Functor.h"
#include "Log.h"


template <class T>
class Tile
{
   public:
      template <class U>
      friend  Tile<U>* TileFactory(Tile<U> const& that);


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


      void operator=(Tile<T> const& that)
      {
         if (this != &that) copy(that);
      }


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
      // Does not reallocate if bound.
      Tile<T>& alloc() 
      {
         if (m_data) return *this;
         m_nData   = numData();
         m_data    = new T[m_nData];
         m_ownData = true;
         return *this;
      }


      // Returns if data is actually delelted.
      bool dealloc() 
      {
         bool dataDeleted(m_ownData && m_data);

         if (dataDeleted) delete [] m_data;

         m_data    = 0;
         m_nData   = 0;
         m_ownData = false;

         return dataDeleted;
      }


      bool isBound() const 
      { 
         return m_data != 0; 
      }


      // This 'borrows' the memory (shared pointer?)
      Tile& bind(T* data)
      {
         dealloc();
         m_ownData = false;
         m_nData   = numData();
         m_data    = data;
         return *this;
      }


      virtual void fill()
      {
         alloc();
         memset(m_data, 0, m_nData*sizeof(T));
      }


      // The following functions are inefficient and are for convenience only.
      void fill(Functor<T> const& functor)
      {
         alloc();

         for (size_t j = 0; j < m_nCols; ++j) {
             for (size_t i = 0; i < m_nRows; ++i) {
                 set(i,j,functor(i,j));
             }
         }
      }


      T operator()(size_t const i, size_t const j) const
      {
          size_t idx(indexOf(i,j));
          return idx < m_nData ? m_data[idx] : T();
      }


      T get(size_t const i, size_t const j) const
      {
          size_t idx(indexOf(i,j));
          return idx < m_nData ? m_data[idx] : T();
      }


      void set(size_t const i, size_t const j, T const value)
      {
          size_t idx(indexOf(i,j));
          if (idx < m_nData) m_data[idx] = value;
      }


      // This is general, but very inefficient.  The copy() functions
      // should be better when Tile classes and data types match.
      template <class U>
      void from(Tile<U> const& that)
      {
         resize(that.nRows(),that.nCols());

         if (that.isBound()) {
            alloc();

            for (size_t j = 0; j < m_nCols; ++j) {
                for (size_t i = 0; i < m_nRows; ++i) {
                    set(i,j,that(i,j));
                }
            }
         }
      }


      // Adds the value t to the diagonals
      void addToDiag(T const t)
      {
         unsigned k(std::min(m_nCols, m_nRows));
         T val;
         for (unsigned i = 0; i < k; ++i) {
             val = (*this)(i,i);
             set(i,i,val+t);
         }
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
            os << "Print called on unbound Tile";
         }
      }


      virtual void info(const char* msg = 0, std::ostream& os = std::cout) const
      {
         if (msg) std::cout << msg << std::endl;
         os << "Storage:     " << toString(storage()) << std::endl;
         os << "Num data:    " << m_nData << std::endl;
         os << "Type size:   " << sizeof(T) << std::endl;
         os << "Dimensions:  " << m_nRows << "x" << m_nCols << std::endl;
         os << "Own data:    " << m_ownData << std::endl;
      }


      virtual Tile<T>& operator+=(Tile<T> const& that) = 0;

      virtual Tile<T>& operator-=(Tile<T> const& that) = 0;

      virtual Tile<T>& scale(T const) = 0;

      virtual size_t numData() const = 0;

      virtual size_t indexOf(size_t const i, size_t const j) const = 0;

      virtual StorageT storage() const = 0;

      virtual std::string id() const = 0;

      virtual double norm2() const;

      virtual void invert() 
      {
          Log::error("ERROR: invert() called on invalid Tile");
      }


   protected:
      size_t m_nRows;       
      size_t m_nCols;       
      size_t m_nData;       
      bool   m_ownData;
      T*     m_data;


      virtual void copy(Tile<T> const& that) = 0;


      virtual void resize(size_t nRows, size_t nCols) 
      {
         if (m_nRows == nRows && m_nCols == nCols) return;
         m_nRows = nRows;
         m_nCols = nCols;
         if (dealloc()) alloc();
      }
};

#endif
