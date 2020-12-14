#ifndef FUNCTOR_H
#define FUNCTOR_H
/******************************************************************************
 * 
 *  Functor definitions for conveniently populating Tiles
 *
 *****************************************************************************/

#include <cmath>
#include "Types.h"



// Function prototype used for evaluating elements of a given sub matrix 

template <class T>
class Functor
{
   public:
      virtual T operator()(unsigned const i, unsigned const j) const = 0;
};


template <class T>
class ZeroFunctor : public Functor<T>
{
   public:
      T operator()(unsigned const i, unsigned const j) const { return T(); } 
};


template <class T>
class DiagonalFunctor : public Functor<T>
{
   public:
      DiagonalFunctor(T const val = T(1.0)) : m_val(val) { }
      T operator()(unsigned const i, unsigned const j) const 
      { 
         //return i == j ? 1.0*(i+1) : 0.0; 
         return i == j ? m_val : T(0.0); 
      } 

   private:
      T m_val;
};


template <class T>
class DebugFunctor : public Functor<T>
{
   public:
      T operator()(unsigned const i, unsigned const j) const 
      { 
         return T((i+1) + 0.01*(j+1)); 
      }
};



template <class T>
class StencilFunctor : public Functor<T>
{
   public:
      StencilFunctor(T const scale = 1.0) : m_scale(scale) 
      { }

      T operator()(unsigned const i, unsigned const j) const 
      { 
         int d(std::abs((int)i-(int)j));
         T val(0.0);
         switch (d) {
            case 0: val = 6.0;         break;
            case 1: val = 3.0*m_scale; break;
            case 2: val = 1.0*m_scale; break;
            case 3: val = 1.0*m_scale; break;
         }
         return val;
      }

   private:
      T m_scale;
};


class TestFunctor : public Functor<double>
{
   public:
      TestFunctor(unsigned const rowOffset = 0, unsigned const colOffset = 0) :
         m_rowOffset(rowOffset), m_colOffset(colOffset) { }

      double operator()(unsigned const i, unsigned const j) const 
      { 
         unsigned offi(i+m_rowOffset);
         unsigned offj(j+m_colOffset);
         return (offi == offj) ? 1.0 + offi : 0.2*std::sin(1.0*(offi+offj));
      }
   private:
      unsigned m_rowOffset;
      unsigned m_colOffset;
};

class ComplexTestFunctor : public Functor<complex>
{
   public:
      ComplexTestFunctor(unsigned const rowOffset = 0, unsigned const colOffset = 0) :
         m_rowOffset(rowOffset), m_colOffset(colOffset) { }

      complex operator()(unsigned const i, unsigned const j) const 
      { 
         unsigned offi(i+m_rowOffset);
         unsigned offj(j+m_colOffset);
         return (offi == offj) ? complex(1.0 + offi,0.0) : 
             complex(0.1*std::cos(1.0*(offi+offj)), 0.1*std::sin(1.0*(offi+offj)));
      }
   private:
      unsigned m_rowOffset;
      unsigned m_colOffset;
};


#endif
