#ifndef FUNCTOR_H
#define FUNCTOR_H
/******************************************************************************
 * 
 *  Class declarations for managing block matrices.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero matrices.
 * 
 *  The data in the matrix is not stored explicitly, evaluate() must be 
 *  called to load the give array with  
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


class DebugFunctor : public Functor<double>
{
   public:
      double operator()(unsigned const i, unsigned const j) const 
      { 
         return (i+1) + 0.01*(j+1); 
      }
};


class ComplexDebugFunctor : public Functor<complex>
{
   public:
      complex operator()(unsigned const i, unsigned const j) const 
      { 
         return complex(1.0*(i+1)+ 0.01*(j+1),0.0);
      }
};


class StencilFunctor : public Functor<double>
{
   public:
      double operator()(unsigned const i, unsigned const j) const 
      { 
         int d(std::abs((int)i-(int)j));
         double val(0.0);
         switch (d) {
            case 0: val = 6.0; break;
            case 1: val = 3.0; break;
            case 2: val = 1.0; break;
            case 3: val = 1.0; break;
         }
         return val;
      }
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
         return (offi == offj) ? 1.0 + offi : 0.1*std::sin(1.0*(offi+offj));
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
