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


class DiagonalFunctor : public Functor<double>
{
   public:
      double operator()(unsigned const i, unsigned const j) const 
      { 
         return i == j ? 1.0*(i+1) : 0.0; 
      } 
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
         return complex((i+1), 0.01*(j+1)); 
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

#endif
