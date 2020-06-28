/******************************************************************************
 * 
 *  Class declarations for managing block matrices.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero matrices.
 * 
 *  The data in the matrix is not stored explicitly, evaluate() must be 
 *  called to load the give array with  
 *
 *****************************************************************************/

// Function prototype used for evaluating elements of a given sub matrix 

class Functor
{
   public:
      virtual double operator()(unsigned const i, unsigned const j) const = 0;
};

class ZeroFunctor : public Functor
{
   public:
      double operator()(unsigned const i, unsigned const j) const { return 0.0; } 
};


class DiagonalFunctor : public Functor
{
   public:
      double operator()(unsigned const i, unsigned const j) const 
      { 
         return i == j ? 1.0*(i+1) : 0.0; 
      } 
};



class DebugFunctor : public Functor
{
   public:
      double operator()(unsigned const i, unsigned const j) const 
      { 
         return (i+1) + 0.01*(j+1); 
      }
};
