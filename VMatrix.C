#include "VMatrix.h"
#include "Types.h"


VMatrix<complex>& VMatrix<complex>::fromDouble(VMatrix<double> const& that)
{
   m_nRows = that.m_nRows;
   m_nCols = that.m_nCols;

   m_storage = that.m_storage;
   m_layout  = that.m_layout;
   m_stripes = that.m_stripes;
   m_data    = 0;
   
   if (that.isBound()) {
      m_nData = that.m_nData;
      m_data = new complex[m_nData];
      for (unsigned i = 0; i < m_nData; ++i) {    
          m_data[i] = complex(that.m_data[i],0.0);
      }
   }
   return *this;
}


template <>
void VMatrix<double>::invert()
{
    if (m_nRows != m_nCols || m_storage != Dense || !isBound()) {
       std::cerr << "invert() called on nvalid matrix (" 
                 << m_nRows << "," << m_nCols << ") -> " << m_storage << std::endl;
       return;
    }

    int n(m_nRows);
    int *ipiv = new int[n+1];
    int lwork = n*n;
    int info;
    double* work = new double[lwork];

#ifdef __INTEL_COMPILER
    dgetrf(&n,&n,m_data,&n,ipiv,&info);
    dgetri(&n,m_data,&n,ipiv,work,&lwork,&info);
#else
    dgetrf_(&n,&n,m_data,&n,ipiv,&info);
    dgetri_(&n,m_data,&n,ipiv,work,&lwork,&info);
#endif

    delete ipiv;
    delete work;
}

template <>
void VMatrix<complex>::invert()
{
    if (m_nRows != m_nCols || m_storage != Dense || !isBound()) {
       std::cerr << "invert() called on invalid matrix (" 
                 << m_nRows << "," << m_nCols << ") -> " << m_storage << std::endl;
       return;
    }

    int n(m_nRows);
    int *ipiv = new int[n+1];
    int lwork = n*n;
    int info;
    complex* work = new complex[lwork];

#ifdef __INTEL_COMPILER
    zgetrf(&n,&n,m_data,&n,ipiv,&info);
    zgetri(&n,m_data,&n,ipiv,work,&lwork,&info);
#else
    zgetrf_(&n,&n,m_data,&n,ipiv,&info);
    zgetri_(&n,m_data,&n,ipiv,work,&lwork,&info);
#endif

    delete ipiv;
    delete work;
}


template <>
double VMatrix<double>::norm2() const
{
   double norm(0.0);
   for (unsigned i = 0; i < m_nData; ++i) {
       norm += m_data[i] * m_data[i];
   }

   return norm;
}

template <>
double VMatrix<complex>::norm2() const
{
   double r(0.0), n;
   for (unsigned i = 0; i < m_nData; ++i) {
       n  = std::norm(m_data[i]);
       r += n*n;
   }

   return r;
}


// The following should be able to be pushed to a single template in the header
template <>
void VMatrix<double>::print(const char* msg) const
{
   if (msg) {
      std::cout << msg << std::endl;
   }
   std::cout << std::fixed << std::showpoint << std::setprecision(2);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j) {
    //std::cout << (*this)(i,j) << "  ";
           std::cout << std::setw(5) << (*this)(i,j) << " ";
       }
       std::cout << std::endl;
   }
   std::cout << std::endl;
}


template <>
void VMatrix<complex>::print(const char* msg) const
{
   if (msg) {
      std::cout << msg << std::endl;
   }
   std::cout << std::fixed << std::showpoint << std::setprecision(2);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j) {
    //std::cout << (*this)(i,j) << "  ";
           std::cout << std::setw(5) << (*this)(i,j) << " ";
       }
       std::cout << std::endl;
   }
   std::cout << std::endl;
}


