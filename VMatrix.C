#include "VMatrix.h"
#include "Types.h"
#include <mkl.h>

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
   double r(0.0);
   for (unsigned i = 0; i < m_nData; ++i) {
       r += std::norm(m_data[i]);
   }

   return r;
}
