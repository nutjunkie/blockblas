#include "Tile.h"
 

template <>
double Tile<double>::norm2() const
{
   double norm(0.0);
   for (unsigned i = 0; i < m_nData; ++i) {
       norm += m_data[i] * m_data[i];
   }

   return norm;
}


template <>
double Tile<complex>::norm2() const
{
   double r(0.0), n;
   for (unsigned i = 0; i < m_nData; ++i) {
       n  = std::norm(m_data[i]);
       r += n*n;
   }

   return r;
}
