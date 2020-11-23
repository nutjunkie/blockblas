/******************************************************************************
 * 
 *  Member function definitions involving type template specialization.
 *
 *****************************************************************************/


template<>
LayoutT VMatrix<TYPE,RowMajor>::layout() const { return RowMajor; }

template<>
LayoutT VMatrix<TYPE,ColumnMajor>::layout() const { return ColumnMajor; }



template<>
TYPE VMatrix<TYPE,RowMajor>::operator()(unsigned const i, unsigned const j) const
{
   TYPE value = TYPE();

   switch (m_storage) {
      case Zero:
         break;

      case Diagonal:
         if (i==j) {
            value = m_data[i];
         }
         break;

      case Banded: {
         int kl(m_stripes[0]);
         int ku(m_stripes[1]);
//       std::cout << "Access element: (" << i << "," << j << ") -> ";
         if (std::max(0,(int)i-kl) <= j && j <= std::min((int)m_nCols,(int)i+ku)) {
            int k(j-i+kl+i*(kl+ku+1));
            value = m_data[k];
//          std::cout <<  k << " = ";
         }
//       std::cout <<  value << std::endl;
      }  break;

      case Striped: {
         int stripe(j-i);
         std::vector<int>::const_iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);

         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            unsigned m(std::min(m_nRows,m_nCols));
            unsigned index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            value = m_data[ij + index*m];
         }

      }  break;

      case Dense:
         value = m_data[i*m_nCols + j];
         break;
   }

   return value;
}
 
template<>
TYPE VMatrix<TYPE,ColumnMajor>::operator()(unsigned const i, unsigned const j) const
{
   TYPE value = TYPE();

   switch (m_storage) {
      case Zero:
         break;

      case Diagonal:
         if (i==j) {
            value = m_data[i];
         }
         break;

      case Banded: {
         int kl(m_stripes[0]);
         int ku(m_stripes[1]);
//       std::cout << "Access element: (" << i << "," << j << ") -> ";
         if (std::max(0,(int)j-ku) <= i && i <= std::min((int)m_nRows,(int)j+kl)) {
            int k(i-j+ku+j*(kl+ku+1));
            value = m_data[k];
//          std::cout <<  k << " = ";
         }
//       std::cout <<  value << std::endl;
      }  break;

      case Striped: {
         int stripe(j-i);
         std::vector<int>::const_iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);

         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            unsigned m(std::min(m_nRows,m_nCols));
            unsigned index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            value = m_data[ij + index*m];
         }
      }  break;

      case Dense:
         value = m_data[i + j*m_nRows];
         break;
   }

   return value;
}
 


template<>
void VMatrix<TYPE,RowMajor>::fillDense(Functor<TYPE> const& functor)
{
   unsigned k(0);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j, ++k) {
           m_data[k] = functor(i,j); 
       }
   }
}

template<>
void VMatrix<TYPE,ColumnMajor>::fillDense(Functor<TYPE> const& functor)
{
   unsigned k(0);
   for (unsigned j = 0; j < m_nCols; ++j) {
       for (unsigned i = 0; i < m_nRows; ++i, ++k) {
           m_data[k] = functor(i,j); 
       }
   }
}



template<>
void VMatrix<TYPE,RowMajor>::fillBanded(Functor<TYPE> const& functor)
{
   int kl(m_stripes[0]);
   int ku(m_stripes[1]);
   int k;

   for (int i = 0; i < m_nRows; ++i) {
       int jmin = std::max(0,i-kl);
       int jmax = std::min((int)m_nCols,i+ku+1);
//     std::cout << "j range for i = " << i << " -> (" << jmin << "..." << jmax << ")" <<std::endl;
       for (int j = jmin ; j < jmax; ++j) {
           k = j-i+kl+i*(kl+ku+1);
           m_data[k] = functor(i,j);
       }
   }
}

template<>
void VMatrix<TYPE,ColumnMajor>::fillBanded(Functor<TYPE> const& functor)
{
   int kl(m_stripes[0]);
   int ku(m_stripes[1]);
   int k;

   for (int j = 0; j < m_nCols; ++j) {
       int imin = std::max(0, j-ku);
       int imax = std::min((int)m_nRows, j+kl+1);
       //std::cout << "i range for j = " << j << " -> (" << imin << "..." << imax << ")" <<std::endl;
       for (int i = imin ; i < imax; ++i) {
           k = i-j+ku+j*(kl+ku+1);
           m_data[k] = functor(i,j);
       }
   }
}


template<>
void VMatrix<TYPE,ColumnMajor>::fillStriped(Functor<TYPE> const& functor)
{
   unsigned nStripes(m_stripes.size());
   unsigned m(std::min(m_nRows,m_nCols));

   for (unsigned k = 0; k < nStripes; ++k) {
       int offset(m_stripes[k]);
       if (offset < 0) {
          unsigned max(std::min(m_nRows + offset,m_nCols));
//        std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned j = 0; j < max; ++j) {
//            std::cout << "  setting data = " << j + k*m<< " to " 
//                      << functor(j-offset,j) << std::endl;
              m_data[j + k*m] = functor(j-offset,j);
          }
       }else {
          unsigned max(std::min(m_nRows, m_nCols-offset));
//        std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned i = 0; i < max; ++i) {
//            std::cout << "  setting data = " << i + k*m << " to " 
//                      << functor(i,i+offset) << std::endl;
              m_data[i + k*m] = functor(i,i+offset);
          }
       } 
   }
/*
   for (int i = 0; i < nStripes*m; ++i) {
       std::cout << "striped data " << i << " "<< m_data[i]<< std::endl;
   }
*/
}

template<>
void VMatrix<TYPE,RowMajor>::fillStriped(Functor<TYPE> const& functor)
{
   unsigned nStripes(m_stripes.size());
   unsigned m(std::min(m_nRows,m_nCols));

   for (unsigned k = 0; k < nStripes; ++k) {
       int offset(m_stripes[k]);
       if (offset < 0) {
          unsigned max(std::min(m_nRows + offset,m_nCols));
//        std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned j = 0; j < max; ++j) {
//            std::cout << "  setting data = " << j + k*m<< " to " 
//                      << functor(j-offset,j) << std::endl;
              m_data[j + k*m] = functor(j-offset,j);
          }
       }else {
          unsigned max(std::min(m_nRows, m_nCols-offset));
//        std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned i = 0; i < max; ++i) {
//            std::cout << "  setting data = " << i + k*m << " to " 
//                      << functor(i,i+offset) << std::endl;
              m_data[i + k*m] = functor(i,i+offset);
          }
       } 
   }
/*
   for (int i = 0; i < nStripes*m; ++i) {
       std::cout << "striped data " << i << " "<< m_data[i]<< std::endl;
   }
*/
}

template<>
void VMatrix<TYPE,RowMajor>::set(unsigned const i, unsigned const j, TYPE const value)
{
   switch (m_storage) {
      case Zero:
         break;
 
      case Diagonal:
         if (i==j) {
            m_data[i] = value;
         }
         break;
 
      case Banded:
         std::cerr << "VMatrix::set NYI for Banded matrices" << std::endl;
         break;
 
      case Striped: 
         int stripe(j-i);
         std::vector<int>::iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);

         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            unsigned m(std::min(m_nRows,m_nCols));
            unsigned index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            m_data[ij + index*m] = value;
         }
         break;
 
      case Dense:
         m_data[i*m_nCols + j] = value;
         break;
   }
}



template<>
void VMatrix<TYPE,ColumnMajor>::set(unsigned const i, unsigned const j, TYPE const value)
{
   switch (m_storage) {
      case Zero:
         break;
 
      case Diagonal:
         if (i==j) {
            m_data[i] = value;
         }
         break;
 
      case Banded:
         std::cerr << "VMatrix::set NYI for Banded matrices" << std::endl;
         break;
 
      case Striped: 
         int stripe(j-i);
         std::vector<int>::iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);

         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            unsigned m(std::min(m_nRows,m_nCols));
            unsigned index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            m_data[ij + index*m] = value;
         }
         break;
 
      case Dense:
         m_data[i + j*m_nRows] = value;
         break;
   }
}



template<>
void VMatrix<TYPE,RowMajor>::toDense()
{
   size_t n(m_nRows*m_nCols);
   TYPE* data = new TYPE[n];
   unsigned k(0);

   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j, ++k) {
           data[k] = (*this)(i,j);
       }
   }  

   release();
   m_nData = n;
   m_data = data;
   m_storage = Dense;
   m_stripes.clear();
}

template<>
void VMatrix<TYPE,ColumnMajor>::toDense()
{
   if (m_storage == Dense) return;

   size_t n(m_nRows*m_nCols);
   TYPE* data = new TYPE[n];
   unsigned k(0);

   for (unsigned j = 0; j < m_nCols; ++j) {
       for (unsigned i = 0; i < m_nRows; ++i, ++k) {
           data[k] = (*this)(i,j);
       }
   }  

   release();
   m_nData = n;
   m_data = data;
   m_storage = Dense;
   m_stripes.clear();
}
