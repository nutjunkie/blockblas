/******************************************************************************
 * 
 *  Member function definitions involving type template specialization.
 *
 *****************************************************************************/


template<>
BlockMatrix<complex,LAYOUT>& 
   BlockMatrix<complex,LAYOUT>::fromDouble(BlockMatrix<double,LAYOUT> const& that)
{
   m_nRowBlocks = that.m_nRowBlocks;
   m_nColBlocks = that.m_nColBlocks;

   if (m_blocks) delete [] m_blocks;
   m_blocks = new VMatrix<complex,LAYOUT>[m_nRowBlocks*m_nColBlocks];

   for (unsigned row = 0; row < m_nRowBlocks; ++row) {
       for (unsigned col = 0; col < m_nColBlocks; ++col) {
           m_blocks[row*m_nColBlocks + col].fromDouble(that(row,col));
      }
   }
   return *this;
}
