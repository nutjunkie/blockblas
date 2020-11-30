/******************************************************************************
 * 
 *  Member function definitions involving type template specialization.
 *
 *****************************************************************************/



template<>
void BlockMatrix<TYPE,ColumnMajor>::bind(TYPE const* data)
{
   unsigned ld(nRows());

   for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
       TYPE const* d0(data+colOffset(bj)*ld);
       for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
           (*this)(bi,bj).bind(d0,ld); 
           d0 += (*this)(bi,bj).nRows();
       }
   }
}


template<>
void BlockMatrix<TYPE,ColumnMajor>::unbind(TYPE* data)
{
   unsigned ld(nRows());

   for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
       TYPE* d0(data+colOffset(bj)*ld);
       for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
           (*this)(bi,bj).unbind(d0,ld); 
           d0 += (*this)(bi,bj).nRows();
       }
   }
}
