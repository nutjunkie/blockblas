#include "BlockMatrix.h"
#include "Types.h"


#define TYPE double
#include "BlockMatrix_templateT.cpp"
#undef TYPE

#define TYPE complex
#include "BlockMatrix_templateT.cpp"
#undef TYPE

#define LAYOUT RowMajor
#include "BlockMatrix_templateL.cpp"
#undef LAYOUT

#define LAYOUT ColumnMajor
#include "BlockMatrix_templateL.cpp"
#undef LAYOUT
