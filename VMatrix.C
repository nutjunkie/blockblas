#include "VMatrix.h"
#include "Types.h"


#define TYPE double
#include "VMatrix_templateT.cpp"
#undef TYPE

#define TYPE complex
#include "VMatrix_templateT.cpp"
#undef TYPE

#define LAYOUT RowMajor
#include "VMatrix_templateL.cpp"
#undef LAYOUT

#define LAYOUT ColumnMajor
#include "VMatrix_templateL.cpp"
#undef LAYOUT
