add_library(blockblas 
   CMTile.C
   DiagonalTile.C
   EigenSolver.C
   JacobiSolver.C
   Log.C
   Tile.C
   TileProduct.C
   Diagonalize.C
)

target_include_directories( blockblas PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
