CXX = mpiicc -DMYMPI
#CXX = icc
CXXFLAGS = -std=c++11 -g -pg  -O2 -qopenmp -funroll-loops -ffast-math
LIBS = -mkl 


HEADERS = Tile.h ZeroTile.h DiagonalTile.h StripedTile.h CMTile.h util.h \
          JacobiSolver.h Functor.h TileProduct.h EigenSolver.h ConjugateSolver.h \
          TileArray.h
OBJECTS = Tile.o CMTile.o TileProduct.o JacobiSolver.o EigenSolver.o



%.o : %.C 
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@

feast.o: feast.C $(HEADERS)
	$(CXX) -c $(LIBS) $(CXXFLAGS) feast.C -o feast.o

tile_test.o: tile_test.C $(HEADERS)
	$(CXX) -c $(LIBS) $(CXXFLAGS) tile_test.C -o tile_test.o


tile_test: tile_test.o $(HEADERS) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o tile_test $(LIBS) $(OBJECTS) tile_test.o
	

feast: feast.o $(OBJECTS) 
	$(CXX) $(CXXFLAGS) -o feast $(LIBS) $(OBJECTS) feast.o

timing: timing.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o timing $(LIBS) $(OBJECTS) timing.o


blas: blas.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o blas $(LIBS) $(OBJECTS) blas.o

clean:
	rm -f $(OBJECTS)  timing.o blas.o feast.o timing blas feast
