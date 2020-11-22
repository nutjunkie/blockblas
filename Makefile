CXX = icc
CXXFLAGS = -std=c++11 -g -pg  -O2 -fopenmp -funroll-loops -ffast-math
LIBS = -mkl

HEADERS = util.h Timer.h VMatrix.h BlockMatrix.h JacobiSolver.h Functor.h
OBJECTS =  MatMult.o VMatrix.o

%.o : %.C %.h
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@


feast: feast.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o feast $(LIBS) $(OBJECTS) feast.o

timing: timing.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o timing $(LIBS) $(OBJECTS) timing.o

test: test.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o test $(LIBS) $(OBJECTS) test.o

blas: blas.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o blas $(LIBS) $(OBJECTS) blas.o

clean:
	rm -f $(OBJECTS) test.o test timing.o timing blas.o blas

VMatrix.o: VMatrix_templateT.cpp VMatrix_templateL.cpp
