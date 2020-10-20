CXX = icc
CXXFLAGS = -std=c++11 -g -pg  -O2 -fopenmp -funroll-loops -ffast-math
LIBS = -mkl

HEADERS = util.h Timer.h VMatrix.h BlockMatrix.h JacobiSolver.h Functor.h
OBJECTS =  MatMult.o VMatrix.o

%.o : %.C %.h
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@

timing: $(OBJECTS) timing.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o timing $(LIBS) $(OBJECTS) timing.o

test: $(OBJECTS) test.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o test $(LIBS) $(OBJECTS) test.o

blas: blas.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o blas $(LIBS) $(OBJECTS) blas.o

clean:
	rm -f $(OBJECTS) test.o test timing.o timing blas.o blas
