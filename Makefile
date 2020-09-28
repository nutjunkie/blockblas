CXX = g++ 
CXXFLAGS = -flax-vector-conversions -std=c++11 -g -pg  -O3 -fopenmp -funroll-loops -ffast-math
LIBS = -framework Accelerate

HEADERS = util.h Timer.h
OBJECTS = VMatrix.o BlockMatrix.o MatMult.o JacobiSolver.o

%.o : %.C %.h
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@

timing: $(OBJECTS) timing.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o timing $(LIBS) $(OBJECTS) timing.o

test: $(OBJECTS) test.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o test $(LIBS) $(OBJECTS) test.o

blas: blas.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o blas $(LIBS) $(OBJECTS) blas.o

clean:
	rm -f $(OBJECTS) spherium.o test.o test timing timing.o test.o
