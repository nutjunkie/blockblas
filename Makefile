CXX = g++ 
CXXFLAGS = -flax-vector-conversions -std=c++11 -g -pg  -O2
LIBS = -framework Accelerate

HEADERS = util.h Timer.h
OBJECTS = VMatrix.o BlockMatrix.o MatMult.o

%.o : %.C %.h
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@

timing: $(OBJECTS) Timing.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o timing $(LIBS) $(OBJECTS) Timing.o

test: $(OBJECTS) Test.o $(HEADERS)
	$(CXX) -o test $(LIBS) $(OBJECTS) Test.o

clean:
	rm -f $(OBJECTS) spherium.o test.o test timing
