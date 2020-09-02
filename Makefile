CXX = g++ -O2
CXXFLAGS = -flax-vector-conversions -std=c++11
LIBS = -framework Accelerate



OBJECTS = VMatrix.o BlockMatrix.o MatMult.o

%.o : %.C %.h
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@

timing: $(OBJECTS) Timing.o
	$(CXX) -o timing $(LIBS) $(OBJECTS) Timing.o

test: $(OBJECTS) Test.o
	$(CXX) -o test $(LIBS) $(OBJECTS) Test.o

spherium: $(OBJECTS) spherium.o
	$(CXX) -o spherium $(LIBS) $(OBJECTS) spherium.o

clean:
	rm -f $(OBJECTS) spherium.o test.o spherium test
