CXX = g++
CXXFLAGS = -flax-vector-conversions
LIBS = -framework Accelerate



OBJECTS = VMatrix.o BlockMatrix.o 

%.o : %.C
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@

test: $(OBJECTS) Test.o
	$(CXX) -o test $(LIBS) $(OBJECTS) Test.o

spherium: $(OBJECTS) spherium.o
	$(CXX) -o spherium $(LIBS) $(OBJECTS) spherium.o

clean:
	rm -f $(OBJECTS) spherium.o test.o
