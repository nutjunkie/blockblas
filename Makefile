CXX = g++
CXXFLAGS = -flax-vector-conversions
LIBS = -framework Accelerate



OBJECTS = spherium.o VMatrix.o BlockMatrix.o

%.o : %.C
	$(CXX) -c $(LIBS) $(CXXFLAGS) $< -o $@


spherium: $(OBJECTS)
	$(CXX) -o spherium $(LIBS) $(OBJECTS) 

clean:
	rm -f $(OBJECTS)
