CXX = g++
CXXFLAGS = 

OBJECTS = spherium.o Matrix.o BlockMatrix.o

%.o : %.C
	$(CXX) -c $(CXXFLAGS) $< -o $@


spherium: $(OBJECTS)
	$(CXX) -o spherium $(OBJECTS) 

clean:
	rm -f $(OBJECTS)
