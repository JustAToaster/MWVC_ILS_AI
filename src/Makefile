CXX = g++
CXXFLAGS = -Wall -O3

OBJECTS = main.o Graph.o optimization.o utilities.o

vertex_cover_ils: $(OBJECTS)
		$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJECTS): Graph.h optimization.h utilities.h