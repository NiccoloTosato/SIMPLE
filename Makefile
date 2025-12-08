# Makefile for SIMPLE fluid dynamics simulation
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -march=native -fopenmp -I.
LDFLAGS = -lhdf5

TARGET = simple_gr
TEST_TARGET = compare_hdf5
SOURCES = simple_gr.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

$(TEST_TARGET): compare_hdf5.cpp
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) compare_hdf5.cpp $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: $(TEST_TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGET) $(OBJECTS) *.bin *.out

# Optional: run target
run: $(TARGET)
	./$(TARGET)
