# Makefile for SIMPLE fluid dynamics simulation
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -march=native -fopenmp -I.
LDFLAGS = -lhdf5

TARGET = simple_gr
SOURCES = simple_gr.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS) *.bin *.out *.h5

# Optional: run target
run: $(TARGET)
	./$(TARGET)
