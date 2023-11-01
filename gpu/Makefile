# makefile
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -fopenmp
NVCC = nvcc
NVCCFLAGS = -O3

LDFLAGS = -lcudart

OBJS = main.o poisson.o

all: program

program: $(OBJS)
	$(CXX) $(CXXFLAGS) -L${CUDA_PATH}/lib64 poisson.o main.o -o program.out $(LDFLAGS) 

main.o: main.cpp wrapper.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

poisson.o: poisson.cu wrapper.hpp
	$(NVCC) $(NVCCFLAGS) -c poisson.cu 
clean:
	rm *.out *.o
