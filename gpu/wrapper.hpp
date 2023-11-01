#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#include "cuda_runtime.h"
#include <stdio.h>
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void poisson_wrapper(double* p, double* p1, double* div, int nx, int ny,double a,double b,double c);
