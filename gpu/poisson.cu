#include"wrapper.hpp"
#include<stdio.h>

__global__ void poisson(double* p, double* p1, double* div, int nx, int ny,double a,double b,double c){
  //considero le ghost cell nella dimensione del dominio, la pressione e' [nx+2,ny+2], la divergenza invece [nx,ny]
   int nrows=ny+2;
   int ncols=nx+2;
   int nrows_div=ny;
   int ncols_div=nx;
   int idx_x = threadIdx.x + blockIdx.x * blockDim.x;
   int idx_y = threadIdx.y + blockIdx.y * blockDim.y;
   //printf("%d %d\n",idx_x,idx_y);

   //relax factor di jacobi, mettiamo 1 
   double rel=1.0;
   if((idx_x>0) && (idx_x <nx+1) && (idx_y>0) && (idx_y < ny+1)) {
     double al=a*(p[ncols*(idx_y+1)+idx_x]+p[ncols*(idx_y-1)+idx_x])+b*(p[ncols*(idx_y)+idx_x+1]+p[ncols*(idx_y)+idx_x-1])-c*div[ncols_div*(idx_y)+idx_x];
     p1[ncols*(idx_y)+idx_x]=al*rel+(1-rel)*p[ncols*(idx_y)+idx_x];
   }
   //forziamo le condizioni al contorno                                                                                                                                                   
   if ((idx_x==0) && (idx_y < ny+1)) 
     p1[ncols*(idx_y)+0]=p1[ncols*(idx_y)+1]; //parete sx
   if(idx_x==(nx+1) && (idx_y < ny+1))	
     p1[ncols*(idx_y)+nx+1]=p1[ncols*(idx_y)+nx]; //parete dx
   if(idx_y==0 && (idx_x < nx+1))
     p1[ncols*(0)+idx_x]=p1[ncols*(1)+idx_x]; //parete top
   if(idx_y==(nx+1) && (idx_x < nx+1))
     p1[ncols*(ny+1)+idx_x]=p1[ncols*(ny)+idx_x] ; //parete bottom

   //condizioni contorno agli spigoli
   if((idx_x==0) && (idx_y==0)) {
     p1[0]=p1[(ncols)*1+1];
     p1[0+nx+1]=p1[ncols*1+nx];
     p1[ncols*(ny+1)+nx+1]=p1[ncols*ny+nx];
     p1[ncols*(ny+1)+0]=p1[ncols*ny+1];
   }

   
}

void poisson_wrapper(double* p, double* p1, double* div, int nx, int ny,double a,double b,double c){
  int it=0;
  while(it<3000) {
    it++;
     

    poisson<<<dim3(ceil((nx+2)/16.0),ceil((nx+2)/16.0),1),dim3(16,16,1)>>>(p, p1,div, nx, ny, a,b,c);
    gpuErrchk( cudaPeekAtLastError() );
    double *tmp;
    tmp=p;
    p=p1;
    p1=tmp;
  
  }
}

