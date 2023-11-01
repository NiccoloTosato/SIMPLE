#include <iostream>
#include <cmath>
#include <utility>
#include <memory>
#include <fstream>
#include <omp.h>
#include "wrapper.hpp"
#include "cuda_runtime.h"

template <typename T=double,typename index_type=unsigned int>
class Grid {
private:

  std::unique_ptr<T[]> data;
  index_type nrows;
  index_type ncols;
  
public:
  Grid(const index_type rows,const index_type cols) : data{new T[ rows* cols]()},  nrows{ rows }, ncols{ cols} {}

  T& operator()(index_type row, index_type col) {
    return data[row*ncols + col];
  }

  const T& operator()(index_type row, index_type col) const {
    return data[row*ncols + col];
  }
  
  void swap(Grid& a) {
    //assume that 2 grid are compatible ! 
    data.swap(a.data);
    return ;
  }

  void to_device(void* device_grid) {
    cudaMemcpy(device_grid,data.get(), nrows*ncols*sizeof(T), cudaMemcpyHostToDevice);
    return;
  }
  void from_device(void* device_grid) {
    cudaMemcpy(data.get(),device_grid, nrows*ncols*sizeof(T), cudaMemcpyDeviceToHost);
    return;
  }

};

int main() {

  const unsigned int nx {64};
  const unsigned int ny {64};

  const double lx {1};
  const double ly {1};
  
  const double dx { lx / nx };
  const double dy { ly / ny };
  const double dx2 { dx * dx };
  const double dy2 { dy * dy };
  
  const double i_dx { 1 / dx };
  const double i_dy { 1 / dy };
  const double i_dx2 { 1 / dx2 };
  const double i_dy2 { 1 / dy2 };
  
  const double dt { 0.00001 };
  const double pr { 1 };
  const double i_pr { 1/pr };
  const double gr { 500 };
  const double i_gr { 1 / gr };
  const double rel { 1.1 };
  
  const double eps_p { 1.0E-11 };
  const double eps { 1.0E-11 };
  

  Grid<double,unsigned short int> vx{ny+2,nx+1};
  Grid<double,unsigned short int> vx1{ny+2,nx+1};
  
  Grid<double,unsigned short int> vy{ny+1,nx+2};
  Grid<double,unsigned short int> vy1{ny+1,nx+2};
 
  Grid<double,unsigned short int> p{ny+2,nx+2};
  Grid<double,unsigned short int> p1{ny+2,nx+2};
  Grid<double,unsigned short int> p2{ny+2,nx+2};

  Grid<double,unsigned short int> div{ny,nx};
  
  Grid<double,unsigned short int> T{ny+2,nx+2};
  Grid<double,unsigned short int> T1{ny+2,nx+2};

  double *p_device,*p1_device,*p2_device;

  cudaMalloc(&p_device, (nx+2)*(ny+2)*sizeof(double));  
  cudaMalloc(&p1_device, (nx+2)*(ny+2)*sizeof(double));
  cudaMalloc(&p2_device, (nx+2)*(ny+2)*sizeof(double));
  cudaMemset(p_device, 0, (nx+2)*(ny+2)*sizeof(double));
  cudaMemset(p1_device, 0, (nx+2)*(ny+2)*sizeof(double));
  cudaMemset(p2_device, 0, (nx+2)*(ny+2)*sizeof(double));
  double *div_device;
  cudaMalloc(&div_device, nx*ny*sizeof(double));  
  cudaMemset(div_device, 0, (nx)*(ny)*sizeof(double));



  
  
  std::cout<<"Fine init\n";
  std::size_t iteration=0;
  while(iteration < 1000) {
    std::cout<<"Iter: "<<iteration<<std::endl;
    #pragma omp parallel for
    //componente x velocità
    for( std::size_t j=1; j < ny; ++j)  //new row 
      for( std::size_t i=1; i < nx-1; ++i) { //next col
	const double vt=0.5*(vy(j-1,i)+vy(j-1,i+1));
	const double vb=0.5*(vy(j,i)+vy(j,i+1));

	const double ut=0.5*(vx(j-1,i)+vx(j,i));
	const double ub=0.5*(vx(j+1,i)+vx(j,i));
	const double ur=0.5*(vx(j,i)+vx(j,i+1));
	const double ul=0.5*(vx(j,i)+vx(j,i-1));
      
	const double a=-(ur*ur-ul*ul)*i_dx-(ub*vb-ut*vt)*i_dy; //convettivo
	const double b=(vx(j,i+1)-2*vx(j,i)+vx(j,i-1))*i_dx2; //diffusivo 1
	const double c=(vx(j+1,i)-2*vx(j,i)+vx(j-1,i))*i_dy2; //diffusivo 2
	const double AA=-a+sqrt(i_gr)*(b+c);
	vx1(j,i)=vx(j,i)+dt*(AA-(p(j,i+1)-p(j,i))*i_dx); //velocita vx al passo n+1    
      } //innermost for


    
    //componente y velocità
    #pragma omp parallel for
    for( std::size_t j=1; j < ny - 1; ++j) //new row
      for( std::size_t i=1; i < nx; ++i) { //next col
	const double ur=0.5*(vx(j,i)+vx(j+1,i));
	const double ul=0.5*(vx(j,i-1)+vx(j+1,i-1));
	const double vr=0.5*(vy(j,i+1)+vy(j,i));
	const double vl=0.5*(vy(j,i-1)+vy(j,i));
	const double vt=0.5*(vy(j-1,i)+vy(j,i));
	const double vb=0.5*(vy(j+1,i)+vy(j,i));

	const double a=-(vb*vb-vt*vt)*i_dy-(vr*ur-vl*ul)*i_dx; //convettivo
	const double b=(vy(j,i+1)-2*vy(j,i)+vy(j,i-1))*i_dx2; //diffusivo 1
	const double c=(vy(j+1,i)-2*vy(j,i)+vy(j-1,i))*i_dy2; //diffusivo 2
	const double temp=-(T1(j+1,i)+T1(j,i))*0.5; //termine galleggiamento
	const double BB=-a-temp+sqrt(i_gr)*(b+c); //approssimo gr>>1 gr/re^2
	vy1(j,i)=vy(j,i)+dt*(BB-(p(j+1,i)-p(j,i))*i_dy); //velocita vy al passo n+1
      } //innermost for

    //condizioni al contorno vx
    //parete left right impermeabilità
    for(std::size_t i=1; i<ny; ++i) {
      vx1(i, 0)=0;
      vx1(i, nx)=0;
    }
    
    //parete top e bottom noslip
    for(std::size_t i=0; i < nx+1 ; ++i) {
      vx1(0,i)=2*0-vx1(1,i);
      vx1(ny+1,i)=2*0-vx1(ny,i);
    }
    
    //condizioni al contorno vy
    //parete left right noslip
    for(std::size_t i=0; i<ny+1;++i) {
      vy1(i,0)=2*0-vy1(i,1);
      vy1(i,nx+1)=2*0-vy1(i,nx);
    }

    //parete top e bottom impermiabilità
    for(std::size_t i=1; i< nx;++i) {
      vy1(0,i)=0;
      vy1(ny,i)=0;
    }

    #pragma omp parallel for
    //calcolo la divergenza, sarà source di poisson 
    for(std::size_t j=1; j< ny;++j)
      for(std::size_t i=1; i< nx;++i)
	div(j,i)=(+(vx1(j,i)-vx1(j,i-1))*i_dx+(vy1(j,i)-vy1(j-1,i))*i_dy)*dx2*dy2/dt;

    //prepare poisson solver
    const double a=dx2*0.5/(dy2+dx2);
    const double b=dy2*0.5/(dy2+dx2);
    const double c=0.5/(dx2+dy2);


    div.to_device(div_device);
    poisson_wrapper(p1_device,p2_device,div_device,nx,ny,a,b,c);
    gpuErrchk( cudaPeekAtLastError() );
    p1.from_device(p1_device);

    /*    
    //start poisson solver, now
    std::size_t it=0;
    while(it<1500) {
      it++;
      //      p2.swap(p1);

      #pragma omp parallel for
      //laplace operator
      for (std::size_t j=1;j<ny;++j)
	for (std::size_t i=1;i<ny;++i) {
	  double al=a*(p1(j+1,i)+p1(j-1,i))+b*(p1(j,i+1)+p1(j,i-1))-c*div(j,i);
	  p1(j,i)=al*rel+(1-rel)*p1(j,i);
	}
      
      //forziamo le condizioni al contorno
      for(std::size_t i=0;i<ny+1;++i) {
	p1(i,0)=p1(i,1); //parete sx
	p1(i,nx+1)=p1(i,nx); //parete dx
      }
      
      for(std::size_t i=0;i<nx+1;++i) {
       	p1(0,i)=p1(1,i); //parete top
       	p1(ny+1,i)=p1(ny,i); //parete bottom
      }

      //condizioni contorno agli spigoli
      p1(0,0)=p1(1,1);
      p1(0,nx+1)=p1(1,nx);
      p1(ny+1,nx+1)=p1(ny,nx);
      p1(ny+1,0)=p1(ny,1);

      //TO DO: CHECK CONVERGENCE
      
    }
    */
    #pragma omp parallel for
    //correggere velocita x
    for(std::size_t j=1;j<ny;++j)
      for(std::size_t i=1;i<nx-1;++i) 
	vx1(j,i)=vx1(j,i)-(p1(j,i+1)-p1(j,i))*i_dx*dt;
    #pragma omp parallel for
    //correggere velocita y
    for(std::size_t j=1;j<ny-1;++j)
      for(std::size_t i=1;i<nx  ;++i) 
	vy1(j,i)=vy1(j,i)-(p1(j+1,i)-p1(j,i))*i_dy*dt;

    #pragma omp parallel for
    for(std::size_t j=1;j<ny+1;++j)
      for(std::size_t i=1;i<nx+1;++i) {
	
	const double u=0.5*(vx1(j,i)+vx1(j,i-1));
	const double v=0.5*(vy1(j,i)+vy1(j-1,i));

	const double a=0.5*u*i_dx*(T(j,i+1)-T(j,i-1));
	const double b=0.5*v*i_dy*(T(j+1,i)-T(j-1,i));
      	const double c=sqrt(i_gr)*i_pr*((T(j,i+1)-2*T(j,i)+T(j,i-1))*i_dx2+(T(j+1,i)-2*T(j,i)+T(j-1,i))*i_dy2);
	T1(j,i)=T(j,i)+dt*(-a-b+c);
      }


    //condizioni al contorno temperatura
    for (std::size_t i=1;i<ny+1;++i) {
      T1(i,0)=2*1-T1(i,1); //parete sx fredda
      T1(i,nx+1)=2*0-T1(i,nx); //parete dx sorgente calda
    }
    for (std::size_t i=1;i<nx+1;++i) {
      T1(0,i)=T1(1,i); //parete top adiabatica
      T1(ny+1,i)=T1(ny,i); //parete bottom adiabatica
    }

    T.swap(T1);
    vx.swap(vx1);
    vy.swap(vy1);
    for(std::size_t j=0;j<ny+1;++j)
      for(std::size_t i=0;i<nx+1;++i){
	p(j,i)=p1(j,i)+p(j,i);
	p1(j,i)=0;
      }
    
    iteration++;
   
  }

  Grid<double,unsigned short int> vxs{ny,nx};
  Grid<double,unsigned short int> vys{ny,nx};
  #pragma omp parallel for
  for(std::size_t j=0;j<ny;++j)
    for(std::size_t i=0;i<nx;++i) {
      vys(j,i)=0.5*(vy1(j,i+1)+vy1(j+1,i+1));
      vxs(j,i)=0.5*(vx1(j+1,i+1)+vx1(j+1,i));
    }



  //////////////////////////ALTAMENTE SPERIMENTALE E SCHIFOSO///////////////////////////////
  std::ofstream vx_file;
  vx_file.open("vex.bin", std::ios::binary | std::ios::out);
  vx_file.write((const char*)&vxs(0,0), sizeof(double)*nx*ny);
  vx_file.close();

  std::ofstream vy_file;
  vy_file.open("vey.bin", std::ios::binary | std::ios::out);
  vy_file.write((const char*)&vys(0,0), sizeof(double)*nx*ny);
  vy_file.close();

  std::ofstream temp_file;
  temp_file.open("temp.bin", std::ios::binary | std::ios::out);
  temp_file.write((const char*)&T(0,0), sizeof(double)*(nx+2)*(ny+2));
  temp_file.close();
  //////////////////////////################################///////////////////////////////

}
