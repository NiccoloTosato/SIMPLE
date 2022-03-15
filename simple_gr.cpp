#include <iostream>
#include <cmath>
#include <utility>
#include <memory>

template <typename T=double,typename index_type=unsigned int>

class Grid {
private:

  std::unique_ptr<T[]> data;
  index_type nrows;
  index_type ncols;
  
public:
  Grid(const std::size_t rows,const std::size_t cols) : data{new T[ rows* cols]},  nrows{ rows }, ncols{ cols} {}

  T& operator()(index_type row, index_type col) {
    return data[row*ncols + col];
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
  
  const double dt { 0.00005 };
  const double pr { 1 };
  const double  gr { 10 };
  const double i_gr { 1 / gr };
  const double rel { 1 };
  
  const double eps_p { 1.0E-10 };
  const double eps { 1.0E-10 };
  

  Grid<double,unsigned short int> vx{ny+1,nx};
  Grid<double,unsigned short int> vx1{ny+1,nx};
  
  Grid<double,unsigned short int> vy{ny,nx+1};
  Grid<double,unsigned short int> vy1{ny,nx+1};
 
  Grid<double,unsigned short int> p{ny+1,nx+1};
  Grid<double,unsigned short int> p1{ny+1,nx+1};
  Grid<double,unsigned short int> p2{ny+1,nx+1};

  Grid<double,unsigned short int> div{ny,nx};
  
  Grid<double,unsigned short int> T{ny+1,nx+1};
  Grid<double,unsigned short int> T1{ny+1,nx+1};
  while(true) {
    //componente x velocità
    for( std::size_t j=1; j < ny; ++j) //new row
      for( std::size_t i=1; i < nx-1; ++i) { //next col
	const double vt=0.5*(vy(j-1,i)+vy(j-1,i+1));
	const double vb=0.5*(vy(j,i)+vy(j,i+1));
	const double ut=0.5*(vx(j-1,i)+vx(j,i));
	const double ub=0.5*(vx(j+1,i)+vx(j,i));
	const double ur=0.5*(vx(j,i)+vx(j,i+1));
	const double ul=0.5*(vx(j,i)+vx(j,i-1));
  
	const double a=-(ur*ur-ul*ul)*i_dx-(ub*vb-ut*vt)*i_dy; //convettivo
	const double b=(vx(j,i+1)-2*vx(j,i)+vx(j,i-1))*i_dx2; //diffusivo 1	// !b=(vy(j,i+1)-2*vy(j,i)+vy(j,i)-1)*i_dx2; !----ERRORE--------ERRORE--------ERRORE--------ERRORE--------ERRORE----
	const double c=(vx(j+1,i)-2*vx(j,i)+vx(j-1,i))*i_dy2; //diffusivo 2
	const double AA=-a+sqrt(i_gr)*(b+c);
	vx1(j,i)=vx(j,i)+dt*(AA-(p(j,i+1)-p(j,i))*i_dx); //velocita vx al passo n+1    
      }
    
    //componente y velocità
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
    }

    //condizioni al contorno vx
    //parete left right impermeabilità
    for(std::size_t i=0; i<ny+1; ++i) {
      vx1(i, 0)=0;
      vx1(i, nx-1)=0;
    }
    
    //parete top e bottom noslip
    for(std::size_t i=0; i < nx ; ++i) {
      vx1(0,i)=2*0-vx1(1,i);
      vx1(ny,i)=2*0-vx1(ny-1,i);
    }
    
    //condizioni al contorno vy
    //parete left right noslip
    for(std::size_t i=0; i<ny;++i) {
      vy1(i,0)=2*0-vy1(i,1);
      vy1(i,nx)=2*0-vy1(i,nx-1);
    }

    //parete top e bottom impermiabilità
    for(std::size_t i=0; i< nx+1;++i) {
      vy1(1,i)=0;
      vy1(ny-1,i)=0;
    }

    //calcolo la divergenza, sarà source di poisson 
    for(std::size_t j=1; j< ny;++j)
      for(std::size_t i=1; i< nx;++i)
	div(j,i)=(+(vx1(j,i)-vx1(j,i-1))*i_dx+(vy1(j,i)-vy1(j-1,i))*i_dy)*dx2*dy2/dt;

    //prepare poisson solver
    const double a=dx2*0.5/(dy2+dx2);
    const double b=dy2*0.5/(dy2+dx2);
    const double c=0.5/(dx2+dy2);

    //
    
    
    break;
  }
}
