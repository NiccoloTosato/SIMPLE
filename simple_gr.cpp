#include <iostream>
#include <cmath>
#include <utility>
#include <memory>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <hdf5.h>
#include "libs/json.hpp"

using json = nlohmann::json;

// Configuration structure
struct Config {
  // Grid parameters
  unsigned int nx;
  unsigned int ny;
  double lx;
  double ly;
  
  // Time parameters
  double dt;
  std::size_t max_iterations;
  
  // Physical parameters
  double pr;  // Prandtl number
  double gr;  // Grashof number
  
  // Solver parameters
  double rel;  // Relaxation factor
  double eps_p;  // Pressure tolerance
  double eps;  // Convergence tolerance
  std::size_t poisson_max_iterations;
  
  // Boundary conditions
  double T_left;
  double T_right;
};

// Load configuration from JSON file
Config load_config(const std::string& filename) {
  Config config;
  
  try {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open configuration file: " + filename);
    }
    
    json j = json::parse(file);
    
    // Load grid parameters
    config.nx = j["grid"]["nx"].get<unsigned int>();
    config.ny = j["grid"]["ny"].get<unsigned int>();
    config.lx = j["grid"]["lx"].get<double>();
    config.ly = j["grid"]["ly"].get<double>();
    
    // Load time parameters
    config.dt = j["time"]["dt"].get<double>();
    config.max_iterations = j["time"]["max_iterations"].get<std::size_t>();
    
    // Load physical parameters
    config.pr = j["physical"]["prandtl_number"].get<double>();
    config.gr = j["physical"]["grashof_number"].get<double>();
    
    // Load solver parameters
    config.rel = j["solver"]["relaxation_factor"].get<double>();
    config.eps_p = j["solver"]["pressure_tolerance"].get<double>();
    config.eps = j["solver"]["convergence_tolerance"].get<double>();
    config.poisson_max_iterations = j["solver"]["poisson_max_iterations"].get<std::size_t>();
    
    // Load boundary conditions
    config.T_left = j["boundary_conditions"]["temperature_left"].get<double>();
    config.T_right = j["boundary_conditions"]["temperature_right"].get<double>();
    
    // Validate parameters
    if (config.nx == 0 || config.ny == 0) {
      throw std::runtime_error("Grid dimensions must be positive");
    }
    if (config.lx <= 0.0 || config.ly <= 0.0) {
      throw std::runtime_error("Domain lengths must be positive");
    }
    if (config.dt <= 0.0) {
      throw std::runtime_error("Time step must be positive");
    }
    if (config.pr <= 0.0) {
      throw std::runtime_error("Prandtl number must be positive");
    }
    if (config.gr <= 0.0) {
      throw std::runtime_error("Grashof number must be positive");
    }
    if (config.rel <= 0.0 || config.rel >= 2.0) {
      throw std::runtime_error("Relaxation factor must be in (0, 2)");
    }
    
    std::cout << "Configuration loaded from: " << filename << std::endl;
    std::cout << "Grid: " << config.nx << "x" << config.ny << std::endl;
    std::cout << "Grashof number: " << config.gr << ", Prandtl number: " << config.pr << std::endl;
    
  } catch (const json::exception& e) {
    throw std::runtime_error("JSON parsing error: " + std::string(e.what()));
  }
  
  return config;
}

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

};

int main(int argc, char* argv[]) {

  // Load configuration
  std::string config_file = "config.json";
  if (argc > 1) {
    config_file = argv[1];
  }
  
  Config cfg = load_config(config_file);
  
  // Print simulation configuration summary
  std::cout << "\n========================================" << std::endl;
  std::cout << "   SIMULATION CONFIGURATION SUMMARY" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Grid size:        " << cfg.nx << " x " << cfg.ny << std::endl;
  std::cout << "Domain size:      " << cfg.lx << " x " << cfg.ly << std::endl;
  std::cout << "Time step (dt):   " << cfg.dt << std::endl;
  std::cout << "Max iterations:   " << cfg.max_iterations << std::endl;
  std::cout << "Prandtl number:   " << cfg.pr << std::endl;
  std::cout << "Grashof number:   " << cfg.gr << std::endl;
  std::cout << "Relaxation:       " << cfg.rel << std::endl;
  std::cout << "Poisson max iter: " << cfg.poisson_max_iterations << std::endl;
  std::cout << "T_left (cold):    " << cfg.T_left << std::endl;
  std::cout << "T_right (hot):    " << cfg.T_right << std::endl;
  std::cout << "========================================\n" << std::endl;
  
  // Base parameters from config
  const unsigned int nx = cfg.nx;
  const unsigned int ny = cfg.ny;
  const double lx = cfg.lx;
  const double ly = cfg.ly;
  const double dt = cfg.dt;
  const double pr = cfg.pr;
  const double gr = cfg.gr;
  const double rel = cfg.rel;
  // Note: eps_p and eps are available in cfg but not used in current implementation
  
  // Compute derived parameters automatically
  const double dx = lx / nx;
  const double dy = ly / ny;
  const double dx2 = dx * dx;
  const double dy2 = dy * dy;
  
  const double i_dx = 1.0 / dx;
  const double i_dy = 1.0 / dy;
  const double i_dx2 = 1.0 / dx2;
  const double i_dy2 = 1.0 / dy2;
  
  const double i_pr = 1.0 / pr;
  const double i_gr = 1.0 / gr;
  

  Grid<double,unsigned int> vx{ny+2, nx+1};
  Grid<double,unsigned int> vx1{ny+2, nx+1};
  
  Grid<double,unsigned int> vy{ny+1, nx+2};
  Grid<double,unsigned int> vy1{ny+1, nx+2};
 
  Grid<double,unsigned int> p{ny+2, nx+2};
  Grid<double,unsigned int> p1{ny+2, nx+2};
  Grid<double,unsigned int> p2{ny+2, nx+2};

  Grid<double,unsigned int> div{ny, nx};
  
  Grid<double,unsigned int> T{ny+2, nx+2};
  Grid<double,unsigned int> T1{ny+2, nx+2};
  
  std::cout<<"Fine init\n";
  std::size_t iteration=0;
  auto start_time = std::chrono::steady_clock::now();
  while(iteration < cfg.max_iterations) {
    std::cout<<"\rIter: "<<iteration<<std::flush;
    #pragma omp parallel for schedule(static, 8)
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
    #pragma omp parallel for schedule(static, 8)
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

    #pragma omp parallel for schedule(static, 8)
    //calcolo la divergenza, sarà source di poisson 
    for(std::size_t j=1; j< ny;++j)
      for(std::size_t i=1; i< nx;++i)
	div(j,i)=(+(vx1(j,i)-vx1(j,i-1))*i_dx+(vy1(j,i)-vy1(j-1,i))*i_dy)*dx2*dy2/dt;

    //prepare poisson solver
    const double a=dx2*0.5/(dy2+dx2);
    const double b=dy2*0.5/(dy2+dx2);
    const double c=0.5/(dx2+dy2);

    //start poisson solver, now
    std::size_t it=0;
    while(it<cfg.poisson_max_iterations) {
      it++;
      p2.swap(p1);

      #pragma omp parallel for schedule(static, 8)
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

    #pragma omp parallel for schedule(static, 8)
    //correggere velocita x
    for(std::size_t j=1;j<ny;++j)
      for(std::size_t i=1;i<nx-1;++i) 
	vx1(j,i)=vx1(j,i)-(p1(j,i+1)-p1(j,i))*i_dx*dt;
    #pragma omp parallel for schedule(static, 8)
    //correggere velocita y
    for(std::size_t j=1;j<ny-1;++j)
      for(std::size_t i=1;i<nx  ;++i) 
	vy1(j,i)=vy1(j,i)-(p1(j+1,i)-p1(j,i))*i_dy*dt;

    #pragma omp parallel for schedule(static, 8)
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
      T1(i,0)=2*cfg.T_left-T1(i,1); //parete sx fredda
      T1(i,nx+1)=2*cfg.T_right-T1(i,nx); //parete dx sorgente calda
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
    
    // Calculate and display progress
    auto current_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
    double elapsed_sec = elapsed_ms / 1000.0;
    double iter_per_sec = (elapsed_sec > 0) ? iteration / elapsed_sec : 0.0;
    double remaining_iters = cfg.max_iterations - iteration;
    double eta_seconds = (iter_per_sec > 0) ? remaining_iters / iter_per_sec : 0.0;
    
    int eta_hours = static_cast<int>(eta_seconds / 3600);
    int eta_minutes = static_cast<int>((eta_seconds - eta_hours * 3600) / 60);
    int eta_secs = static_cast<int>(eta_seconds) % 60;
    
    std::cout << "\rIter: " << iteration << "/" << cfg.max_iterations 
              << " [" << std::fixed << std::setprecision(2) << iter_per_sec << " it/s]" 
              << " ETA: " << std::setfill('0') << std::setw(2) << eta_hours << ":"
              << std::setw(2) << eta_minutes << ":" << std::setw(2) << eta_secs
              << "    " << std::flush;
   
  }

  Grid<double,unsigned int> vxs{ny, nx};
  Grid<double,unsigned int> vys{ny, nx};
  #pragma omp parallel for schedule(static, 8)
  for(std::size_t j=0;j<ny;++j)
    for(std::size_t i=0;i<nx;++i) {
      vys(j,i)=0.5*(vy1(j,i+1)+vy1(j+1,i+1));
      vxs(j,i)=0.5*(vx1(j+1,i+1)+vx1(j+1,i));
    }



  // Save results to HDF5 file with metadata
  std::cout << "\nSaving results to HDF5 file..." << std::endl;
  
  hid_t file_id = H5Fcreate("simulation_results.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  // Save velocity fields
  hsize_t dims_vel[2] = {ny, nx};
  hid_t dataspace_vel = H5Screate_simple(2, dims_vel, NULL);
  
  hid_t dataset_vx = H5Dcreate2(file_id, "/velocity_x", H5T_NATIVE_DOUBLE, dataspace_vel,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_vx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vxs(0,0));
  H5Dclose(dataset_vx);
  
  hid_t dataset_vy = H5Dcreate2(file_id, "/velocity_y", H5T_NATIVE_DOUBLE, dataspace_vel,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_vy, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vys(0,0));
  H5Dclose(dataset_vy);
  H5Sclose(dataspace_vel);
  
  // Save temperature field
  hsize_t dims_temp[2] = {ny+2, nx+2};
  hid_t dataspace_temp = H5Screate_simple(2, dims_temp, NULL);
  hid_t dataset_temp = H5Dcreate2(file_id, "/temperature", H5T_NATIVE_DOUBLE, dataspace_temp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_temp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T(0,0));
  H5Dclose(dataset_temp);
  H5Sclose(dataspace_temp);
  
  // Save metadata as attributes
  hid_t attr_space = H5Screate(H5S_SCALAR);
  
  // Grid parameters
  hid_t attr_nx = H5Acreate2(file_id, "nx", H5T_NATIVE_UINT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_nx, H5T_NATIVE_UINT, &nx);
  H5Aclose(attr_nx);
  
  hid_t attr_ny = H5Acreate2(file_id, "ny", H5T_NATIVE_UINT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_ny, H5T_NATIVE_UINT, &ny);
  H5Aclose(attr_ny);
  
  hid_t attr_lx = H5Acreate2(file_id, "lx", H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_lx, H5T_NATIVE_DOUBLE, &lx);
  H5Aclose(attr_lx);
  
  hid_t attr_ly = H5Acreate2(file_id, "ly", H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_ly, H5T_NATIVE_DOUBLE, &ly);
  H5Aclose(attr_ly);
  
  // Physical parameters
  hid_t attr_gr = H5Acreate2(file_id, "grashof_number", H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_gr, H5T_NATIVE_DOUBLE, &gr);
  H5Aclose(attr_gr);
  
  hid_t attr_pr = H5Acreate2(file_id, "prandtl_number", H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_pr, H5T_NATIVE_DOUBLE, &pr);
  H5Aclose(attr_pr);
  
  hid_t attr_dt = H5Acreate2(file_id, "dt", H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_dt, H5T_NATIVE_DOUBLE, &dt);
  H5Aclose(attr_dt);
  
  std::size_t iter = iteration;
  hid_t attr_iter = H5Acreate2(file_id, "iterations", H5T_NATIVE_ULLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_iter, H5T_NATIVE_ULLONG, &iter);
  H5Aclose(attr_iter);
  
  H5Sclose(attr_space);
  H5Fclose(file_id);
  
  std::cout << "Results saved to simulation_results.h5" << std::endl;

}
