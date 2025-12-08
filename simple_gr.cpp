#include <iostream>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <hdf5.h>
#include "config.hpp"
#include "simulation.hpp"

int main(int argc, char *argv[])
{

  // Load configuration
  std::string config_file = "config.json";
  if (argc > 1)
  {
    config_file = argv[1];
  } else {
    std::cout << "No configuration file provided. Using default: config.json" << std::endl;
    return 1;
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
  std::cout << "========================================\n"
            << std::endl;

  // Create simulation parameters structure
  SimulationParams params(cfg);
  
  // Base parameters from config (for compatibility)
  const unsigned int nx = params.nx;
  const unsigned int ny = params.ny;
  const double lx = cfg.lx;
  const double ly = cfg.ly;
  const double dt = params.dt;
  const double pr = cfg.pr;
  const double gr = cfg.gr;
  const double rel = cfg.rel;
  const double dx2 = params.dx2;
  const double dy2 = params.dy2;
  const double i_dx = params.i_dx;
  const double i_dy = params.i_dy;
  const double i_dx2 = params.i_dx2;
  const double i_dy2 = params.i_dy2;
  const double i_pr = params.i_pr;
  const double i_gr = params.i_gr;

  // Initialize field structures
  VelocityField<double, unsigned int> velocity(ny, nx);
  PressureField<double, unsigned int> pressure_field(ny, nx);
  Grid<double, unsigned int> div{ny, nx};
  Grid<double, unsigned int> T{ny + 2, nx + 2};
  Grid<double, unsigned int> T1{ny + 2, nx + 2};

  init_grids(velocity, pressure_field, T1, params);
  std::cout << "Fine init\n";
  std::size_t iteration_final = 0;
  std::size_t iteration = 0;

  #pragma omp parallel firstprivate(iteration)
  {
  auto start_time = std::chrono::steady_clock::now();
  while (iteration < cfg.max_iterations)
  {
    // Update velocity components using new functions
    update_velocity_x(velocity, pressure_field, params);
    update_velocity_y(velocity, pressure_field, T1, params);
    #pragma omp barrier
    
    #pragma omp single 
    {
    std::cout << "\rIter: " << iteration << std::flush;
    

    // condizioni al contorno vx
    // parete left right impermeabilità
    for (std::size_t i = 1; i < ny; ++i)
    {
      velocity.x1(i, 0) = 0;
      velocity.x1(i, nx) = 0;
    }

    // parete top e bottom noslip
    for (std::size_t i = 0; i < nx + 1; ++i)
    {
      velocity.x1(0, i) = 2 * 0 - velocity.x1(1, i);
      velocity.x1(ny + 1, i) = 2 * 0 - velocity.x1(ny, i);
    }

    // condizioni al contorno vy
    // parete left right noslip
    for (std::size_t i = 0; i < ny + 1; ++i)
    {
      velocity.y1(i, 0) = 2 * 0 - velocity.y1(i, 1);
      velocity.y1(i, nx + 1) = 2 * 0 - velocity.y1(i, nx);
    }

    // parete top e bottom impermiabilità
    for (std::size_t i = 1; i < nx; ++i)
    {
      velocity.y1(0, i) = 0;
      velocity.y1(ny, i) = 0;
    }
  }// end single

    // calcolo la divergenza, sarà source di poisson
    calculate_divergence(velocity, div, params);
    #pragma omp barrier
    // start poisson solver, now
    std::size_t it = 0;
    while (it < cfg.poisson_max_iterations)
    {
      it++;
      #pragma omp single
      {
      pressure_field.swap();
      }
      // laplace operator
      apply_laplace_operator(pressure_field, div, params, rel);
      #pragma omp single
      {
      // forziamo le condizioni al contorno
      for (std::size_t i = 0; i < ny + 1; ++i)
      {
        pressure_field.correction(i, 0) = pressure_field.correction(i, 1);       // parete sx
        pressure_field.correction(i, nx + 1) = pressure_field.correction(i, nx); // parete dx
      }

      for (std::size_t i = 0; i < nx + 1; ++i)
      {
        pressure_field.correction(0, i) = pressure_field.correction(1, i);       // parete top
        pressure_field.correction(ny + 1, i) = pressure_field.correction(ny, i); // parete bottom
      }

      // condizioni contorno agli spigoli
      pressure_field.correction(0, 0) = pressure_field.correction(1, 1);
      pressure_field.correction(0, nx + 1) = pressure_field.correction(1, nx);
      pressure_field.correction(ny + 1, nx + 1) = pressure_field.correction(ny, nx);
      pressure_field.correction(ny + 1, 0) = pressure_field.correction(ny, 1);
      }
      // TO DO: CHECK CONVERGENCE
    }

    #pragma omp single 
    {

    // correggere velocita x
    for (std::size_t j = 1; j < ny; ++j)
      for (std::size_t i = 1; i < nx - 1; ++i)
        velocity.x1(j, i) = velocity.x1(j, i) - (pressure_field.correction(j, i + 1) - pressure_field.correction(j, i)) * i_dx * dt;

    // correggere velocita y
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t i = 1; i < nx; ++i)
        velocity.y1(j, i) = velocity.y1(j, i) - (pressure_field.correction(j + 1, i) - pressure_field.correction(j, i)) * i_dy * dt;

    for (std::size_t j = 1; j < ny + 1; ++j)
      for (std::size_t i = 1; i < nx + 1; ++i)
      {

        const double u = 0.5 * (velocity.x1(j, i) + velocity.x1(j, i - 1));
        const double v = 0.5 * (velocity.y1(j, i) + velocity.y1(j - 1, i));

        const double a = 0.5 * u * i_dx * (T(j, i + 1) - T(j, i - 1));
        const double b = 0.5 * v * i_dy * (T(j + 1, i) - T(j - 1, i));
        const double c = sqrt(i_gr) * i_pr * ((T(j, i + 1) - 2 * T(j, i) + T(j, i - 1)) * i_dx2 + (T(j + 1, i) - 2 * T(j, i) + T(j - 1, i)) * i_dy2);
        T1(j, i) = T(j, i) + dt * (-a - b + c);
      }

    // condizioni al contorno temperatura
    for (std::size_t i = 1; i < ny + 1; ++i)
    {
      T1(i, 0) = 2 * cfg.T_left - T1(i, 1);        // parete sx fredda
      T1(i, nx + 1) = 2 * cfg.T_right - T1(i, nx); // parete dx sorgente calda
    }
    for (std::size_t i = 1; i < nx + 1; ++i)
    {
      T1(0, i) = T1(1, i);       // parete top adiabatica
      T1(ny + 1, i) = T1(ny, i); // parete bottom adiabatica
    }

    T.swap(T1);
    velocity.swap();
    for (std::size_t j = 0; j < ny + 1; ++j)
      for (std::size_t i = 0; i < nx + 1; ++i)
      {
        pressure_field.pressure(j, i) = pressure_field.correction(j, i) + pressure_field.pressure(j, i);
        pressure_field.correction(j, i) = 0;
      }


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
    if(iteration % 25 == 0)
    std::cout << "\rIter: " << iteration << "/" << cfg.max_iterations
              << " [" << std::fixed << std::setprecision(2) << iter_per_sec << " it/s]"
              << " ETA: " << std::setfill('0') << std::setw(2) << eta_hours << ":"
              << std::setw(2) << eta_minutes << ":" << std::setw(2) << eta_secs
              << "    " << std::flush;
    
  } //end single
  iteration++;
    } //end while
    #pragma omp single
    iteration_final = iteration;
  } //end parallel

Grid<double, unsigned int> vxs{ny, nx};
Grid<double, unsigned int> vys{ny, nx};
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i)
    {
      vys(j, i) = 0.5 * (velocity.y1(j, i + 1) + velocity.y1(j + 1, i + 1));
      vxs(j, i) = 0.5 * (velocity.x1(j + 1, i + 1) + velocity.x1(j + 1, i));
    }

  // Save results to HDF5 file with metadata
  std::cout << "\nSaving results to HDF5 file..." << std::endl;

  hid_t file_id = H5Fcreate("simulation_results.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Save velocity fields
  hsize_t dims_vel[2] = {ny, nx};
  hid_t dataspace_vel = H5Screate_simple(2, dims_vel, NULL);

  hid_t dataset_vx = H5Dcreate2(file_id, "/velocity_x", H5T_NATIVE_DOUBLE, dataspace_vel,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_vx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vxs(0, 0));
  H5Dclose(dataset_vx);

  hid_t dataset_vy = H5Dcreate2(file_id, "/velocity_y", H5T_NATIVE_DOUBLE, dataspace_vel,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_vy, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vys(0, 0));
  H5Dclose(dataset_vy);
  H5Sclose(dataspace_vel);

  // Save temperature field
  hsize_t dims_temp[2] = {ny + 2, nx + 2};
  hid_t dataspace_temp = H5Screate_simple(2, dims_temp, NULL);
  hid_t dataset_temp = H5Dcreate2(file_id, "/temperature", H5T_NATIVE_DOUBLE, dataspace_temp,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_temp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &T(0, 0));
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

  hid_t attr_iter = H5Acreate2(file_id, "iterations", H5T_NATIVE_ULLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_iter, H5T_NATIVE_ULLONG, &iteration_final);
  H5Aclose(attr_iter);

  H5Sclose(attr_space);
  H5Fclose(file_id);

  std::cout << "Results saved to simulation_results.h5" << std::endl;
}
