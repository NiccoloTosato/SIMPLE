#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
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

#endif // CONFIG_HPP
