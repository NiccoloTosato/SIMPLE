#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cmath>
#include <utility>
#include <memory>
#include "config.hpp"

// Grid template class for storing 2D data
template <typename T=double, typename index_type=unsigned int>
class Grid {
private:
  std::unique_ptr<T[]> data;
  index_type nrows;
  index_type ncols;
  
public:
  Grid(const index_type rows, const index_type cols) 
    : data{new T[rows * cols]()}, nrows{rows}, ncols{cols} {}

  T& operator()(index_type row, index_type col) {
    return data[row * ncols + col];
  }

  const T& operator()(index_type row, index_type col) const {
    return data[row * ncols + col];
  }
  
  void swap(Grid& a) {
    // Assume that 2 grids are compatible
    data.swap(a.data);
    return;
  }
};

// Velocity field structure with automatic swap
template <typename T=double, typename index_type=unsigned int>
struct VelocityField {
  Grid<T, index_type> x;      // vx component
  Grid<T, index_type> x1;     // vx temporary/next step
  Grid<T, index_type> y;      // vy component
  Grid<T, index_type> y1;     // vy temporary/next step
  
  // Constructor
  VelocityField(index_type ny, index_type nx)
    : x{ny+2, nx+1},   // vx dimensions
      x1{ny+2, nx+1},
      y{ny+1, nx+2},   // vy dimensions
      y1{ny+1, nx+2} {}
  
  // Swap current and next time step for both components
  void swap() {
    x.swap(x1);
    y.swap(y1);
  }
};

// Pressure field structure with automatic swap
template <typename T=double, typename index_type=unsigned int>
struct PressureField {
  Grid<T, index_type> pressure;      // Actual pressure field
  Grid<T, index_type> correction;    // Pressure correction (p1)
  Grid<T, index_type> correction_old; // Previous correction (p2)
  
  // Constructor
  PressureField(index_type ny, index_type nx)
    : pressure{ny+2, nx+2},
      correction{ny+2, nx+2},
      correction_old{ny+2, nx+2} {}
  
  // Swap for Poisson solver iteration (swap correction fields)
  void swap() {
    correction_old.swap(correction);
  }
};

// Structure for simulation parameters (grid size and derived/computed parameters)
struct SimulationParams {
  unsigned int nx, ny;    // Grid size
  double dx, dy;          // Grid spacing
  double dx2, dy2;        // Grid spacing squared
  double i_dx, i_dy;      // Inverse grid spacing
  double i_dx2, i_dy2;    // Inverse grid spacing squared
  double i_pr, i_gr;      // Inverse Prandtl and Grashof numbers
  double dt;              // Time step
  
  // Constructor from Config
  SimulationParams(const Config& cfg)
    : nx(cfg.nx),
      ny(cfg.ny),
      dx(cfg.lx / cfg.nx),
      dy(cfg.ly / cfg.ny),
      dx2(dx * dx),
      dy2(dy * dy),
      i_dx(1.0 / dx),
      i_dy(1.0 / dy),
      i_dx2(1.0 / dx2),
      i_dy2(1.0 / dy2),
      i_pr(1.0 / cfg.pr),
      i_gr(1.0 / cfg.gr),
      dt(cfg.dt) {}
};

// Update x component of velocity
template <typename T, typename index_type>
void update_velocity_x(
    VelocityField<T, index_type>& velocity,
    PressureField<T, index_type>& pressure_field,
    const SimulationParams& params)
{
  #pragma omp parallel for schedule(static, 8)
  for (std::size_t j = 1; j < params.ny; ++j)
    for (std::size_t i = 1; i < params.nx - 1; ++i)
    {
      const double vt = 0.5 * (velocity.y(j - 1, i) + velocity.y(j - 1, i + 1));
      const double vb = 0.5 * (velocity.y(j, i) + velocity.y(j, i + 1));

      const double ut = 0.5 * (velocity.x(j - 1, i) + velocity.x(j, i));
      const double ub = 0.5 * (velocity.x(j + 1, i) + velocity.x(j, i));
      const double ur = 0.5 * (velocity.x(j, i) + velocity.x(j, i + 1));
      const double ul = 0.5 * (velocity.x(j, i) + velocity.x(j, i - 1));

      const double a = -(ur * ur - ul * ul) * params.i_dx - (ub * vb - ut * vt) * params.i_dy;
      const double b = (velocity.x(j, i + 1) - 2 * velocity.x(j, i) + velocity.x(j, i - 1)) * params.i_dx2;
      const double c = (velocity.x(j + 1, i) - 2 * velocity.x(j, i) + velocity.x(j - 1, i)) * params.i_dy2;
      const double AA = -a + sqrt(params.i_gr) * (b + c);
      velocity.x1(j, i) = velocity.x(j, i) + params.dt * (AA - (pressure_field.pressure(j, i + 1) - pressure_field.pressure(j, i)) * params.i_dx);
    }
}

// Update y component of velocity
template <typename T, typename index_type>
void update_velocity_y(
    VelocityField<T, index_type>& velocity,
    PressureField<T, index_type>& pressure_field,
    const Grid<T, index_type>& T1,
    const SimulationParams& params)
{
  #pragma omp parallel for schedule(static, 8)
  for (std::size_t j = 1; j < params.ny - 1; ++j)
    for (std::size_t i = 1; i < params.nx; ++i)
    {
      const double ur = 0.5 * (velocity.x(j, i) + velocity.x(j + 1, i));
      const double ul = 0.5 * (velocity.x(j, i - 1) + velocity.x(j + 1, i - 1));
      const double vr = 0.5 * (velocity.y(j, i + 1) + velocity.y(j, i));
      const double vl = 0.5 * (velocity.y(j, i - 1) + velocity.y(j, i));
      const double vt = 0.5 * (velocity.y(j - 1, i) + velocity.y(j, i));
      const double vb = 0.5 * (velocity.y(j + 1, i) + velocity.y(j, i));

      const double a = -(vb * vb - vt * vt) * params.i_dy - (vr * ur - vl * ul) * params.i_dx;
      const double b = (velocity.y(j, i + 1) - 2 * velocity.y(j, i) + velocity.y(j, i - 1)) * params.i_dx2;
      const double c = (velocity.y(j + 1, i) - 2 * velocity.y(j, i) + velocity.y(j - 1, i)) * params.i_dy2;
      const double temp = -(T1(j + 1, i) + T1(j, i)) * 0.5;
      const double BB = -a - temp + sqrt(params.i_gr) * (b + c);
      velocity.y1(j, i) = velocity.y(j, i) + params.dt * (BB - (pressure_field.pressure(j + 1, i) - pressure_field.pressure(j, i)) * params.i_dy);
    }
}

// Calculate divergence of velocity field
template <typename T, typename index_type>
void calculate_divergence(
    const VelocityField<T, index_type>& velocity,
    Grid<T, index_type>& div,
    const SimulationParams& params)
{
  #pragma omp parallel for schedule(static, 8)
  for (std::size_t j = 1; j < params.ny; ++j)
    for (std::size_t i = 1; i < params.nx; ++i)
      div(j, i) = (+(velocity.x1(j, i) - velocity.x1(j, i - 1)) * params.i_dx + 
                   (velocity.y1(j, i) - velocity.y1(j - 1, i)) * params.i_dy) * 
                   params.dx2 * params.dy2 / params.dt;
}

// Apply Laplace operator for Poisson solver
template <typename T, typename index_type>
void apply_laplace_operator(
    PressureField<T, index_type>& pressure_field,
    const Grid<T, index_type>& div,
    const SimulationParams& params,
    double relaxation)
{
  const double a = params.dx2 * 0.5 / (params.dy2 + params.dx2);
  const double b = params.dy2 * 0.5 / (params.dy2 + params.dx2);
  const double c = 0.5 / (params.dx2 + params.dy2);
  
  #pragma omp parallel for schedule(static, 8)
  for (std::size_t j = 1; j < params.ny; ++j)
    for (std::size_t i = 1; i < params.nx; ++i)
    {
      double al = a * (pressure_field.correction(j + 1, i) + pressure_field.correction(j - 1, i)) + 
                  b * (pressure_field.correction(j, i + 1) + pressure_field.correction(j, i - 1)) - 
                  c * div(j, i);
      pressure_field.correction(j, i) = al * relaxation + (1 - relaxation) * pressure_field.correction(j, i);
    }
}

#endif // SIMULATION_HPP
