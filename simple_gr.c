#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main() {

  const unsigned int nx = 64;
  const unsigned int ny = 64;

  const double lx = 1;
  const double ly = 1;

  const double dx = lx / nx;
  const double dy = ly / ny;
  const double dx2 = dx * dx;
  const double dy2 = dy * dy;

  const double i_dx = 1 / dx;
  const double i_dy = 1 / dy;
  const double i_dx2 = 1 / dx2;
  const double i_dy2 = 1 / dy2;

  const double dt = 0.00005;
  const double pr = 1;
  const double gr = 10;
  const double i_gr = 1 / gr;
  const double rel = 1;

  const double eps_p = 1.0E-10;
  const double eps = 1.0E-10;

  double *vx = calloc((ny + 1) * (nx), sizeof(double));
  double *vx1 = calloc((ny + 1) * (nx), sizeof(double));

  double *vy = calloc((ny) * (nx + 1), sizeof(double));
  double *vy1 = calloc((ny) * (nx + 1), sizeof(double));

  double *p = calloc((ny + 1) * (nx + 1), sizeof(double));
  double *p1 = calloc((ny + 1) * (nx + 1), sizeof(double));
  double *p2 = calloc((ny + 1) * (nx + 1), sizeof(double));

  double *div = calloc((ny) * (nx), sizeof(double));

  double *T = calloc((ny + 1) * (nx + 1), sizeof(double));
  double *T1 = calloc((ny + 1) * (nx + 1), sizeof(double));
  while (1) {
    // x direction 
    for (unsigned int j = 1; j < ny ; ++j)
      for (unsigned int i = 1; i < nx -1; ++i) {

        double vt = 0.5 * (vy[j - 1 + i * (nx + 1)] + vy[j - 1 + (i + 1) * (nx + 1)]);
        double vb = 0.5 * (vy[j + i * (nx + 1)] + vy[j + (i + 1) * (nx + 1)]);

	double ut = 0.5 * (vx[j - 1 + i * nx] + vx[j + i * nx]);
        double ub = 0.5 * (vx[j + 1 + i * nx] + vx[j + i * nx]);
        double ur = 0.5 * (vx[j + i * nx] + vx[j + (i + 1) * nx]);
        double ul = 0.5 * (vx[j + i * nx] + vx[j + (i - 1)] * nx);

        double a = -(ur * ur - ul * ul) * i_dx - (ub * vb - ut * vt) * i_dy;
        double b = (vx[j + (i + 1) * nx] - 2 * vx[j + i * nx] + vx[j + (i - 1) * nx]) * i_dx2;
        double c = (vx[j + 1 + i*nx] - 2 * vx[j +  i*nx] + vx[j - 1 +  i*nx]) * i_dy2;

        double AA = -a + sqrt(i_gr) * (b + c);
        vx1[j + i*nx] = vx[j + i*nx] + dt * (AA - (p[j + (i + 1)*(nx+1)] - p[j+  i*nx]) * i_dx);
      }
    for (unsigned int j = 1;
	  
    break;
  }
}
