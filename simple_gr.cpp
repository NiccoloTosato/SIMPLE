#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
    for (unsigned int j = 1; j < ny; ++j)
      for (unsigned int i = 1; i < nx - 1; ++i) {

        double vt =
            0.5 * (vy[j - 1 + i * (nx + 1)] + vy[j - 1 + (i + 1) * (nx + 1)]);
        double vb = 0.5 * (vy[j + i * (nx + 1)] + vy[j + (i + 1) * (nx + 1)]);

        double ut = 0.5 * (vx[j - 1 + i * nx] + vx[j + i * nx]);
        double ub = 0.5 * (vx[j + 1 + i * nx] + vx[j + i * nx]);
        double ur = 0.5 * (vx[j + i * nx] + vx[j + (i + 1) * nx]);
        double ul = 0.5 * (vx[j + i * nx] + vx[j + (i - 1)] * nx);

      }
    
    // y direction

    for (unsigned int j = 1; j < ny - 1; ++j)
      for (unsigned int i = 1; i < nx; ++i) {

        double ur = 0.5 * (vx[j + i * nx] + vx[j + 1 + i * nx]);
        double ul = 0.5 * (vx[j + (i - 1) * nx] + vx[j + 1 + (i-1) * nx]);
        double vr = 0.5 * (vy[j + i*(nx+1) + 1] + vy[j+ i*(nx+1)]);
	double vl =0.5 * (vy[j + (i - 1)*(nx+1)] + vy[j + i*(nx+1)]) ;
	double vt =0.5 * (vy[j - 1 + i*(nx+1)] + vy[j + i*(nx+1)]);
	double vb =   0.5 * (vy[j + 1 + i*(nx+1)] + vy[j + i*(nx+1)]);
	
	double a=-(vb*vb-vt*vt)*i_dy-(vr*ur-vl*ul)*i_dx; //convettivo
	double b=(vy[j + (i+1) * (nx+1)]-2*vy[j + i * (nx+1)]+vy[ j + (i-1) * (nx+1) ])*i_dx2; //diffusivo 1
	double c=(vy[j+1+i*(nx+1)]-2*vy[j + i*(nx+1)]+vy[j-1+i*(nx+1)])*i_dy2; //diffusivo 2
	double temp=-(T1[j+1+i*(nx+1)]+T1[j+i*(nx+1)])*0.5; //termine galleggiamento
	double BB=-a-temp+sqrt(i_gr)*(b+c); //approssimo gr>>1 gr/re^2
	vy1[j + i*(nx+1)]=vy[j + i*(nx+1)]+dt*(BB-(p[j+1+i*(nx+1)]-p[j + i*(nx+1)])*i_dy); //velocita vy al passo n+1
      }

    break;
  }
}
