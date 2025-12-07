import numpy as np
import matplotlib.pyplot as plt
import h5py

# Load HDF5 file
with h5py.File("simulation_results.h5", "r") as f:
    # Read datasets
    vex = f["/velocity_x"][:]
    vey = f["/velocity_y"][:]
    temp = f["/temperature"][:]
    
    # Read metadata attributes
    nx = f.attrs["nx"]
    ny = f.attrs["ny"]
    lx = f.attrs["lx"]
    ly = f.attrs["ly"]
    gr = f.attrs["grashof_number"]
    pr = f.attrs["prandtl_number"]
    dt = f.attrs["dt"]
    iterations = f.attrs["iterations"]
    
    print("=" * 50)
    print("Simulation Parameters from HDF5:")
    print("=" * 50)
    print(f"Grid size:         {nx} x {ny}")
    print(f"Domain size:       {lx} x {ly}")
    print(f"Grashof number:    {gr}")
    print(f"Prandtl number:    {pr}")
    print(f"Time step:         {dt}")
    print(f"Iterations:        {iterations}")
    print("=" * 50)

# Compute grid spacing
dx = lx / nx
dy = ly / ny

# Create mesh grid
Y, X = np.mgrid[0:ly:dy, 0:lx:dx]

# Compute velocity magnitude
velocity_magnitude = np.sqrt(vex**2 + vey**2)

# Figure 1: Streamlines with velocity magnitude colormap
plt.figure(1, figsize=(10, 8))
speed = np.sqrt(vex**2 + vey**2)
strm = plt.streamplot(X, Y, vex, vey, 
                      density=[4.5, 4.5], linewidth=0.9, arrowsize=0.5,
                      color=speed, cmap='viridis')
plt.colorbar(strm.lines, label='Velocity Magnitude')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title(f'Velocity Field - Streamlines (Gr={gr:.0f}, Pr={pr:.1f})', 
          fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3, linestyle='--')
plt.axis('equal')
plt.tight_layout()

# Figure 2: Temperature with heat map + isolines
plt.figure(2, figsize=(10, 8))
# Remove ghost cells for temperature
temp_interior = temp[1:-1, 1:-1]

# Heat map with plasma colormap (better for scientific visualization)
heat_map = plt.contourf(X, Y, temp_interior, 100, cmap='plasma', alpha=0.95)
cbar = plt.colorbar(heat_map, label='Temperature', pad=0.02)
cbar.ax.tick_params(labelsize=10)

# Overlay isolines (contour lines)
contour_lines = plt.contour(X, Y, temp_interior, 15, colors='black', linewidths=0.8, alpha=0.6)
plt.clabel(contour_lines, inline=True, fontsize=9, fmt='%.3f')

plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title(f'Temperature Field - Heat Map with Isolines (Gr={gr:.0f}, Pr={pr:.1f})', 
          fontsize=14, fontweight='bold')
plt.axis('equal')
plt.xlim(0, lx)
plt.ylim(0, ly)
plt.tight_layout()

plt.show()
