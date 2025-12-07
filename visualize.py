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

# Figure 2: Temperature contours
plt.figure(2, figsize=(10, 8))
# Remove ghost cells for temperature
temp_interior = temp[1:-1, 1:-1]
contours = plt.contour(X, Y, temp_interior, 30, cmap='RdBu_r')
plt.clabel(contours, inline=True, fontsize=8, fmt='%.3f')
plt.colorbar(contours, label='Temperature')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title(f'Temperature Field - Contour Plot (Gr={gr:.0f}, Pr={pr:.1f})', 
          fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3, linestyle='--')
plt.axis('equal')
plt.tight_layout()

# Figure 3: Combined plot with temperature and velocity
plt.figure(3, figsize=(12, 8))
plt.contourf(X, Y, temp_interior, 50, cmap='RdBu_r', alpha=0.7)
plt.colorbar(label='Temperature')
skip = max(1, nx // 20)  # Adjust vector density based on grid size
plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
           vex[::skip, ::skip], vey[::skip, ::skip],
           scale=None, alpha=0.8, width=0.003)
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title(f'Temperature Field with Velocity Vectors (Gr={gr:.0f}, Pr={pr:.1f})', 
          fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3, linestyle='--', color='white')
plt.axis('equal')
plt.tight_layout()

plt.show()
