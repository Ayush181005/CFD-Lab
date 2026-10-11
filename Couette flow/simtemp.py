import numpy as np
import matplotlib.pyplot as plt

ny = 100
H = 1
U = 1
T = 1
Re = 100
dy = H / ny
non_dimensional = (1/Re)*(T*U/H)  # nu * T / H^2

# Target physical time to compare all schemes
target_time = 0.2

F_values = [0.25, 0.45, 0.55]
colors = ['blue', 'green', 'red']

plt.figure(figsize=(8, 6))

for F, color in zip(F_values, colors):
    dt = F * dy**2 / non_dimensional
    nit = int(target_time / dt)
    
    u = np.zeros(ny)
    u_old = u.copy()
    
    # Boundary conditions
    def apply_boundary(u):
        u[0] = 0
        u[-1] = U

    for n in range(nit):
        apply_boundary(u)
        d2u_dy2 = (u_old[2:] + u_old[0:-2] - 2*u_old[1:-1]) / (dy**2)
        u[1:-1] = u_old[1:-1] + dt * non_dimensional * d2u_dy2
        u_old = u.copy()
        
        # Break early if it blows up to infinity to prevent math errors
        if np.max(np.abs(u)) > 10 * U:
            print(f"Fo={F} blew up at iteration {n}")
            break
            
    apply_boundary(u)
    
    # Plotting
    label = f'Fo = {F}'
    if F > 0.5:
        label += ' (Unstable)'
    plt.plot(u, np.linspace(0, H, ny), label=label, color=color)

plt.xlabel('Velocity (u)')
plt.ylabel('y / H')
plt.title(f'Explicit Solution Stability at t = {target_time}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('Figure_2.png')
plt.show()