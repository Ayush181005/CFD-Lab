# SIMPLE, all u, v, p stored on the same (ny, nx) mesh - lid driven cavity

import numpy as np
import matplotlib.pyplot as plt

lx = 1.0
ly = 1.0
nx = 60
ny = 60
dx = lx/(nx-1)
dy = ly/(ny-1)
nt = 20000
u_up = 1.0
rho = 1.0
Re = 100
nu = (u_up*lx)/Re
Co = 0.5
dt = Co*dx/u_up
alpha_p = 0.3
tol = 1e-6

x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
X, Y = np.meshgrid(x, y)

def A(P):
    return np.maximum(0, (1 - 0.1*np.abs(P))**5) # Power Law Scheme

def solve_v(u_old, v_old, u, v, p):
    # u_old, v_old : previous time level | u, v, p : current guess (u*, v*, p*)
    u_old_p = u_old[1:-1, 1:-1]
    v_old_p = v_old[1:-1, 1:-1]

    u_P = u[1:-1, 1:-1]
    u_E = u[1:-1, 2:]
    u_W = u[1:-1, 0:-2]
    u_N = u[2:, 1:-1]
    u_S = u[0:-2, 1:-1]
    v_P = v[1:-1, 1:-1]
    v_E = v[1:-1, 2:]
    v_W = v[1:-1, 0:-2]
    v_N = v[2:, 1:-1]
    v_S = v[0:-2, 1:-1]

    p_E = p[1:-1, 2:]
    p_W = p[1:-1, 0:-2]
    p_N = p[2:, 1:-1]
    p_S = p[0:-2, 1:-1]

    # faces (uniform grid)
    u_e = (u_E + u_P)/2
    u_w = (u_W + u_P)/2
    v_n = (v_N + v_P)/2
    v_s = (v_S + v_P)/2

    A_x = dy
    A_y = dx

    D_x = nu*dy/dx
    D_y = nu*dx/dy
    F_e = (rho*u_e)*dy
    F_w = (rho*u_w)*dy
    F_n = (rho*v_n)*dx
    F_s = (rho*v_s)*dx
    P_e = F_e/D_x
    P_w = F_w/D_x
    P_n = F_n/D_y
    P_s = F_s/D_y
    a_e = D_x*A(P_e) + np.maximum(-F_e, 0)
    a_w = D_x*A(P_w) + np.maximum(F_w, 0)
    a_n = D_y*A(P_n) + np.maximum(-F_n, 0)
    a_s = D_y*A(P_s) + np.maximum(F_s, 0)
    a_P0 = rho*dx*dy/dt
    b_u = a_P0*u_old_p
    b_v = a_P0*v_old_p
    a_P = a_e + a_w + a_n + a_s + a_P0

    u_star = u.copy()
    v_star = v.copy()
    u_star[1:-1, 1:-1] = (a_e*u_E + a_w*u_W + a_n*u_N + a_s*u_S + b_u + (p_W - p_E)/2*A_x)/a_P
    v_star[1:-1, 1:-1] = (a_e*v_E + a_w*v_W + a_n*v_N + a_s*v_S + b_v + (p_S - p_N)/2*A_y)/a_P

    return u_star, v_star, a_P

def pressure_correction(u_star, v_star, a_P):
    A_x = dy
    A_y = dx

    u_star_e = u_star[1:-1, 2:]
    u_star_w = u_star[1:-1, 0:-2]
    v_star_n = v_star[2:, 1:-1]
    v_star_s = v_star[0:-2, 1:-1]

    d_e = A_x/a_P
    d_w = A_x/a_P
    d_n = A_y/a_P
    d_s = A_y/a_P
    a_E = rho*d_e*dy
    a_W = rho*d_w*dy
    a_N = rho*d_n*dx
    a_S = rho*d_s*dx
    a_P_p = a_E + a_W + a_N + a_S
    b = (rho*(u_star_w - u_star_e))*dy/2 + (rho*(v_star_s - v_star_n))*dx/2   # /2 : central difference over 2 cells

    p_correction = np.zeros((ny, nx))
    # for k in range(50): # internal loop for time relaxation
    p_correction_E = p_correction[1:-1, 2:]
    p_correction_W = p_correction[1:-1, 0:-2]
    p_correction_N = p_correction[2:, 1:-1]
    p_correction_S = p_correction[0:-2, 1:-1]
    p_correction[1:-1, 1:-1] = (a_E*p_correction_E + a_W*p_correction_W + a_N*p_correction_N + a_S*p_correction_S + b)/a_P_p
    # zero gradient at walls
    p_correction[:, 0] = p_correction[:, 1]
    p_correction[:, -1] = p_correction[:, -2]
    p_correction[0, :] = p_correction[1, :]
    p_correction[-1, :] = p_correction[-2, :]

    residual = np.abs(b).max()
    return p_correction, residual

def velocity_correction(u_star, v_star, p_correction, a_P):
    A_x = dy
    A_y = dx
    d_e = A_x/a_P
    d_n = A_y/a_P

    u = u_star.copy()
    v = v_star.copy()
    u[1:-1, 1:-1] = u_star[1:-1, 1:-1] + d_e*(p_correction[1:-1, 0:-2] - p_correction[1:-1, 2:])/2
    v[1:-1, 1:-1] = v_star[1:-1, 1:-1] + d_n*(p_correction[0:-2, 1:-1] - p_correction[2:, 1:-1])/2
    return u, v

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx)) # p*
u[-1, :] = u_up

for n in range(nt):
    u_old = u.copy()
    v_old = v.copy()
    # u*, v*
    u_star, v_star, a_P = solve_v(u_old, v_old, u, v, p)
    # boundary conditions
    u_star[:, 0] = 0; u_star[:, -1] = 0; u_star[0, :] = 0; u_star[-1, :] = u_up
    v_star[:, 0] = 0; v_star[:, -1] = 0; v_star[0, :] = 0; v_star[-1, :] = 0
    # p'
    p_correction, residual = pressure_correction(u_star, v_star, a_P)
    # p = p* + alpha_p p'
    p = p + alpha_p*p_correction
    p[:, 0] = p[:, 1]; p[:, -1] = p[:, -2]; p[0, :] = p[1, :]; p[-1, :] = p[-2, :]
    # u = u* + u', v = v* + v'
    u, v = velocity_correction(u_star, v_star, p_correction, a_P)
    # boundary conditions
    u[:, 0] = 0; u[:, -1] = 0; u[0, :] = 0; u[-1, :] = u_up
    v[:, 0] = 0; v[:, -1] = 0; v[0, :] = 0; v[-1, :] = 0

    if n % 200 == 0:
        print(f"Step {n}, Mass Residual: {residual:.3e}")
    if residual < tol:
        print(f"Converged at step {n}")
        break

# --- POST PROCESSING ---
velocity_mag = np.sqrt(u**2 + v**2)

plt.figure(figsize=(8, 6))
plt.contourf(X, Y, velocity_mag, levels=20)
plt.colorbar(label='Velocity Magnitude')
plt.streamplot(x, y, u, v, color='white', density=1.2)
plt.title('Lid-Driven Cavity Flow - SIMPLE')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

# Validation with Ghia et al.
# at y=0.5
x_test = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9375, 0.9063, 0.8594, 0.8125, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0625, 0.0000]
# v_test = v[ny//2,x_test[:]*nx//lx]
v_ghia_100 = [0.00000, -0.05906, -0.07391, -0.08814, -0.10150, -0.11401, -0.15661, -0.20371, -0.23112, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.09266, 0.00000]
v_ghia_400 = [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.26253, -0.23827, -0.44993, -0.42447, 0.00906, 0.30174, 0.30203, 0.28299, 0.22217, 0.19837, 0.16945, 0.00000]
v_ghia_1000 = [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.43764, -0.51750, -0.42665, -0.27404, 0.02526, 0.33230, 0.33075, 0.29012, 0.21191, 0.18262, 0.14830, 0.00000]
v_ghia_10000 = [0.00000, -0.54302, -0.52987, -0.49008, -0.41341, -0.31688, -0.00412, 0.17540, 0.23186, 0.00831, -0.12143, -0.10515, 0.07889, 0.27224, 0.33157, 0.40020, 0.00000]

plt.plot(x, v[nx//2,:])
plt.scatter(x_test, v_ghia_100)
plt.title('Comparison with Ghia et al. (1982) for Re=100, v values at y=0.5')
plt.xlabel('x')
plt.ylabel('v')
plt.show()