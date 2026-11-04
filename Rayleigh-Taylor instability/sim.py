"""
Rayleigh-Taylor instability - SIMPLE + VOF
rho = rho_l + (rho_h-rho_l)*alpha (alpha = 1 : heavy fluid, on top)
"""

import numpy as np
import matplotlib.pyplot as plt

lx = 1
ly = 2
rho_l = 1
rho_h = 1.2
rho_ref = 0.5*(rho_l + rho_h) # reference density
At = (rho_h-rho_l)/(rho_h+rho_l) # Atwood number
g = 1.0
mu = 2e-3
nx = 40
ny = 80
dx = lx/(nx-1)
dy = ly/(ny-1)
Co = 0.5
dt = 0.004
nt = 1500
n_inner = 3 # SIMPLE iterations per time step
alpha_p = 0.3
tol = 1e-6

CFL_max  = 0.3  # convective limit
visc_max = 0.45 # viscous limit
dt_max   = 0.01 # cap, needed because u = v = 0 at t = 0
t_end    = 6.0
nu_max   = mu/rho_l # kinematic viscosity is largest in the light fluid

def compute_dt(u, v):
    conv = np.abs(u).max()/dx + np.abs(v).max()/dy
    dt_conv = CFL_max/conv if conv > 1e-12 else np.inf
    dt_visc = visc_max/(2*nu_max*(1/dx**2 + 1/dy**2))
    return min(dt_conv, dt_visc, dt_max)

x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
X, Y = np.meshgrid(x, y)

def A(P):
    return np.maximum(0, (1 - 0.1*np.abs(P))**5) # Power Law Scheme

def solve_v(u_old, v_old, u, v, p, rho):
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

    rho_P = rho[1:-1, 1:-1]
    rho_E = rho[1:-1, 2:]
    rho_W = rho[1:-1, 0:-2]
    rho_N = rho[2:, 1:-1]
    rho_S = rho[0:-2, 1:-1]

    # faces (uniform grid)
    u_e = (u_E + u_P)/2
    u_w = (u_W + u_P)/2
    v_n = (v_N + v_P)/2
    v_s = (v_S + v_P)/2
    rho_e = (rho_E + rho_P)/2
    rho_w = (rho_W + rho_P)/2
    rho_n = (rho_N + rho_P)/2
    rho_s = (rho_S + rho_P)/2

    A_x = dy
    A_y = dx

    D_x = mu*dy/dx
    D_y = mu*dx/dy
    F_e = (rho_e*u_e)*dy
    F_w = (rho_w*u_w)*dy
    F_n = (rho_n*v_n)*dx
    F_s = (rho_s*v_s)*dx
    P_e = F_e/D_x
    P_w = F_w/D_x
    P_n = F_n/D_y
    P_s = F_s/D_y
    a_e = D_x*A(P_e) + np.maximum(-F_e, 0)
    a_w = D_x*A(P_w) + np.maximum(F_w, 0)
    a_n = D_y*A(P_n) + np.maximum(-F_n, 0)
    a_s = D_y*A(P_s) + np.maximum(F_s, 0)
    a_P0 = rho_P*dx*dy/dt
    b_u = a_P0*u_old_p
    b_v = a_P0*v_old_p - (rho_P - rho_ref)*g*dx*dy     # buoyancy (gravity in -y)
    a_P = a_e + a_w + a_n + a_s + a_P0

    u_star = u.copy()
    v_star = v.copy()
    u_star[1:-1, 1:-1] = (a_e*u_E + a_w*u_W + a_n*u_N + a_s*u_S + b_u + (p_W - p_E)/2*A_x)/a_P
    v_star[1:-1, 1:-1] = (a_e*v_E + a_w*v_W + a_n*v_N + a_s*v_S + b_v + (p_S - p_N)/2*A_y)/a_P

    return u_star, v_star, a_P

def pressure_correction(u_star, v_star, a_P, rho):
    A_x = dy
    A_y = dx
    rho_P = rho[1:-1, 1:-1]

    u_star_e = u_star[1:-1, 2:]
    u_star_w = u_star[1:-1, 0:-2]
    v_star_n = v_star[2:, 1:-1]
    v_star_s = v_star[0:-2, 1:-1]

    d_e = A_x/a_P
    d_w = A_x/a_P
    d_n = A_y/a_P
    d_s = A_y/a_P
    a_E = rho_P*d_e*dy
    a_W = rho_P*d_w*dy
    a_N = rho_P*d_n*dx
    a_S = rho_P*d_s*dx
    a_P_p = a_E + a_W + a_N + a_S
    b = (rho_P*(u_star_w - u_star_e))*dy/2 + (rho_P*(v_star_s - v_star_n))*dx/2   # /2 : central difference over 2 cells

    p_correction = np.zeros((ny, nx))
    for k in range(50): # internal loop (Jacobi) for p'
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

def solve_vof(alpha_old, u, v):
    # d(alpha)/dt + u.grad(alpha) = 0, first-order upwind (bounded)
    a = alpha_old
    aP = a[1:-1, 1:-1]
    uP = u[1:-1, 1:-1]
    vP = v[1:-1, 1:-1]
    dadx = np.where(uP > 0, aP - a[1:-1, 0:-2], a[1:-1, 2:] - aP)/dx
    dady = np.where(vP > 0, aP - a[0:-2, 1:-1], a[2:, 1:-1] - aP)/dy
    alpha_new = a.copy()
    alpha_new[1:-1, 1:-1] = aP - dt*(uP*dadx + vP*dady)
    alpha_new = np.clip(alpha_new, 0, 1)
    # zero gradient at walls
    alpha_new[:, 0] = alpha_new[:, 1]; alpha_new[:, -1] = alpha_new[:, -2]
    alpha_new[0, :] = alpha_new[1, :]; alpha_new[-1, :] = alpha_new[-2, :]
    return alpha_new

def apply_bc(u, v):
    u[:, 0] = 0; u[:, -1] = 0; u[0, :] = 0; u[-1, :] = 0
    v[:, 0] = 0; v[:, -1] = 0; v[0, :] = 0; v[-1, :] = 0

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx)) # p*

# initial interface: heavy fluid on top, single-mode perturbation
y_int = ly/2 + 0.05*lx*np.cos(2*np.pi*x/lx)
alpha = (Y > y_int[None, :]).astype(float)
rho = rho_l + (rho_h-rho_l)*alpha

t = 0.0
n = 0
snap_times = [0, 1.5, 3.0, 4.5, 6.0]
snapshots = {}
next_snap = 0

while t<t_end:
    dt = compute_dt(u, v) # updates the global dt used by solve_v / solve_vof
    dt = min(dt, t_end - t)

    u_old = u.copy()
    v_old = v.copy()
    alpha_old = alpha.copy()

    for it in range(n_inner):
        # u*, v*
        u_star, v_star, a_P = solve_v(u_old, v_old, u, v, p, rho)
        # boundary conditions
        apply_bc(u_star, v_star)
        # p'
        p_correction, residual = pressure_correction(u_star, v_star, a_P, rho)
        # p = p* + alpha_p p'
        p = p + alpha_p*p_correction
        p[:, 0] = p[:, 1]; p[:, -1] = p[:, -2]; p[0, :] = p[1, :]; p[-1, :] = p[-2, :]
        # u = u* + u', v = v* + v'
        u, v = velocity_correction(u_star, v_star, p_correction, a_P)
        apply_bc(u, v)

    # VOF + density update
    alpha = solve_vof(alpha_old, u, v)
    rho = rho_l + (rho_h-rho_l)*alpha

    t += dt
    n += 1
    if next_snap < len(snap_times) and t >= snap_times[next_snap] - 1e-12:
        snapshots[t] = (alpha.copy(), u.copy(), v.copy())
        next_snap += 1
    if n % 100 == 0:
        print(f"Step {n}, t = {t:.3f}, dt = {dt:.3e}, residual = {residual:.3e}")


# --- POST PROCESSING ---
fig, axs = plt.subplots(1, len(snapshots), figsize=(3*len(snapshots), 6), sharey=True)
fig_num = 0
for ax, (n, (a, uu, vv)) in zip(np.atleast_1d(axs), snapshots.items()):
    ax.contourf(X, Y, a, levels=np.linspace(0, 1, 21), cmap='coolwarm')
    ax.contour(X, Y, a, levels=[0.5], colors='k', linewidths=1)
    ax.set_title(f't = {snap_times[fig_num]:.2f}')
    ax.set_xlabel('X')
    ax.set_aspect('equal')
    fig_num += 1
np.atleast_1d(axs)[0].set_ylabel('Y')
plt.suptitle(f'Rayleigh-Taylor Instability (At = {At:.3f}) - alpha field')
plt.tight_layout()
plt.show()

velocity_mag = np.sqrt(u**2 + v**2)
plt.figure(figsize=(5, 8))
plt.contourf(X, Y, velocity_mag, levels=20)
plt.colorbar(label='Velocity Magnitude')
plt.streamplot(x, y, u, v, color='white', density=1.2)
plt.contour(X, Y, alpha, levels=[0.5], colors='r')
plt.title('Rayleigh-Taylor - velocity at final time')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()