"""
1D wave propagation (linear shallow water) - FVM on the characteristic variables r and l
========================================================================================

d(eta)/dt + H d(u)/dx = 0,   d(u)/dt + g d(eta)/dx = 0,   c = sqrt(g*H)

Characteristic variables (they decouple the two equations):
    r = eta + (H/c) u   moves right:  r_t + c r_x = 0
    l = eta - (H/c) u   moves left :  l_t - c l_x = 0
Recover:  eta = (r + l)/2,   u = c/(2H) (r - l)

FVM, first-order upwind (take the value from the side the wave comes from):
    r_i : from the left cell  (i-1)
    l_i : from the right cell (i+1)
  EXPLICIT : r_i^{n+1} = r_i^n - C (r_i^n - r_{i-1}^n),   l_i^{n+1} = l_i^n - C (l_i^n - l_{i+1}^n)
             C = c dt/dx, stable for C <= 1. Only ONE time level is needed (new = old + dt*()).
  IMPLICIT : backward Euler,  (1 + C) r_i^{n+1} - C r_{i-1}^{n+1} = r_i^n   (same for l)
             unconditionally stable; solved here with Jacobi sweeps.
Wall (u = 0)  ->  r = l at the wall.
Initial condition: hump at rest (u = 0), so r = l = eta.

Upwind is first-order: numerical diffusion ~ c dx/2 (1 - C) explicit, c dx/2 (1 + C) implicit,
so the pulse is smeared (peak drops), but it never produces wiggles.
Compare with wave_1d_fvm_slice.py (second-order eta scheme): keeps the height, can disperse.

Time step: dt = CFL*dx/c. The number of steps comes from t_end (nt = ceil(t_end/dt)),
and the last step is shortened so the run ends exactly at t_end (needed to compare with exact).
"""
import numpy as np
from matplotlib import pyplot as plt

g = 9.81
H = 1.0
c = np.sqrt(g*H)
L = 10.0
N = 200
dx = L/N
x = (np.arange(N) + 0.5)*dx          # N cell centres
t_end = 0.8
x0 = L/2
w = 0.4
A = 0.1

eta0 = A*np.exp(-((x - x0)/w)**2)
eta_exact = 0.5*A*(np.exp(-((x - x0 - c*t_end)/w)**2) + np.exp(-((x - x0 + c*t_end)/w)**2))

# explicit
CFL = 0.8
dt = CFL*dx/c
nt = int(np.ceil(t_end/dt))
print(f'Explicit: dt = {dt:.4f}, CFL = {c*dt/dx:.3f}, steps = {nt}')

r = eta0.copy()                      # u = 0 initially -> r = l = eta
l = eta0.copy()
for n in range(nt):
    dt_n = min(dt, t_end - n*dt)     # last step shortened to land on t_end
    r_old = r.copy()
    l_old = l.copy()
    r[1:] = r_old[1:] - c*dt_n/dx*(r_old[1:] - r_old[0:-1])        # r moves right: upwind = left cell
    l[:-1] = l_old[:-1] - c*dt_n/dx*(l_old[:-1] - l_old[1:])       # l moves left : upwind = right cell
    # Boundary (wall, u = 0 -> r = l)
    r[0] = l[0]
    l[-1] = r[-1]
eta_explicit = 0.5*(r + l)

# implicit
C = 4
dt = C*dx/c
nt = int(np.ceil(t_end/dt))
tol = 1e-10
max_it = 5000
print(f'Implicit: dt = {dt:.4f}, CFL = {C:.3f}, steps = {nt}')

r = eta0.copy()
l = eta0.copy()
total_it = 0
for n in range(nt):
    dt_n = min(dt, t_end - n*dt)     # last step shortened to land on t_end
    C_n = c*dt_n/dx
    r_old = r.copy()
    l_old = l.copy()
    for it in range(max_it):
        r_it = r.copy()
        l_it = l.copy()
        r_it[1:] = (r_old[1:] + C_n*r[0:-1])/(1 + C_n)
        l_it[:-1] = (l_old[:-1] + C_n*l[1:])/(1 + C_n)
        # Boundary
        r_it[0] = l_it[0]
        l_it[-1] = r_it[-1]
        err = max(np.max(np.abs(r_it - r)), np.max(np.abs(l_it - l)))
        r = r_it.copy()
        l = l_it.copy()
        if err < tol:
            break
    total_it += it + 1
eta_implicit = 0.5*(r + l)
print(f'Implicit: total Jacobi sweeps = {total_it}')

# ------------------------------------------------------------------ check
m = x < x0
print(f'Peak (left pulse): exact = {eta_exact[m].max():.4f}, '
      f'explicit = {eta_explicit[m].max():.4f}, implicit = {eta_implicit[m].max():.4f}')
for name, e in [('Explicit', eta_explicit), ('Implicit', eta_implicit)]:
    print(f'{name}: RMS error = {np.sqrt(np.mean((e - eta_exact)**2)):.3e}')

plt.figure(figsize=(9, 4.5))
plt.plot(x, eta0, 'k:', label='initial')
plt.plot(x, eta_exact, 'k-', label=f'exact at t={t_end} s')
plt.plot(x, eta_explicit, label=f'explicit, CFL={CFL}')
plt.plot(x, eta_implicit, label=f'implicit, CFL={C}')
plt.xlabel('x')
plt.ylabel('eta')
plt.title(f'1D wave, characteristic upwind FVM, c = {c:.2f}')
plt.legend()
plt.tight_layout()
plt.show()