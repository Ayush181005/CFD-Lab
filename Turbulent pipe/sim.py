"""
Turbulent pipe flow - mean velocity with algebraic (zero-equation) RANS closures, simple FVM
============================================================================================

WHAT IS SOLVED
  Fully developed, steady, axisymmetric TIME-AVERAGED flow in a smooth pipe (radius R).
  Only the axial mean velocity u(r) exists. Force balance on the mean flow:
      0 = G + (1/r) d/dr [ r (nu + nu_t) du/dr ]          G = -(1/rho) dp/dx = 2 u_tau^2 / R
  nu   : molecular viscosity (fluid property)
  nu_t : eddy viscosity (flow property): stands in for the stress of the turbulent
         fluctuations, the Reynolds stress  -<u'v'> = nu_t du/dr   (Boussinesq idea).
  Integrating once gives the exact shear stress (nu + nu_t) du/dr = -G r/2: linear in r,
  zero on the axis, rho*u_tau^2 at the wall. The model only decides how that stress is
  split between viscous (nu du/dr) and turbulent (nu_t du/dr) parts.

CLOSURES (nu_t = l_m^2 |du/dr|, y = R - r is the distance from the wall)
  Prandtl mixing length   : l_m = kappa*y
  van Driest damped       : l_m = kappa*y*(1 - exp(-y+/A+)),  y+ = y u_tau/nu  (viscous layer
                            damps the eddies near the wall)
  The pressure gradient is FIXED through u_tau, so the bulk velocity Ub is an OUTPUT.

BOUNDARY CONDITIONS : r = 0 axis: du/dr = 0 (symmetry, no flux) ; r = R wall: u = 0 (no slip)

FVM : N concentric rings of width h = R/N, velocity at ring centres. Integrate over a ring
  (volume V_i, face area A_f = 2 pi r_f per unit length):
      (D_w + D_e) u_i - D_w u_{i-1} - D_e u_{i+1} = G V_i ,   D_f = A_f nu_eff,f / h
  Axis face has A = 0. Wall face is only h/2 away, and nu_t = 0 there: D_wall = A_R nu/(h/2).
  The system is tridiagonal; solved here as a small dense matrix (a TDMA would be used for big N).

NONLINEAR -> PICARD : nu_t depends on du/dr, i.e. on the unknown. Freeze nu_eff from the last
  u, solve the linear system, under-relax (u = 0.25 u_new + 0.75 u), repeat until converged.

CHECKS : (1) nu_t = 0 must give laminar Poiseuille u = G/(4 nu)(R^2 - r^2).
         (2) f_D = 8 (u_tau/Ub)^2 compared with Blasius f = 0.3164 Re^-1/4 (smooth pipe, Re > 4000).
         Agreement with Blasius alone does NOT validate a turbulence model.
"""
import numpy as np
from matplotlib import pyplot as plt

R = 0.025
nu = 1.5e-5
rho = 1.2
u_tau = 0.18
kappa = 0.41
A_plus = 26.0
G = 2*u_tau**2/R
N = 120
h = R/N
max_it = 1000
tol = 1e-10
relax = 0.25

r = (np.arange(N) + 0.5)*h           # ring centres
rf = np.arange(N + 1)*h              # faces: rf[0] = axis, rf[N] = wall
V = np.pi*(rf[1:]**2 - rf[:-1]**2)   # ring volume per unit length
Af = 2*np.pi*rf                      # face area per unit length
y = R - rf                           # wall distance at the faces
print(f'N = {N}, R = {R}, nu = {nu}, u_tau = {u_tau}, Re_tau = {u_tau*R/nu:.1f}, G = {G:.4f}')

results = {}
for model in ['Laminar', 'Prandtl', 'van Driest']:
    u = np.zeros(N)
    for count in range(1, max_it + 1):
        # gradient du/dr at the faces
        slope = np.zeros(N + 1)                      # axis face: 0 (symmetry)
        slope[1:-1] = (u[1:] - u[:-1])/h             # internal faces: centre to centre
        slope[-1] = -u[-1]/(0.5*h)                   # wall face: half-cell, u_wall = 0

        # mixing length and effective viscosity at the faces
        mixing = kappa*y
        if model == 'van Driest':
            mixing = mixing*(1 - np.exp(-y*u_tau/(A_plus*nu)))
        if model == 'Laminar':
            mixing = np.zeros(N + 1)
        eff = nu + mixing**2*np.abs(slope)           # nu_eff = nu + nu_t  (nu_t = 0 at the wall)

        # tridiagonal system  (D_w + D_e) u_i - D_w u_(i-1) - D_e u_(i+1) = G V_i
        cond = Af[1:-1]*eff[1:-1]/h                  # conductance of internal faces 1..N-1
        cond_wall = Af[-1]*nu/(0.5*h)                # wall face
        diag = np.zeros(N)
        diag[:-1] += cond                            # east face of rings 0..N-2
        diag[1:] += cond                             # west face of rings 1..N-1
        diag[-1] += cond_wall
        M = np.diag(diag) - np.diag(cond, 1) - np.diag(cond, -1)
        u_new = np.linalg.solve(M, G*V)

        diff = np.max(np.abs(u_new - u))
        if count == 1:
            u = u_new.copy()                         # first pass = laminar solution (u = 0 -> nu_t = 0)
        else:
            u = relax*u_new + (1 - relax)*u
        if count > 1 and diff < tol:
            u = u_new.copy()
            break
    if count == max_it:
        raise RuntimeError(model + ' Picard iteration did not converge')

    Ub = np.sum(V*u)/(np.pi*R**2)                    # bulk velocity
    Re = 2*R*Ub/nu
    f = 8*(u_tau/Ub)**2                              # Darcy friction factor
    results[model] = (u.copy(), Ub, Re, f, count - 1)

# checks
u_lam_exact = G/(4*nu)*(R**2 - r**2)
print(f'Laminar check: max |FVM - exact| = {np.max(np.abs(results["Laminar"][0] - u_lam_exact)):.3e} m/s '
      f'(centre velocity {u_lam_exact[0]:.2f} m/s)')
for model in ['Prandtl', 'van Driest']:
    u, Ub, Re, f, count = results[model]
    f_B = 0.3164*Re**(-0.25)
    print(f'{model:11s}: Ub = {Ub:.4f} m/s, Re_D = {Re:.1f}, f_D = {f:.5f}, Blasius = {f_B:.5f}, '
          f'Picard iterations = {count}')

# plots
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for model, style in [('Laminar', 'k--'), ('Prandtl', 'C0-'), ('van Driest', 'C1-')]:
    u, Ub, Re, f, count = results[model]
    ax[0].plot(u/Ub, r/R, style, label=model)
ax[0].set_xlabel('u / Ubulk')
ax[0].set_ylabel('r / R')
ax[0].set_title('Mean velocity profiles')
ax[0].legend()

yc = R - r                                         # wall distance of ring centres
yp = yc*u_tau/nu
yl = np.logspace(0, np.log10(300), 100)
for model, style in [('Prandtl', 'C0-'), ('van Driest', 'C1-')]:
    ax[1].semilogx(yp, results[model][0]/u_tau, style, label=model)
ax[1].semilogx(yl[yl < 12], yl[yl < 12], 'k:', label='u+ = y+')
ax[1].semilogx(yl[yl > 20], np.log(yl[yl > 20])/kappa + 5.0, 'k--', label='u+ = ln(y+)/0.41 + 5')
ax[1].set_xlabel('y+')
ax[1].set_ylabel('u+')
ax[1].set_title('Law of the wall')
ax[1].legend()

re_grid = np.linspace(4000, 14000, 250)
ax[2].plot(re_grid, 0.3164*re_grid**(-0.25), 'C2-', label='Blasius')
for model, col in [('Prandtl', 'C0'), ('van Driest', 'C1')]:
    ax[2].plot(results[model][2], results[model][3], 'o', color=col, label=model)
ax[2].set_xlabel('Re (bulk)')
ax[2].set_ylabel('Darcy friction factor')
ax[2].set_title('Friction factor')
ax[2].legend()
for a in ax:
    a.grid(alpha=0.2)
plt.tight_layout()
plt.savefig('turbulent_pipe_fvm.png', dpi=120)
plt.show()