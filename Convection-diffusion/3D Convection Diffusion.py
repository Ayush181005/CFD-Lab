import numpy as np
from matplotlib import pyplot as plt, cm
# For animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation

# Variables
l = 2
nx = 21 # nx = ny = nz
dx = l/(nx-1) # dx = dy = dz
nu = 0.3
courant = 0.1
dt = courant * dx**2 / nu
tolerance = 1e-5 # threshold for convergence

# Generating the mesh
x = np.linspace(0, l, nx)
y = np.linspace(0, l, nx)
z = np.linspace(0, l, nx)

zvalue = 0.5 # the part of z where we want to see the contour plot

# Initialization
u = np.ones((nx, nx, nx))
v = np.ones((nx, nx, nx))
w = np.ones((nx, nx, nx))
un = np.ones((nx, nx, nx))
vn = np.ones((nx, nx, nx))
wn = np.ones((nx, nx, nx))

# Intital and Boundary Condition
u[int(.5/dx):int(1/dx+1),int(.5/dx):int(1/dx+1),int(.5/dx):int(1/dx+1)]=2
v[int(.5/dx):int(1/dx+1),int(.5/dx):int(1/dx+1),int(.5/dx):int(1/dx+1)]=2
w[int(.5/dx):int(1/dx+1),int(.5/dx):int(1/dx+1),int(.5/dx):int(1/dx+1)]=2

# Function to update the field
def conv_diffuse (u, v, w):
    global dt, dx, nu

    un = u.copy()
    vn = v.copy()
    wn = w.copy()

    u[1:-1,1:-1,1:-1] = un[1:-1,1:-1,1:-1] + (dt/dx**2)*(nu*(un[2:,1:-1,1:-1] + un[0:-2,1:-1,1:-1] + un[1:-1,2:,1:-1] + un[1:-1,0:-2,1:-1] + un[1:-1,1:-1,2:] + un[1:-1,1:-1,0:-2] - (6*un[1:-1,1:-1,1:-1])) - dx*(un[1:-1,1:-1,1:-1]*(un[2:,1:-1,1:-1] - un[1:-1,1:-1,1:-1]) + vn[1:-1,1:-1,1:-1]*(un[1:-1,2:,1:-1] - un[1:-1,1:-1,1:-1]) + wn[1:-1,1:-1,1:-1]*(un[1:-1,1:-1,2:] - un[1:-1,1:-1,1:-1])))

    v[1:-1,1:-1,1:-1] = vn[1:-1,1:-1,1:-1] + (dt/dx**2)*(nu*(vn[2:,1:-1,1:-1] + vn[0:-2,1:-1,1:-1] + vn[1:-1,2:,1:-1] + vn[1:-1,0:-2,1:-1] + vn[1:-1,1:-1,2:] + vn[1:-1,1:-1,0:-2] - (6*vn[1:-1,1:-1,1:-1])) - dx*(un[1:-1,1:-1,1:-1]*(vn[2:,1:-1,1:-1] - vn[1:-1,1:-1,1:-1]) + vn[1:-1,1:-1,1:-1]*(vn[1:-1,2:,1:-1] - vn[1:-1,1:-1,1:-1]) + wn[1:-1,1:-1,1:-1]*(vn[1:-1,1:-1,2:] - vn[1:-1,1:-1,1:-1])))

    w[1:-1,1:-1,1:-1] = wn[1:-1,1:-1,1:-1] + (dt/dx**2)*(nu*(wn[2:,1:-1,1:-1] + wn[0:-2,1:-1,1:-1] + wn[1:-1,2:,1:-1] + wn[1:-1,0:-2,1:-1] + wn[1:-1,1:-1,2:] + wn[1:-1,1:-1,0:-2] - (6*wn[1:-1,1:-1,1:-1])) - dx*(un[1:-1,1:-1,1:-1]*(wn[2:,1:-1,1:-1] - wn[1:-1,1:-1,1:-1]) + vn[1:-1,1:-1,1:-1]*(wn[1:-1,2:,1:-1] - wn[1:-1,1:-1,1:-1]) + wn[1:-1,1:-1,1:-1]*(wn[1:-1,1:-1,2:] - wn[1:-1,1:-1,1:-1])))

    # Boundary Conditions
    u[:,0,:] = 1 # Bottom face
    u[:,-1,:] = 1 # Top face
    u[0,:,:] = 1 # Left face
    u[-1,:,:] = 1 # Right face
    u[:,:,0] = 1 # Back face
    u[:,:,-1] = 1 # Front face

    v[:,0,:] = 1 # Bottom face
    v[:,-1,:] = 1 # Top face
    v[0,:,:] = 1 # Left face
    v[-1,:,:] = 1 # Right face
    v[:,:,0] = 1 # Back face
    v[:,:,-1] = 1 # Front face

    w[:,0,:] = 1 # Bottom face
    w[:,-1,:] = 1 # Top face
    w[0,:,:] = 1 # Left face
    w[-1,:,:] = 1 # Right face
    w[:,:,0] = 1 # Back face
    w[:,:,-1] = 1 # Front face

    return u,v,w

# Post Processing
velocity = np.sqrt(u**2 + v**2 + w**2)

# Plot the contour AT GIVEN Z LEVEL
# fig, ax = plt.subplots(figsize=(10, 8))
X, Y,Z = np.meshgrid(x, y,z)
# contour = ax.contourf(X, Y, velocity[:, :, int(nx*zvalue)], alpha=0.75, cmap=cm.viridis)
# ax.set_xlabel('$x$')
# ax.set_ylabel('$y$')
# ax.set_title(f'2D Contour Plot of 3D Diffused Field at z = 0.5 times length {l}')
# plt.colorbar(contour)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.scatter(X,Y,Z, c=velocity, cmap=cm.viridis, alpha=0.01)
ax.set_xlabel('$u$')
ax.set_ylabel('$v$')
ax.set_zlabel('$w$')
ax.set_title(f'3D Contour Plot')
# plt.show()

def update(frame):
    global u, v, w, un, vn, tolerance, ani

    un = u.copy()
    vn = v.copy()
    wn = w.copy()

    u, v, w = conv_diffuse(u, v, w)

    velocity = np.sqrt(u**2 + v**2 + w**2)

    # ax.clear()
    # contour = ax.contourf(X, Y, velocity[:, :, int(nx*zvalue)], alpha=0.75, cmap=cm.viridis) # Intital Contour Plot
    # ax.set_xlabel('$x$')
    # ax.set_ylabel('$y$')
    # ax.set_title(f'2D Contour Plot of 3D Diffused Field at z = {zvalue} times length {l}')

    # fig = plt.figure()
    plt.cla()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.scatter(X,Y,Z, c=velocity, cmap=cm.viridis, alpha=0.1)
    ax.set_xlabel('$u$')
    ax.set_ylabel('$v$')
    ax.set_zlabel('$w$')
    ax.set_title(f'3D Scatter Plot of 3D Velocity Field (frame {frame})')
    # print(frame)

    diff = np.linalg.norm(u-un) + np.linalg.norm(v-vn) + np.linalg.norm(w-wn)
    if diff<tolerance:
        ani.event_source.stop()
    return ax.collections

# create the animation
ani = FuncAnimation(fig, update, frames=200, interval = 100, blit=False)
# ani.save('3D_Convection_Diffusion.gif', writer='imagemagick', fps=60)
plt.show()