
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, animation


L = 2*np.pi
N = 256
DX = L/N
GRID = np.meshgrid(np.arange(-L/2,L/2,DX), np.arange(-L/2,L/2,DX))
GRID_X, GRID_Y = GRID
DT = 0.01
mu = np.meshgrid(np.fft.fftfreq(N, d=DX),np.fft.fftfreq(N, d=DX))
mu_sq = mu[0]**2 + mu[1]**2


def shi_init(x):
    """[summary]
    Args:
        x is type of np.meshgrid
    Returns:
        returns the initial condition of the wave function
    """
    return np.exp(-(x[0]**2+x[1]**2)/2)/np.sqrt(np.pi)


def potential(x):
    """[summary]
    V(x, y) = (x**2 + y**2) / 2
    Args:
        x is type of np.meshgrid
    Returns:
        returns the potential 
    """
    return (np.conjugate(x[0])*x[0] + np.conjugate(x[1])*x[1])/ 2


def v_shi(x, shi):
    """[summary]
    Args:
        x is type of np.meshgrid
        psi is the wave function
    Returns:
        returns the v(x)*shi
    """
    return potential(x)*shi


def shi_cube(shi):
    """[summary]
    Args:
        psi is the wave function
    Returns:
        returns the shi**3
    """
    return shi*np.conjugate(shi)*shi


psi = shi_init(GRID)
steps = 1000

fig = plt.figure(figsize=(14,8))

fig.suptitle('Time evolution of the wavefunction-2D',fontsize=16)

ax = fig.add_subplot(221, projection='3d')
ax1 = fig.add_subplot(222, projection='3d')
ax2 = fig.add_subplot(223, projection='3d')
ax3 = fig.add_subplot(224)

def update(frame):
    global psi
    
    #2-D integration of psi
    #print(np.sum(np.abs(psi)**2)*DX**2)
    
    psi  = np.fft.fft2(psi*np.exp(-0.5j*DT*potential(GRID)))
    psi = np.exp(-0.5j*DT*mu_sq)*psi
    psi = np.fft.ifft2(psi)
    psi = psi*np.exp(-0.5j*DT*(np.absolute(psi)**2 + potential(GRID)))
    ax.clear(), ax1.clear(), ax2.clear(), ax3.clear()
    
    ax.set_title('Time evolution real part of the wavefunction at t={}'.format(np.round(frame*DT, 2)))
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    ax.set_zlabel('real($\psi(x,y)$)')
    
    ax1.set_title('Time evolution imaginary part of the wavefunction at t={}'.format(np.round(frame*DT, 2)))
    ax1.set_ylabel('y')
    ax1.set_xlabel('x')
    ax1.set_zlabel('imag($\psi(x,y)$)')
    
    ax2.set_title('Time evolution of the absolute of the wavefunction at t={}'.format(np.round(frame*DT, 2)))
    ax2.set_ylabel('y')
    ax2.set_xlabel('x')
    ax2.set_zlabel('$\psi(x,y)$')
    
    ax3.set_title('HEAT PLOT: Time evolution of the wavefunction at t={}'.format(np.round(frame*DT, 2)))
    ax3.set_ylabel('y')
    ax3.set_xlabel('x')
    
    ax.plot_surface(GRID_X, GRID_Y, np.real(psi), cmap=cm.jet, linewidth=0, antialiased=False)
    ax1.plot_surface(GRID_X, GRID_Y, np.imag(psi), cmap=cm.jet, linewidth=0, antialiased=False)
    ax2.plot_surface(GRID_X, GRID_Y, np.absolute(psi), cmap=cm.jet, linewidth=0, antialiased=False)
    ax3.imshow(np.absolute(psi), cmap=cm.jet)

ani = animation.FuncAnimation(fig, update, interval=10, frames=steps, repeat=False)
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
#ani.save('animation-2D.gif', fps=30,writer='ffmpeg')
plt.show()



