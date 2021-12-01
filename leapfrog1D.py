"""
The License of this code is reserved to the author of the code,
Author: Rohit Ranjan
This code is not for commercial use. and it is taken from the github repository
of any other author.
All Right Reserved @Rohit Ranjan.
"""


import spectral1D as sp1
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
psi0 = sp1.psi0
steps = 1000
dt = 0.01

psi_pp = psi0
psi_p = psi0

#psi_p = 

L = sp1.L
GRID = sp1.GRID

def update(frame):
    global psi_pp, psi_p
    v_t = np.fft.fft(sp1.v_shi(sp1.GRID, psi_p))
    shi_cube_t = np.fft.fft(sp1.shi_cube(psi_p))
    mul = np.exp(-0.5j*dt*sp1.alpha*sp1.mu_sq)
    psi = np.fft.fft(psi_pp)*(mul**2) - 2j*dt*(v_t + sp1.kappa* shi_cube_t)*mul
    psi = np.fft.ifft(psi)
    psi_pp = psi_p
    psi_p = psi
    
    ax.clear()
    ax2.clear()
    ax.set_xlabel('x')
    ax.set_ylabel('$\psi(x)$')
    ax2.set_xlabel('x')
    ax2.set_ylabel('$\psi(x)$')
    ax.set_title('Time evolution of the wavefunction at time = {}s'.format(np.round(frame*dt, 1)))
    ax2.set_title('Heat plot for the absolute of wavefunction at time = {}s'.format(np.round(frame*dt, 1)))
    ax.set_xlim(-L/2,L/2)
    ax.set_ylim(-1.0,1.0)
    ax.plot(np.real(GRID), np.real(psi), 'r-')
    ax.plot(np.real(GRID), np.imag(psi), 'b-')
    ax.plot(np.real(GRID), np.abs(psi), 'g-')
    ax.plot(np.real(GRID), np.real(psi0), 'r--')
    ax.legend(['Real', 'Imaginary', 'Absolute', 'Initial'])
    # ax2.imshow(np.real(psi)[np.newaxis, :], cmap='plasma', interpolation='nearest', extent=[-L/2, L/2, 0, 1])
    # ax2.imshow(np.imag(psi)[np.newaxis, :], cmap='plasma', interpolation='nearest', extent=[-L/2, L/2, 0, 1])
    ax2.imshow(np.abs(psi)[np.newaxis, :], cmap='plasma', interpolation='nearest', extent=[-L/2, L/2, 0, 1])

if __name__ == '__main__':
    fig, (ax, ax2) = plt.subplots(nrows=1, ncols =2, figsize=(12,5))
    fig.suptitle('Time evolution of the wavefunction-1D using Leapfrog scheme',fontsize=16)
    ani  = animation.FuncAnimation(fig, update, frames=steps, interval=10)
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
    #ani.save('leapfrog-1D.mp4', writer='ffmpeg', fps=30)
    #save ani as gif
    #ani.save('leapfrog-1D.gif', writer='ffmpeg', fps=30)
    plt.show()

def leapfrog1D():
    global psi_pp, psi_p
    v_t = np.fft.fft(sp1.v_shi(sp1.GRID, psi_p))
    shi_cube_t = np.fft.fft(sp1.shi_cube(psi_p))
    mul = np.exp(-0.5j*dt*sp1.alpha*sp1.mu_sq)
    psi = np.fft.fft(psi_pp)*(mul**2) - 2j*dt*(v_t + sp1.kappa* shi_cube_t)*mul
    psi = np.fft.ifft(psi)
    psi_pp = psi_p
    psi_p = psi
