from spectral2D import *
import numpy as np

L = 2*np.pi
N = 256
DX = L/N
GRID = np.meshgrid(np.arange(-L/2,L/2,DX), np.arange(-L/2,L/2,DX))
GRID_X, GRID_Y = GRID
DT = 0.01
alpha = 1
kappa = 1
mu = np.meshgrid(2*np.pi*np.fft.fftfreq(N, d=DX),2*np.pi*np.fft.fftfreq(N, d=DX))
mu_sq = mu[0]**2 + mu[1]**2

def update(frame):
    global psi
    psi  = np.fft.fft2(psi*np.exp(-0.5j*DT*potential(GRID)))
    psi = np.exp(-0.5j*DT*alpha*mu_sq)*psi
    psi = np.fft.ifft2(psi)
    psi = psi*np.exp(-0.5j*DT*(kappa*np.absolute(psi)**2 + potential(GRID)))
    ax.clear(), ax1.clear(), ax2.clear(), ax3.clear()