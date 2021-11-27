
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from scipy.integrate import quad

L = 2*np.pi #Length of the domain
N = 500 #Number of Grid points
DX = L/N #Grid spacing
GRID = np.arange(-L/2,L/2,DX, dtype=complex) #Grid points
DT = 0.01 #Time step
alpha = 1
kappa = 1
mu_sq  = (2*np.pi*np.fft.fftfreq(N, d=DX))**2 #Fourier frequencies

def shi_init(x):
    """This function initializes the wavefunction with a gaussian distribution
        centered at x=0 and with a width of 0.5
        shi(x,0) = exp(-x^2/2)/sqrt(pi)
    Args:
        x (complex): Grid points

    Returns:
        complex: Initial wavefunction
    """
    return np.exp(-(x**2)/2)/(np.pi**0.25)

def potential(x):
    """This function defines the potential
        V(x) = 0.5*x^2
    Args:
        x (complex): Grid points

    Returns:
        complex: Potential
    """
    return np.conjugate(x)*x/2


def shi_cube(shi):
    """This function defines the cube of the wavefunction
    Args:
        shi (complex): Wavefunction values

    Returns:
        complex: Cube of the wavefunction
    """
    return np.conjugate(shi)*shi*shi


def v_shi(x, shi):
    """This function defines the potential*wavefunction
    Args:
        x (complex): Grid points

    Returns:
        complex: Potential*wavefunction of the wavefunction
    """
    return shi*potential(x)


psi = shi_init(GRID)
steps = 1000

psi0 = psi

def update(frame):
    global psi
    #print(np.sum(np.abs(psi)**2)*DX)
    psi  = np.fft.fft(psi*np.exp(-0.5j*DT*potential(GRID)))
    psi = np.exp(-0.5j*DT*alpha*mu_sq)*psi
    psi = np.fft.ifft(psi)
    psi = psi*np.exp(-0.5j*DT*(kappa*np.absolute(psi)**2 + potential(GRID)))
    ax.clear()
    ax2.clear()
    ax.set_xlabel('x')
    ax.set_ylabel('$\psi(x)$')
    ax2.set_xlabel('x')
    ax2.set_ylabel('$\psi(x)$')
    ax.set_title('Time evolution of the wavefunction at time = {}s'.format(np.round(frame*DT, 1)))
    ax2.set_title('Heat plot for the absolute of wavefunction at time = {}s'.format(np.round(frame*DT, 1)))
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
    fig.suptitle('Time evolution of the wavefunction-1D using TSSP',fontsize=16)
    ani = animation.FuncAnimation(fig, update, frames=steps,  interval=50, repeat=False)
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
    ani.save('spectral-1D.mp4', writer='ffmpeg', fps=30)
    #save ani as gif
    #ani.save('spectral-1D.gif', writer='ffmpeg', fps=30)
    #plt.show()
    

def spectral1D():
    global psi
    #print(np.sum(np.abs(psi)**2)*DX)
    psi  = np.fft.fft(psi*np.exp(-0.5j*DT*potential(GRID)))
    psi = np.exp(-0.5j*DT*alpha*mu_sq)*psi
    psi = np.fft.ifft(psi)
    psi = psi*np.exp(-0.5j*DT*(kappa*np.absolute(psi)**2 + potential(GRID)))
    



