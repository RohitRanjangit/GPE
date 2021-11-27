import spectral2D as sp2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, cm

psi0 = sp2.psi
steps = 9385
dt = 0.001

psi_pp = psi0
psi_p = psi0



def update(frame):
    global psi_pp, psi_p
    v_t = np.fft.fft2(sp2.v_shi(sp2.GRID, psi_p))
    shi_cube_t = np.fft.fft2(sp2.shi_cube(psi_p))
    mul = np.exp(-0.5j*dt*sp2.alpha*sp2.mu_sq)
    psi = np.fft.fft2(psi_pp)*(mul**2) - 2j*dt*(v_t + sp2.kappa* shi_cube_t)*mul
    psi = np.fft.ifft2(psi)
    psi_pp = psi_p
    psi_p = psi
    
    ax.clear(), ax1.clear(), ax2.clear(), ax3.clear()
    ax.set_title('Time evolution real part of the wavefunction at t={}'.format(np.round(frame*dt, 2)))
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    ax.set_zlabel('real($\psi(x,y)$)')
    
    ax1.set_title('Time evolution imaginary part of the wavefunction at t={}'.format(np.round(frame*dt, 2)))
    ax1.set_ylabel('y')
    ax1.set_xlabel('x')
    ax1.set_zlabel('imag($\psi(x,y)$)')
    
    ax2.set_title('Time evolution of the absolute of the wavefunction at t={}'.format(np.round(frame*dt, 2)))
    ax2.set_ylabel('y')
    ax2.set_xlabel('x')
    ax2.set_zlabel('$\psi(x,y)$')
    
    ax3.set_title('HEAT PLOT: Time evolution of the wavefunction at t={}'.format(np.round(frame*dt, 2)))
    ax3.set_ylabel('y')
    ax3.set_xlabel('x')
    
    ax.plot_surface(sp2.GRID_X, sp2.GRID_Y, np.real(psi), cmap=cm.jet, linewidth=0, antialiased=False)
    ax1.plot_surface(sp2.GRID_X, sp2.GRID_Y, np.imag(psi), cmap=cm.jet, linewidth=0, antialiased=False)
    ax2.plot_surface(sp2.GRID_X, sp2.GRID_Y, np.absolute(psi), cmap=cm.jet, linewidth=0, antialiased=False)
    ax3.imshow(np.absolute(psi), cmap=cm.jet)

if __name__ == '__main__':
    fig = plt.figure(figsize=(14,8))
    fig.suptitle('Time evolution of the wavefunction-2D using leapfrog scheme',fontsize=16)

    ax = fig.add_subplot(221, projection='3d')
    ax1 = fig.add_subplot(222, projection='3d')
    ax2 = fig.add_subplot(223, projection='3d')
    ax3 = fig.add_subplot(224)
    ani = animation.FuncAnimation(fig, update, interval=10, frames=steps, repeat=False)
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
    ani.save('leapfrog-2D.gif', fps=30,writer='ffmpeg')
    #plt.show()
    
def leapfrog2D():
    global psi_pp, psi_p
    v_t = np.fft.fft2(sp2.v_shi(sp2.GRID, psi_p))
    shi_cube_t = np.fft.fft2(sp2.shi_cube(psi_p))
    mul = np.exp(-0.5j*dt*sp2.alpha*sp2.mu_sq)
    psi = np.fft.fft2(psi_pp)*(mul**2) - 2j*dt*(v_t + sp2.kappa* shi_cube_t)*mul
    psi = np.fft.ifft2(psi)
    psi_pp = psi_p
    psi_p = psi
