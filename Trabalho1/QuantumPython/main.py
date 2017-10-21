import numpy as np
from matplotlib import pyplot as pl, animation
from scipy.fftpack import fft,ifft

def time_step():
    global t, psi_mod_x, psi_mod_k, t_max
    psi_mod_x *= x_evolve_half
    for i in range(N_steps):
        psi_mod_k = fft(psi_mod_x)
        psi_mod_k *= k_evolve
        psi_mod_x = ifft(psi_mod_k)
        psi_mod_x *= x_evolve_half*x_evolve_half
    psi_mod_k = fft(psi_mod_x)
    t += dt * N_steps
    t_max = t_max - 1
#print(t_max)

def get_psi_x():
    return (psi_mod_x * np.exp(1j * k[0] * x) * np.sqrt(2 * np.pi) / dx)

def get_psi_k():
    return psi_mod_k * np.exp(-1j * x[0] * dk * np.arange(Nx))

def heaviside(x):
    x,y = np.asarray(x), np.zeros(x.shape)
    y[x > 0] = 1.0
    return y

def f1(x):
    return 3.*np.sqrt(-(x/7.)**2.+1.)

def f2(x):
    return -f1(x)

def f3(x):
    return np.abs(x/2.)-((3.*np.sqrt(33.)-7.)*x**2.)/112.+np.sqrt(1.-(np.abs(np.abs(x)-2.)-1.)**2.)-3.

def f4(x):
    return 9.-8.*np.abs(x)

def f5(x):
    return 3.*np.abs(x)+.75

def f6(x):
    return 2.25

def f7(x):
    return 1.5 - .5*np.abs(x)-6.*np.sqrt(10.)*(np.sqrt(3.-x**2.+2.*np.abs(x))-2.)/14.

def top(x):
    if(np.abs(x) > 7.):
        return 0.
    if(np.abs(x) > 3.):
        return f1(x)
    if(np.abs(x) > .75 and abs(x) < 1.):
        return f4(x)
    if(np.abs(x) < .5):
        return f6(x)
    if(np.abs(x) > 1.):
        return f7(x)

def bottom(x):
    if(np.abs(x) > 7.):
        return 0.
    if(np.abs(x) < 4):
        return f3(x)
    if(np.abs(x) >= 4):
        return f2(x)

def batman_barrier(x,width):
    x,y = np.asarray(x), np.zeros(x.shape)
    for i in range(len(x)):
        #if(i%2==0):
        #    y[i] = bottom(x[i])
        #else:
        #    y[i] = top(x[i])
        y[i] = bottom(x[i])
        y[i] += 3
    y[x < -width] = 0.
    y[x > width] = 0.
    return y

def triangular_barrier(x, width, height):
    x,y = np.asarray(x), np.zeros(x.shape)
    y[x < 0] = (height/width)*x[x < 0] + height
    y[x >= 0] = -(height/width)*x[x >= 0] + height
    y[x < -width] = 0.
    y[x > width] = 0.
    return y

def square_barrier(x, width, height):
    return height * (heaviside(x) - heaviside(x - width))

def gauss_x(x, d, x0, k0):
    return ((d * np.sqrt(np.pi)) ** (-0.5)
            * np.exp(-0.5 * ((x - x0) * 1. / d) ** 2 + 1j * x * k0))

dt, N_steps, t_max = 0.005, 50, 460
hbar, m, N, dx = 1.0, 1.0, 2 ** 11, 0.1
x = dx * (np.arange(N) - 0.5 * N)
V0 = 1.5
L = hbar / np.sqrt(2 * m * V0)
a, x0 = 3 * L, -60 * L
V_x = square_barrier(x, a, V0)
#V_x = triangular_barrier(x, 3.5, V0)
#V_x = batman_barrier(x,7)
#print(*V_x)
V_x[x < -100] = 1E6
V_x[x > 100] = 1E6
#changed 0.2 to 0.4
p0 = np.sqrt(2 * m * 0.4 * V0)
d = hbar / np.sqrt(2 * p0 ** 2. * 1./80)
psi_x0 = gauss_x(x,d,x0,p0/hbar)
k0, t = -28, 0.0
Nx = len(x)
x, psi_x0 = map(np.asarray, (x, psi_x0))
dx = x[1] - x[0]
dk = 2 * np.pi / (Nx * dx)
x_evolve_half = np.exp(-0.5 * 1j * V_x / hbar * dt )
k = k0 + dk * np.arange(Nx)
k_evolve = np.exp(-0.5 * 1j * hbar / m * (k * k) * dt)
psi_x = np.asarray(psi_x0)
psi_mod_x = (psi_x * np.exp(-1j * k[0] * x) * dx / np.sqrt(2 * np.pi))
psi_mod_k = fft(psi_mod_x)

def writeToFile():
    fileX = open("QuantumPythonX.txt", "w")
    fileK = open("QuantumPythonK.txt", "w")
    fileV = open("QuantumPythonV.txt", "w")
    fileC = open("QuantumPythonC.txt", "w")
    for j in range(Nx):
        fileV.write("{0:.4f} {1:.4f}\n".format(x[j], V_x[j]))
    for i in range(t_max):
        time_step()
        center = x0 + t * p0 / m
        fileC.write("{0:.4f} {1:.4f}\n".format(center, 0))
        fileC.write("{0:.4f} {1:.4f}\n".format(center, 10))
        psi_xf = 4 * abs(get_psi_x())
        psi_kf = abs(get_psi_k())
        for j in range(Nx):
            fileX.write("{0:.4f} {1:.4f}\n".format(x[j], psi_xf[j]))
            fileK.write("{0:.4f} {1:.4f}\n".format(k[j], psi_kf[j]))
    fileX.close()
    fileK.close()
    fileV.close()
    fileC.close()

#writeToFile();exit()
############
# PLOTTING #
############

fig = pl.figure()

# plotting limits
xlim = (-100, 100)
klim = (-5, 5)

# top axes show the x-space data
ymin = 0
ymax = V0
ax1 = fig.add_subplot(211, xlim=xlim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_x_line, = ax1.plot([], [], c='r', label=r'$|\psi(x)|$')
V_x_line, = ax1.plot([], [], c='k', label=r'$V(x)$')
center_line = ax1.axvline(0, c='k', ls=':',
                          label = r"$x_0 + v_0t$")

title = ax1.set_title("")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')

# bottom axes show the k-space data
ymin = abs(get_psi_k()).min()
ymax = abs(get_psi_k()).max()
ax2 = fig.add_subplot(212, xlim=klim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_k_line, = ax2.plot([], [], c='r', label=r'$|\psi(k)|$')

p0_line1 = ax2.axvline(-p0 / hbar, c='k', ls=':', label=r'$\pm p_0$')
p0_line2 = ax2.axvline(p0 / hbar, c='k', ls=':')
mV_line = ax2.axvline(np.sqrt(2 * V0) / hbar, c='k', ls='--',
                      label=r'$\sqrt{2mV_0}$')
ax2.legend(prop=dict(size=12))
ax2.set_xlabel('$k$')
ax2.set_ylabel(r'$|\psi(k)|$')

V_x_line.set_data(x, V_x)

#############
# ANIMATION #
#############
frames = int(t_max / float(N_steps * dt))

def init():
    psi_x_line.set_data([], [])
    V_x_line.set_data([], [])
    center_line.set_data([], [])

    psi_k_line.set_data([], [])
    return (psi_x_line, V_x_line, center_line, psi_k_line, title)

def animate(i):
    time_step()
    psi_x_line.set_data(x, 4 * abs(get_psi_x()))
    V_x_line.set_data(x, V_x)
    center_line.set_data(2 * [x0 + t * p0 / m], [0, 1])
    psi_k_line.set_data(k, abs(get_psi_k()))
    return (psi_x_line, V_x_line, center_line, psi_k_line, title)

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frames, interval=30, blit=True)
pl.show()
