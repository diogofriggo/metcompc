import numpy as np
from matplotlib import pyplot as pl
from matplotlib import animation
from scipy.fftpack import fft,ifft


class Schrodinger(object):
    def __init__(self, x, psi_x0, V_x, k0 = None, hbar=1, m=1, t0=0.0):
        self.x, psi_x0, self.V_x = map(np.asarray, (x, psi_x0, V_x))
        N = self.x.size
        self.hbar = hbar
        self.m = m
        self.t = t0
        self.dt_ = None
        self.N = len(x)
        self.dx = self.x[1] - self.x[0]
        self.dk = 2 * np.pi / (self.N * self.dx)
        self.k0 = k0
        self.k = self.k0 + self.dk * np.arange(self.N)
        
        self.psi_x = psi_x0
        self.psi_mod_k = fft(self.psi_mod_x)
        
        # variables which hold steps in evolution of the
        self.x_evolve_half = None
        self.x_evolve = None
        self.k_evolve = None
        
        # attributes used for dynamic plotting
        self.psi_x_line = None
        self.psi_k_line = None
        self.V_x_line = None
    
    def _set_psi_x(self, psi_x):
        print("_set_psi_x")
        self.psi_mod_x = (psi_x * np.exp(-1j * self.k[0] * self.x)
                          * self.dx / np.sqrt(2 * np.pi))
    
    def _get_psi_x(self):
        return (self.psi_mod_x * np.exp(1j * self.k[0] * self.x)
                * np.sqrt(2 * np.pi) / self.dx)
    
    def _set_psi_k(self, psi_k):
        print("_set_psi_k")
        self.psi_mod_k = psi_k * np.exp(1j * self.x[0]
                                        * self.dk * np.arange(self.N))
    
    def _get_psi_k(self):
        return self.psi_mod_k * np.exp(-1j * self.x[0] *
                                       self.dk * np.arange(self.N))
    
    def _get_dt(self):
        return self.dt_
    
    def _set_dt(self, dt):
        if dt != self.dt_:
            self.dt_ = dt
            self.x_evolve_half = np.exp(-0.5 * 1j * self.V_x
                                         / self.hbar * dt )
            self.x_evolve = self.x_evolve_half * self.x_evolve_half
            self.k_evolve = np.exp(-0.5 * 1j * self.hbar /
                                    self.m * (self.k * self.k) * dt)
    
    psi_x = property(_get_psi_x, _set_psi_x)
    psi_k = property(_get_psi_k, _set_psi_k)
    dt = property(_get_dt, _set_dt)
    
    def time_step(self, dt, Nsteps = 1):
        self.dt = dt
        self.psi_mod_x *= self.x_evolve_half
        for i in range(Nsteps):
            self.psi_mod_k = fft(self.psi_mod_x)
            self.psi_mod_k *= self.k_evolve
            self.psi_mod_x = ifft(self.psi_mod_k)
            self.psi_mod_x *= self.x_evolve
        self.psi_mod_k = fft(self.psi_mod_x)
        self.t += dt * Nsteps

def gauss_x(x, a, x0, k0):
    return ((a * np.sqrt(np.pi)) ** (-0.5)
            * np.exp(-0.5 * ((x - x0) * 1. / a) ** 2 + 1j * x * k0))

def theta(x):
    y = np.zeros(x.shape)
    y[x > 0] = 1.0
    return y

def square_barrier(x, width, height):
    return height * (theta(x) - theta(x - width))

#definicoes de fora de Schrodinger
my_dt = 0.005
my_N_steps = 50
my_t_max = 120
my_frames = int(my_t_max / float(my_N_steps * my_dt))
my_hbar = 1.0
my_m = 1.0
my_N = 2 ** 11
my_dx = 0.1
my_x = my_dx * (np.arange(my_N) - 0.5 * my_N)
my_V0 = 1.5
my_L = my_hbar / np.sqrt(2 * my_m * my_V0)
my_a = 3 * my_L
my_x0 = -60 * my_L
#<this replaces the square barrier function>
my_V_x = my_V0 * (theta(my_x) - theta(my_x - my_a))
#</this replaces the square barrier function>
my_V_x[my_x < -98] = 1E6
my_V_x[my_x > 98] = 1E6
my_p0 = np.sqrt(2 * my_m * 0.2 * my_V0)
my_dp2 = my_p0 * my_p0 * 1./80
my_d = my_hbar / np.sqrt(2 * my_dp2)
my_k0 = my_p0 / my_hbar
my_v0 = my_p0 / my_m
#<this replaces the gauss function>
my_psi_x0 = ((my_d * np.sqrt(np.pi)) ** (-0.5)
             * np.exp(-0.5 * ((my_x - my_x0) * 1. / my_d) ** 2 + 1j * my_x * my_k0))
#</this replaces the gauss function>
my_t0 = 0.0
#definicoes de dentro de Schrodinger
my_sx, my_spsi_x0, my_sV_x = map(np.asarray, (my_x, my_psi_x0, my_V_x))
my_sN = my_sx.size
my_st = my_t0
my_sN = len(my_sx)
my_sdx = my_sx[1] - my_sx[0]
my_sdk = 2 * np.pi / (my_sN * my_sdx)
my_sk0 = my_k0
my_sk = my_sk0 + my_sdk * np.arange(my_sN)
my_spsi_x = my_spsi_x0
#<this in particular comes from a set method: _set_psi_x>
my_spsi_mod_x = (my_spsi_x * np.exp(-1j * my_sk[0] * my_sx)* my_sdx / np.sqrt(2 * np.pi))
#</this in particular comes from a set method: _set_psi_x>
my_spsi_mod_k = fft(my_spsi_mod_x)
my_sx_evolve_half = None
my_sx_evolve = None
my_sk_evolve = None
my_spsi_x_line = None
my_spsi_k_line = None
my_sV_x_line = None
#<this in particular comes from a set method: _set_dt>
my_sx_evolve_half = np.exp(-0.5 * 1j * my_sV_x / my_hbar * my_dt )
my_sx_evolve = my_sx_evolve_half * my_sx_evolve_half
my_sk_evolve = np.exp(-0.5 * 1j * my_hbar / my_m * (my_sk * my_sk) * my_dt)
#<this in particular comes from a set method: _set_dt>

def get_psi_x():
    return (my_spsi_mod_x * np.exp(1j * my_sk[0] * my_sx)
            * np.sqrt(2 * np.pi) / my_sdx)

def get_psi_k():
    return my_spsi_mod_k * np.exp(-1j * my_sx[0] *
                                  my_sdk * np.arange(my_sN))

def my_time_step():
    global my_spsi_mod_x, my_spsi_mod_k, my_st
    my_spsi_mod_x *= my_sx_evolve_half
    for i in range(my_N_steps):
        my_spsi_mod_k = fft(my_spsi_mod_x)
        my_spsi_mod_k *= my_sk_evolve
        my_spsi_mod_x = ifft(my_spsi_mod_k)
        my_spsi_mod_x *= my_sx_evolve
    my_spsi_mod_k = fft(my_spsi_mod_x)
    my_st += my_dt * my_N_steps

######################################################################
# Create the animation

# specify time steps and duration
dt = 0.005
N_steps = 50
t_max = 120
frames = int(t_max / float(N_steps * dt))

# specify constants
hbar = 1.0   # planck's constant
m = 1.0      # particle mass

# specify range in x coordinate
N = 2 ** 11
dx = 0.1
x = dx * (np.arange(N) - 0.5 * N)

# specify potential
V0 = 1.5
L = hbar / np.sqrt(2 * m * V0)
a = 3 * L
x0 = -60 * L
V_x = square_barrier(x, a, V0)
V_x[x < -98] = 1E6
V_x[x > 98] = 1E6

# specify initial momentum and quantities derived from it
p0 = np.sqrt(2 * m * 0.2 * V0)
dp2 = p0 * p0 * 1./80
d = hbar / np.sqrt(2 * dp2)

k0 = p0 / hbar
v0 = p0 / m
psi_x0 = gauss_x(x, d, x0, k0)

# define the Schrodinger object which performs the calculations
S = Schrodinger(x=x,
                psi_x0=psi_x0,
                V_x=V_x,
                hbar=hbar,
                m=m,
                k0=-28)
######################################################################
# Set up plot
fig = pl.figure()

# plotting limits
xlim = (-100, 100)
klim = (-5, 5)

# top axes show the x-space data
ymin = 0
ymax = my_V0#ymax = V0
ax1 = fig.add_subplot(211, xlim=xlim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_x_line, = ax1.plot([], [], c='r', label=r'$|\psi(x)|$')
V_x_line, = ax1.plot([], [], c='k', label=r'$V(x)$')
center_line = ax1.axvline(0, c='k', ls=':', label = r"$x_0 + v_0t$")

title = ax1.set_title("")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')

# bottom axes show the k-space data
#ymin = abs(S.psi_k).min()
#ymax = abs(S.psi_k).max()
ymin = abs(get_psi_k()).min()
ymax = abs(get_psi_k()).max()
ax2 = fig.add_subplot(212, xlim=klim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_k_line, = ax2.plot([], [], c='r', label=r'$|\psi(k)|$')

#p0_line1 = ax2.axvline(-p0 / hbar, c='k', ls=':', label=r'$\pm p_0$')
#p0_line2 = ax2.axvline(p0 / hbar, c='k', ls=':')
#mV_line = ax2.axvline(np.sqrt(2 * V0) / hbar, c='k', ls='--', label=r'$\sqrt{2mV_0}$')
p0_line1 = ax2.axvline(-my_p0 / my_hbar, c='k', ls=':', label=r'$\pm p_0$')
p0_line2 = ax2.axvline(my_p0 / my_hbar, c='k', ls=':')
mV_line = ax2.axvline(np.sqrt(2 * my_V0) / my_hbar, c='k', ls='--', label=r'$\sqrt{2mV_0}$')

ax2.legend(prop=dict(size=12))
ax2.set_xlabel('$k$')
ax2.set_ylabel(r'$|\psi(k)|$')

#V_x_line.set_data(S.x, S.V_x)
V_x_line.set_data(my_sx, my_sV_x)

######################################################################
# Animate plot
def init():
    psi_x_line.set_data([], [])
    V_x_line.set_data([], [])
    center_line.set_data([], [])
    
    psi_k_line.set_data([], [])
    return (psi_x_line, V_x_line, center_line, psi_k_line, title)

def animate(i):
    #    S.time_step(dt, N_steps)
    #    psi_x_line.set_data(S.x, 4 * abs(S.psi_x))
    #    V_x_line.set_data(S.x, S.V_x)
    #    center_line.set_data(2 * [x0 + S.t * p0 / m], [0, 1])
    #    psi_k_line.set_data(S.k, abs(S.psi_k))
    my_time_step()
    psi_x_line.set_data(my_sx, 4 * abs(get_psi_x()))
    V_x_line.set_data(my_sx, my_sV_x)
    center_line.set_data(2 * [my_x0 + my_st * my_p0 / my_m], [0, 1])
    psi_k_line.set_data(my_sk, abs(get_psi_k()))
    return (psi_x_line, V_x_line, center_line, psi_k_line, title)

# call the animator.  blit=True means only re-draw the parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=30, blit=True)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=my_frames, interval=30, blit=True)


# uncomment the following line to save the video in mp4 format.  This
# requires either mencoder or ffmpeg to be installed on your system

#anim.save('schrodinger_barrier.mp4', fps=15, extra_args=['-vcodec', 'libx264'])

pl.show()