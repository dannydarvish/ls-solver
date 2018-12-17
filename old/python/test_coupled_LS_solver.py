import numpy as np
from math import pi
from coupled_LS_solver import TransitionMatrixSolver
import cmath
import matplotlib.pyplot as plt

# units: MeV
g_sigma_pi_pi = 2.0
c_sigma_pi_pi = 0.6722
G_pi_pi_pi_pi = 2.4998
d_pi_pi = 0.2440 #fm
g_sigma_k_k = 0.6451
c_sigma_k_k = 1.0398
G_k_k_k_k = 0.020
d_k_k = 0.10
G_pi_pi_k_k = 0.350
m_pi = 139.0
m_k = 497.0
m_sigma = [700.0]
m_alpha = [(m_pi,m_pi),(m_k,m_k)]
beta = 0
kp = 100.0


# Let's try calculating the transition matrix
# with incoming momentum ranging from 56 MeV to
# 740 MeV.
def g_00(k):
    return g_sigma_pi_pi/np.sqrt(np.pi)/(1+(c_sigma_pi_pi*k)**2)
def g_01(k):
    return g_sigma_k_k/np.sqrt(np.pi)/(1+(c_sigma_k_k*k)**2)
def g(k):
    return np.matrix([[g_00(k), g_01(k)]])

def v_00(k, kp):
    return G_pi_pi_pi_pi/m_pi**2 * 1.0/(1+(d_pi_pi*k)**2)**2 * \
        1.0/(1+(d_pi_pi*kp)**2)**2
def v_01(k, kp):
    return G_pi_pi_k_k/m_pi**2 * 1.0/(1+(d_pi_pi*k)**2)**2 * \
        1.0/(1+(d_k_k*kp)**2)**2
def v_10(k, kp):
    return G_pi_pi_k_k/m_pi**2 * 1.0/(1+(d_k_k*k)**2)**2 * \
        1.0/(1+(d_pi_pi*kp)**2)**2
def v_11(k, kp):
    return G_k_k_k_k/m_pi**2 * 1.0/(1+(d_k_k*k)**2)**2 * \
        1.0/(1+(d_k_k*kp)**2)**2
def v(k, kp):
    return np.matrix([[v_00(k, kp), v_01(k, kp)], [v_10(k, kp), v_11(k, kp)]])

tms = TransitionMatrixSolver(beta, kp, g, v, m_sigma, m_alpha)

t = tms.t_matrix

# just the pi-pi channel portion of the t-matrix
t_pipi = np.array(t.T)[0,:len(t)/2]
k = tms.k_vector
E_pipi = np.sqrt(k**2+m_pi**2)

delta_pipi = 1.0/2.0*np.angle(1.0 - 1j*pi*k*np.sqrt(k**2 + m_pi**2)*t_pipi,deg=False)

plt.plot(E_pipi, np.real(delta_pipi))
plt.savefig('test_plot.png')
