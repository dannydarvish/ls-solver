import numpy as np
from math import pi
from coupled_LS_solver import TransitionMatrixSolver
import cmath
import matplotlib.pyplot as plt

# units: MeV
g_sigma_pi_pi = 1.6380
c_sigma_pi_pi = 1.0200
G_pi_pi_pi_pi = 0.5560
d_pi_pi = 0.5140 #fm
m_pi = 139.0
m_sigma = [700.0]
m_alpha = [(m_pi,m_pi)]
beta = 0
fm_times_MeV = 0.005067731
# fm_times_MeV = 1.0
f = plt.figure(1)
for E in 100.0*np.arange(3,6):
    kp = np.sqrt((E/2)**2 - m_pi**2)
    def g_00(k):
        return g_sigma_pi_pi/np.sqrt(np.pi)/(1+(c_sigma_pi_pi*k*fm_times_MeV)**2)
    def g(k):
        return np.array([[g_00(k)]])

    def v_00(k, kp):
        return G_pi_pi_pi_pi/m_pi**2 * 1.0/(1+(d_pi_pi*k*fm_times_MeV)**2)**2 * \
            1.0/(1+(d_pi_pi*kp*fm_times_MeV)**2)**2
    def v(k, kp):
        return np.array([[v_00(k, kp)]])

    tms = TransitionMatrixSolver(beta, kp, g, v, m_sigma, m_alpha)

    t = tms.t_matrix

    k = tms.k_vector
    E_pipi = np.sqrt(k[100]**2+m_pi**2)

    delta_pipi_re = np.real(1.0/2.0*np.angle(1.0 - 1j*pi*k[100]*np.sqrt(k[100]**2 + m_pi**2)*t[100],deg=False))
    print E, delta_pipi_re
    plt.figure(1)
    plt.plot(E_pipi, delta_pipi_re,'.k')
plt.figure(1)
plt.show()
