import numpy as np
from numpy import random
# We ultimately want to find the t-matrix that is a fixed point of a known operator, F.
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

n_chan = 1
def g_00(k):
    return g_sigma_pi_pi/np.sqrt(np.pi)/(1+(c_sigma_pi_pi*k*fm_times_MeV)**2)
def g(k):
    return np.array([[g_00(k)]])

def v_00(k, kp):
    return G_pi_pi_pi_pi/m_pi**2 * 1.0/(1+(d_pi_pi*k*fm_times_MeV)**2)**2 * \
        1.0/(1+(d_pi_pi*kp*fm_times_MeV)**2)**2
def v(k, kp):
    return np.array([[v_00(k, kp)]])

def V(k,kp):
    g_k = g(k)
    g_kp = g(kp)
    v_k_kp = v(k,kp)
    return np.array([[sum([np.asmatrix(g_k).H[alpha, i] * 1.0/(2*np.sqrt(m_pi**2+kp**2) - m_sigma[i]) *\
             g_kp[i, beta] for i in range(g_k.shape[0])]) +\
             v_k_kp[alpha, beta]
             for beta in range(n_chan)]
             for alpha in range(n_chan)])

def norm(t):
    return np.max(np.abs(t))

# pick a kp.
kp = np.sqrt(500.0**2 - m_pi**2)

k_max = 1e2 * kp
dk = kp * 1e-2

k_vec = np.arange(0,k_max,dk)
t1 = 1e-2*random.rand(int(k_max/dk))
t2 = 1e-2*random.rand(int(k_max/dk))

ep = 1e-4*kp
def F(t):
    # V_k_kp = np.array([V(k,kp) for k in k_vec])
    # return V_k_kp + np.array([sum([k_vec[j]**2*\
    #     (1.0-0.5*int(abs(k_vec[j]-k_vec[0])<dk*1e-3 or abs(k_vec[j]-k_vec[-1])<dk*1e-3))*V(k_vec[i],k_vec[j])*\
    #     1.0/(2*np.sqrt(m_pi**2+kp**2)-2*np.sqrt(m_pi**2+k_vec[j]**2)+1j*ep)*t[j] for j in range(len(k_vec))])\
    #     for i in range(len(k_vec))])

    F_t = [0 for i in t]
    for i in range(len(F_t)):
        rest = 0
        for j in range(len(F_t)):
            rest += k_vec[j]**2 * (1.0-0.5*int(abs(k_vec[j]-k_vec[0])<dk*1e-3 or abs(k_vec[j]-k_vec[-1])<dk*1e-3))*\
                V(k_vec[i],k_vec[j])[0][0]*\
                1.0/(2*np.sqrt(m_pi**2+kp**2)-2*np.sqrt(m_pi**2+k_vec[j]**2)+1j*ep)*\
                t[j]
        F_t[i] = V(k_vec[i],kp)[0,0] + rest
        print i
    return F_t

print F(t1)
