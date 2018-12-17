import numpy as np
from math import sqrt
# for debugging only
import matplotlib.pyplot as plt

def delta(i,j):
    return int(i == j)

class TransitionMatrixSolver:
    '''
    TransitionMatrixSolver solves the coupled-channel Lippman-Schwinger
    equations to determine the transition matrix in the center-of-momentum frame.

    A TransitionMatrixSolver instance is defined for a given incoming
    channel and momentum, and solves for the matrix elements for each
    outgoing channel and momentum. These matrix elements are stored
    as a vector that is indexed by momentum.
    '''
    def __init__(self, beta, kp, g, v, m_sigma, m_alpha,k_max_tol=1e-3, ep_tol=1e-3):
        '''
        Create a TransitionMatrixSolver.

        Parameters
        ----------

        beta : int
            index of incoming channel
        kp : float
            incoming momentum
        g : function
            Function that returns an n-by-N array of 1-to-2 interactions
            (function of momentum), where n is the number of single-particle
            states and N is the number of two-particle decay channels.
            E.g. g(k)[i, alpha]
        v : function
            Function that returns an N by N (numpy) array of 2-to-2
            interactions (function of momentum) between the incoming channel
            beta and the outgoing channel alpha.
            E.g. v(k,kp)[alpha,beta].
        m_sigma : list of floats
            list of single-particle masses
        m_alpha : list of tuples of floats
            list of tuples containing the particle masses of the two-particle
            decay channels
        '''
        self._beta = beta
        self._kp = kp
        if callable(g):
            g_test = g(1.0)
        else:
            raise TypeError('g should be a function.')
        if isinstance(g_test, np.ndarray) and g_test.ndim == 2:
            self._g = g
        else:
            raise TypeError('g should be a function that returns a numpy \
array of dimension 2.')
        if callable(v):
            v_test = v(1.0,1.0)
        else:
            raise TypeError('v should be a function.')
        if isinstance(v_test, np.ndarray) and v_test.ndim == 2 and v_test.shape[0] == v_test.shape[1]:
            self._v = v
        else:
            raise TypeError('v is not a function that returns a square numpy \
array of dimension 2.')
        if g_test.shape[1] != v_test.shape[0]:
            raise TypeError('g should be n-by-N and v should be N-by-N.')
        del(g_test, v_test)

        self._m_sigma = m_sigma
        self._m_alpha = m_alpha
        self._E = sqrt(m_alpha[beta][0]**2 + kp**2) +\
            sqrt(m_alpha[beta][1]**2 + kp**2)
        self._n_bare = len(m_sigma)
        self._n_chan = len(m_alpha)
        self._k_vector = None
        self._t_matrix = None
        self._k_max_tol = k_max_tol
        self._ep_tol = ep_tol

    # Set up the matrix of coupled channel potentials
    def _V(self, k, kp):
        '''
        Returns a matrix of functions V(k,kp)[alpha,beta] that are the coupled-
        channel potentials.

        Parameters
        ---------

        k : float
            incoming momentum
        '''
        g_k = self._g(k)
        g_kp = self._g(kp)
        v = self._v(k,kp)
        # print self._m_alpha
        # print k, kp, v, np.array([[sum([np.asmatrix(g_k).H[alpha, i] * 1.0/(self._E - self._m_sigma[i]) *\
        #          g_kp[i, beta] for i in range(g_k.shape[0])]) +\
        #          v[alpha, beta]
        #          for beta in range(self._n_chan)]
        #          for alpha in range(self._n_chan)])
        # print self._E, self._m_sigma, g_k, g_kp
        # print '\n'
        return np.array([[sum([np.asmatrix(g_k).H[alpha, i] * 1.0/(self._E - self._m_sigma[i]) *\
                 g_kp[i, beta] for i in range(g_k.shape[0])]) +\
                 v[alpha, beta]
                 for beta in range(self._n_chan)]
                 for alpha in range(self._n_chan)])
        # return  self._g(k).H * 1.0/(self._E - self._m_sigma)*\
        #         self._g(kp) + self._v(k, kp)

    @property
    def t_matrix(self):
        '''
        Returns a vector of transition matrix elements between the in-channel & in-momentum, and
        various out-channels & out-momenta. The order of the two indices in the 1-D vector is
        channel first, followed by momentum.

        When called for the first time, evaluates the entire t matrix (and so will take some time).

        Parameters
        ---------

        None.
        '''
        if self._t_matrix is None:
            print 'First call of TransitionMatrixSolver.t_matrix. Solving for \
the transition matrix.'
            # We should also adaptively refine dk, but we don't right now because
            # we are lazy
            dk = self._kp / 100.0

            # For a given epsilon and a given k_max, define the matrix M and the vector b
            # in the linear system Mv = b, where v is the vector we are solving for
            # (the t matrix between kp and all other k in our discretization)
            def construct_linear_system(ep, k_max):
                self._k_vector = np.arange(0, k_max, dk)
                dim = len(self._k_vector) * self._n_chan
                M = np.zeros((dim, dim),dtype=complex)
                i = 0
                # In 'flattening' our two-indexed kets (channel space and momentum space),
                # channel is the first index and momentum is the second index
                for alpha in range(self._n_chan):
                    for k in self._k_vector:
                        j = 0
                        for gamma in range(self._n_chan):
                            for kpp in self._k_vector:
                                # The 0.5 makes this a trapozoidal integration.
                                M[i, j] = delta(alpha, gamma)*delta(k, kpp) -\
                                    (1.0-0.5*int(kpp==self._k_vector[0] or kpp==self._k_vector[-1]))*\
                                    dk * kpp**2 * self._V(k, kpp)[alpha, gamma] *\
                                    1.0 / (self._E - sqrt(self._m_alpha[gamma][0]**2 + kpp**2) \
                                                   - sqrt(self._m_alpha[gamma][1]**2 + kpp**2) \
                                                   + 1j * ep)
                                j += 1
                        i += 1
                b = np.zeros(dim)
                i = 0
                for alpha in range(self._n_chan):
                    for k in self._k_vector:
                        b[i] = self._V(k, self._kp)[alpha, self._beta]
                        i += 1
                return M, b

            # Adaptively find a sufficiently large k_max to approximate the integral as
            # k_max -> inf, starting with k_max = kp. Double until the fractional difference
            # between two successive solutions is below som tolerance.

            # N = K_max/k_p, should be even so we can easily sample t later in our comparison.
            N = 4
            assert(N % 2 == 0)
            k_max = N*self._kp
            # Epsilon will be adaptively adjusted later, but for now we hold it constant
            # while adjusting k_max. Is this okay to do?
            ep = 1e-3*self._kp
            max_iter = 10
            n_iter = 1
            M, b = construct_linear_system(ep, k_max)
            new = np.linalg.solve(M, b)
            old_k_vec = None
            # while n_iter == 1 or\
            #     (max(np.abs((new - old)/old)) >= k_max_tol and n_iter <= max_iter):
            while n_iter == 1 or\
                (max(np.abs((new[:len(new)/self._n_chan/4] -\
                    old[:len(old)/self._n_chan/2])/old[:len(old)/self._n_chan/2])) >=\
                    self._k_max_tol and n_iter <= max_iter):
                print 'k_max iteration: %d' % n_iter
                old = new
                old_k_vec = self._k_vector
                k_max *= 2.0
                M, b = construct_linear_system(ep, k_max)
                # make sure we sample the same points as before
                # new = np.linalg.solve(M, b)[[2*i for i in range(len(old))], 0]
                new = np.linalg.solve(M, b)
                n_iter += 1
            if n_iter > max_iter:
                print 'Did not converge as k_max -> inf with max_iter = %d.' % max_iter
            k_max /= 2.0
            self._k_vector = old_k_vec

            # Now take the limit as epsilon -> 0^+ in the same way we took the limit
            # in k_max.
            max_iter = 10
            n_iter = 1
            M, b = construct_linear_system(ep, k_max)
            new = np.linalg.solve(M, b)
            # plt.figure(2)
            while n_iter == 1 or\
                (max(np.abs((new[:len(new)/self._n_chan/2] -\
                    old[:len(old)/self._n_chan/2])/old[:len(old)/self._n_chan/2])) >=\
                    self._ep_tol and n_iter <= max_iter):
                # plt.plot(np.sqrt(self._k_vector**2+self._m_alpha[0][0]**2),new[:len(new)/self._n_chan],'.k')
                # plt.show()
                # plt.savefig('debug_plt%d.png'%n_iter,dpi=500)
                print 'epsilon teration: %d' % n_iter
                ep /= 100.0
                M, b = construct_linear_system(ep, k_max)
                old = new
                new = np.linalg.solve(M, b)
                n_iter += 1

            if n_iter > max_iter:
                print 'Did not converge as i\\epsilon -> 0 with max_iter = %d.' % max_iter

            self._t_matrix = new
        return self._t_matrix

    @property
    def k_vector(self):
        if self._k_vector is None:
            raise RuntimeError('TransitionMatrixSolver._k_vector not assigned\
                                yet.')
        else:
            return self._k_vector
