#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <armadillo>
#include <numeric>
#include <assert.h>
#include "coupled_LS_solver.h"

// debugging
#include <iomanip>
//

using namespace std;
using namespace arma;

inline int delta(int i, int j){return int(i==j);}
inline double delta(double x, double y){return int(x==y);}
template <typename T>
Col<T> every_other(Col<T> v);

// Set up the matrix of coupled channel potentials
vector<vector<complex<double> > > TransitionMatrixSolver::V(double k, double kp)
{
    vector<vector<complex<double> > > V(n_chan, vector<complex<double> >(n_chan));
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (int beta = 0; beta < n_chan; beta++)
        {
            vector<complex<double> > sum_me (n_bare);
            for (int i = 0; i < n_bare; i++)
            {
                sum_me[i] = conj(g(k)[i][alpha]) * 1.0/(E - m_sigma[i]) * g(kp)[i][beta];
            }
            complex<double> sum = accumulate(sum_me.begin(), sum_me.end(), complex<double>(0.0));
            V[alpha][beta] = sum + v(k,kp)[alpha][beta];
        }
    }
    return V;
}

// Helper function for get_t_matrix
//
// For a given epsilon and a given k_max, define the matrix M and the vector b
// in the linear system Mv = b, where v is the vector we are solving for
// (the t matrix between kp and all other k in our discretization)
void TransitionMatrixSolver::construct_linear_system(Mat<complex<double> > & M,
    Col<complex<double> > & b, double & ep, double & k_max, double & dk)
{
    k_vector.clear();
    for (double k = 0; k < k_max && abs(k-k_max) > 1e-6*dk; k += dk)
        k_vector.push_back(k);

    int dim = k_vector.size() * n_chan;
    M.resize(dim, dim);
    b.resize(dim);
    int i = 0;
    // In 'flattening' our two-indexed kets (channel space and momentum space),
    // channel is the first index and momentum is the second index
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (double k : k_vector)
        {
            int j = 0;
            for (int gamma = 0; gamma < n_chan; gamma++)
            {
                for (double kpp : k_vector)
                {
                    // The 0.5 makes this a trapozoidal integration.
                    M(i,j) = delta(alpha, gamma)*delta(k, kpp) -
                    (1.0 - 0.5*int(kpp==k_vector.front() || kpp==k_vector.back()))*
                    dk * sqr(kpp) * V(k,kpp)[alpha][gamma] *
                    1.0 / (E - sqrt(sqr(m_alpha[gamma].first)  + sqr(kpp))
                             - sqrt(sqr(m_alpha[gamma].second) + sqr(kpp))
                             + complex<double>(1i) * ep);
                    j++;
                }
            }
            i++;
        }
    }
    cout << endl;

    i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (auto k : k_vector)
        {
            b(i) = V(k, kp)[alpha][beta];
            i++;
        }
    }
}

// Returns a vector of transition matrix elements between the in-channel & in-momentum, and
// various out-channels & out-momenta. The order of the two indices in the 1-D vector is
// channel first, followed by momentum.
//
// When called for the first time, evaluates the entire t matrix (and so will take some time).
vector<complex<double> > & TransitionMatrixSolver::get_t_matrix()
{
    if (!started_t_matrix)
    {
        cout << "First call of get_t_matrix. Solving for the transition matrix." << endl;
        started_t_matrix = true;
    }
    // k_max will be adaptively adjusted later.
    // N = K_max/k_p, should be even. The idea is that we will only sample half of the domain.
    // Change this soon.
    int N = 4;
    double k_max = N * kp;

    // Epsilon will be adaptively adjusted later.
    double ep = 1e-3 * kp;

    // We will start with adaptively adjustking dk.
    double dk = kp / 10.0;

    Mat<complex<double> > M;
    Col<complex<double> > b;
    construct_linear_system(M, b, ep, k_max, dk);
    Col<complex<double> > new_sol = solve(M, b);
    Col<complex<double> > old_sol;
    Col<complex<double> > new_sol_every_other;
    vector<double> old_k_vec;

    int max_iter = 10;
    int n_iter = 1;
    do
    {
        cout << "dk iteration: " << n_iter << endl;
        dk /= 2.0;
        construct_linear_system(M, b, ep, k_max, dk);
        old_sol = new_sol;
        old_k_vec = k_vector;
        new_sol = solve(M, b);
        // In order to do our comparison, our k_vector for the new solution will be sampled
        // at twice the rate of the old solution. Therefore, we need to compare every other
        // point of new_sol to old_sol.
        new_sol_every_other = every_other(new_sol);
        n_iter++;
        assert(old_sol.n_elem == new_sol_every_other.n_elem);
        cout << "Fractional difference: " << max(abs((new_sol(span(0,new_sol_every_other.n_elem/n_chan/2 - 1)) - old_sol(span(0,old_sol.n_elem/n_chan/2 - 1))) /
            old_sol(span(0,old_sol.n_elem/n_chan/2 - 1)))) << endl;
        cout << "k_vector.size(): " << k_vector.size() << endl <<
            "dk/k_vector.back(): " << dk/k_vector.back() << endl << endl;

    }
    while (max(abs((new_sol(span(0,new_sol_every_other.n_elem/n_chan/2 - 1)) - old_sol(span(0,old_sol.n_elem/n_chan/2 - 1))) /
            old_sol(span(0,old_sol.n_elem/n_chan/2 - 1)))) >= dk_tol && n_iter <= max_iter);

    if (n_iter > max_iter)
        cout << "Did not converge as dk -> 0 with max_iter = << " << max_iter << endl;
    // Since the one before was good enough
    new_sol = old_sol;
    k_vector = old_k_vec;

    // Adaptively find a sufficiently large k_max to approximate the integral as
    // k_max -> inf, starting with k_max = kp. Double until the fractional difference
    // between two successive solutions is below some tolerance.
    n_iter = 1;
    do
    {
        cout << "k_max iteration: " << n_iter << endl;
        k_max *= 2.0;
        construct_linear_system(M, b, ep, k_max, dk);
        old_sol = new_sol;
        old_k_vec = k_vector;
        new_sol = solve(M, b);
        n_iter ++;
    }
    // The /4 accounts for the fact that k_max in new_sol is twice that as in
    // old_sol, and we can only compare the points they have in common.
    while (max(abs((new_sol(span(0,new_sol.n_elem/n_chan/4 - 1)) - old_sol(span(0,old_sol.n_elem/n_chan/2 - 1))) /
            old_sol(span(0,old_sol.n_elem/n_chan/2 - 1)))) >= k_max_tol && n_iter <= max_iter);

    if (n_iter > max_iter)
        cout << "Did not converge as k_max -> inf with max_iter = " << max_iter << endl;
    // Since the one before was good enough
    k_max /= 2.0;
    k_vector = old_k_vec;
    new_sol = old_sol;

    // Now take the limit as epsilon -> 0^+ in the same way we took the other limits
    n_iter = 1;
    do
    {
        cout << "epsilon iteration: " << n_iter << endl;
        ep /= 100.0;
        construct_linear_system(M, b, ep, k_max, dk);
        old_sol = new_sol;
        new_sol = solve(M, b);
        n_iter++;
    }
    while (max(abs((new_sol(span(0,new_sol.n_elem/n_chan/2 - 1)) - old_sol(span(0,old_sol.n_elem/n_chan/2 - 1))) /
            old_sol(span(0,old_sol.n_elem/n_chan/2 - 1)))) >= ep_tol && n_iter <= max_iter);
    if (n_iter > k_max)
        cout << "Did not converge as i\\epsilon -> 0 with max_iter: " << max_iter << endl;

    // use the old solution, since it was good enough.
    ep *= 100.0;
    t_matrix = conv_to<vector<complex<double> > >::from(old_sol);
    return t_matrix;
}
vector<double> & TransitionMatrixSolver::get_k_vector()
{
    if (!started_t_matrix)
        throw runtime_error("TransitionMatrixSolver._k_vector not assigned yet.");
    return k_vector;
}

template <typename T>
Col<T> every_other(Col<T> v)
{
    vector<T > v1 = conv_to<vector<T> >::from(v);
    assert(v1.size() % 2 == 0);
    vector<T> v2(v1.size()/2);
    for (int i = 0; i < v1.size()/2; i++)
        v2[i] = v1[2*i];
    Col<T> c(v2);
    return c;
}
