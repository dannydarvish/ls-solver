#include "coupled_LS_solver.h"
#include <iostream>
#include <random>
#include <algorithm>
using namespace std;

// Define the norm as an L^2-like norm, where we sum over the L^2 norms of the different channels
double TransitionMatrixSolver::t_norm(cmat t)
{
    vector<double> channel_norms(n_chan);
    for (int c : range(n_chan))
    {
        for (int i : range(k_vec.size()))
        {
            // No idea why std::complex calls the magnitude-squared the norm, but
            // norm(z) = z^* z.
            channel_norms[c] += norm(t[c][i]);
        }
        channel_norms[c] = sqrt(channel_norms[c]);
    }
    return accumulate(channel_norms.begin(),channel_norms.end(),0.0);
}

cmat V(double k, double kp)
{
    cmat g_k = g(k);
    cmat g_kp = g(kp);
    cmat v_k_kp = v(k,kp);
    cmat V(n_chan, cvec(n_chan));
    for (int alpha : range(n_chan))
    {
        for (int beta : range(n_chan))
        {
            cvec sum_me (n_bare);
            for (int i : range(n_bare))
            {
                // This is totally wrong. E is fixed, and "kp" doesn't necessarily mean
                // kp as we globally define it.
                // sum_me[i] = conj(g_k[i][alpha]) * 1.0/(sqrt(sqr(m_pi)+sqr(kp)) - m_sigma[i]) *
                //                                                                     g_kp[i][beta];
            }
            complex<double> sum = accumulate(sum_me.begin(), sum_me.end(), complex<double>(0.0));
            V[alpha][beta] = sum + v_k_kp[alpha][beta];
        }
    }
    return V;
}

cmat TransitionMatrixSolver::F(cmat const & t)
{
    cmat F_t(t.size());
    complex<double> rest;
    // Trapezoidal integration. To do: improve.
    for (int a : range(n_chan))
    {
        for (int i : range(k_vec.size()))
        {
            rest = 0;
            for (int g : range(k_vec.size()))
            {
                for (j : range(k_vec.size()))
                {
                    rest += dk * (1.0-0.5*int(j == 0 || j == k_vec.size()-1) )*
                    sqr(k_vec[j]) * V(k_vec[i],k_vec[j])[a][g]*
                    1.0/(sqrt(sqr(m_sigma[beta].first())+sqr(kp))+sqrt(sqr(m_sigma[beta].second()))-
                    sqrt(sqr(m_sigma[g].first())+sqr(k_vec[j]))-sqrt(sqr(m_sigma[g].second())+sqr(k_vec[j]))
                    +complex<double>(1i)*ep)*t[g][j];
                }
            }
            F_t[a][i] = V(k_vec[i],kp)[a][i] + rest;
        }
    }
}

cvec & TransitionMatrixSolver::get_t()
{
    if (!get_t_called)
    {
        cout << "First call of TransitionMatrixSolver::get_t. Calculating the t." << endl << endl;
        get_t_called = true;
    }

    // Pick a random t
    mt19937 rng(time(0));
    rng.seed(random_device()());
    uniform_int_distribution<mt19937::result_type> dist(0,1e6);
    t.resize(n_chan);
    for (int i : range(n_chan))
    {
        t[i].resize(k_vec.size());
        for (int j : range(k_vec.size()))
            t[i][j] = 100.0*complex<double>(double(dist(rng))/1e6);
    }

    cmat t_prev, t_prev_prev;
    double frac_diff;
    int n_iter = -1;
    do
    {   
        if (n_iter > 0)
            cout << "Iteration: " << n_iter << end;
        t_prev_prev = t_prev;
        t_prev = t;
        t = F(t);
        double t_prev_dist = t_norm(t_prev-t_prev_prev);
        double t_dist = t_norm(t - t_prev);
        frac_diff = abs((t_dist-t_prev_dist)/t_dist);
        if (n_iter > 0)
            cout  << "Fractional difference: " << frac_diff << endl << endl;
    }
    while ((frac_diff > tol && n_iter <= max_iter) || n_iter < 1);
}