#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <armadillo>
#include <numeric>
#include <algorithm>
#include <assert.h>
#include "coupled_LS_solver.h"

// debugging
#include <iomanip>

using namespace std;
using namespace arma;

inline bool fequals(double a, double b, double tol){return abs(a-b)<tol;}
inline int delta(int i, int j){return int(i==j);}
inline double delta(double x, double y, double tol){return int(fequals(x,y,tol));}
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

// Calculate k as a function of omega
inline double k_from_om(double om, double m1, double m2)
{
    return sqrt(sqr(1.0/(2*om)*(sqr(om)+sqr(m1)-sqr(m2)))-sqr(m2));
}

// Helper function to compare two solutions
// double TransitionMatrixSolver::max_frac_diff(Col<complex<double> > old_sol,
//                      Col<complex<double> > new_sol,
//                      bool skip = false)
// {
//     double max_so_far = 0;
//     double comp;
//     // for (auto a : old_sol)
//     //     cout << a << endl;
//     // cout << endl << endl << endl;
//     // for (auto a : new_sol)
//     //     cout << a << endl;
//     for (int chan = 0; chan < n_chan; chan++)
//     {
//         for (int i = 0; i < old_sol.n_elem/n_chan; i++)
//         {
//             // If we need to sample at twice the rate (i.e. if dk has been halved)
//             if (skip)
//             {
//                 // assert(old_sol.n_elem - 1 == new_sol.n_elem/2);
//                 comp = abs(old_sol[i+chan*old_sol.n_elem/n_chan] -
//                             new_sol[2*i + chan*new_sol.n_elem/n_chan]);
//             }
//             else
//             {
//                 assert(old_sol.n_elem == new_sol.n_elem);
//                 comp = abs(old_sol[i+chan*old_sol.n_elem/n_chan] - 
//                             new_sol[i + chan*new_sol.n_elem/n_chan]);
//             }
//             if (comp > max_so_far)
//                 max_so_far = comp;
//         }
//         // // compare endpoints for each channel
//         // comp = abs(old_sol[chan*old_sol.n_elem/n_chan] -
//         //             new_sol[chan*new_sol.n_elem/n_chan]);
//         // if (comp > max_so_far)
//         //     max_so_far = comp;
//     }
//     // cout << "max: " << max_so_far << endl << "abs(old_sol).max(): " << abs(old_sol).max() << endl;
//     return max_so_far / abs(old_sol).max();
// }
double TransitionMatrixSolver::max_frac_diff(Col<complex<double> > old_sol,
                     Col<complex<double> > new_sol, vector<vector<double> > & old_om_vec,
                     vector<vector<double> >& new_om_vec, double tol)
{
    double max_so_far = 0;
    double comp;
    for (int chan = 0; chan < n_chan; chan++)
    {
        for (int i = 0; i < old_om_vec[chan].size(); i++)
        {
            int count = 0;
            for (int j = 0; j < new_om_vec[chan].size(); j++)
                if (fequals(old_om_vec[chan][i], new_om_vec[chan][j],tol))
                {
                    int to_add_old = 0;
                    int to_add_new = 0;
                    for (int a = 0; a < chan; a++)
                    {
                        to_add_old += old_om_vec[a].size();
                        to_add_new += new_om_vec[a].size();
                    }
                    comp = abs((old_sol[i + to_add_old] - new_sol[j + to_add_new]) / 
                                old_sol[i + to_add_old]);
                    if (comp > max_so_far)
                        max_so_far = comp;
            
                    count++;
                }
            assert(count == 1);
        }
    }
    return max_so_far;
}

double TransitionMatrixSolver::mean_frac_diff(Col<complex<double> > old_sol,
                     Col<complex<double> > new_sol,
                     bool skip = false)
{
    double sum = 0;
    for (int chan = 0; chan < n_chan; chan++)
    {
        for (int i = 0; i < old_sol.n_elem/n_chan; i++)
        {
            if (skip)
            {
                // assert(old_sol.n_elem - 1 == new_sol.n_elem/2);
                sum += abs((old_sol[i+chan*old_sol.n_elem/n_chan] -
                            new_sol[2*i + chan*new_sol.n_elem/n_chan])/
                                old_sol[i+chan*old_sol.n_elem/n_chan]);
                // cout << abs((old_sol[i+chan*old_sol.n_elem/n_chan] -
                //             new_sol[2*i + chan*new_sol.n_elem/n_chan])/
                //                 old_sol[i+chan*old_sol.n_elem/n_chan]) << endl;
            }
            else
            {
                assert(old_sol.n_elem == new_sol.n_elem);
                sum += abs((old_sol[i+chan*old_sol.n_elem/n_chan] - 
                            new_sol[i + chan*new_sol.n_elem/n_chan])/
                                old_sol[i+chan*old_sol.n_elem/n_chan]);
            }
        }
    }
    // cout << "sum: " << sum << endl << "n_elem: " << old_sol.n_elem << endl << endl;
    return (sum / old_sol.n_elem);
}

// Helper function for contruct_linear_system
complex<double> TransitionMatrixSolver::integrand_numerator(double om, double om_pp, int alpha, int gamma)
{
    double m1g, m2g, q1;
    m1g = m_alpha[gamma].first;
    m2g = m_alpha[gamma].second;
    q1 = sqr(om_pp) + sqr(m1g) - sqr(m2g);

    return sqrt(sqr(1.0/(2*om_pp)*q1) - sqr(m1g)) / (4*pow(om_pp,3.0)) * 
                (sqr(om_pp)+(sqr(m1g)-sqr(m2g))) * (sqr(om_pp)-(sqr(m1g)-sqr(m2g))) *
                V(k_from_om(om,m_alpha[alpha].first,m_alpha[alpha].second),
                    k_from_om(om_pp,m1g,m2g))[alpha][gamma];
}

// Helper function for get_t_matrix
//
// For a given d_omega and max value of omega - omega(0), define the matrix M and the vector b
// in the linear system Mv = b, where v is the vector we are solving for
// (the t matrix between omega' and all other omega in our discretization)
void TransitionMatrixSolver::construct_linear_system(Mat<complex<double> > & M,
    Col<complex<double> > & b, double & om_minus_om0_max, vector<double> & dom_vec, vector<double> & Delta_vec)
{
    cout << "Delta: " << Delta_vec[0] << ", d_om: " << dom_vec[0] << endl;
    for (int c = 0; c < n_chan; c++)
    {
        // Delta should be some integer multiple of d\omega
        assert(fmod(Delta_vec[c], dom_vec[c]) < 1e-6 * dom_vec[c]);
        assert(dom_vec[c] > 0);
        assert(om_minus_om0_max > dom_vec[c]);
    }
    for (int c = 0; c < n_chan; c++)
    {
        double dom = dom_vec[c];
        om_vec[c].clear();
        for (double om = m_alpha[c].first + m_alpha[c].second;
                om < om_minus_om0_max+m_alpha[c].first+m_alpha[c].second ||
                fequals(om,om_minus_om0_max+m_alpha[c].first+m_alpha[c].second,1e-6*dom);
                                                                    om += dom)
            om_vec[c].push_back(om);
    }
    cout << "om_vec[0][0] = " << om_vec[0][0] << ". om_vec[0][1] = " << om_vec[0][1] <<
            ". om_vec[0].back() = " << om_vec[0].back() << ". om_max = " <<
            om_minus_om0_max + m_alpha[0].first + m_alpha[0].second << endl;
    // sanity check
    // ************
    bool E_in_om_vec = false;
    for (int c = 0; c < n_chan; c++)
    {
        vector<double> & ov = om_vec[c];
        for (double om : ov)
            if (fequals(E,om,1e-6*dom_vec[c])) E_in_om_vec = true;
        assert(E_in_om_vec);
    }
    // ************
    for (int c = 0; c < n_chan; c++)
    {
        // assert(fequals(*max_element(om_vec[c].begin(),om_vec[c].end()),
        //                om_minus_om0_max + m_alpha[c].first + m_alpha[c].second, 1e-6*dom_vec[c]));
        assert(om_vec[0].size() == om_vec[c].size());
    }
    int dim = om_vec[0].size() * n_chan;
    M.resize(dim, dim);
    b.resize(dim);
    cout << "dimension: " << dim << " x " << dim << endl;
    int i = 0;
    // In 'flattening' our two-indexed kets (channel space and momentum space),
    // channel is the first index and momentum is the second index
    double q1, q2, m1g, m2g, d_omega, dom, Delta;
    double om_p = E;

    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (int om_ind = 0; om_ind < om_vec[alpha].size(); om_ind++)
        {
            double om = om_vec[alpha][om_ind];
            int j = 0;
            for (int gamma = 0; gamma < n_chan; gamma++)
            {
                Delta = Delta_vec[gamma];
                int delta_accounted_for = 0;
                // for (double om_pp : om_vec[gamma])
                bool E_accounted_for_check = 0;
                for (int om_pp_ind = 0; om_pp_ind < om_vec[gamma].size(); om_pp_ind++)
                {
                    double om_pp = om_vec[gamma][om_pp_ind];
                    if (om_pp_ind == 0 || om_pp_ind == om_vec[gamma].size() - 1)
                        dom = 0.5 * (om_vec[gamma][om_pp_ind + 1] - om_vec[gamma][om_pp_ind]);
                    // dom will be changed later for om_pp = E +- Delta
                    else 
                        dom = 0.5 * (om_vec[gamma][om_pp_ind + 1] - om_vec[gamma][om_pp_ind - 1]);
                    // Useful quantities
                    m1g = m_alpha[gamma].first;
                    m2g = m_alpha[gamma].second;
                    q1 = sqr(om_pp) + sqr(m1g) - sqr(m2g);
                    q2 = sqr(E) + sqr(m1g) - sqr(m2g);
                    double om_max = *max_element(om_vec[gamma].begin(),om_vec[gamma].end());

                    // Trapezoidal integration
                    if (abs(om_pp - E) > Delta && !fequals(om_pp,E-Delta,1e-6*dom) &&
                                                  !fequals(om_pp,E+Delta,1e-6*dom))
                    {
                        M(i,j) = delta(alpha,gamma) * delta(om,om_pp,1e-6*dom) - dom *
                                    integrand_numerator(om, om_pp, alpha, gamma) / (E - om_pp);
                    }
                    else if (abs(om_pp - E) < Delta && !fequals(om_pp,E-Delta,1e-6*dom) &&
                                                       !fequals(om_pp,E+Delta,1e-6*dom))
                    {
                        if (fequals(om_pp,E,1e-6*dom))
                        {
                            E_accounted_for_check += 1;

                            M(i,j) = delta(alpha,gamma)*delta(om, om_pp, 1e-6*dom) - dom *
                                        0.5*(integrand_numerator(om, om_pp+dom, alpha, gamma)+
                                             integrand_numerator(om, om_pp-dom, alpha, gamma)) +
                                    complex<double>(1i) * M_PI * integrand_numerator(om, E, alpha, gamma);
                        }
                        else
                        {
                            M(i,j) = delta(alpha,gamma) * delta(om,om_pp,1e-6*dom) - 
                                (integrand_numerator(om,om_pp,alpha,gamma) - integrand_numerator(om,E,alpha,gamma)) / 
                                    (E - om_pp);
                        }
                    }
                    else if (fequals(om_pp,E - Delta,1e-6*dom))
                    {
                        delta_accounted_for++;
                        M(i,j) = delta(alpha,gamma)*delta(om,om_pp,1e-6*dom) - 
                                0.5*(om_vec[gamma][om_pp_ind]-om_vec[gamma][om_pp_ind-1]) * 
                                    integrand_numerator(om,om_pp,alpha,gamma) / (E - om_pp) -
                                0.5*(om_vec[gamma][om_pp_ind+1]-om_vec[gamma][om_pp_ind]) *
                                    (integrand_numerator(om,om_pp,alpha,gamma)-integrand_numerator(om,E,alpha,gamma)) /
                                        (E - om_pp);
                    }
                    else if (fequals(om_pp,E + Delta,1e-6*dom))
                    {
                        delta_accounted_for++;
                        M(i,j) = delta(alpha,gamma)*delta(om,om_pp,1e-6*dom) - 
                                0.5*(om_vec[gamma][om_pp_ind+1]-om_vec[gamma][om_pp_ind]) * 
                                    integrand_numerator(om,om_pp,alpha,gamma) / (E - om_pp) -
                                0.5*(om_vec[gamma][om_pp_ind]-om_vec[gamma][om_pp_ind-1]) *
                                    (integrand_numerator(om,om_pp,alpha,gamma)-integrand_numerator(om,E,alpha,gamma)) /
                                        (E - om_pp);
                    }
                    else
                        assert(false);
                    j++;
                }
                assert(delta_accounted_for == 2);
                assert(E_accounted_for_check == 1);
            }   
            i++;
        }
    }
    i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (auto om : om_vec[alpha])
        {
            b(i) = V(k_from_om(om,m_alpha[alpha].first,m_alpha[alpha].second),
                     k_from_om(om_p,m_alpha[beta].first,m_alpha[beta].second))[alpha][beta];
            i++;
        }
    }
    // cout << "Det(M) = " << det(M) << endl;
}

// Returns a vector of transition matrix elements between the in-channel & in-momentum, and
// various out-channels & out-momenta. The order of the two indices in the 1-D vector is
// channel first, followed by momentum.
//
// When called for the first time, evaluates the entire t matrix (and so will take some time).
const vector<complex<double> > & TransitionMatrixSolver::get_t_matrix()
{
    if (!started_t_mat)
    {
        cout << "First call of get_t_matrix. Solving for the transition matrix." << endl;
        started_t_mat = true;
    }
    // om_minus_om0_max will be adaptively adjusted later.
    const int N = 4;
    double om_minus_om0_max = N * om_p;

    // We will start with adaptively adjusting dom.
    vector<double> dom_vec(n_chan);
    vector<double> Delta_vec(n_chan);
    const int N_SLICES_BETWEEN_E_AND_M1_PLUS_M2 = 5;
    const int DELTA_OVER_DOM = 2;
    for (int c = 0; c < n_chan; c++)
    {
        dom_vec[c] = (E - (m_alpha[c].first + m_alpha[c].second))/N_SLICES_BETWEEN_E_AND_M1_PLUS_M2;
        Delta_vec[c] = dom_vec[c] * DELTA_OVER_DOM;
        assert(E - Delta_vec[c] > m_alpha[c].first + m_alpha[c].second);
    }

    Mat<complex<double> > M;
    Col<complex<double> > b;
    construct_linear_system(M, b, om_minus_om0_max, dom_vec, Delta_vec);
    Col<complex<double> > new_sol = solve(M, b);
    Col<complex<double> > old_sol;
    Col<complex<double> > new_sol_every_other;
    vector<vector<double> > old_om_vec;

    int max_iter = 10;
    int n_iter = 1;
    double frac_diff;
    // Don't change this, our comparison function depends on the old solution being
    // sampled at twice the rate as the new solution.
    const int DOM_DIVISION = 2;
    do
    {
        cout << "dom iteration: " << n_iter << endl;
        for (double & dom : dom_vec)
            dom /= DOM_DIVISION;
        old_sol = new_sol;
        old_om_vec = om_vec;
        construct_linear_system(M, b, om_minus_om0_max, dom_vec, Delta_vec);
        // old_om_vec = om_vec;
        new_sol = solve(M, b);
        // assert(new_sol.n_elem == 2*old_sol.n_elem-1);
        frac_diff = max_frac_diff(old_sol, new_sol, old_om_vec, om_vec, 1e-6*dom_vec[0]);
        n_iter++;
        cout << "Fractional difference: " << frac_diff << endl;
    }
    while (frac_diff >= dom_tol && n_iter <= max_iter);
    // The last solution was good enough
    new_sol = old_sol;
    om_vec = old_om_vec;
    for (double & dom : dom_vec)
        dom *= DOM_DIVISION;

    if (n_iter > max_iter)
        cout << "Did not converge as dom -> 0 with max_iter = << " << max_iter << endl;

    // Adaptively find a sufficiently large om_max to approximate the integral as
    // om_max -> inf. Double until the fractional difference
    // between two successive solutions is below some tolerance.
    // n_iter = 1;
    // do
    // {
    //     cout << "om_max iteration: " << n_iter << endl;
    //     om_minus_om0_max *= 2.0;
    //     old_sol = new_sol;
    //     old_om_vec = om_vec;
    //     construct_linear_system(M, b, om_minus_om0_max, dom_vec, Delta_vec);
    //     new_sol = solve(M, b);
    //     // assert(new_sol.n_elem/n_chan == 2*old_sol.n_elem/n_chan-1);
    //     frac_diff = mean_frac_diff(old_sol(span(0,old_sol.n_elem/n_chan-1)),
    //                               new_sol(span(0,(new_sol.n_elem/n_chan-1)/2)));
    //     n_iter++;
    //     cout << "Fractional difference: " << frac_diff << endl;
    // }
    // while (frac_diff >= om_max_tol && n_iter <= max_iter);
    // // The last solution was good enough
    // new_sol = old_sol;
    // om_vec = old_om_vec;
    // om_minus_om0_max /= 2.0;
    
    // if (n_iter > max_iter)
    //     cout << "Did not converge as om_max -> inf with max_iter = " << max_iter << endl;

    t_mat = conv_to<vector<complex<double> > >::from(new_sol);
    return t_mat;
}
// vector<double> & TransitionMatrixSolver::get_k_vector()
// {
//     if (!started_t_mat)
//         throw runtime_error("TransitionMatrixSolver._k_vec not assigned yet.");
    
// }
const vector<vector<double> > & TransitionMatrixSolver::get_om_vector()
{
    if (!started_t_mat)
        throw runtime_error("TransitionMatrixSolver._om_vec not assigned yet.");
    return om_vec;
}