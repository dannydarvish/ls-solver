#include "coupled_LS_solver.h"

#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using namespace std;
using namespace arma;

cdouble TransitionMatrixSolver::V(double k, double kp,
                                  int alpha, int beta)
{
    cdouble ret;
    vector<complex<double> > sum_me (n_bare);
    for (int i = 0; i < n_bare; i++)
    {
        sum_me[i] = conj(g[i][alpha](k)) * 1.0/(E - m_bare[i]) * g[i][beta](kp);
    }
    complex<double> sum = accumulate(sum_me.begin(), sum_me.end(), complex<double>(0.0));
    ret = sum + v[alpha][beta](k, kp);
    return ret;
}

cvec kernel_moment(double b, double msq, cdouble E_plus_iep)
{
    cvec moments(4);

    cdouble q1 = sqrt(4*msq-sqr(E_plus_iep));
    cdouble q2 = sqrt(msq + sqr(b));
    moments[0] = -0.5*1.0/q1 * (q1 * log(q2+b) +
                                    E_plus_iep*atan(E_plus_iep*b/(q1*q2)) +
                                    E_plus_iep*atan(2*b/q1));
    moments[1] = -0.25 * E_plus_iep * log(2.0*q2 - E_plus_iep) -
                                    0.5 * q2;
    moments[2] = 1.0/(8.0*q1) *
                (q1*((2*msq-sqr(E_plus_iep))*log(q2+b) -
                 2*b*(q2+E_plus_iep)) -
                 E_plus_iep*(sqr(E_plus_iep)-4*msq)*
                 (atan(E_plus_iep*b/(q1*q2)) + atan(2*b/q1)));
    moments[3] = 1.0/48.0 *
                     (-3.0*(pow(E_plus_iep,3.0)-4*msq*E_plus_iep)*log(2.0*q2-E_plus_iep) - 
                       6.0*(sqr(E_plus_iep)-4*msq)*q2 -
                       6.0*E_plus_iep*sqr(q2) - 8.0*pow(q2,3.0));
    return moments;
}

void TransitionMatrixSolver::construct_linear_system(Mat<cdouble> & M, Col<cdouble> & b)
{
    // For this implementation to work, we also need the two-particle channel
    // masses to be equal
    for (const auto & chan : m_alpha)
            assert(chan.first == chan.second);
    int dim = k_vec.size() * n_chan;
    M.resize(dim, dim);
    b.resize(dim);
    int i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
        for (int k_ind = 0; k_ind < k_vec.size(); k_ind++)
        {
            double k = k_vec[k_ind];
            int j = 0;
            for (int gamma = 0; gamma < n_chan; gamma++)
            {
                // dk''
                double dkpp = k_max/n_steps;
                double dkpp_inv = 1.0/dkpp;
                cvec weights(k_vec.size(), 0.0);
                cvec w_lower(4), w_upper(4), w(4);
                w_upper = kernel_moment(k_vec[0], sqr(m_alpha[gamma].first), E + cdouble(1i) * ep);
                double B = k_vec[0];
                for (int kpp_ind = 0; kpp_ind < k_vec.size() - 3; kpp_ind++)
                {
                    double A = B;
                    w_lower = w_upper;
                    B += dkpp;
                    double kpp_ind_d = kpp_ind;
                    // Stitch together a Gaussian quadrature using a singular kernel.
                    if (kpp_ind == k_vec.size() - 4)
                        B = dkpp * n_steps;
                    w_upper = kernel_moment(B, sqr(m_alpha[gamma].first), E + cdouble(1i) * ep);
                    for (int ii = 0; ii < 4; ii++)
                        w[ii] = (w_upper[ii] - w_lower[ii]) * pow(dkpp_inv, ii);
                    weights[kpp_ind] += (((kpp_ind_d+1.0)*(kpp_ind_d+2.0)*(kpp_ind_d+3.0)*w[0] -
                                          (11.0+kpp_ind_d*(12.0+kpp_ind_d*3.0))*w[1] +
                                          3.0*(kpp_ind_d+2.0)*w[2] -
                                          w[3]
                                         )/6.0
                                        );
                    weights[kpp_ind+1] += ((-kpp_ind_d*(kpp_ind_d+2.0)*(kpp_ind_d+3.0)*w[0] +
                                           (6.0+kpp_ind_d*(10.0+kpp_ind_d*3.0))*w[1] - 
                                           (3.0*kpp_ind_d+5.0)*w[2] +
                                           w[3])*0.5
                                          );
                    weights[kpp_ind+2] += ((kpp_ind_d*(kpp_ind_d+1.0)*(kpp_ind_d+3.0)*w[0] -
                                           (3.0+kpp_ind_d*(8.0+kpp_ind_d*3.0))*w[1] +
                                           (3.0*kpp_ind_d+4.0)*w[2] -
                                           w[3])*0.5
                                          );
                    weights[kpp_ind+3] += ((-kpp_ind_d*(kpp_ind_d+1.0)*(kpp_ind_d+2.0)*w[0] +
                                           (2.0+kpp_ind_d*(6.0+kpp_ind_d*3.0))*w[1] -
                                           3.0*(kpp_ind_d+1.0)*w[2] +
                                           w[3]
                                           )/6.0
                                          );
                }
                for (int kpp_ind = 0; kpp_ind < k_vec.size(); kpp_ind++)
                {
                    M(i,j) = delta(alpha, gamma) * delta(k, k_vec[kpp_ind], 1e-6 * dkpp) - 
                            sqr(k_vec[kpp_ind]) * V(k, k_vec[kpp_ind], alpha, gamma) *
                            weights[kpp_ind];
                    j++;
                }
            }
            i++;
        }
    i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
        for (double k : k_vec)
        {
            b(i) = V(k, kp, alpha, beta);
            i++;
        }
}

const cvec & TransitionMatrixSolver::get_t_matrix()
{
    if (!started_calc)
    {        
        cout << "Solving for transition matrix. Linear system has dimension " + 
                to_string(n_chan*k_vec.size()) << "." << endl;
        started_calc = true;
        if (k_max <= kp)
            throw runtime_error("k_max is not > kp. k_max = " + to_string(k_max));
        Mat<complex<double> > M;
        Col<complex<double> > b;
        construct_linear_system(M, b);
        t = conv_to<cvec>::from(solve(M, b));
    }
    return t;
}

const vector<double> & TransitionMatrixSolver::get_k_vec()
{   
    if (!started_calc)
        throw runtime_error("Can't return k_vec without first running get_t_matrix.");
    return k_vec;
}