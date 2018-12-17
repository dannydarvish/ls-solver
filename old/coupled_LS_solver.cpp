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

ostream & operator<< (ostream & out, vector<double> v)
{
    for (auto a : v)
        out << a << " ";
    cout << endl;   
    return out;
}

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
//                      Col<complex<double> > new_sol, vector<vector<double> > & old_om_vec,
//                      vector<vector<double> >& new_om_vec, double tol)
// {
//     double max_so_far = 0;
//     double comp;
//     for (int chan = 0; chan < n_chan; chan++)
//     {
//         for (int i = 0; i < old_om_vec[chan].size(); i++)
//         {
//             int count = 0;
//             for (int j = 0; j < new_om_vec[chan].size(); j++)
//                 if (fequals(old_om_vec[chan][i], new_om_vec[chan][j],tol))
//                 {
//                     int to_add_old = 0;
//                     int to_add_new = 0;
//                     for (int a = 0; a < chan; a++)
//                     {
//                         to_add_old += old_om_vec[a].size();
//                         to_add_new += new_om_vec[a].size();
//                     }
//                     comp = abs((real(old_sol[i + to_add_old] - new_sol[j + to_add_new])) / 
//                                 real(old_sol[i + to_add_old]));
//                     if (comp > max_so_far)
//                         max_so_far = comp;
//                     comp = abs((imag(old_sol[i + to_add_old] - new_sol[j + to_add_new])) / 
//                                 imag(old_sol[i + to_add_old]));
//                     if (comp > max_so_far)
//                         max_so_far = comp;
            
//                     count++;
//                 }
//             assert(count == 1);
//         }
//     }
//     return max_so_far;
// }
double TransitionMatrixSolver::max_frac_diff(Col<complex<double> > & old_sol,
                     Col<complex<double> > & new_sol, vector<vector<double> > & old_om_vec,
                     vector<vector<double> >& new_om_vec, int gamma, int i)
{
    double max_so_far = 0;
    double comp;
    for (int j = 0; j < old_om_vec[gamma].size(); j++)
    {
        int jp;
        assert((old_om_vec[gamma].size() - 1) % n_subint_steps == 0);
        if (double(j)/n_subint_steps < i)
            jp = j;
        else if (i <= double(j)/n_subint_steps && double(j)/n_subint_steps < i + 1)
            jp = 2*j - i * n_subint_steps;
        else if (double(j)/n_subint_steps >= i + 1)
            jp = j + n_subint_steps;
        else
            assert(false);
        int to_add_old = 0;
        int to_add_new = 0;
        for (int a = 0; a < gamma; a++)
        {
            to_add_old += old_om_vec[a].size();
            to_add_new += new_om_vec[a].size();
        }
        assert(fequals(old_om_vec[gamma][j], om_vec[gamma][jp], 1e-6*(om_vec[gamma][jp] -
                                                                      om_vec[gamma][jp-1])));
        comp = abs(real(old_sol[j + to_add_old] - new_sol[jp + to_add_new]) /
                        real(old_sol[j + to_add_old]));
        if (comp > max_so_far && !isinf(comp))
            max_so_far = comp;

        comp = abs(imag(old_sol[j + to_add_old] - new_sol[jp + to_add_new]) /
                        imag(old_sol[j + to_add_old]));
        if (comp > max_so_far && !isinf(comp))
            max_so_far = comp;
    }
    return max_so_far;
}

// Helper function for contruct_linear_system
complex<double> TransitionMatrixSolver::integrand_numerator(double om, double om_pp, int alpha, int gamma)
{
    double m1g, m2g, q1;
    m1g = m_alpha[gamma].first;
    m2g = m_alpha[gamma].second;
    q1 = sqr(om_pp) + sqr(m1g) - sqr(m2g);

    return -1*sqrt(sqr(1.0/(2*om_pp)*q1) - sqr(m1g)) / (4*pow(om_pp,3.0)) * 
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
    Col<complex<double> > & b, double & om_minus_om0_max, vector<vector<double> > & subint_pts,
     vector<double> & Delta_vec)
{
    assert(n_subint_steps > 1); // want this for trapezoidal rule to work and so
                                // we can preperly deal with the integrand at E

    // keep track of which index, for each channel, points to om_pp = E.
    vector<int> E_ind(n_chan, 0.0);

    for (int c = 0; c < n_chan; c++)
    {
        om_vec[c].clear();
        for (int i = 0; i < subint_pts[c].size() - 1; i++)
        {
            double dom = (subint_pts[c][i+1] - subint_pts[c][i])/n_subint_steps;
            for (int j = 0; j < n_subint_steps; j++)
            {   
                om_vec[c].push_back(subint_pts[c][i] + j * dom);
            }
        }

        om_vec[c].push_back(subint_pts[c].back());
        if (!bool(count(om_vec[c].begin(), om_vec[c].end(), E)))
            om_vec[c].push_back(E);
        sort(om_vec[c].begin(), om_vec[c].end());

        // determine the index of E in om_vec
        for (int i =0; i < om_vec[c].size(); i++)
            if (fequals(om_vec[c][i],E,1e-6*(om_vec[c][i+1]-om_vec[c][i])))
                E_ind[c] = i;
    }
  
    int dim = 0;
    for (int c = 0; c < n_chan; c++)
        dim += om_vec[c].size();
    M.resize(dim, dim);
    b.resize(dim);

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
                    
                    // BUG HERE
                    assert(false);
                    if (om_pp_ind == 0 || om_pp_ind == om_vec[gamma].size() - 1)
                        dom = 0.5 * (om_vec[gamma][om_pp_ind + 1] - om_vec[gamma][om_pp_ind]);
                    // dom will be changed later for om_pp = E +- Delta
                    else 
                        dom = 0.5 * (om_vec[gamma][om_pp_ind + 1] - om_vec[gamma][om_pp_ind - 1]);
                    
                    // Trapezoidal integration
                    if (abs(om_pp - E) > Delta && !fequals(om_pp,E-Delta,1e-6*dom) &&
                                                  !fequals(om_pp,E+Delta,1e-6*dom))
                    {
                        M(i,j) = delta(alpha,gamma) * delta(om,om_pp,1e-6*dom) - dom *
                                    integrand_numerator(om, om_pp, alpha, gamma) / (om_pp - E);
                    }
                    else if (abs(om_pp - E) < Delta && !fequals(om_pp,E-Delta,1e-6*dom) &&
                                                       !fequals(om_pp,E+Delta,1e-6*dom))
                    {
                        if (fequals(om_pp,E,1e-6*dom))
                        {
                            E_accounted_for_check += 1;

                            // M(i,j) = delta(alpha,gamma)*delta(om, om_pp, 1e-6*dom) - dom *
                            //             0.5*(integrand_numerator(om, om_pp+dom, alpha, gamma) -
                            //                  integrand_numerator(om, om_pp-dom, alpha, gamma)) +
                            //         complex<double>(1i) * M_PI * integrand_numerator(om, E, alpha, gamma);
                            M(i,j) = delta(alpha,gamma)*delta(om, om_pp, 1e-6*dom) +
                                    complex<double>(1i) * M_PI * integrand_numerator(om, E, alpha, gamma);
                        }
                        else if (om_pp_ind == E_ind[gamma] + 1)
                        {
                            M(i,j) = delta(alpha,gamma) * delta(om,om_pp,1e-6*dom) -
                                dom*(integrand_numerator(om,om_pp,alpha,gamma) - integrand_numerator(om,E,alpha,gamma)) / 
                                        (om_pp - E) -
                                0.5*((integrand_numerator(om,om_vec[gamma][om_pp_ind+1],alpha,gamma) - integrand_numerator(om,E,alpha,gamma)) / 
                                        (om_vec[gamma][om_pp_ind+1] - E));
                        }
                        else if (om_pp_ind == E_ind[gamma] - 1)
                        {
                            M(i,j) = delta(alpha,gamma) * delta(om,om_pp,1e-6*dom) -
                                dom*(integrand_numerator(om,om_pp,alpha,gamma) - integrand_numerator(om,E,alpha,gamma)) / 
                                        (om_pp - E) +
                                0.5*((integrand_numerator(om,om_vec[gamma][om_pp_ind-1],alpha,gamma) - integrand_numerator(om,E,alpha,gamma)) / 
                                        (om_vec[gamma][om_pp_ind-1] - E));
                        }
                        else
                        {
                            M(i,j) = delta(alpha,gamma) * delta(om,om_pp,1e-6*dom) - 
                                dom*(integrand_numerator(om,om_pp,alpha,gamma) - integrand_numerator(om,E,alpha,gamma)) / 
                                    (om_pp - E);
                        }
                    }
                    else if (fequals(om_pp,E - Delta,1e-6*dom))
                    {
                        delta_accounted_for++;
                        M(i,j) = delta(alpha,gamma)*delta(om,om_pp,1e-6*dom) - 
                                0.5*(om_vec[gamma][om_pp_ind]-om_vec[gamma][om_pp_ind-1]) * 
                                    integrand_numerator(om,om_pp,alpha,gamma) / (om_pp - E) -
                                0.5*(om_vec[gamma][om_pp_ind+1]-om_vec[gamma][om_pp_ind]) *
                                    (integrand_numerator(om,om_pp,alpha,gamma)-integrand_numerator(om,E,alpha,gamma)) /
                                        (om_pp - E);
                    }
                    else if (fequals(om_pp,E + Delta,1e-6*dom))
                    {
                        delta_accounted_for++;
                        M(i,j) = delta(alpha,gamma)*delta(om,om_pp,1e-6*dom) - 
                                0.5*(om_vec[gamma][om_pp_ind+1]-om_vec[gamma][om_pp_ind]) * 
                                    integrand_numerator(om,om_pp,alpha,gamma) / (om_pp - E) -
                                0.5*(om_vec[gamma][om_pp_ind]-om_vec[gamma][om_pp_ind-1]) *
                                    (integrand_numerator(om,om_pp,alpha,gamma)-integrand_numerator(om,E,alpha,gamma)) /
                                        (om_pp - E);
                    }
                    else
                        assert(false);
                    j++;
                }
                if (delta_accounted_for != 2)
                    cout << delta_accounted_for << endl << subint_pts[0] << endl;
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
        // cout << "First call of get_t_matrix. Solving for the transition matrix." << endl;
        started_t_mat = true;
    }
    vector<vector<double> > subint_pts(n_chan); 
    vector<double> Delta_vec(n_chan);
    for (int c = 0; c < n_chan; c++)
    {
        Delta_vec[c] = 1*(E-m_alpha[c].first-m_alpha[c].second);
        // subint_pts[c] = {m_alpha[c].first + m_alpha[c].second, E - Delta_vec[c], E + Delta_vec[c],
        //                  om_minus_om0_max + m_alpha[c].first + m_alpha[c].second};
        subint_pts[c] = {E - Delta_vec[c], E + Delta_vec[c],
                         om_minus_om0_max + m_alpha[c].first + m_alpha[c].second};
    }
    

    Mat<complex<double> > M;
    Col<complex<double> > b;
    Col<complex<double> > new_sol;
    Col<complex<double> > old_sol;
    vector<vector<double> > old_om_vec;
    double frac_diff;
    vector<vector<double> > pts_to_add(n_chan);
    vector<vector<double> > temp_subint_pts(n_chan);
    double pt_to_add;
    // cout << endl << endl;
    // cout << subint_pts[0] << endl;
    while(true)
    {
        for (int c = 0; c < n_chan ; c++)
            for (int i = 0; i < subint_pts[c].size() - 1; i++)
                assert(subint_pts[c][i+1] > subint_pts[c][i]);
        construct_linear_system(M, b, om_minus_om0_max, subint_pts, Delta_vec);
        cout << "Dimension: " << om_vec[0].size() << endl;
        old_sol = solve(M, b);
        old_om_vec = om_vec;

        // CAN BE EASILY OPTIMIZED
        assert(false);
        for (int gamma = 0; gamma < n_chan; gamma++)
            for (int i = 0; i < subint_pts[gamma].size() - 1; i++)
            {
                temp_subint_pts = subint_pts;
                pt_to_add = 0.5*(subint_pts[gamma][i+1] + subint_pts[gamma][i]);
                temp_subint_pts[gamma].push_back(pt_to_add);
                sort(temp_subint_pts[gamma].begin(), temp_subint_pts[gamma].end());
                construct_linear_system(M, b, om_minus_om0_max, temp_subint_pts, Delta_vec);
                new_sol = solve(M, b);

                //////////////////////////////////////////////////////////
                // figure out tolerance for comparing if two points in the
                // om_vec[gamma] are the same
                double min = 1e100;
                for (int j = 0; j < temp_subint_pts[gamma].size() - 1; j++)
                    if (temp_subint_pts[gamma][j+1] - temp_subint_pts[gamma][j] < min)
                        min = temp_subint_pts[gamma][j+1] - temp_subint_pts[gamma][j];
                //////////////////////////////////////////////////////////
                // frac_diff = max_frac_diff(old_sol, new_sol, old_om_vec, om_vec, 1e-2 * min / n_subint_steps);
                // cout << "Starting max frac diff" << endl;
                frac_diff = max_frac_diff(old_sol, new_sol, old_om_vec, om_vec, gamma, i);
                // cout << "max frac diff done" << endl;
                if (frac_diff > dom_tol)
                {  
                    // cout << pt_to_add << " ";
                    pts_to_add[gamma].push_back(pt_to_add);
                }
            }
        // cout << endl;
        int total_num_pts_to_add = 0;
        for (int gamma = 0; gamma < n_chan; gamma++)
        {
            total_num_pts_to_add += pts_to_add[gamma].size();
            for (auto pt : pts_to_add[gamma])
                subint_pts[gamma].push_back(pt);
            sort(subint_pts[gamma].begin(), subint_pts[gamma].end());
            pts_to_add[gamma].clear();
        }
        if (total_num_pts_to_add == 0)
            break;
        else
            cout << total_num_pts_to_add << " points added." << endl;
    }
    cout << endl << endl;
    // cout << subint_pts[0] << endl;
    t_mat = conv_to<vector<complex<double> > >::from(new_sol);
    // cout << endl << endl;
    return t_mat;
}

const vector<vector<double> > & TransitionMatrixSolver::get_om_vector()
{
    if (!started_t_mat)
        throw runtime_error("TransitionMatrixSolver._om_vec not assigned yet.");
    return om_vec;
}