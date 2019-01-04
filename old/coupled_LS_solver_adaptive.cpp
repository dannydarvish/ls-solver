#include "coupled_LS_solver_adaptive.h"

#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using namespace std;
using namespace arma;

// To do: change assertions to exceptions

// debugging
///////////////////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <fstream>
ostream & operator<< (ostream & out, vector<double> v)
{
    for (auto a : v)
        out << a << endl;
    cout << endl;   
    return out;
}
//////////////////////////////////////////////////////////////////////////////////////////

cdouble TransitionMatrixSolver::V(double k, double kp,
                                  int alpha, int beta)
{
    cdouble ret;
    vector<complex<double> > sum_me (n_bare);
    for (int i = 0; i < n_bare; i++)
    {
        sum_me[i] = conj(g[i][alpha](k)) * 1.0/(E - m_sigma[i]) * g[i][beta](kp);
    }
    complex<double> sum = accumulate(sum_me.begin(), sum_me.end(), complex<double>(0.0));
    ret = sum + v[alpha][beta](k, kp);
    return ret;
}

double TransitionMatrixSolver::max_frac_diff(Col<cdouble > & old_sol,
                     Col<cdouble > & new_sol, vector<vector<double> > & old_k_vec,
                     vector<vector<double> > & new_k_vec, int gamma, int i)
{
    double max_so_far = 0;
    double comp;
    for (int j = 0; j < old_k_vec[gamma].size(); j++)
    {
        int jp;
        if (!((old_k_vec[gamma].size() - 1) % n_subint_steps == 0))
            cout << k_max << endl << endl << n_subint_steps << endl << endl <<  old_k_vec[gamma] << endl;

        assert((old_k_vec[gamma].size() - 1) % n_subint_steps == 0);
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
            to_add_old += old_k_vec[a].size();
            to_add_new += new_k_vec[a].size();
        }

        assert(fequals(old_k_vec[gamma][j], k_vec[gamma][jp], 1e-6*(k_vec[gamma][1] -
                                                                    k_vec[gamma][0])));
        comp = abs(real(old_sol[j + to_add_old] - new_sol[jp + to_add_new]) /
                        real(old_sol[j + to_add_old]));
        if (comp > max_so_far && !std::isinf(comp))
            max_so_far = comp;

        comp = abs(imag(old_sol[j + to_add_old] - new_sol[jp + to_add_new]) /
                        imag(old_sol[j + to_add_old]));
        if (comp > max_so_far && !std::isinf(comp))
            max_so_far = comp;
    }
    return max_so_far;
}

double TransitionMatrixSolver::kp_frac_diff(Col<cdouble> & old_sol, Col<cdouble> & new_sol,
                        vector<vector<double> > & old_k_vec, vector<vector<double> > & new_k_vec,
                        int gamma, int i, int old_kp_ind, int new_kp_ind)
{
    int j = old_kp_ind;
    int jp = new_kp_ind;
    if (!((old_k_vec[gamma].size() - 1) % n_subint_steps == 0))
        cout << k_max << endl << endl << n_subint_steps << endl << endl <<  old_k_vec[gamma] << endl;

    assert((old_k_vec[gamma].size() - 1) % n_subint_steps == 0);
    
    int to_add_old = 0;
    int to_add_new = 0;
    for (int a = 0; a < gamma; a++)
    {
        to_add_old += old_k_vec[a].size();
        to_add_new += new_k_vec[a].size();
    }
    if (fequals(old_k_vec[gamma][j],480.291,1e-5)) cout << old_k_vec[gamma][j] << " " << k_vec[gamma][jp] << endl;
    if (!fequals(old_k_vec[gamma][j], k_vec[gamma][jp], 1e-6*(k_vec[gamma][1] -
                                                                k_vec[gamma][0])))
        cout << old_k_vec[gamma][j] << " " << k_vec[gamma][jp] << endl;

    assert(fequals(old_k_vec[gamma][j], k_vec[gamma][jp], 1e-6*(k_vec[gamma][1] -
                                                                k_vec[gamma][0])));
    double re_comp= abs(real(old_sol[j + to_add_old] - new_sol[jp + to_add_new]) /
                    real(old_sol[j + to_add_old]));

    double im_comp = abs(imag(old_sol[j + to_add_old] - new_sol[jp + to_add_new]) /
                    imag(old_sol[j + to_add_old]));
    
    return (re_comp > im_comp) ? re_comp : im_comp;
}

cvec kermom(double b, double msq, cdouble E_plus_iep)
{
    cvec moments(4);

    cdouble q1 = sqrt(4*msq-sqr(E_plus_iep));
    cdouble q2 = sqrt(msq + sqr(b));
    assert(!isnan(real(q1)) && !isnan(imag(q1)));
    assert(!isnan(real(q2)) && !isnan(imag(q2)));
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

void TransitionMatrixSolver::construct_linear_system(Mat<cdouble> & M, Col<cdouble> & b,
                                                double k_max, vector<vector<double> > & subint_pts)
{
    // // For this implemented quadrature to work
    // assert(n_subint_steps == 3);
    assert(n_subint_steps % 3 == 0);

    // // for the uniform solver method
    // for (auto & pts : subint_pts)
    //     assert(pts.size() == 2);

    // For this implementation to work, we also need the two-particle channel
    // masses to be equal
    for (const auto & chan : m_alpha)
            assert(chan.first == chan.second);

    for (int c = 0; c < n_chan; c++)
    {   
        //DEBUG 
        #warning "kp assertion turned off"
        // assert(count(subint_pts[c].begin(), subint_pts[c].end(), kp) == 1);
        k_vec[c].clear();
        for (int i = 0; i < subint_pts[c].size() - 1; i++)
        {
            double dk = (subint_pts[c][i+1] - subint_pts[c][i])/n_subint_steps;
            for (int j = 0; j < n_subint_steps; j++)
            {
                k_vec[c].push_back(subint_pts[c][i] + j * dk);
            }
        }
        k_vec[c].push_back(subint_pts[c].back());
        sort(k_vec[c].begin(), k_vec[c].end());
        // determine the index of kp in k_vec
        #warning "kp_ind determination turned off"
        // for (int i = 0; i < k_vec[c].size(); i++)
        //     if (fequals(k_vec[c][i], kp, 1e-6*(k_vec[c][i+1]-k_vec[c][i])))
        //         kp_ind[c] = i;
        // assert(kp_ind[c] >= 0);
    }
    int dim = 0;
    for (int c = 0; c < n_chan; c++)
        dim += k_vec[c].size();
    if (dim > max_dim)
        throw runtime_error("Max dimension is " + to_string(max_dim) + ". Dimension is " + 
                                to_string(dim));
    
    M.zeros(dim, dim);
    b.zeros(dim);

    cdouble w0, w1, w2, w3;
    int i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (int k_ind = 0; k_ind < k_vec[alpha].size(); k_ind++)
        {
            double k = k_vec[alpha][k_ind];
            int j = 0;
            for (int gamma = 0; gamma < n_chan; gamma++)
            {
                double dk;
                for (int ii = 0; ii < subint_pts[gamma].size() - 1; ii++)
                    for (int jj = 0; jj < 4; jj++)
                    {
                        int kpp_ind = n_subint_steps * ii + jj;
                        if (jj == 0)
                        {
                            dk = (subint_pts[gamma][ii+1]-subint_pts[gamma][ii])/n_subint_steps;
                            cdouble E_plus_iep = E + cdouble(1i)*ep;
                            double A = subint_pts[gamma][i];
                            double B = subint_pts[gamma][i+1];
                            double msq = sqr(m_alpha[gamma].first);
                            cdouble q1 = sqrt(4*msq-sqr(E_plus_iep));
                            cdouble q2 = sqrt(msq + sqr(B));
                            cdouble q3 = sqrt(msq + sqr(A));
                            assert(!isnan(real(q1)) && !isnan(imag(q1)));
                            assert(!isnan(real(q2)) && !isnan(imag(q2)));
                            assert(!isnan(real(q3)) && !isnan(imag(q3)));
                            w0 = -0.5*1.0/q1 * (q1 * log(q2+B) +
                                    E_plus_iep*atan(E_plus_iep*B/(q1*q2)) +
                                    E_plus_iep*atan(2*B/q1))
                                    +
                                    0.5*1.0/q1 * (q1 * log(q3+A) +
                                    E_plus_iep*atan(E_plus_iep*A/(q1*q3)) +
                                    E_plus_iep*atan(2*A/q1));
                            assert(!isnan(real(w0)) && !isnan(imag(w0)));
                            w1 = 1.0/dk*(
                                -0.25 * E_plus_iep * log(2.0*q2 - E_plus_iep) -
                                    0.5 * q2
                                +0.25 * E_plus_iep * log(2.0*q3 - E_plus_iep) +
                                    0.5 * q3);
                            assert(!isnan(real(w1)) && !isnan(imag(w1)));
                            w2 = 1.0/pow(dk,2.0)*1.0/(8.0*q1) * (
                                    (q1*((2*msq-sqr(E_plus_iep))*log(q2+B) -
                                        2*B*(q2+E_plus_iep)) -
                                    E_plus_iep*(sqr(E_plus_iep)-4*msq)*
                                        (atan(E_plus_iep*B/(q1*q2)) + atan(2*B/q1)))
                                -
                                (q1*((2*msq-sqr(E_plus_iep))*log(q3+A) -
                                        2*A*(q3+E_plus_iep)) -
                                    E_plus_iep*(sqr(E_plus_iep)-4*msq)*
                                        (atan(E_plus_iep*A/(q1*q3)) + atan(2*A/q1))));
                            assert(!isnan(real(w2)) && !isnan(imag(w2)));
                            w3 = 1.0/pow(dk,3.0)*1.0/48.0 * (
                                (-3.0*(pow(E_plus_iep,3.0)-4*msq*E_plus_iep)*log(2.0*q2-E_plus_iep) - 
                                                6.0*(sqr(E_plus_iep)-4*msq)*q2 -
                                                6.0*E_plus_iep*sqr(q2) - 8.0*pow(q2,3.0))

                                                -
                                (-3.0*(pow(E_plus_iep,3.0)-4*msq*E_plus_iep)*log(2.0*q3-E_plus_iep) - 
                                                6.0*(sqr(E_plus_iep)-4*msq)*q3 -
                                                6.0*E_plus_iep*sqr(q3) - 8.0*pow(q3,3.0)));
                            assert(!isnan(real(w3)) && !isnan(imag(w3)));
                        }

                        cdouble weight;
                        switch(jj)
                        {
                            case 0:
                                weight = 1.0/6.0 * (
                                        double((kpp_ind+1)*(kpp_ind+2)*(kpp_ind+3))*w0 -
                                        double((3*sqr(kpp_ind)+12*kpp_ind+11))*w1 +
                                        double(3*(kpp_ind+2))*w2 -
                                        w3);
                                break;
                            case 1:
                                weight = 1.0/2.0 * (
                                        double(-kpp_ind*(kpp_ind+2)*(kpp_ind+3))*w0 +
                                        double(3*sqr(kpp_ind)+10*kpp_ind+6)*w1 -
                                        double(3*kpp_ind+5)*w2 +
                                        w3);
                                break;
                            case 2:
                                weight = 1.0/2.0 * (
                                        double(kpp_ind*(kpp_ind+1)*(kpp_ind+3))*w0 -
                                        double(3*sqr(kpp_ind)+8*kpp_ind+3)*w1 +
                                        double(3*kpp_ind+4)*w2 -
                                        w3);
                                break;
                            case 3:
                                weight = 1.0/6.0 * (
                                        double(-kpp_ind*(kpp_ind+1)*(kpp_ind+2))*w0 +
                                        double(3*sqr(kpp_ind)+6*kpp_ind+2)*w1 -
                                        double(3*(kpp_ind+1))*w2 +
                                        w3);
                                break;
                            default:
                                throw runtime_error("Error in assigning quadrature weights.");
                        }
                        double kpp = k_vec[gamma][kpp_ind];
                        if (kpp_ind == 0 || jj != 0 || kpp_ind == k_vec[gamma].size() - 1)
                            M(i,j) = delta(alpha, gamma) * delta(k, kpp, 1e-6 * dk) - 
                                    sqr(kpp) * V(k, kpp, alpha, gamma) *
                                    weight;
                        else
                            M(i,j) -= sqr(kpp) * V(k, kpp, alpha, gamma) *
                                    weight; 
                        // if (abs(kpp)>1e-9 && abs(real(M(i,j))) < 1e-9)
                        // {
                        //     cout << kpp << endl << M(i,j) << endl;
                        // }
                        if (jj != 3)
                            j++;
                //     }
                }
            }
            // cout << "incrementing i" << endl;
            i++;
        }
    }
    i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (double k : k_vec[alpha])
        {
            assert (i < dim);
            b(i) = V(k, kp, alpha, beta);
            i++;
        }
    }
}

const cvec & TransitionMatrixSolver::get_t_matrix()
{
    if (!started_t_mat)
    {        
        cout << "Solving for transition matrix." << endl;
        started_t_mat = true;
        vector<vector<double> > subint_pts(n_chan);
        if (k_max <= kp)
            throw runtime_error("k_max is not > kp. k_max = " + to_string(k_max));
        
        // Should I change this so instead of kp, it picks the k that gives \omega'' = E for each channel?
        // Currently it only does that for g = beta.
        // for (int g = 0; g < n_chan; g++)
        //      subint_pts[g] = {0.0, k_max};
        for (int g = 0; g < n_chan; g++)
        {   
            subint_pts[g].resize(300);
            for (int i = 0; i < 300; i++)
                subint_pts[g][i] = k_max/299.0 * i;
        }
        

        Mat<complex<double> > M;
        Col<complex<double> > b;
        Col<complex<double> > new_sol;
        Col<complex<double> > old_sol;
        vector<vector<double> > old_k_vec;
        vector<int> old_kp_ind(n_chan);
        double frac_diff;
        vector<vector<double> > pts_to_add(n_chan);
        vector<vector<double> > temp_subint_pts(n_chan);
        double pt_to_add;

        // an optimization where we ignore an interval that has already been checked.
        vector<vector<double> > subint_pts_to_skip(n_chan);
        // We will use adaptive refinement to determine how k_vec is discretized.

        // TO DO: Optimize such that we don't recompute interval endpoints every time.
        while(true)
        {
            // for (int c = 0; c < n_chan; c++)
            //     cout << subint_pts[c] << endl << endl;
            for (int c = 0; c < n_chan; c++)
                for (int i = 0; i < subint_pts[c].size() - 1; i++)
                    assert(subint_pts[c][i+1] > subint_pts[c][i]);
            construct_linear_system(M, b, k_max, subint_pts);

            int temp = 0;
            for (int c = 0; c < n_chan; c++)
                temp += k_vec[c].size();
            cout << "Dimension: " << temp << endl;

            old_sol = solve(M, b);
            
            // debugging
            t = conv_to<cvec>::from(old_sol);
            return t;
            ////////////
            old_k_vec = k_vec;
            old_kp_ind = kp_ind;
            int total_num_pts_to_add = 0;
            for (int gam = 0; gam < n_chan; gam++)
            {
                for (int i = 0; i < subint_pts[gam].size()-1; i++)
                {
                    if (bool(count(subint_pts_to_skip[gam].begin(),
                                   subint_pts_to_skip[gam].end(),
                                   subint_pts[gam][i])))
                        continue;
                        
                    temp_subint_pts = subint_pts;
                    pt_to_add = 0.5*(subint_pts[gam][i+1] + subint_pts[gam][i]);
                    temp_subint_pts[gam].push_back(pt_to_add);
                    sort(temp_subint_pts[gam].begin(), temp_subint_pts[gam].end());
                    construct_linear_system(M, b, k_max, temp_subint_pts);
                    new_sol = solve(M, b);
                    // cout << M << endl << endl << endl;
                    // frac_diff = max_frac_diff(old_sol, new_sol, old_k_vec, k_vec, gam, i);
                    frac_diff = kp_frac_diff(old_sol, new_sol, old_k_vec, k_vec, gam, i,
                                             old_kp_ind[gam], kp_ind[gam]);
                    if (frac_diff > dk_tol)
                        pts_to_add[gam].push_back(pt_to_add);
                    else
                        subint_pts_to_skip[gam].push_back(subint_pts[gam][i]);
                }
                total_num_pts_to_add += pts_to_add[gam].size();
                for (double pt : pts_to_add[gam])
                    subint_pts[gam].push_back(pt);
                sort(subint_pts[gam].begin(), subint_pts[gam].end());
                pts_to_add[gam].clear();
            }
            if (total_num_pts_to_add == 0)
                break;
        }
        k_vec = old_k_vec;
        t = conv_to<cvec>::from(old_sol);
    }
    return t;
}

const vector<vector<double> > & TransitionMatrixSolver::get_k_vec()
{
    return k_vec;
}