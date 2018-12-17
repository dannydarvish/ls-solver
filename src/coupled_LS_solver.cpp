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

void TransitionMatrixSolver::construct_linear_system(Mat<cdouble> & M, Col<cdouble> & b,
                                                double k_max, vector<vector<double> > & subint_pts)
{
    // The integral in question is along the real line and goes from k = 0 to k = \infty.
    // However, 
    // Need this for trapezoidal rule to work.
    assert(n_subint_steps > 1);
    // keep track of which index, for each channel, points to kpp = kp.
    vector<int> kp_ind(n_chan, 0.0);
    for (int c = 0; c < n_chan; c++)
    {   
        assert(count(subint_pts[c].begin(), subint_pts[c].end(), kp) == 1);
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
        for (int i = 0; i < k_vec[c].size(); i++)
            if (fequals(k_vec[c][i], kp, 1e-6*(k_vec[c][i+1]-k_vec[c][i])))
                kp_ind[c] = i;
    }
    int dim = 0;
    for (int c = 0; c < n_chan; c++)
        dim += k_vec[c].size();
    assert(dim <= max_dim);
    M.resize(dim, dim);
    b.resize(dim);

    double dk;
    int i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (int k_ind = 0; k_ind < k_vec[alpha].size(); k_ind++)
        {
            double k = k_vec[alpha][k_ind];
            int j = 0;
            for (int gamma = 0; gamma < n_chan; gamma++)
            {
                for (int kpp_ind = 0; kpp_ind < k_vec[gamma].size(); kpp_ind++)
                {
                    double kpp = k_vec[gamma][kpp_ind];
                    if (kpp_ind == 0 )
                        dk = 0.5 * (k_vec[gamma][kpp_ind + 1] - k_vec[gamma][kpp_ind]);
                    else if (kpp_ind == k_vec[gamma].size() - 1)
                        dk = 0.5 * (k_vec[gamma][kpp_ind] - k_vec[gamma][kpp_ind - 1]);
                    else
                        dk = 0.5 * (k_vec[gamma][kpp_ind + 1] - k_vec[gamma][kpp_ind - 1]);

                    M(i,j) = delta(alpha, gamma) * delta(k, kpp, 1e-6 * dk) - 
                            dk * sqr(kpp) * V(k, kpp, alpha, gamma) * 
                            1.0 / (E - sqrt(sqr(m_alpha[gamma].first) + sqr(kpp)) - 
                                        sqrt(sqr(m_alpha[gamma].second)+ sqr(kpp)) +
                                        cdouble(1i)*ep);
                    j++;
                }
            }
            i++;
        }
    }
    i = 0;
    for (int alpha = 0; alpha < n_chan; alpha++)
    {
        for (double k : k_vec[alpha])
        {
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
        assert(k_max > kp);
        for (int g = 0; g < n_chan; g++)
            subint_pts[g] = {0.0, kp, k_max};
        Mat<complex<double> > M;
        Col<complex<double> > b;
        Col<complex<double> > new_sol;
        Col<complex<double> > old_sol;
        vector<vector<double> > old_k_vec;
        double frac_diff;
        vector<vector<double> > pts_to_add(n_chan);
        vector<vector<double> > temp_subint_pts(n_chan);
        double pt_to_add;

        // an optimization where we ignore an interval that has already been checked.
        vector<vector<double> > subint_pts_to_skip(n_chan);
        // We will use adaptive refinement to determine how k_vec is discretized.
        while(true)
        {
            for (int c = 0; c < n_chan; c++)
                for (int i = 0; i < subint_pts[c].size() - 1; i++)
                    assert(subint_pts[c][i+1] > subint_pts[c][i]);
            construct_linear_system(M, b, k_max, subint_pts);
            
            int temp = 0;
            for (int c = 0; c < n_chan; c++)
                temp += k_vec[c].size();
            cout << "Dimension: " << temp << endl;

            old_sol = solve(M, b);
            old_k_vec = k_vec;
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

                    frac_diff = max_frac_diff(old_sol, new_sol, old_k_vec, k_vec, gam, i);
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