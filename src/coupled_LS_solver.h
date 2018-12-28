#ifndef COUPLED_LS_SOLVER
#define COUPLED_LS_SOLVER
#include <vector>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

typedef complex<double> cdouble;
typedef vector<complex<double> > cvec;

inline double sqr(double x){return x*x;}
inline bool fequals(double a, double b, double tol){return abs(a-b)<tol;}
inline int delta(int i, int j){return int(i==j);}
inline double delta(double x, double y, double tol){return int(fequals(x,y,tol));}

class TransitionMatrixSolver
{
private:
    double kp, k_max, E, dk_tol, ep;
    vector<int> kp_ind; // keep track of which index, for each channel, points to kp
    int beta, n_bare, n_chan, n_subint_steps;
    bool started_t_mat;
    int max_dim;
    vector<double> m_sigma;
    // We discretize k differently for each channel,
    // since the integrands are different for each 
    // channel, and so the refinements will in general
    // be different.
    vector<vector<double> > k_vec;
    vector<pair<double, double> > m_alpha;
    vector<cdouble> t;
    vector<vector<cdouble (*)(double)> > g;
    vector<vector<cdouble (*)(double, double)> > v;
    cdouble V(double k, double kp, int alpha, int beta);
    double max_frac_diff(Col<cdouble > & old_sol,
                     Col<cdouble > & new_sol, vector<vector<double> > & old_k_vec,
                     vector<vector<double> > & new_k_vec, int gamma, int i);
    double kp_frac_diff(Col<cdouble> & old_sol, Col<cdouble> & new_sol,
                        vector<vector<double> > & old_k_vec, vector<vector<double> > & new_k_vec,
                        int gamma, int i, int old_kp_ind, int new_kp_ind);
    void construct_linear_system(Mat<cdouble> & M, Col<cdouble> & b, double k_max,
                                 vector<vector<double> > & subint_pts);
public:
    TransitionMatrixSolver(double kp, int beta, vector<vector<cdouble (*)(double)> > g,
        vector<vector<cdouble (*)(double, double)> > v, vector<double> m_sigma,
        vector<pair<double, double> > m_alpha, double k_max, double dk_tol,
        int n_subint_steps, double ep, int max_dim = 10000):
            kp(kp), beta(beta), g(g), v(v), m_sigma(m_sigma), m_alpha(m_alpha),
            k_max(k_max), dk_tol(dk_tol), n_subint_steps(n_subint_steps), ep(ep),
            started_t_mat(false), n_chan(m_alpha.size()), n_bare(m_sigma.size()),
            max_dim(max_dim)
    {
        E = sqrt(sqr(m_alpha[beta].first) + sqr(kp)) + sqrt(sqr(m_alpha[beta].second) + sqr(kp));
        k_vec.resize(n_chan);
        kp_ind.resize(n_chan, -1);
    }
    const cvec & get_t_matrix();
    const vector<vector<double> > & get_k_vec();
    void set_kp(double kp)
    {
        started_t_mat = false;
        this->kp = kp;
    }
};
#endif