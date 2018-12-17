#ifndef COUPLED_LS_SOLVER
#define COUPLED_LS_SOLVER
#include <vector>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

typedef complex<double> cdouble;
typedef vector<complex<double> > cvec;
typedef Mat<cdouble> cmat;

inline double sqr(double x){return x*x;}

class TransitionMatrixSolver
{
private:
    double om_p, E, om_max_tol, dom_tol, om_minus_om0_max;
    int beta, n_bare, n_chan, n_subint_steps;
    vector<double> m_sigma;
    vector<vector<double> > om_vec;
    vector<pair<double, double> > m_alpha;
    vector<complex<double> > t_mat;
    vector<vector<complex<double> > >  (*g) (double);
    vector<vector<complex<double> > > (*v) (double, double);
    vector<vector<complex<double> > > V(double k, double kp);
    double max_frac_diff(Col<complex<double> > & old_sol,
                     Col<complex<double> > & new_sol, vector<vector<double> > & old_om_vec,
                     vector<vector<double> >& new_om_vec, int gamma, int i);
    // double mean_frac_diff(Col<cdouble> old_sol, Col<cdouble> new_sol, bool skip);
    bool started_t_mat;
    // om refers to omega (energy).
    void construct_linear_system(Mat<complex<double> > & M,
    Col<complex<double> > & b, double & om_minus_om0_max, vector<vector<double> > & subint_pts,
     vector<double> & Delta_vec);
    complex<double> integrand_numerator(double om, double om_pp, int alpha, int gamma);

public:
    TransitionMatrixSolver(int beta, double om_p, vector<vector<complex<double> > > (*g)(double),
        vector<vector<complex<double> > > (*v)(double, double), vector<double> m_sigma,
        vector<pair<double, double> > m_alpha, double dom_tol,
        double om_minus_om0_max, int n_subint_steps):
            beta(beta), om_p(om_p), E(om_p), g(g), v(v), m_sigma(m_sigma), m_alpha(m_alpha),
            dom_tol(dom_tol), started_t_mat(false), n_bare(m_sigma.size()), n_chan(m_alpha.size()),
            om_minus_om0_max(om_minus_om0_max), n_subint_steps(n_subint_steps)
    {
        om_vec.resize(n_chan);
    }
    const vector<complex<double> > & get_t_matrix();
    // vector<double> & get_k_vector(int alpha);
    const vector<vector<double> > & get_om_vector();
    void set_om_p(double om_p)
    {
        started_t_mat = false;
        this->om_p = om_p;
    }
};
#endif
