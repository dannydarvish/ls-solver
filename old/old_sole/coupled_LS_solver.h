#ifndef COUPLED_LS_SOLVER
#define COUPLED_LS_SOLVER
#include <vector>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

inline double sqr(double x){return x*x;}
class TransitionMatrixSolver
{
private:
    double kp, E, k_max_tol, ep_tol, dk_tol;
    int beta, n_bare, n_chan;
    vector<double> m_sigma, k_vector;
    vector<pair<double, double> > m_alpha;
    vector<complex<double> > t_matrix;
    vector<vector<complex<double> > >  (*g) (double);
    vector<vector<complex<double> > > (*v) (double, double);
    vector<vector<complex<double> > > V(double k, double kp);
    bool started_t_matrix;
    void construct_linear_system(Mat<complex<double> > & M,
        Col<complex<double> > & b, double & ep, double & k_max, double & dk);

public:
    TransitionMatrixSolver(int beta, double kp, vector<vector<complex<double> > > (*g)(double),
        vector<vector<complex<double> > > (*v)(double, double), vector<double> m_sigma,
        vector<pair<double, double> > m_alpha, double k_max_tol=1e-2, double ep_tol=1e-2,
        double dk_tol=1e-2):
            beta(beta), kp(kp), g(g), v(v), m_sigma(m_sigma), m_alpha(m_alpha), k_max_tol(k_max_tol),
            ep_tol(ep_tol), started_t_matrix(false), n_bare(m_sigma.size()), n_chan(m_alpha.size()),
            E(sqrt(sqr(m_alpha[beta].first)+sqr(kp))+sqrt(sqr(m_alpha[beta].second)+sqr(kp))) {}

    vector<complex<double> > & get_t_matrix();
    vector<double> & get_k_vector();
};
#endif
