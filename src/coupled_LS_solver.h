#ifndef COUPLED_LS_SOLVER
#define COUPLED_LS_SOLVER
#include <vector>
#include <armadillo>
#include <cmath>
#include <assert.h>

using namespace std;
using namespace arma;

typedef complex<double> cdouble;
typedef vector<complex<double> > cvec;

inline double sqr(double x){return x*x;}
inline cdouble sqr(cdouble x){return x*x;}
inline int sqr(int x){return x*x;}
inline bool fequals(double a, double b, double tol){return abs(a-b)<tol;}
inline int delta(int i, int j){return int(i==j);}
inline double delta(double x, double y, double tol){return int(fequals(x,y,tol));}

class TransitionMatrixSolver
{
private:
    double kp, k_max, E, ep;
    int beta, n_bare, n_chan, n_steps;
    bool started_t_mat;
    vector<double> m_sigma;
    vector<double> k_vec;
    vector<pair<double, double> > m_alpha;
    vector<cdouble> t;
    vector<vector<cdouble (*)(double)> > g;
    vector<vector<cdouble (*)(double, double)> > v;
    cdouble V(double k, double kp, int alpha, int beta);
    void construct_linear_system(Mat<cdouble> & M, Col<cdouble> & b);
public:
    TransitionMatrixSolver(double kp, int beta, vector<vector<cdouble (*)(double)> > g,
        vector<vector<cdouble (*)(double, double)> > v, vector<double> m_sigma,
        vector<pair<double, double> > m_alpha, double k_max, int n_steps,
        double ep):
            kp(kp), beta(beta), g(g), v(v), m_sigma(m_sigma), m_alpha(m_alpha),
            k_max(k_max), n_steps(n_steps), ep(ep),
            started_t_mat(false), n_chan(m_alpha.size()), n_bare(m_sigma.size())
    {
        if (n_steps < 4)
            throw runtime_error("n_steps must be at least 4 (and should probably much larger).");

        E = sqrt(sqr(m_alpha[beta].first) + sqr(kp)) + sqrt(sqr(m_alpha[beta].second) + sqr(kp));
        k_vec.resize(n_steps + 1);
        for (int i = 0; i < n_steps + 1; i++)
            k_vec[i] = k_max/n_steps * i;
        // In this branch, we only support two-particle channels with identical particles.
        // Integrating the quadrature weightm, w(k) = 1/(E-sqrt(m1^2+k^2)-sqrt(m2^2+k^2)+i\ep),
        // analytically is much more difficult without this constraint
        for (const auto & chan : m_alpha)
            assert(chan.first == chan.second);
    }
    const cvec & get_t_matrix();
    const vector<double> & get_k_vec();
    void set_kp(double kp)
    {
        started_t_mat = false;
        this->kp = kp;
        E = sqrt(sqr(m_alpha[beta].first) + sqr(kp)) + sqrt(sqr(m_alpha[beta].second) + sqr(kp));
    }
};
#endif