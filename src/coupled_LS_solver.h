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

// Notation here follows the Wu et. al paper
// Unless otherwise stated, units are in MeV and fm.
class TransitionMatrixSolver
{
public:
    TransitionMatrixSolver(double kp, int beta, vector<vector<cdouble (*)(double)> > g,
        vector<vector<cdouble (*)(double, double)> > v, const vector<double> & m_bare,
        const vector<pair<double, double> > & m_tp, double k_max, int n_steps,
        double ep):
            kp(kp), beta(beta), g(g), v(v), m_bare(m_bare), m_tp(m_tp),
            k_max(k_max), n_steps(n_steps), ep(ep),
            started_calc(false), n_chan(m_tp.size()), n_bare(m_bare.size())
    {
        if (n_steps < 4)
            throw runtime_error("n_steps must be at least 4 (and should probably much larger).");

        E = sqrt(sqr(m_tp[beta].first) + sqr(kp)) + sqrt(sqr(m_tp[beta].second) + sqr(kp));
        k_vec.resize(n_steps + 1);
        for (int i = 0; i < n_steps + 1; i++)
            k_vec[i] = k_max/n_steps * i;
        // In this branch, we only support two-particle channels with identical particles.
        // Integrating the quadrature weightm, w(k) = 1/(E-sqrt(m1^2+k^2)-sqrt(m2^2+k^2)+i\ep),
        // analytically is much more difficult without this constraint
        for (const auto & chan : m_tp)
            assert(chan.first == chan.second);
    }
    ~TransitionMatrixSolver() {}
    const cvec & get_t_matrix();
    const vector<double> & get_k_vec();
    void set_kp(double kp)
    {
        started_calc = false;
        this->kp = kp;
        E = sqrt(sqr(m_tp[beta].first) + sqr(kp)) + sqrt(sqr(m_tp[beta].second) + sqr(kp));
    }
private:
    // kp: k'
    // k_max: upper limit in the integral
    // E: total energy
    // ep: \epsilon in the integral
    // beta: the two-particle incoming channel number
    // n_bare: number of bare particles
    // n_chan: number of two-particle channels
    // n_steps: the number of integration steps
    // m_bare: vector of single particle masses
    // m_tp: vector of pairs of two-particle channel masses
    // k_vec: discretized integration domain
    // t: t_{\alpha,\beta,k,k^\prime}
    //    \beta and k^\prime are fixed. \alpha and k are arranged row-major
    // g: User supplied 1-to-2 coupling
    // v: User supplied 2-to-2 coupling
    // V: The total potential in the Lippmann-Schwinger eq.
    double kp, k_max, E, ep;
    int beta, n_bare, n_chan, n_steps;
    bool started_calc;
    vector<double> m_bare;
    vector<double> k_vec;
    vector<pair<double, double> > m_tp;
    vector<cdouble> t;
    vector<vector<cdouble (*)(double)> > g;
    vector<vector<cdouble (*)(double, double)> > v;
    cdouble V(double k, double kp, int alpha, int beta);
    // This constructs the linear system to solve, M t = b
    void construct_linear_system(Mat<cdouble> & M, Col<cdouble> & b);
};
#endif