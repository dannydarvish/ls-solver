#ifndef COUPLED_LS_SOLVER
#define COUPLED_LS_SOLVER
#include <cmath>
#include <vector>
#include <complex>
#include <assert.h>

using namespace std;
inline double sqr(double in){return in*in;}
inline complex<double> sqr(complex<double> in){return in*in;}
typedef vector<vector<complex<double> > > cmat;
typedef vector<complex<double> > cvec;
typedef vector<double> dvec;

class TransitionMatrixSolver
{
private:
    double kp, E, ep, tol, dk;
    int beta, n_bare, n_chan, max_iter;
    dvec m_sigma, k_vec;
    vector<pair<double, double> > m_alpha;
    cmat t;
    bool get_t_called;
    cmat  (*g) (double);
    cmat (*v) (double, double);
    cmat V(double k, double kp);
    cmat F(cmat const & t);
    double t_norm(cvec v);
public:
    TransitionMatrixSolver(double kp, int beta, dvec m_sigma, vector<pair<double,double> > m_alpha,
                           double k_max, double dk, double ep, double tol, cmat (*g) (double),
                           cmat (*v) (double, double)) :
        kp(kp), E(E), beta(beta), m_sigma(m_sigma), m_alpha(m_alpha), k_vec(range(0.0,k_max,dk)),
        ep(ep), tol(tol), n_bare(m_sigma.size()), n_chan(m_alpha.size()), max_iter(max_iter),
        get_t_called(false), g(g), v(v) {}
    
    // Getter for the t-matrix. Also calculates the t-matrix the first time it's called.
    cvec & get_t();
    cvec & get_k();
};

// Helper functions
vector<int> range(int len)
{
    vector<int> ret(len);
    for (int i = 0; i < len; i++)
        ret[i] = i;
    return ret;
}
vector<int> range(int start, int end_m1)
{
    vector<int> ret(end_m1 - start);
    for (int i = start; i < end_m1; i++)
        ret[i] = i;
    return ret;
}
vector<double> range(double start, double end, double dt)
{
    int len = (end-start)/dt;
    vector<double> ret (len);
    for (int i : range(len))
    {
        ret[i] = dt * i;
    }
    return ret;
}

cmat operator- (cmat m1, cmat m2)
{
    assert(m1.size()==m2.size());
    for (int i = 0; i < m1.size(); i++)
        assert(m1[i].size() == m2[i].size());
    cmat m3(m1.size());
    for (int i = 0; i < m1.size(); i++)
        for (int j = 0; j < m1[i].size(); j++)
            m3[i][j] = m1[i][j] - m2[i][j];
    return m3;
}
#endif
