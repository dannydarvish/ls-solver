#include <iostream>
#include <fstream>
#include <cmath>
#include "coupled_LS_solver.h"

using namespace std;

double g_sigma_pi_pi = 1.6380;
double c_sigma_pi_pi = 1.0200;
double G_pi_pi_pi_pi = 0.5560;
double d_pi_pi = 0.5140; //fm
double m_pi = 139.0;
vector<double> m_sigma = {700.0};
vector<pair<double,double> > m_alpha = {{m_pi,m_pi}};
int beta = 0;
double fm_times_MeV = 0.005067731;
// double fm_times_MeV = 1.0;
vector<double> E_vec = {300.0,400.0,500.0};

complex<double> g_00(double k)
{
    return g_sigma_pi_pi/sqrt(M_PI)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
}
vector<vector<complex<double> > > g(double k)
{
    vector<vector<complex<double > > > foo {{g_00(k)}};
    return foo;
}

complex<double> v_00(double k, double kp)
{
    return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}

vector<vector<complex<double> > > v(double k, double kp)
{
    vector<vector<complex<double > > > foo({{v_00(k, kp)}});
    return foo;
}

int main()
{
    for (double E : E_vec)
    {
        double kp = sqrt(sqr(E/2.0) - sqr(m_pi));

        TransitionMatrixSolver tms(beta, kp, g, v, m_sigma, m_alpha);

        vector<complex<double> > t = tms.get_t_matrix();
        vector<double> k = tms.get_k_vector();

        cout << "kp = " << kp << ", k[20] = " << k[20] << endl << endl;
        double E_pi = sqrt(sqr(k[20])+sqr(m_pi));
        double delta_pipi_re = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k[20]*
                                                sqrt(sqr(k[20])+sqr(m_pi))*t[20]));
        cout << "E: " << E << ", delta = " << delta_pipi_re << "." << endl;
    }
    return 0;
}
