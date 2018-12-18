#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include "coupled_LS_solver.h"

using namespace std;

inline double k_from_om(double om, double m1, double m2)
{
    return sqrt(sqr(1.0/(2*om)*(sqr(om)+sqr(m1)-sqr(m2)))-sqr(m2));
}

double g_sigma_pi_pi = 1.6380;
double c_sigma_pi_pi = 1.0200; // fm
double G_pi_pi_pi_pi = 0.5560;
double d_pi_pi = 0.5140; //fm
double m_pi = 139.0;
vector<double> m_sigma = {700.0};
vector<pair<double,double> > m_alpha = {{m_pi,m_pi}};
int beta = 0;
double fm_times_MeV = 0.005067731;

vector<double> E_vec;

complex<double> g_00(double k)
{
    return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
}
// vector<vector<complex<double> > > g(double k)
// {
//     vector<vector<complex<double > > > ret {{g_00(k)}};
//     return ret;
// }
vector<vector<complex<double> (*)(double)> > g = {{g_00}};

complex<double> v_00(double k, double kp)
{
    return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}

// vector<vector<complex<double> > > v(double k, double kp)
// {
//     vector<vector<complex<double > > > ret({{v_00(k, kp)}});
//     return ret;
// }

vector<vector<complex<double> (*)(double, double)> > v = {{v_00}};


int main(int argc, const char * argv[])
{
    double E =stod(argv[1]);
    double kp = sqrt(sqr(E/2.0) - sqr(m_pi));
    double delta_pipi_re;
    double ep = 1e-1;
    double k_max = 100*kp;
    TransitionMatrixSolver tms(kp, beta, g, v, m_sigma, m_alpha, k_max, 1e-2, 10, 1.0);
    const vector<complex<double> > t = tms.get_t_matrix();
    const vector<vector<double> > k_vec = tms.get_k_vec();
    for (int i = 0; i < k_vec[0].size(); i++)
    {
        double k = k_vec[0][i];
        if (abs(k - kp) < 1e-6*(k_vec[0][1]-k_vec[0][0]))
        {
            delta_pipi_re = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
                                            sqrt(sqr(k)+sqr(m_pi))*t[i]));
            if (delta_pipi_re < 0)
                delta_pipi_re += M_PI;
            cout << "E: " << E << 
                ", k_max: " << k_max << ", ep : " << ep << ", t = " << t[i] <<
                    ", delta = " << delta_pipi_re << "." << endl;
        }
    }
    return 0;
}
