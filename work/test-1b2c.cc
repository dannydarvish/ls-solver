#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <fstream>
#include "coupled_LS_solver.h"

using namespace std;

inline double k_from_om(double om, double m1, double m2)
{
    return sqrt(sqr(1.0/(2*om)*(sqr(om)+sqr(m1)-sqr(m2)))-sqr(m2));
}

double g_sigma_pi_pi = 2.0000;
double c_sigma_pi_pi = 0.6722; // fm
double G_pi_pi_pi_pi = 2.4998;
double d_pi_pi = 0.2440; //fm
double m_pi = 134.976;
double m_k = 497.648;
double g_sigma_k_kbar = 0.6451;
double c_sigma_k_kbar = 1.0398;
double G_k_kbar_k_kbar = 0.0200;
double d_k_kbar = 0.1000;
double G_pi_pi_k_kbar = 0.3500;

vector<double> m_sigma = {700.0};
vector<pair<double,double> > m_alpha = {{m_pi,m_pi},
                                        {m_k, m_k}};
double fm_times_MeV = 0.005067731;

vector<double> E_vec;

complex<double> g_00(double k)
{
    return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
}
complex<double> g_01(double k)
{
    return g_sigma_k_kbar/sqrt(m_pi)/(1+sqr(c_sigma_k_kbar*k*fm_times_MeV));
}
vector<vector<complex<double> (*)(double)> > g = {{g_00, g_01}};

complex<double> v_00(double k, double kp)
{
    return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}
complex<double> v_01(double k, double kp)
{
    return G_pi_pi_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_k_kbar*kp*fm_times_MeV));
}
complex<double> v_10(double k, double kp)
{
    return G_pi_pi_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_k_kbar*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}
complex<double> v_11(double k, double kp)
{
    return G_k_kbar_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_k_kbar*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_k_kbar*kp*fm_times_MeV));
}
vector<vector<complex<double> (*)(double, double)> > v = {{v_00, v_01},
                                                          {v_10, v_11}};

int main(int argc, const char * argv[])
{
    
    double E =stod(argv[1]);
    int beta = 0;
    double kp = sqrt(sqr(E/2.0) - sqr(m_pi));
    double delta_pipi_re;
    double ep = 1e-2;
    double k_max = 20000.0;
    int n_steps = 2000;
    
    TransitionMatrixSolver tms(kp, beta, g, v, m_sigma, m_alpha, k_max, n_steps, ep);
    vector<complex<double> > t = tms.get_t_matrix();
    vector<double> k_vec = tms.get_k_vec();
    double min = 1e6;
    int argmin;
    double submin = 1e6;
    int argsubmin;
    for (int i = 0; i < k_vec.size(); i++)
    {
        double k = k_vec[i];
        if (abs(k-kp) < min)
        {
            submin = min;
            argsubmin = argmin;
            min = abs(k-kp);
            argmin = i;
        }
    }
    double k = k_vec[argmin];
    delta_pipi_re = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
                                            sqrt(sqr(k)+sqr(m_pi))*0.5*(t[argmin]+t[argsubmin])));
            if (delta_pipi_re < 0)
                delta_pipi_re += M_PI;
            cout << "E: " << E << 
                ", k_max: " << k_max << ", ep : " << ep << ", t = " << 0.5*(t[argmin]+t[argsubmin]) <<
                    ", delta = " << delta_pipi_re << "." << endl;
    cout << argmin << " " << k << " " << kp << " " << t[argmin] << " " << t[argmin-1] << " " << t[argmin+1] << endl; 
    return 0;
}
