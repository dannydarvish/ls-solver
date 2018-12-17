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

double g_sigma_pi_pi = 2.0000;
double c_sigma_pi_pi = 0.6722; // fm
double G_pi_pi_pi_pi = 2.4998;
double d_pi_pi = 0.2440; //fm
double m_pi = 139.0;
double m_k = 497.0;
double g_sigma_k_kbar = 0.6451;
double c_sigma_k_kbar = 1.0398;
double G_k_kbar_k_kbar = 0.0200;
double d_k_kbar = 0.1000;
double G_pi_pi_k_kbar = 0.3500;

vector<double> m_sigma = {700.0};
vector<pair<double,double> > m_alpha = {{m_pi,m_pi},
                                        {m_k, m_k}};
int beta = 0;
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
    
    // for (double E = 300.0; E <= 1500.0; E += 100.0)
    //     E_vec.push_back(E);
    E_vec.push_back(stod(argv[1]));
    ofstream fout("delta_pipi_re_1b2c_E_" + string(argv[1]) + ".txt");
    for (double E : E_vec)
    {
        double kp = sqrt(sqr(E/2.0) - sqr(m_pi));
        double delta_pipi_re;
        double ep = 1e-1;
        double k_max = 4*kp;
        TransitionMatrixSolver tms(kp, beta, g, v, m_sigma, m_alpha, k_max, 1e-2, 100, 10.0);
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
                    delta_pipi_re += 2*M_PI;
                // delta_pipi_re = fmod(delta_pipi_re, M_PI);
                cout << "E: " << E << 
                    ", k_max: " << k_max << ", ep : " << ep << ", t = " << t[i] <<
                        ", delta = " << delta_pipi_re << "." << endl;
                fout << E << " " << delta_pipi_re * 180.0 / M_PI << endl;
            }
        }
    }
    return 0;
}
