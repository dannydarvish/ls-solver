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
    double k_max = 20000.0;

    double delta_pipi_re, delta_pipi_re_ep100, delta_pipi_re_ep75, delta_pipi_re_ep50;

    double ep = 37.5;
    TransitionMatrixSolver tms(kp, beta, g, v, m_sigma, m_alpha, k_max, 1e-3, 100, ep, 15000);
    vector<complex<double> > t = tms.get_t_matrix();
    vector<vector<double> > k_vec = tms.get_k_vec();
    for (int i = 0; i < k_vec[0].size(); i++)
    {
        double k = k_vec[0][i];
        if (abs(k - kp) < 1e-6*(k_vec[0][1]-k_vec[0][0]))
        {
            delta_pipi_re = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
                                            sqrt(sqr(k)+sqr(m_pi))*t[i]));
            if (delta_pipi_re < 0)
                delta_pipi_re += M_PI;
            // if (E > 2*m_pi)
            //     delta_pipi_re += M_PI;
            cout << "E: " << E << 
                ", k_max: " << k_max << ", ep : " << ep << ", t = " << t[i+k_vec[0].size()] <<
                    ", delta_pipi_re = " << delta_pipi_re << "." << endl;
        }
    }

    // double ep = 100;
    // TransitionMatrixSolver tms(kp, beta, g, v, m_sigma, m_alpha, k_max, 1e-3, 100, ep, 15000);
    // vector<complex<double> > t = tms.get_t_matrix();
    // vector<vector<double> > k_vec = tms.get_k_vec();
    // for (int i = 0; i < k_vec[0].size(); i++)
    // {
    //     double k = k_vec[0][i];
    //     if (abs(k - kp) < 1e-6*(k_vec[0][1]-k_vec[0][0]))
    //     {
    //         delta_pipi_re_ep100 = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
    //                                         sqrt(sqr(k)+sqr(m_pi))*t[i]));
    //         if (delta_pipi_re_ep100 < 0)
    //             delta_pipi_re_ep100 += M_PI;
    //         // if (E > 2*m_pi)
    //         //     delta_pipi_re += M_PI;
    //         cout << "E: " << E << 
    //             ", k_max: " << k_max << ", ep : " << ep << ", t = " << t[i+k_vec[0].size()] <<
    //                 ", delta_pipi_re_ep100 = " << delta_pipi_re_ep100 << "." << endl;
    //     }
    // }
    // ep = 75;
    // tms = TransitionMatrixSolver(kp, beta, g, v, m_sigma, m_alpha, k_max, 1e-3, 100, ep, 15000);
    // t = tms.get_t_matrix();
    // k_vec = tms.get_k_vec();
    // for (int i = 0; i < k_vec[0].size(); i++)
    // {
    //     double k = k_vec[0][i];
    //     if (abs(k - kp) < 1e-6*(k_vec[0][1]-k_vec[0][0]))
    //     {
    //         delta_pipi_re_ep75 = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
    //                                         sqrt(sqr(k)+sqr(m_pi))*t[i]));
    //         if (delta_pipi_re_ep75 < 0)
    //             delta_pipi_re_ep75 += M_PI;
    //         // if (E > 2*m_pi)
    //         //     delta_pipi_re += M_PI;
    //         cout << "E: " << E << 
    //             ", k_max: " << k_max << ", ep : " << ep << ", t = " << t[i+k_vec[0].size()] <<
    //                 ", delta_pipi_re_ep75 = " << delta_pipi_re_ep75 << "." << endl;
    //     }
    // }
    // ep = 50;
    // tms = TransitionMatrixSolver(kp, beta, g, v, m_sigma, m_alpha, k_max, 1e-3, 100, ep, 15000);
    // t = tms.get_t_matrix();
    // k_vec = tms.get_k_vec();
    // for (int i = 0; i < k_vec[0].size(); i++)
    // {
    //     double k = k_vec[0][i];
    //     if (abs(k - kp) < 1e-6*(k_vec[0][1]-k_vec[0][0]))
    //     {
    //         delta_pipi_re_ep50 = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
    //                                         sqrt(sqr(k)+sqr(m_pi))*t[i]));
    //         if (delta_pipi_re_ep50 < 0)
    //             delta_pipi_re_ep50 += M_PI;
    //         // if (E > 2*m_pi)
    //         //     delta_pipi_re += M_PI;
    //         cout << "E: " << E << 
    //             ", k_max: " << k_max << ", ep : " << ep << ", t = " << t[i+k_vec[0].size()] <<
    //                 ", delta_pipi_re_ep50 = " << delta_pipi_re_ep50 << "." << endl;
    //     }
    // }
    // cout << "E: " << E << 
    //             ", k_max: " << k_max << ", ep : " << ep <<
    //                 ", delta_pipi_re = " << 
    //                 delta_pipi_re_ep75 + 75.0 * (delta_pipi_re_ep50 - delta_pipi_re_ep100) / 50.0 + 
    //                 0.5 * (delta_pipi_re_ep50 + delta_pipi_re_ep100 - 2*delta_pipi_re_ep75) / sqr(25.0) * sqr(75.0)
    //                 << "." << endl;

    return 0;
}
