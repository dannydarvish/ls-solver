#include <iostream>
#include <fstream>
#include <cmath>
#include "coupled_LS_solver.h"

using namespace std;

inline double k_from_om(double om, double m1, double m2)
{
    return sqrt(sqr(1.0/(2*om)*(sqr(om)+sqr(m1)-sqr(m2)))-sqr(m2));
}

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
vector<double> E_vec = {800.0};

complex<double> g_00(double k)
{
    return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
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
    for (int mult = 10; mult < pow(10,10); mult*=10){
    for (double E : E_vec)
    {
        TransitionMatrixSolver tms(beta, E, g, v, m_sigma, m_alpha, 1e-2, .1);

        const vector<complex<double> > t = tms.get_t_matrix();
        
        vector<vector<double> > om_vec = tms.get_om_vector();
        for (int i = 0; i < om_vec[0].size(); i++)
        {
            double om = om_vec[0][i];
            if (abs(E - om) < 1e-6*(om_vec[0][1]-om_vec[0][0]))
            {
                double k = k_from_om(om,m_pi,m_pi);
                double delta_pipi_re = (0.5 * arg(1.0 - complex<double>(1i)*M_PI*k*
                                                sqrt(sqr(k)+sqr(m_pi))*t[i]));
                cout << "E: " << E << 
                    ", om_minus_om0_max: " << mult*E << ", delta = " << delta_pipi_re << "." << endl;
            }
        }
    }
    }
    return 0;
}
