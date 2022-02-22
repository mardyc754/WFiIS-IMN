#include <iostream>
#include <fstream>
#include <cmath>

int main()
{
    
    double y0 = 1;
    double lambda = -1;
    double dt[3] = {0.01, 0.1, 1.0};
    int n[3] = {5.0/dt[0], 5.0/dt[1], 5.0/dt[2]};

    // metoda Eulera
    std::ofstream euler("euler.dat");
    std::ofstream euler_err("euler_err.dat");
    for(int j=0; j<3; j++)
    {
        double y1 = y0;
        double t = 0;
        double delta;
        for(int i=0; i<=n[j]; i++)
        {
            delta = fabs(y1 - exp(lambda * t));
            euler << t << " " << y1 << std::endl;
            euler_err << t << " " << delta << std::endl;
            y1 = y1 + dt[j] * lambda * y1;  
            t += dt[j];
        }
        euler << std::endl << std::endl;
        euler_err << std::endl << std::endl;
    }
    euler.close();
    euler_err.close();

    // metoda RK2
    std::ofstream rk2("rk2.dat");
    std::ofstream rk2_err("rk2_err.dat");
    for(int j=0; j<3; j++)
    {
        double y1 = y0;
        double t = 0;
        double delta;
        double k1, k2;
        for(int i=0; i<=n[j]; i++)
        {
            k1 = lambda * y1;
            k2 = lambda * (y1 + dt[j] * k1);
            delta = fabs(y1 - exp(lambda * t));
            rk2 << t << " " << y1 << std::endl;
            rk2_err << t << " " << delta << std::endl;
            y1 = y1 + dt[j]/2 * (k1 + k2);  
            t += dt[j];
        }
        rk2 << std::endl << std::endl;
        rk2_err << std::endl << std::endl;
    }
    rk2.close();
    rk2_err.close();

    // metoda RK4
    std::ofstream rk4("rk4.dat");
    std::ofstream rk4_err("rk4_err.dat");
    for(int j=0; j<3; j++)
    {
        double y1 = y0;
        double t = 0;
        double delta;
        double k1, k2, k3, k4;
        for(int i=0; i<=n[j]; i++)
        {
            k1 = lambda * y1;
            k2 = lambda * (y1 + dt[j]/2 * k1);
            k3 = lambda * (y1 + dt[j]/2 * k2);
            k4 = lambda * (y1 + dt[j] * k3);
            delta = fabs(y1 - exp(lambda * t));
            rk4 << t << " " << y1 << std::endl;
            rk4_err << t << " " << delta << std::endl;
            y1 = y1 + dt[j]/6 * (k1 + 2*k2 + 2*k3 + k4);  
            t += dt[j];
        }
        rk4 << std::endl << std::endl;
        rk4_err << std::endl << std::endl;
    }
    rk4.close();
    rk4_err.close();

    // RRZ 2 rzedu - obwod RLC
    double delta_t = 1e-4;
    double R = 100;
    double L = 0.1;
    double C = 0.001;
    double w_0 = 1/sqrt(L * C);
    double w_V[4] = {w_0*0.5, w_0*0.8, w_0, w_0*1.2}; 
    double T_0 = 2*M_PI/w_0;
    int m = 4*T_0 / delta_t;

    std::ofstream rlc_Q("rlc_q.dat");
    std::ofstream rlc_I("rlc_i.dat");
    for(int j=0; j<4; j++)
    {
        double Q_n = 0;
        double I_n = 0;
        double t = 0;
        double k1_Q, k2_Q, k3_Q, k4_Q;
        double k1_I, k2_I, k3_I, k4_I;
        for(int i=0; i<=m; i++)
        {
            
            double V = 10 * sin(w_V[j] * t);
            double V_hdt = 10 * sin(w_V[j] * (t + delta_t * 0.5));
            double V_dt = 10 * sin(w_V[j] * (t + delta_t));
            
            rlc_Q << t << " " << Q_n << std::endl;
            rlc_I << t << " " << I_n << std::endl;

            k1_Q = I_n; 
            k1_I = V/L - 1/(L*C) * Q_n - R/L * I_n;
            k2_Q = I_n + delta_t/2 * k1_I;
            k2_I = V_hdt/L - 1/(L*C) * (Q_n + delta_t/2 * k1_Q) - R/L * (I_n + delta_t/2 * k1_I);
            k3_Q = I_n + delta_t/2 * k2_I;
            k3_I = V_hdt/L - 1/(L*C) * (Q_n + delta_t/2 * k2_Q) - R/L * (I_n + delta_t/2 * k2_I);
            k4_Q = I_n + delta_t * k3_I;
            k4_I = V_dt/L - 1/(L*C) * (Q_n + delta_t * k3_Q) - R/L * (I_n + delta_t * k3_I);

            Q_n += delta_t/6 * (k1_Q + 2*k2_Q + 2*k3_Q + k4_Q);
            I_n += delta_t/6 * (k1_I + 2*k2_I + 2*k3_I + k4_I);
            t += delta_t;
        }
        rlc_Q << std::endl << std::endl;
        rlc_I << std::endl << std::endl;
    }
    rlc_Q.close();
    rlc_I.close();
}