#include <iostream>
#include <fstream>
#include <cmath>

double fun(double beta, double N, double gamma, double u)
{
    return (beta*N - gamma)*u - beta*u*u;
}


int main()
{
    double beta = 0.001;
    double N = 500;
    double gamma = 0.1;
    double t_max = 100;
    double dt = 0.1;
    double TOL = 1e-6;
    int mu_max = 20;
    double u_0 = 1;
    double alpha = beta * N - gamma;

    // iteracja Picarda
    std::ofstream picard("picard.dat");

    for(double t = 0; t < t_max; t += dt)
    {
        double u_n = u_0;
        picard << t << " " << u_n << " " << N - u_n << std::endl;
        for(int mu=0; mu<=mu_max; mu++)
        {
            double temp = u_n;
            u_n = u_0 + dt/2 * ( (alpha * u_0 - beta*u_0*u_0) + (alpha*u_n - beta*u_n*u_n));
            if(fabs(u_n - temp) < TOL){ break; }
        }
        u_0 = u_n;
    }
    picard.close();


    u_0 = 1;
    // iteracja Newtona
    std::ofstream newton("newton.dat");

    for(double t = 0; t < t_max; t += dt)
    {
        double u_n = u_0;
        newton << t << " " << u_n << " " << N - u_n << std::endl;
        for(int mu=0; mu<=mu_max; mu++)
        {
            double temp = u_n;
            double numerator = u_n - u_0 - dt/2 * ((alpha*u_0 - beta*u_0*u_0) + (alpha*u_n - beta*u_n*u_n));
            double denominator = 1 - dt/2 * (alpha - 2*beta*u_n);
            u_n -= numerator / denominator;
            if(fabs(u_n - temp) < TOL){ break; }
        }
        u_0 = u_n;
    }
    newton.close();


    u_0 = 1;
    
    // metoda nrk2
    std::ofstream nrk2("nrk2.dat");

    double a_1_1 = 0.25;
    double a_1_2 = 0.25 - sqrt(3)/6;
    double a_2_1 = 0.25 + sqrt(3)/6;
    double a_2_2 = 0.25;
    double b_1 = 0.5;
    double b_2 = 0.5;
    double c_1 = 0.5 - sqrt(3)/6;
    double c_2 = 0.5 + sqrt(3)/6;

    for(double t = 0; t < t_max; t += dt)
    {
        double u_n = u_0;
        nrk2 << t << " " << u_n << " " << N - u_n << std::endl;
        double U_1, U_2;
        U_1 = U_2 = u_n;

        for(int mu=0; mu<=mu_max; mu++)
        {
            double temp = u_n;
            u_n = u_0 + dt*(b_1*fun(beta, N, gamma, U_1) + b_2*fun(beta, N, gamma, U_2));  
            if(fabs(u_n - temp) < TOL){ break; }

            double m_1_1, m_1_2, m_2_1, m_2_2;
            m_1_1 = 1 - dt * a_1_1 * (alpha - 2*beta*U_1);
            m_1_2 = -dt * a_1_2 * (alpha - 2*beta*U_2);
            m_2_1 = -dt * a_2_1 * (alpha - 2*beta*U_1);
            m_2_2 = 1 - dt * a_2_2 * (alpha - 2*beta*U_2);
            double F_1, F_2;
            F_1 = U_1 - u_n - dt*(a_1_1 * (alpha*U_1 - beta*U_1*U_1) + a_1_2 *(alpha*U_2 - beta*U_2*U_2));
            F_2 = U_2 - u_n - dt*(a_2_1 * (alpha*U_1 - beta*U_1*U_1) + a_2_2 *(alpha*U_2 - beta*U_2*U_2));

            double dU_1, dU_2;
            dU_1 = (F_2*m_1_2 - F_1*m_2_2)/(m_1_1*m_2_2 - m_1_2*m_2_1);
            dU_2 = (F_1*m_2_1 - F_2*m_1_1)/(m_1_1*m_2_2 - m_1_2*m_2_1);
            U_1 += dU_1;
            U_2 += dU_2;
        }
        u_0 = u_n;
    }
    nrk2.close();
    
}
