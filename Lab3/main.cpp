#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>


double f(double x, double v)
{
    return v;
}

double g(double x, double v, double alpha)
{
    return alpha * (1 - x*x) * v - x;
}

std::vector<double> metoda_trapezow(double x_n, double v_n, double dt, double alpha)
{
    double x_new = x_n;
    double v_new = v_n;

    double dx, dv;
    do
    {
        double F = x_new - x_n - dt/2 * (f(x_n, v_n) + f(x_new, v_new));
        double G = v_new - v_n - dt/2 * (g(x_n, v_n, alpha) + g(x_new, v_new, alpha));

        double a_1_1 = 1;
        double a_1_2 = -dt/2;
        double a_2_1 = -dt/2 * (-2*alpha*x_new*v_new - 1);
        double a_2_2 = 1 - dt/2 * alpha * (1 - x_new*x_new);

        dx = (-F * a_2_2 + G * a_1_2)/(a_1_1*a_2_2 - a_1_2*a_2_1);
        dv = (-G * a_1_1 + F * a_2_1)/(a_1_1*a_2_2 - a_1_2*a_2_1);

        x_new += dx;
        v_new += dv;
    } while (fabs(dx) < 1e-10 && fabs(dv) < 1e-10);
    
    return std::vector<double>{x_new, v_new};

}

std::vector<double> rk2(double x_n, double v_n, double dt, double alpha)
{
    
    double k1x = v_n;
    double k1v = alpha * (1 - x_n*x_n) * v_n - x_n;

    double k2x = v_n + dt * k1v;
    double k2v = alpha * (1 - pow(x_n + dt * k1x, 2)) * (v_n + dt * k1v) - (x_n + dt*k1x);

    double x_new = x_n + dt/2 * (k1x + k2x);
    double v_new = v_n + dt/2 * (k1v + k2v);
    return std::vector<double>{x_new, v_new};
}



void control_t(std::vector<double> (*schemat_numeryczny)(double, double, double, double), double tol, std::string filename)
{
    std::ofstream file(filename);
    file.precision(8);
    double x_0 = 0.01;
    double v_0 = 0;
    double dt_0 = 1;
    double S = 0.75;
    double p = 2;
    double t_max = 40;
    double alpha = 5;

    double t = 0;
    double dt = dt_0;
    double x_n = x_0;
    double v_n = v_0;

    do
    {
        std::vector<double> xv2_1 = schemat_numeryczny(x_n, v_n, dt, alpha);
        std::vector<double> xv2_2 = schemat_numeryczny(xv2_1[0], xv2_1[1], dt, alpha);
        
        std::vector<double> xv1_2 = schemat_numeryczny(x_n, v_n, 2*dt, alpha);

        double Ex = (xv2_2[0] - xv1_2[0])/(pow(2,p) - 1);
        double Ev = (xv2_2[1] - xv1_2[1])/(pow(2,p) - 1);
        if(std::max(fabs(Ex), fabs(Ev)) < tol)
        {
            t += 2 * dt;
            x_n = xv2_2[0];
            v_n = xv2_2[1];
            file << t << "\t" << dt << "\t" << x_n << "\t" << v_n << std::endl;
        }
        dt *= pow((S*tol/std::max(fabs(Ex), fabs(Ev))), 1.0/(p+1.0));

    } while (t<t_max);
    
    file.close();
}


int main()
{

    control_t(metoda_trapezow, 1e-2, "metoda_trapezow_10_2.dat");
    control_t(metoda_trapezow, 1e-5, "metoda_trapezow_10_5.dat");
    
    control_t(rk2, 1e-2, "rk2_10_2.dat");
    control_t(rk2, 1e-5, "rk2_10_5.dat");
}
