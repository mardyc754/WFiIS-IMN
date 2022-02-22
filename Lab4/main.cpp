// uruchamiac z flaga -O3 (dotyczy kompilatora g++)
#include <iostream>
#include <cmath>

const double eps = 1.0;
const double delta = 0.1;
const int n_x = 150;
const int n_y = 100;
const double V_1 = 10;
const double V_2 = 0;
const double x_max = delta * n_x;
const double y_max = delta * n_y;
const double TOL = 1e-8;

// funkcja gestosci
double get_rho(double x, double y)
{
    double sigma_x = 0.1 * x_max;
    double sigma_y = 0.1 * y_max;
    double rho_1 = 1.0 * exp(-1.*pow((x-0.35*x_max)/sigma_x, 2) - pow((y-0.5*y_max)/sigma_y,2)); 
    double rho_2 = -1.0 * exp(-1.*pow((x-0.65*x_max)/sigma_x, 2) - pow((y-0.5*y_max)/sigma_y,2));
    return rho_1 + rho_2;
}

// warunek stopu
double S(double (*V)[n_y+1], double (*rho)[n_y+1])
{
    double result = 0;
    for(int i=0; i<n_x; i++)
    {
        for(int j=0; j<n_y; j++)
        {
            result += pow(delta,2) * (0.5*pow((V[i+1][j] - V[i][j])/delta,2) + 0.5*pow((V[i][j+1] - V[i][j])/delta,2) - rho[i][j] * V[i][j]);
        }
    }
    return result;
}


void global_relaxation(double omega_G, FILE *rg_S, FILE *rg_V, FILE *rg_d)
{
    double V_s[n_x+1][n_y+1] = {0.};
    double V_n[n_x+1][n_y+1] = {0.};
    double rho[n_x+1][n_y+1] = {0.};
    
    for(int i=0; i<=n_x; i++)
    {
        for(int j=0; j<=n_y; j++)
        {
            rho[i][j] = get_rho(i*delta, j*delta);
        }
    }

    // warunki brzegowe
    for(int i=0; i<=n_x; i++)
    {
        V_s[i][0] = V_n[i][0] = V_1;
        V_s[i][n_y] = V_n[i][n_y] = V_2;
    }

    int it = 0;
    double S_new = 0.0;

    while(true){
        it++;
        for(int i=1; i<n_x; i++)
        {
            for(int j=1; j<n_y; j++)
            {
                V_n[i][j] = 0.25 * (V_s[i+1][j] + V_s[i-1][j] + V_s[i][j-1] + V_s[i][j+1] + (delta*delta/eps)*rho[i][j]);  
            }
        }

        // warunki brzegowe Neumanna
        for(int j=1; j<n_y; j++)
        {
            V_n[0][j] = V_n[1][j];
            V_n[n_x][j] = V_n[n_x-1][j];
        }

        for(int i=0; i<=n_x; i++)
        {
            for(int j=0; j<=n_y; j++)
            {
                V_s[i][j] = (1. - omega_G) * V_s[i][j] + omega_G * V_n[i][j];  
            }
        }
        
        // Warunek stopu
        double S_old = S_new;
        S_new = S(V_n, rho);

        fprintf(rg_S, "%d %g\n", it, S_new);
        if(fabs((S_new-S_old)/S_old) < TOL) { break; }
    }


    for(int i=0; i<=n_x; i++)
    {
        for(int j=0; j<=n_y; j++)
        {
            fprintf(rg_V, "%g %g %g\n", i*delta, j*delta, V_n[i][j]);
        }
    }

    for(int i=1; i<n_x; i++)
    {
        for(int j=1; j<n_y; j++)
        {
            double nabla_2_V = (V_n[i+1][j] - 2*V_n[i][j] + V_n[i-1][j])/pow(delta, 2) + (V_n[i][j+1] - 2*V_n[i][j] + V_n[i][j-1])/pow(delta,2);
            fprintf(rg_d, "%g %g %g\n", i*delta, j*delta, nabla_2_V + rho[i][j]/eps); 
        }
    }

    fprintf(rg_S, "\n\n");
    fprintf(rg_V, "\n\n");
    fprintf(rg_d, "\n\n");

}

void local_relaxation(double omega_L, FILE *fp)
{
    double V[n_x+1][n_y+1] = {0.};
    double rho[n_x+1][n_y+1] = {0.};
    
    for(int i=0; i<=n_x; i++)
    {
        for(int j=0; j<=n_y; j++)
        {
            rho[i][j] = get_rho(i*delta, j*delta);
        }
    }


    // warunki brzegowe
    for(int i=0; i<=n_x; i++)
    {
        V[i][0] = V_1;
        V[i][n_y] = V_2;
    }

    int it = 0;
    double S_new = 0.0;

    while(true){
        it++;
        for(int i=1; i<n_x; i++)
        {
            for(int j=1; j<n_y; j++)
            {
                V[i][j] = (1.0 - omega_L) * V[i][j] + omega_L/4.0 * (V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + delta*delta/eps * rho[i][j]);
            }
        }

        // warunki brzegowe Neumanna
        for(int j=1; j<n_y; j++)
        {
            V[0][j] = V[1][j];
            V[n_x][j] = V[n_x-1][j];
        }

        // warunek stopu
        double S_old = S_new;
        S_new = S(V, rho);
        fprintf(fp, "%10d %15g\n", it, S_new);
        
        if(fabs((S_new-S_old)/S_old) < TOL) { break; }
    }
    fprintf(fp, "\n\n");
}


int main()
{
    FILE *rg_S = fopen("relaksacja_globalna_S_t.dat", "w");
    FILE *V_x_y = fopen("relaksacja_globalna_V_x_y.dat", "w");
    FILE *d_x_y = fopen("relaksacja_globalna_d_x_y.dat", "w");
    global_relaxation(0.6, rg_S, V_x_y, d_x_y);
    global_relaxation(1.0, rg_S, V_x_y, d_x_y);
    fclose(rg_S);
    fclose(V_x_y);
    fclose(d_x_y);

    FILE *rl_S = fopen("relaksacja_lokalna_S_t.dat", "w");
    local_relaxation(1.0, rl_S);
    local_relaxation(1.4, rl_S);
    local_relaxation(1.8, rl_S);
    local_relaxation(1.9, rl_S);
    fclose(rl_S);
}
