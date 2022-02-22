#include <iostream>
#include <cmath>


const double delta = 0.2;
const int n_x = 128;
const int n_y = 128;
const double x_max = delta * n_x;
const double y_max = delta * n_y;
const double TOL = 1e-8;


double S(double (*V)[n_y+1], int k)
{
    double result = 0.0;
    for(int i=0; i<=n_x-k; i+=k){
        for(int j=0; j<=n_y-k; j+=k){
            result += pow(k*delta,2)/2 * (pow((V[i+k][j]-V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta)  ,2) 
					+ pow((V[i][j+k]-V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta) ,2));
        }
    }
    return result;
}


int main()
{
    
    double V[n_x+1][n_y+1] = {0.};

    //warunki brzegowe Dirichleta
    for(int j=0; j<=n_y; j++){
        V[0][j] = sin(M_PI * (j*delta)/y_max);
        V[j][n_y] = -1. * sin(2*M_PI * (j*delta)/x_max);
        V[n_x][j] = sin(M_PI * (j*delta)/y_max);
        V[j][0] = sin(2*M_PI * (j*delta)/x_max);
    }
    
    int it = 0;

    FILE *V_grid = fopen("V.dat", "w");
    FILE *S_it = fopen("S.dat", "w");

    for(int k=16; k>0; k /= 2)
    {
        double S_old = 0.0;
        double S_new = 0.0;
        while(true)
        {
            it++;
            for(int i=k; i<=n_x-k; i+=k)
                for(int j=k; j<=n_y-k; j+=k)
                    V[i][j] = 0.25 * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);

            S_old = S_new;
            S_new = S(V, k);

            fprintf(S_it, "%15d %15g\n", it, S_new);

            // warunek stopu
            if(fabs((S_new - S_old)/S_new) < TOL){ break; } 
        }
        
        for(int i=0; i<=n_x; i+=k)
            for (int j=0; j<=n_y; j+=k)
                fprintf(V_grid, "%15g %15g %15g\n", i*delta, j*delta, V[i][j]);
            
        // zageszczanie siatki
        if(k > 1)
        {
            for(int i=0; i<=n_x-k; i+=k)
            {
                for(int j=0; j<=n_y-k; j+=k)
                {
                    V[i+k/2][j+k/2] = 0.25 * (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                    if(i!=n_x-k)    
                        V[i+k][j+k/2] = 0.5 * (V[i+k][j] + V[i+k][j+k]);
                    if(j!=n_y-k)
                        V[i+k/2][j+k] = 0.5 * (V[i][j+k] + V[i+k][j+k]);
                    if(j!=0)
                        V[i+k/2][j] = 0.5 * (V[i][j] + V[i+k][j]);
                    if(i!=0)
                        V[i][j+k/2] = 0.5 * (V[i][j] + V[i][j+k]);
                }
            }
        }

        fprintf(S_it, "\n\n");
        fprintf(V_grid, "\n\n");
    }

    fclose(V_grid);
    fclose(S_it);
}


