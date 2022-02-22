#include <iostream>
#include "mgmres.h"
#include <cmath>
#include <fstream>

const double delta = 0.1;

double rho_1(double x, double y, double x_max, double y_max, double sigma){
    return exp( -pow( (x-0.25*x_max)/sigma, 2) - pow((y-0.5*y_max)/sigma,2));
}

double rho_2(double x, double y, double x_max, double y_max, double sigma){
    return -exp( -pow((x-0.75*x_max)/sigma, 2) - pow((y-0.5*y_max)/sigma,2));
}

double eps_l(int l, int n_x, double eps_1, double eps_2){
    int j =l/(n_x+1);  
    return l - (j * (n_x + 1)) <= n_x/2 ? eps_1 : eps_2;
}


void gmres(int n_x, int n_y, double eps_1, double eps_2, double V_1, double V_2, double V_3, double V_4, bool enable_rho, FILE *map_file){
    int N = (n_x+1)*(n_y+1);
    double a[5*N] = {0.}; // niezerowe wartosci el. macierzowych
    int ja[5*N]  = {0}; // przechowuje informacje o numerac kolumn
    int ia[N+1] = {0}; // wskazniki do elementow rozpoczynających dany wiersz

    for(int i=0; i<=N; i++){ ia[i] = -1.0; }
    double b[N] = {0.};
    double V[N] = {0.};

    double x_max = n_x*delta;
    double y_max = n_y*delta;
    double sigma = x_max/10.0;
    // wypelnianie macierzy rzadkich - warunki brzegowe Dirichleta
    
    int k = -1;
    for(int l=0; l<N; l++){

        int j = l/(n_x+1);
        int i = l - j * (n_x+1);

        int brzeg = 0;
        double vb = 0.;

        double epsilon_l = eps_l(l, n_x, eps_1, eps_2);
        double epsilon_l_1 = eps_l(l+1, n_x, eps_1, eps_2);
        double epsilon_l_nx_1 = eps_l(l+n_x+1, n_x, eps_1, eps_2);

        if(i == 0){
            brzeg = 1;
            vb = V_1;
        }

        if(j == n_y){
            brzeg = 1;
            vb = V_2;
        }

        if(i == n_x){
            brzeg = 1;
            vb = V_3;
        }

        if(j==0){
            brzeg = 1;
            vb = V_4;
        }

        b[l] = enable_rho ? -(rho_1(i*delta, j*delta, x_max, y_max, sigma) 
                            + rho_2(i*delta, j*delta, x_max, y_max, sigma)) : 0;

        if(brzeg == 1)
            b[l] = vb;
        
        ia[l] = -1;

        // lewa skrajna przekatna
        if(l - n_x - 1 >= 0 && brzeg == 0){
            k++;
            if(ia[l]<0) ia[l] = k;
            a[k] = epsilon_l/(delta*delta);
            ja[k] = l - n_x - 1;
        }

        // poddiagonala
        if(l - 1 >= 0 && brzeg == 0){
            k++;
            if(ia[l]<0) ia[l]=k;
            a[k] = epsilon_l / (delta*delta);
            ja[k] = l - 1;
        }

        // diagonala
        k++;
        if(ia[l] < 0) ia[l] = k;
        
        if(brzeg==0){
            a[k] = -(2*epsilon_l + epsilon_l_1 + epsilon_l_nx_1)/(delta*delta);
        } else { a[k] = 1; }
        
        ja[k] = l;

        //naddiagonala
        if(l < N && brzeg == 0){
            k++;
            a[k] = epsilon_l_1/(delta*delta);
            ja[k] = l+1;
        }

        // prawa skrajna przekatna
        if(l< N-n_x-1 && brzeg==0){
            k++;
            a[k] = epsilon_l_nx_1/(delta*delta);
            ja[k] = l + n_x + 1;
        }
    }
    int nz_num = k+1; //ilosc niezerowych elementow
    ia[N] = nz_num;

    if(n_x == 4 && n_y == 4){
        std::ofstream a_file("a.dat");
        std::ofstream b_file("b.dat");
        for(int l=0; l<N; l++){
            int j = l/(n_x+1);
            int i = l - j*(n_x+1);
            a_file << l << " " << i << " " << j << " " << a[l] << std::endl;
            b_file << l << " " << i << " " << j << " " << b[l] << std::endl;
        }
        a_file.close();
        b_file.close();
    } else {
        pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, 500, 500, 1e-8, 1e-8);
        for(int l=0; l<N; l++){
            int j = l/(n_x+1);
            int i = l - j * (n_x+1);

            fprintf(map_file, "%15d %15d %15g\n", i, j, V[l]);
        }
        fprintf(map_file, "\n\n");
    }
}


int main(){

    gmres(4,4, 1,1, 10,-10, 10,-10, false, nullptr);
    FILE *zad_5 = fopen("zad_5.dat", "w");
    gmres(50, 50, 1, 1, 10, -10, 10, -10, false, zad_5);
    gmres(100, 100, 1, 1, 10, -10, 10, -10, false, zad_5);
    gmres(200, 200, 1, 1, 10, -10, 10, -10, false, zad_5);
    fclose(zad_5);

    FILE *zad_6 = fopen("zad_6.dat", "w");
    gmres(100, 100, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, true, zad_6);
    gmres(100, 100, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, true, zad_6);
    gmres(100, 100, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, true, zad_6);
    fclose(zad_6);
}