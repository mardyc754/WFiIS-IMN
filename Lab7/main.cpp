#include <iostream>
#include <cmath>
#include <fstream>

const double delta = 0.01;
const double rho = 1;
const double mu = 1;
const int n_x = 200;
const int n_y = 90;
const int i_1 = 50;
const int j_1 = 55;
const int IT_MAX = 20000;

double x(int i){ return i*delta; }
double y(int j){ return j*delta; }

double gamma(double (*psi)[n_y+1], double (*zeta)[n_y+1]){
    double result = 0.0;
    int j_2 = j_1 + 2;
    for(int i=1; i<n_x-1; i++){
        result += (psi[i+1][j_2] + psi[i-1][j_2] + psi[i][j_2+1] + psi[i][j_2-1] - 4*psi[i][j_2] - delta*delta * zeta[i][j_2]);
    }
    return result;
}

void WB_psi(double (*psi)[n_y+1], double Q_we){
    double Q_wy = Q_we * (pow(y(n_y), 3) - pow(y(j_1), 3) - 3*y(j_1)*y(n_y)*y(n_y) + 3*y(j_1)*y(j_1)*y(n_y))/pow(y(n_y), 3);

    // brzeg A (wejscie)
    for(int j=j_1; j<=n_y; j++){
        psi[0][j] = Q_we/(2*mu) * (pow(y(j),3)/3 - pow(y(j), 2)/2 * (y(j_1) + y(n_y)) + y(j)*y(j_1)*y(n_y));
    }

    // brzeg C (wyjscie)
    for(int j=0; j<=n_y; j++){
        psi[n_x][j] = Q_wy/(2*mu) * (pow(y(j), 3)/3 - pow(y(j), 2)/2 * y(n_y)) + (Q_we*y(j_1)*y(j_1) * (-y(j_1) + 3*y(n_y)))/(12*mu);
    }

    // brzeg B
    for(int i=1; i<n_x; i++){
        psi[i][n_y] = psi[0][n_y];
    }

    // brzeg D
    for(int i=i_1; i<n_x; i++){
        psi[i][0] = psi[0][j_1];
    }

    // brzeg E
    for(int j=1; j<=j_1; j++){
        psi[i_1][j] = psi[0][j_1];
    }

    // brzeg F
    for(int i=1; i<=i_1; i++){
        psi[i][j_1] = psi[0][j_1];
    }
}

void WB_zeta(double (*zeta)[n_y+1], double (*psi)[n_y+1], double Q_we){
    double Q_wy = Q_we * (pow(y(n_y), 3) - pow(y(j_1), 3) - 3*y(j_1)*y(n_y)*y(n_y) + 3*y(j_1)*y(j_1)*y(n_y))/pow(y(n_y), 3);
    
    // brzeg A (wejscie)
    for(int j=j_1; j<=n_y; j++){
        zeta[0][j] = Q_we/(2*mu) * (2*y(j) - y(j_1) - y(n_y));
    }

    // brzeg C (wyjscie)
    for(int j=0; j<=n_y; j++){
        zeta[n_x][j] = Q_wy/(2*mu) * (2*y(j) - y(n_y));
    }

    // brzeg B
    for(int i=1; i<n_x; i++){
        zeta[i][n_y] = 2/pow(delta,2) * (psi[i][n_y-1] - psi[i][n_y]);
    }

    // brzeg D
    for(int i=i_1+1; i<n_x; i++){
        zeta[i][0] = 2/pow(delta,2) * (psi[i][1] - psi[i][0]);
    }

    // brzeg E
    for(int j=1; j<j_1; j++){
        zeta[i_1][j] = 2/pow(delta,2) * (psi[i_1+1][j] - psi[i_1][j]);
    }

    // brzeg F
    for(int i=1; i<=i_1; i++){
        zeta[i][j_1] = 2/pow(delta,2) * (psi[i][j_1+1] - psi[i][j_1]);
    }

    // wierzchołek E/F
    zeta[i_1][j_1] = 0.5*(zeta[i_1-1][j_1] + zeta[i_1][j_1-1]);
}


void relaxation_NS(double Q, FILE *psi_file, FILE *zeta_file, FILE *u_file, FILE *v_file){
    double psi[n_x+1][n_y+1] = {0.};
    double zeta[n_x+1][n_y+1] = {0.};

    double u[n_x+1][n_y+1] = {0.};
    double v[n_x+1][n_y+1] = {0.};

    WB_psi(psi, Q);

    for(int it=1; it<=IT_MAX; it++){
        double Omega = it < 2000 ? 0 : 1;

        for(int i=1; i<n_x; i++){
            for(int j=1; j<n_y; j++){
                if(i > i_1 || j > j_1){
                    psi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - delta*delta * zeta[i][j]);
                    zeta[i][j] = 0.25 * (zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]) -
                                Omega * rho/(16*mu) * ( (psi[i][j+1] - psi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j]) -
                                (psi[i+1][j] - psi[i-1][j]) * (zeta[i][j+1] - zeta[i][j-1]));
                    u[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*delta);
                    v[i][j] = -(psi[i+1][j] - psi[i-1][j])/(2*delta);
                } 
            }
        }

        WB_zeta(zeta, psi, Q);
        double Gamma = gamma(psi, zeta);
        if(it%100 == 0){
            std::cout << "Gamma = " << Gamma << std::endl;
        }
    }

    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            double X = x(i);
            double Y = y(j);
            fprintf(psi_file, "%15g %15g %15g\n", X, Y, psi[i][j]);
            if(Q != 4000){
                fprintf(zeta_file, "%15g %15g %15g\n", X, Y, zeta[i][j]);
                fprintf(u_file, "%15g %15g %15g\n", X, Y, u[i][j]);
                fprintf(v_file, "%15g %15g %15g\n", X, Y, v[i][j]);
            }
        }
    }

    std::cout << std::endl << std::endl;

    
    fprintf(psi_file, "\n\n");
    if(Q != 4000){
        fprintf(zeta_file, "\n\n");
        fprintf(u_file, "\n\n");
        fprintf(v_file, "\n\n");
    }
}

int main(){

    FILE *psi_file = fopen("psi.dat", "w");
    FILE *zeta_file = fopen("zeta.dat", "w");
    FILE *u_file = fopen("u.dat", "w");
    FILE *v_file = fopen("v.dat", "w");

    relaxation_NS(-1000, psi_file, zeta_file, u_file, v_file);
    relaxation_NS(-4000, psi_file, zeta_file, u_file, v_file);
    relaxation_NS(4000, psi_file, nullptr, nullptr, nullptr);

    fclose(psi_file);
    fclose(zeta_file);
    fclose(u_file);
    fclose(v_file);
}