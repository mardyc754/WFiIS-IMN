#include <iostream>
#include <cmath>
#include <fstream>

const int n_x = 400;
const int n_y = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 10 * delta;
const double x_A = 0.45;
const double y_A = 0.45;


void load_psi(double (*psi)[n_y+1]){
    
    std::ifstream psi_file("psi.dat");
    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            psi_file >> i >> j >> psi[i][j];
        }
    }
    psi_file.close();
}

double v(double v_x, double v_y){
    return sqrt(v_x*v_x + v_y*v_y);
}

void copy_u(double (*u_0)[n_y+1], double (*u_1)[n_y+1]){
    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            u_1[i][j] = u_0[i][j];
        }
    }
}

double get_v_max(double (*v_x)[n_y+1], double (*v_y)[n_y+1]){
    double v_max = v(v_x[0][0], v_y[0][0]);
    for(int i=0; i<n_x; i++){
        for(int j=0; j<n_y; j++){
            v_max = std::max(v_max, v(v_x[i][j], v_y[i][j]));
        }
    }
    return v_max;
}
double delta_t(double (*v_x)[n_y+1], double (*v_y)[n_y+1]){
    double v_max = get_v_max(v_x, v_y);
    return delta/(4*v_max);   
}

void velocity_field(double (*psi)[n_y+1], double (*v_x)[n_y+1], double (*v_y)[n_y+1]){
    std::ofstream v_file("v.dat");
    

    for(int i=1; i<n_x; i++){
        for(int j=1; j<n_y; j++){
            v_x[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*delta);
            v_y[i][j] = -(psi[i+1][j] - psi[i-1][j])/(2*delta);
        }
    }

    // zastawka
    for(int i=i_1; i<=i_2; i++){
        for(int j=0; j<=j_1; j++){
            v_x[i][j] = v_y[i][j] = 0;
        }
    }

    // dolny i górny brzeg
    for(int i=1; i<n_x; i++){
        v_x[i][0] = v_y[i][n_y] = 0;
    }

    // lewy i prawy brzeg
    for(int j=0; j<=n_y; j++){
        v_x[0][j] = v_x[1][j];
        v_x[n_x][j] = v_x[n_x-1][j];
    }

    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            v_file << i*delta << " " << j*delta << " " << v_x[i][j] << std::endl;
        }
    }
    v_file << std::endl << std::endl;

    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            v_file << i*delta << " " << j*delta << " " << v_y[i][j] << std::endl;
        }
    }
    v_file.close();
}

double c_t(double (*u)[n_y+1]){
    double result = 0;
    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            result += u[i][j];
        }
    }
    return result * delta * delta;
}

double x_sr(double (*u)[n_y+1]){
    double result = 0;
    for(int i=0; i<=n_x; i++){
        double x = i*delta;
        for(int j=0; j<=n_y; j++){
            result += x * u[i][j];
        }
    }
    return result * delta * delta;
}

double u_ij(double (*u_0)[n_y+1], double (*u_1)[n_y+1], int i, int j, double dt, double D, double (*v_x)[n_y+1], double (*v_y)[n_y+1]){
    int left = i-1;
    int right = i+1;
    
    if(i == 0) left = n_x;
    else if(i == n_x)  right = 0;

    return ( 1.0/( 1.0+( (2.0*D*dt) / pow(delta, 2)) ) ) * ( u_0[i][j] - (dt/2.0) * v_x[i][j] *
    ( ( (u_0[right][j] - u_0[left][j])/(2.0*delta) ) + (u_1[right][j] - u_1[left][j])/(2.0*delta) ) - (dt / 2.0) * v_y[i][j] * 
    ( ( u_0[i][j+1] - u_0[i][j-1] )/(2.0*delta) + (u_1[i][j+1] - u_1[i][j-1])/(2.0*delta) ) + (dt/2.0) * D * 
    ( ( u_0[right][j] + u_0[left][j] + u_0[i][j+1] + u_0[i][j-1] - 4*u_0[i][j] )/pow(delta,2) + 
    ( u_1[right][j] + u_1[left][j] + u_1[i][j+1] + u_1[i][j-1] )/pow(delta,2) ));
}


void AD(double D, double (*psi)[n_y+1], double (*v_x)[n_y+1], double (*v_y)[n_y+1], double dt, FILE *c_x_sr_file, std::string u_filename){

    double u_0[n_x+1][n_y+1] = {0.};
    double u_1[n_x+1][n_y+1] = {0.};

    std::ofstream u_file(u_filename);

    for(int i=0; i<=n_x; i++){
        double x = i*delta;
        for(int j=0; j<=n_y; j++){
            double y = j*delta;
            if(i < i_1 || i > i_2 || j > j_1){
                u_0[i][j] = 1/(2*M_PI*pow(sigma, 2)) * exp(- (pow(x-x_A, 2) + pow(y-y_A, 2))/(2*sigma*sigma));
            }
        }
    }

    int k = 0;
    int IT_MAX = 10000;
    for(int it=1; it<=IT_MAX; it++){

        copy_u(u_0, u_1);

        if(k == 0){
            for(int i=0; i<=n_x; i++){
                double x = i * delta;
                for(int j=0; j<=n_y; j++){
                    double y = j * delta;
                    u_file << x << " " <<  y << " " << u_1[i][j] << std::endl;
                }
            }
            std::cout << "k = " << k << std::endl;
            u_file << std::endl << std::endl;
            k++;
        }

        #pragma omp parallel for
        for(int K=1; K<=20; K++){
            for(int i=0; i<=n_x; i++){
                for(int j=1; j<n_y; j++){
                    if (i < i_1 || i > i_2 || j > j_1) {
                        u_1[i][j] = u_ij(u_0, u_1, i, j, dt, D, v_x, v_y);
                    }
                }
            }
        }

        copy_u(u_1, u_0);

        int T = k*IT_MAX/5;
        if(T == it){
            for(int i=0; i<=n_x; i++){
                double x = i * delta;
                for(int j=0; j<=n_y; j++){
                    double y = j * delta;
                    u_file << x << " " <<  y << " " << u_1[i][j] << std::endl;
                }
            }
            std::cout << "k = " << k << std::endl;
            u_file << std::endl << std::endl;
            k++;
        }

        fprintf(c_x_sr_file, "%15g %15g %15g\n", it*dt, x_sr(u_1), c_t(u_1));
    }

    fprintf(c_x_sr_file, "\n\n");
    std::cout << std::endl << std::endl;
    u_file.close();
}


int main(){
    double psi[n_x+1][n_y+1] = {0.};
    double v_x[n_x+1][n_y+1] = {0.};
    double v_y[n_x+1][n_y+1] = {0.};

    load_psi(psi);
    velocity_field(psi, v_x, v_y);
    double dt = delta_t(v_x, v_y);

    FILE *c_x_sr_file = fopen("c_x_sr.dat", "w");
    
    std::cout << "D = 0.0" << std::endl;
    AD(0, psi, v_x, v_y, dt, c_x_sr_file, "u_D_0.dat");

    std::cout << "D = 0.1" << std::endl;
    AD(0.1, psi, v_x, v_y, dt, c_x_sr_file, "u_D_01.dat");
    
    fclose(c_x_sr_file);
}