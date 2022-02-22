#include <iostream>
#include <fstream>

#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"


const int n_x = 40;
const int n_y = 40;
const int N = (n_x + 1) * (n_y + 1);
const double delta = 1.0;
const double dt = 1.0;
const double T_A = 40.0;
const double T_B = 0.0;
const double T_C = 30.0;
const double T_D = 0.0;
const double k_B = 0.1;
const double k_D = 0.6;
const int IT_MAX = 2000;



int main(){
    gsl_matrix *A = gsl_matrix_calloc(N,N);
    gsl_matrix *B = gsl_matrix_calloc(N,N);
	gsl_vector *C = gsl_vector_calloc(N);

    // Wnetrze obszaru
    for(int i=1; i<n_x; i++){
        for(int j=1; j<n_y; j++){
            int l = i + j * (n_x+1);

            gsl_matrix_set(A, l, l-n_x-1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l-1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l+1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l+n_x+1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l, -2.0*dt/(delta*delta)-1);

            gsl_matrix_set(B, l, l-n_x-1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l-1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l+1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l+n_x+1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l, 2.0*dt/(delta*delta)-1.0);

        }
    }	

	// WB Dirichleta (lewy i prawy brzeg)
    for(int i=0; i<=n_x; i+=n_x){
        for(int j=0; j<=n_y; j++){
            int l = i + j * (n_x+1);

            gsl_matrix_set(A, l, l, 1.0);
            gsl_matrix_set(B, l, l, 1.0);
            gsl_vector_set(C, l, 0.0);
        }
    }

    // WB von Neumanna na gornym brzegu dla chwili n+1
    for(int i=1; i<n_x; i++){
        int l = i + n_y * (n_x+1);
        gsl_matrix_set(A, l, l-n_x-1, -1.0/(k_B*delta));
        gsl_matrix_set(A, l, l, 1.0 + 1.0/(k_B*delta));
        gsl_vector_set(C, l, T_B);
        for(int k=0; k<N; k++){
            gsl_matrix_set(B, l, k, 0.0);
        }
    }

    // WB von Neumanna na dolnym brzegu dla chwili n+1
    for(int i=1; i<n_x; i++){
        int l = i;
        gsl_matrix_set(A, l, l+n_x+1, -1.0/(k_D*delta));
        gsl_matrix_set(A, l, l, 1.0 + 1.0/(k_D*delta));
        gsl_vector_set(C, l, T_D);
        for(int k=0; k<N; k++){
            gsl_matrix_set(B, l, k, 0.0);
        }
    }

    gsl_vector *T = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);  
    gsl_permutation *p = gsl_permutation_calloc(N);
    int signum = 0;
    // Warunki poczatkowe dla wektora T
    for(int i=0; i<=n_x; i++){
        for(int j=0; j<=n_y; j++){
            int l = i + j * (n_x+1);

            // lewy brzeg
            if(i == 0){
                gsl_vector_set(T, l, T_A);
            }
            // prawy brzeg
            else if(i == n_x){
                gsl_vector_set(T, l, T_C);
            }
            // pozostaly obszar
            else{
                gsl_vector_set(T, l, 0);
            }
        }
    }

    gsl_linalg_LU_decomp(A, p, &signum);


    std::ofstream T_file("T.dat");
    std::ofstream nabla_2_T_file("nabla_2_T.dat");
    
    // Algorytm CN
    for(int it=1; it<=IT_MAX; it++){

        // T = T_n
        gsl_blas_dgemv(CblasNoTrans, 1, B, T, 0, d);
        gsl_blas_daxpy(1, C, d);

        // T = T_n+1
        gsl_linalg_LU_solve(A, p, d, T);
        if(it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000){
            for(int i=1; i<n_x; i++){
                for(int j=1; j<n_y; j++){
                    int l = i + j * (n_x+1);
                    T_file << i*delta << " " << j*delta << " " << gsl_vector_get(T, l) << std::endl;
                    
                    double nabla_2_T = ((gsl_vector_get(T, l + 1) - 2.0*gsl_vector_get(T, l) + gsl_vector_get(T, l-1))/pow(delta, 2)) + ((gsl_vector_get(T, l + n_x + 1) - 2.0*gsl_vector_get(T, l) + gsl_vector_get(T, l- n_x -1) )/pow(delta, 2));
                    nabla_2_T_file << i*delta << " " << j*delta << " " << nabla_2_T << std::endl;
              }
            }
            T_file << std::endl << std::endl;
            nabla_2_T_file << std::endl << std::endl;
        }
    }

    T_file.close();
    nabla_2_T_file.close();

	gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(C);
    gsl_vector_free(T);
    gsl_vector_free(d);
    gsl_permutation_free(p);
}
