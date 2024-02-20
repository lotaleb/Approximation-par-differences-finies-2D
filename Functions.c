#include <math.h>
#include <stdio.h>
#include <stdlib.h> // Pour malloc et free
#include <time.h>
#include "Functions.h"
#define PI 3.14
// Définition des variables globales
double alpha = 1;
double nu = 1;
int b = 1;
int a = 0;
double test_arret = 1e-10; // Ajustement de la tolérance

// Définition des fonctions
int index(int i, int j, int m) { return i * (m + 2) + j; }

double square(double x) { return x * x; }

double sol_exact(int i, int j, double hx, double hy) {
  double x = a + (double) i * hx;
  double y = a + (double)j * hy;
  return x * y * sin(2 * PI * x) * cos(PI * y / 2);
}


//Méthode de Jacobi
void Jacobi(int m, int n) {
      //coeficients en espace 
      double hx = (b - a) / (m + 1.0);
      double hy = (b - a) / (n + 1.0);
      double lambda = alpha + 2 * nu * ((1 / (hx * hx)) + (1 / (hy * hy)));
      double hx2 = nu / (hx * hx);
      double hy2 = nu / (hy * hy);
      //Tolérance  et test d'arret 
      int max_iter = 10 * (m + 1);
      double sum_diff = 1;
      double sum_uij = 1;
      int iter = 0;
      double test_arret = 1e-10; // tolérance pour l'arrêt de la boucle
      double error;

      // Tableaux
      double *u, *u_prev, *f;

      // Allocation de mémoire
      u = (double *)malloc((m + 2) * (n + 2) * sizeof(double));
      u_prev = (double *)malloc((m + 2) * (n + 2) * sizeof(double));
      f = (double *)malloc((m + 2) * (n + 2) * sizeof(double));

      // Initialisation de u et u_prev à 0
      for (int i = 0; i < m + 2; i++) {
          for (int j = 0; j < n + 2; j++) {
              int ij = index(i, j, m);
              f[ij] = square(17 * PI / 4) * (a + i * hx) * (a + j * hy) *
                          sin(2 * PI * (a + i * hx)) * cos(PI * (a + j * hy) / 2) +
                      PI * (a + i * hx) * sin(2 * PI * (a + i * hx)) *
                          sin(PI * (a + j * hy) / 2) -
                      4 * PI * (a + j * hy) * cos(2 * PI * (a + i * hx)) *
                          cos(PI * (a + j * hy) / 2) +
                      (a + i * hx) * (a + j * hy) * sin(2 * PI * (a + i * hx)) *
                          cos(PI * (a + j * hy) / 2); // Initialisez f correctement
              u[ij] = 0; // Condition initiale : u(x) = 0 sur le bord
              u_prev[ij] = 0;
          }
      }
  // Enregistrer le temps de début
      clock_t start, end;
      start = clock(); 

      while (iter < max_iter) {
          sum_diff = 0;
          sum_uij = 0;

          for (int i = 1; i < m + 1; i++) {
              for (int j = 1; j < n + 1; j++) {
                  int ij = index(i, j, m);
                  int imj = index(i - 1, j, m);
                  int ipj = index(i + 1, j, m);
                  int ijm = index(i, j - 1, m);
                  int ijp = index(i, j + 1, m);

                  double exact_solution = sol_exact(i, j, hx, hy);
                  error = fabs(u[ij] - exact_solution);

                  double u_new = (f[ij] + hx2 * (u_prev[imj] + u_prev[ipj]) +
                                  hy2 * (u_prev[ijm] + u_prev[ijp])) /
                                  lambda;
                  u_prev[ij] = u[ij];
                  u[ij] = u_new;

                  sum_diff += square(u[ij] - u_prev[ij]);
                  sum_uij += square(u[ij]);
              }
          }

          if (sum_diff < test_arret * sum_uij) break; // Arrêter la boucle si la condition de convergence est vérifiée

          iter++;
      }

      end = clock(); // Enregistrer le temps de fin
      double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
      printf("Time taken for Jacobi(%d, %d): %.6f seconds\n", m, n, time_taken);

    // Affichage des résultats
      printf("Nombre d'itérations: %d\n", iter);
      printf("Erreur: %f\n", error);
     
//vider la mémoire 
      free(u);
      free(u_prev);
      free(f);
 
 
 

}
void Gauss(int m, int n) {
  //coeficients en espace
    double hx = (b - a) / (m + 1.0);
    double hy = (b - a) / (n + 1.0);
    double lambda = alpha + 2 * nu * ((1 / (hx * hx)) + (1 / (hy * hy)));
    double hx2 = nu / (hx * hx);
    double hy2 = nu / (hy * hy);
  //Tolérance et test d'arret
    int max_iter = 10 * (m + 1);
    double test_arret = 1e-10; // tolérance pour l'arrêt de la boucle
    double sum_diff, sum_uij;
    int iter = 0;
    double error;
    
    // Tableaux
    double *u, *u_prev, *f;

    // Allocation de mémoire
    u = (double *)malloc((m + 2) * (n + 2) * sizeof(double));
    u_prev = (double *)malloc((m + 2) * (n + 2) * sizeof(double));
    f = (double *)malloc((m + 2) * (n + 2) * sizeof(double));

    // Initialisation de u et u_prev à 0
    for (int i = 0; i < m + 2; i++) {
        for (int j = 0; j < n + 2; j++) {
            int ij = index(i, j, m);
            f[ij] = square(17 * PI / 4) * (a + i * hx) * (a + j * hy) *
                        sin(2 * PI * (a + i * hx)) * cos(PI * (a + j * hy) / 2) +
                    PI * (a + i * hx) * sin(2 * PI * (a + i * hx)) *
                        sin(PI * (a + j * hy) / 2) -
                    4 * PI * (a + j * hy) * cos(2 * PI * (a + i * hx)) *
                        cos(PI * (a + j * hy) / 2) +
                    (a + i * hx) * (a + j * hy) * sin(2 * PI * (a + i * hx)) *
                        cos(PI * (a + j * hy) / 2); // Initialisez f correctement
            u[ij] = 0; // Condition initiale : u(x) = 0 sur le bord
            u_prev[ij] = 0;
        }
    }
  // Enregistrer le temps de début
    clock_t start, end;
    start = clock(); 
    while (iter < max_iter) {
        sum_diff = 0;
        sum_uij = 0;

        for (int i = 1; i < m + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                int ij = index(i, j, m);
                int imj = index(i - 1, j, m);
                int ipj = index(i + 1, j, m);
                int ijm = index(i, j - 1, m);
                int ijp = index(i, j + 1, m);

                double exact_solution = sol_exact(i, j, hx, hy);
                error = fabs(u[ij] - exact_solution);

                double u_new = (f[ij] + hx2 * (u[ipj] + u_prev[imj]) + hy2 * (u[ijp] + u_prev[ijm])) / lambda;
                u_prev[ij] = u[ij];
                u[ij] = u_new;

                sum_diff += square(u[ij] - u_prev[ij]);
                sum_uij += square(u[ij]);
            }
        }

        if (sum_diff < test_arret * sum_uij) break; // Arrêter la boucle si la condition de convergence est vérifiée

        iter++;
    }
    end = clock(); // Enregistrer le temps de fin
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for Gauss-Seidel(%d, %d): %.6f seconds\n", m, n, time_taken);
    // Affichage des résultats
    printf("Nombre d'itérations: %d\n", iter);
    printf("Erreur: %f\n", error);

    // Libération de la mémoire
    free(u);
    free(u_prev);
    free(f);
}

//méthode d'euler implicite 
void Euler_implicite(int M,int N, int use) {
    // Coefficients
    double ht = 1.e-3;
    double alpha = 1/ht, nu = 1.0;
    double a = 3.14, b = 3.14;
    double lambda;
    // Domaine spatial
    double ax = 0.0, bx = 1.0;
    double ay = 0.0, by = 1.0;
    double hx = (bx-ax)/(double)(M+1.0), hy = (by-ay)/(double)(N+1.0);

    // travail
    int iter_Jmax = 10*M*N;
    int i, j, ij, imj, ipj, ijm, ijp;
    int iter_t, iter_j;
    double s, t;
    double hx2, hy2, w0, x, y;
    double Tmax = 1000;

    // Tolerance
    double tol_j = 1.e-5, tol_t = 1.e-8;
    double err_j, err_t;
    // tableaux
    double *u, *uu, *f;
    double *v, *vv, *ff;


    // Allocation tableaux
    u  = (double*)malloc((M+2)*(N+2)*sizeof(double));
    uu = (double*)malloc((M+2)*(N+2)*sizeof(double));
    v  = (double*)malloc((M+2)*(N+2)*sizeof(double));
    vv = (double*)malloc((M+2)*(N+2)*sizeof(double));
    f = (double*)malloc((M+2)*(N+2)*sizeof(double));
    ff = (double*)malloc((M+2)*(N+2)*sizeof(double));

    // Initialisation u et f
    for (i = 0; i < M+2; i++)
    for (j = 0; j < N+2; j++){
        x = ax + (double)i*hx;
        y = ay +(double)j*hy;
        w0 = sin(a*x)*cos(b*y);
        ij = i*(N+2) + j;
        f[ij] = (a*a+b*b)*w0;
        u[ij] = 0.0;
        if ( i == 0 || i == M+1 || j == 0 || j == N+1) u[ij] = w0;
        uu[ij] = 0.0;
    }
    lambda = alpha+2.0/hx/hx+2.0/hy/hy;
    hx2 = nu/hx/hx;
    hy2 = nu/hy/hy;

  // Boucle en temps
  clock_t start, end;
  start = clock(); // Enregistrer le temps de début
  iter_t = 0;
  t = 0.0;
  err_t = 1.0;
  while (err_t > tol_t && t < Tmax){
      iter_t++;
      t = (double)iter_t*ht;

      // Nouveau 2nd membre (Euler implicite)
      for (i = 0; i < M+2; i++)
      for (j = 0; j < N+2; j++){
          ij = i*(N+2)+j;
          ff[ij] = f[ij] + alpha * u[ij];
          vv[ij] = u[ij];
          v[ij] = u[ij];
      }

      iter_j = 0;
      err_j = 1.0;
      // Boucle de Jacobi
    if (use==1) {
        while (err_j > tol_j && iter_j < iter_Jmax){
            iter_j++;
  
            // Mise à jour Jacobi
            for (i = 1; i < M+1; i++)
            for (j = 1; j < N+1; j++){
                ij = i*(N+2)+j;
                imj = (i-1)*(N+2)+j; ipj = (i+1)*(N+2)+j;
                ijm = i*(N+2)+j-1; ijp = i*(N+2)+j+1;
                v[ij] =(ff[ij] + hx2*(vv[imj]+vv[ipj]) + hy2*(vv[ijm]+vv[ijp]))/lambda;
            }

            // Test d'arrêt Jacobi & vv <- v
            err_j = 0;
            s = 0.0;
            for (i = 1; i < M+1; i++)
            for (j = 1; j < N+1; j++){
                ij = i*(N+2)+j;
                err_j += (v[ij]-vv[ij])*(v[ij]-vv[ij]);
                s += v[ij]*v[ij];
                vv[ij] = v[ij];
            }
            err_j = sqrt(err_j/s);
        }
  
        // Récupération de la solution de Jacobi
        for (i = 1; i < M+1; i++)
        for (j = 1; j < N+1; j++){
            ij = i*(N+2)+j;
            u[ij] = v[ij];
        }
    }
    //boucle si le choix de l'utiliateur est la méthode de Gauss-Seidel
    if (use==2) {
      while (err_j > tol_j && iter_j < iter_Jmax){
              iter_j++;

              // Mise à jour Jacobi
              for (i = 1; i < M+1; i++)
              for (j = 1; j < N+1; j++){
                  ij = i*(N+2)+j;
                  imj = (i-1)*(N+2)+j; ipj = (i+1)*(N+2)+j;
                  ijm = i*(N+2)+j-1; ijp = i*(N+2)+j+1;
                  v[ij] =(ff[ij] + hx2*(v[imj] + vv[ipj]) + hy2*(v[ijm] + vv[ijp])) / lambda;
              }

              // Test d'arrêt Jacobi & vv <- v
              err_j = 0;
              s = 0.0;
              for (i = 1; i < M+1; i++)
              for (j = 1; j < N+1; j++){
                  ij = i*(N+2)+j;
                  err_j += (v[ij]-vv[ij])*(v[ij]-vv[ij]);
                  s += v[ij]*v[ij];
                  vv[ij] = v[ij];
              }
              err_j = sqrt(err_j/s);
          }

          // Récupération de la solution de Jacobi
          for (i = 1; i < M+1; i++)
          for (j = 1; j < N+1; j++){
              ij = i*(N+2)+j;
              u[ij] = v[ij];
          }
      }
      
    
      // Test solution stationnaire et uu <- u
      err_t = 0.0;
      s = 0.0;
      for (i = 1; i < M+1; i++)
      for (j = 1; j < N+1; j++){
          ij = i*(N+2)+j;
          err_t += (u[ij]-uu[ij])*(u[ij]-uu[ij]);
          s += u[ij]*u[ij];
          uu[ij] = u[ij];
      }
      err_t = sqrt(err_t/s);
  }
  const char *L[] = {"Jacobi", "Gauss-Seidel"}; // Tableau de chaînes de caractères

  end = clock(); // Enregistrer le temps de fin
  double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Time taken for Euler_implicite(%d, %d) and used %s %.6f seconds\n", M, N, L[use - 1], time_taken);

    printf("iter=%d err=%15.8e \n",iter_t,err_t);
}



void Crank_Nicolson(int M, int N,int use) {
    // Coefficients
    double ht = 1.e-3;
    double alpha = 1.0 / ht, nu = 1.0;
    double a = 3.14, b = 3.14;
    double lambda;
    // Domaine spatial
    double ax = 0.0, bx = 1.0;
    double ay = 0.0, by = 1.0;
    double hx = (bx - ax) / (double)M, hy = (by - ay) / (double)N;

    // travail
    int iter_Jmax = 10 * M * N;
    int i, j, ij, imj, ipj, ijm, ijp;
    int iter_t, iter_j;
    double s, t;
    double hx2, hy2, w0, x, y;
    double Tmax = 1000;

    // Tolerance
    double tol_j = 1.e-5, tol_t = 1.e-8;
    double err_j, err_t;
    // tableaux
    double *u, *uu, *f;
    double *v, *vv, *ff;

    // Allocation tableaux
    u = (double*)malloc((M + 2) * (N + 2) * sizeof(double));
    uu = (double*)malloc((M + 2) * (N + 2) * sizeof(double));
    v = (double*)malloc((M + 2) * (N + 2) * sizeof(double));
    vv = (double*)malloc((M + 2) * (N + 2) * sizeof(double));
    f = (double*)malloc((M + 2) * (N + 2) * sizeof(double));
    ff = (double*)malloc((M + 2) * (N + 2) * sizeof(double));

    // Initialisation u et f
    for (i = 0; i < M + 2; i++)
        for (j = 0; j < N + 2; j++) {
            x = ax + (double)i * hx;
            y = ay + (double)j * hy;
            w0 = sin(a * x) * cos(b * y);
            ij = i * (N + 2) + j;
            f[ij] = (a * a + b * b) * w0;
            u[ij] = 0.0;
            if (i == 0 || i == M + 1 || j == 0 || j == N + 1) u[ij] = w0;
            uu[ij] = 0.0;
        }
    lambda = alpha + 2.0 / hx / hx + 2.0 / hy / hy;
    hx2 = nu / hx / hx;
    hy2 = nu / hy / hy;
    // Boucle en temps
    clock_t start, end;
    start = clock(); // Enregistrer le temps de début
    iter_t = 0;
    t = 0.0;
    err_t = 1.0;
    while (err_t > tol_t && t < Tmax) {
        iter_t++;
        t = (double)iter_t * ht;

        // Nouveau 2nd membre (Crank-Nicolson)
        for (i = 0; i < M + 2; i++)
            for (j = 0; j < N + 2; j++) {
                ij = i * (N + 2) + j;
                ff[ij] = f[ij] + 0.5 * (u[ij] + uu[ij]) / ht;
                vv[ij] = u[ij];
                v[ij] = u[ij];
            }

       
      if (use==1) {
        iter_j = 0;
        err_j = 1.0;
        // Boucle de Jacobi
          while (err_j > tol_j && iter_j < iter_Jmax){
              iter_j++;

              // Mise à jour Jacobi
              for (i = 1; i < M+1; i++)
              for (j = 1; j < N+1; j++){
                  ij = i*(N+2)+j;
                  imj = (i-1)*(N+2)+j; ipj = (i+1)*(N+2)+j;
                  ijm = i*(N+2)+j-1; ijp = i*(N+2)+j+1;
                  v[ij] =(ff[ij] + hx2*(vv[imj]+vv[ipj]) + hy2*(vv[ijm]+vv[ijp]))/lambda;
              }

              // Test d'arrêt Jacobi & vv <- v
              err_j = 0;
              s = 0.0;
              for (i = 1; i < M+1; i++)
              for (j = 1; j < N+1; j++){
                  ij = i*(N+2)+j;
                  err_j += (v[ij]-vv[ij])*(v[ij]-vv[ij]);
                  s += v[ij]*v[ij];
                  vv[ij] = v[ij];
              }
              err_j = sqrt(err_j/s);
          }

          // Récupération de la solution de Jacobi
          for (i = 1; i < M+1; i++)
          for (j = 1; j < N+1; j++){
              ij = i*(N+2)+j;
              u[ij] = v[ij];
          }
      }
      //boucle si le choix de l'utiliateur est la méthode de Gauss-Seidel
      if (use==2) {
        iter_j = 0;
        err_j = 1.0;
        // Boucle de Jacobi
        while (err_j > tol_j && iter_j < iter_Jmax){
                iter_j++;

                // Mise à jour Gauss-Seidel
                for (i = 1; i < M+1; i++)
                for (j = 1; j < N+1; j++){
                    ij = i*(N+2)+j;
                    imj = (i-1)*(N+2)+j; ipj = (i+1)*(N+2)+j;
                    ijm = i*(N+2)+j-1; ijp = i*(N+2)+j+1;
                    v[ij] =(ff[ij] + hx2*(v[imj] + vv[ipj]) + hy2*(v[ijm] + vv[ijp])) / lambda;
                }
          // Test d'arrêt Gauss Seidel  & vv <- v
              err_j = 0;
              s = 0.0;
              for (i = 1; i < M+1; i++)
              for (j = 1; j < N+1; j++){
                  ij = i*(N+2)+j;
                  err_j += (v[ij]-vv[ij])*(v[ij]-vv[ij]);
                  s += v[ij]*v[ij];
                  vv[ij] = v[ij];
              }
              err_j = sqrt(err_j/s);
          }
        // Récupération de la solution de Gauss Seidel 
        for (i = 1; i < M+1; i++)
        for (j = 1; j < N+1; j++){
            ij = i*(N+2)+j;
            u[ij] = v[ij];
        }
      }

        // Test solution stationnaire et uu <- u
        err_t = 0.0;
        s = 0.0;
        for (i = 1; i < M + 1; i++)
            for (j = 1; j < N + 1; j++) {
                ij = i * (N + 2) + j;
                err_t += (u[ij] - uu[ij]) * (u[ij] - uu[ij]);
                s += u[ij] * u[ij];
                uu[ij] = u[ij];
            }
        err_t = sqrt(err_t / s);
    }
    const char *L[] = {"Jacobi", "Gauss-Seidel"}; // Tableau de chaînes de caractères

    end = clock(); // Enregistrer le temps de fin
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for Crank-Nickolson(%d, %d) and used %s %.6f seconds\n", M, N, L[use - 1], time_taken);
    
    printf("iter=%d  err=%15.8e \n", iter_t, err_t);

    // Libération de la mémoire
    free(u);
    free(uu);
    free(v);
    free(vv);
    free(f);
    free(ff);
}


