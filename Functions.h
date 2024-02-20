#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Déclaration des fonctions
int index(int i, int j, int m);
double square(double x);
double sol_exact(int i, int j, double hx, double hy);
void Jacobi(int m,int n);
void Gauss(int m,int n);
void  Euler_implicite(int m,int n,int test);
void Crank_Nicolson(int m,int n,int i);
#endif
