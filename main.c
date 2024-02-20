#include "Functions.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    printf("*********************Jacobi*************************\n");
    Jacobi(15,15);
    Jacobi(31,31);
    Jacobi(63,63);
    printf("**********************Gauss-Seidel*********************\n");
    Gauss(15,15);
    Gauss(31,31);
    Gauss(63,63);
    printf( "*******************Euler_implicite************************\n");
    Euler_implicite(15,15,1);
    Euler_implicite(31,31,2);
    Euler_implicite(63,63,1);
    printf( "******************Cranck-Nicolson************************\n");
    Crank_Nicolson(15,15,1);
    Crank_Nicolson(31,31,2);
    Crank_Nicolson(63,63,1);
   
    
    

    return 0;
}
