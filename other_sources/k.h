//
// Created by user on 19/12/2023.
//

#ifndef UNIFORMITY_TESTING_K_H
#define UNIFORMITY_TESTING_K_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double K(int n,double d) ;                        //Kolmogorov distribution
double D_n(double *pvalues, int n);                 // Kolmogorov


double KS_stat_to_pval(int n,double d, int tails);
double KS(double *pvalues, int n, int tails);

#endif //UNIFORMITY_TESTING_K_H
