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
double KS_tails(double *pvalues, int n, int tails);

double KS_both(double *pvalues, int n);
double KS_left(double *pvalues, int n);
double KS_right(double *pvalues, int n);

#endif //UNIFORMITY_TESTING_K_H
