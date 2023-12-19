//
// Created by user on 19/12/2023.
//

#include <stdio.h>
#include "k.h"

int main()
{   int  n = 10;
    double pvalues[10] = {0.12378337, 0.84600339, 0.6511558 , 0.45466759, 0.29017108,
                          0.74067633, 0.3374257 , 0.16274542, 0.76651674, 0.84278765};
    double d;

    d = D_n(pvalues, n);
    printf("KS statistic %0.10lf =\n",d);
    printf("KS cdf %0.10lf = \n",K(n,d));
    printf("KS left tail p-value  %0.10lf = \n", KS_stat_to_pval(n, d, -1));
    printf("KS right tail p-value %0.10lf = \n", KS_stat_to_pval(n, d, 1));
    printf("KS both tails pvalue %0.10lf = \n", KS_stat_to_pval(n, d, 2));

    printf("KS left tail p-value %0.10lf = \n", KS(pvalues, n, -1));
    printf("KS right tail p-value %0.10lf = \n", KS(pvalues, n, 1));
    printf("KS both tails p-value %0.10lf = \n", KS(pvalues, n, 2));

}