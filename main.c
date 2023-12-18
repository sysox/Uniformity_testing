//
// Created by user on 17/12/2023.
//

#include "utils.h"
#include <stdlib.h>


int main(){
    int i, j, num_pvals, sample_size, bins[11];
    double *pvalues, *sample, second_lvl;
    double bin_edges[12] = {0,0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8,0.9, 1.0, 1.1};

    sample_size = 10000;
    sample = (double *) malloc(sample_size*sizeof(double));
    read_pvals("../uniform_devurand.pval", &pvalues,&num_pvals);
//    hist(num_pvals, pvalues, 12, bin_edges, bins);
//    for(j = 0; j < 12; j++){
//        printf("%i\n", bins[j]);
//    }

    for(i = 0; i < 10; i++){
        random_sample(pvalues, num_pvals, sample, sample_size);
//        second_lvl = KS(); // Tu doplnit prosim

        printf("%lf", second_lvl);
    }

    return 0;
}