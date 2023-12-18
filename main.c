//
// Created by user on 17/12/2023.
//

#include "utils.h"
#include <stdlib.h>
#include "libpvals/libpvals.h"

int main(){
    int i, j, num_pvals, sample_size, bins[11];
    double *pvalues, *sample, second_lvl;
    double bin_edges[12] = {0,0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8,0.9, 1.0, 1.1};

    sample_size = 10000;
    sample = (double *) malloc(sample_size*sizeof(double));
    read_pvals("../uniform_devurand.pval", &pvalues,&num_pvals);

    const double *tst = pvalues;
    unsigned long long tst_num = (unsigned long long)num_pvals;
    printf("Num pvals %llu\n", tst_num);
    printf("DIEHARDER         %.8f\n", dieharder_pvalue(tst, tst_num));
    printf("DIEHARDER_KUIPER  %.8f\n", dieharder_pvalue_kuiper(tst, tst_num));
    printf("NIST              %.8f\n", nist_pvalue(tst, tst_num));
    printf("TESTU01_ClosePair %.8f\n", testu01_pvalue_snpair_ClosePairs(tst, tst_num));
    printf("TESTU01_MaxOft    %.8f\n", testu01_pvalue_sknuth_MaxOft(tst, tst_num));
    printf("TESTU01_KSP       %.8f\n", testu01_pvalue_ksp(tst, tst_num));
    printf("TESTU01_KSM       %.8f\n", testu01_pvalue_ksm(tst, tst_num));
    printf("TESTU01_KS        %.8f\n", testu01_pvalue_ks(tst, tst_num));
    printf("TESTU01_AD        %.8f\n", testu01_pvalue_ad(tst, tst_num));
    printf("TESTU01_CM        %.8f\n", testu01_pvalue_cm(tst, tst_num));
    printf("TESTU01_WG        %.8f\n", testu01_pvalue_wg(tst, tst_num));
    printf("TESTU01_WU        %.8f\n", testu01_pvalue_wu(tst, tst_num));


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