//
// Created by user on 17/12/2023.
//

#include "libpvals/utils.h"
#include <stdlib.h>
#include "libpvals/libpvals.h"

int main(){
    int i, j, num_pvals, sample_size;
    double *pvalues, *sample, GoF_pval;

    sample_size = 10;
    sample = (double *) malloc(sample_size*sizeof(double));
    read_pvals("uniform_pvals_devurand.pval", &pvalues,&num_pvals);
    random_sample(pvalues, num_pvals, sample, sample_size);

    const double *tst = sample;
    unsigned long long tst_num = (unsigned long long)sample_size;

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

    GoF_pvals( "uniform_pvals_devurand.pval", sample_size, 1, 0, &GoF_pval );

    printf("%lf", GoF_pval);


    return 0;
}