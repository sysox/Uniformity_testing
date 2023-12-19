//
// Created by user on 17/12/2023.
//

#include "libpvals/utils.h"
#include <stdlib.h>
#include "libpvals/libpvals.h"
#include "other_sources/k.h"

void test_ks(void)
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

    printf("KS left tail p-value %0.10lf = \n", KS_left(pvalues, n));
    printf("KS right tail p-value %0.10lf = \n", KS_right(pvalues, n));
    printf("KS both tails p-value %0.10lf = \n", KS_both(pvalues, n));

}

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

    printf("%lf\n", GoF_pval);

    test_ks();
    return 0;
}