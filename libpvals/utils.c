//
// Created by user on 17/12/2023.
//

#include "utils.h"

uint32_t xorshift32_state = 1;
void seed_xorshift32(uint32_t seed){
    time_t t;
    if (seed == 0){
        seed = time(NULL);
    }
    xorshift32_state = seed;
}
uint32_t xorshift32()
{
    /* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
    xorshift32_state ^= xorshift32_state << 13;
    xorshift32_state ^= xorshift32_state >> 17;
    xorshift32_state ^= xorshift32_state << 5;
    return xorshift32_state;
}

//file must ends with '\n'
void read_pvals(const char* filename, double** pvalues, int *num_pvalues){
    FILE *fp;
    int size, bytes_read, i, idx, decimal_counter;
    unsigned char *tmp_array, *ptr;
    double tmp, *file_doubles;
    double decimals[20] = {1};

    for(i = 1; i < 20; i++)
    {
        decimals[i] = 10*decimals[i-1];
    }


    //read file to memory
    fp = fopen(filename,"rb");
    if (fp == NULL){
        printf("%s filename not opened", filename);
        return;
    }
    fseek(fp, 0L, SEEK_END);
    size = ftell(fp);
    fseek( fp, 0, SEEK_SET );
    tmp_array = (unsigned char *) malloc(size);
    fread(tmp_array, 1, size, fp);

    //counting number of doubles
    ptr = tmp_array;
    *num_pvalues = 0;
    for(i = 0; i < size; i++) {
        if (tmp_array[i] == '\n') {
            *num_pvalues += 1;
        }
    }

    //parsing doubles
    *pvalues =  (double *) malloc(*num_pvalues*sizeof(double));
    tmp = idx = 0;
    for(i = 0; i < size; i++) {

        if (tmp_array[i] == '.'){
            decimal_counter = 0;
            continue;
        }
        if (tmp_array[i] != '\n'){
            tmp = tmp*10 + (tmp_array[i] - '0');
            decimal_counter += 1;
        }
        else{
            (*pvalues)[idx] = tmp/decimals[decimal_counter];
            idx += 1;
            tmp = 0;
        }
    }
}

void random_sample(const double* pvalues, int num_values, double* sample, int sample_size){
    unsigned int idx;
    int i;

    for (i = 0; i < sample_size; i++) {
        idx = xorshift32() % num_values;
        // random block from
        sample[i] = pvalues[idx];
    }
}


int compare_function(const void *a,const void *b) {
    double *x = (double *) a;
    double *y = (double *) b;
    if (*x < *y) return -1;
    else{
        if (*x > *y) return 1; return 0;
    }
}




void hist(int num_pvals, double* pvalues, double* bin_edges, int* bins){
    int idx_pvals, idx_bin_edges;
//    qsort (pvalues, num_pvals, sizeof(double ), compare_function);
    if(pvalues[num_pvals - 1] > 1){
        printf("Large pvalue %0.10lf", pvalues[num_pvals - 1]);
    }
    idx_bin_edges = 0;
    bins[0] = 0;
    for(idx_pvals = 0; idx_pvals < num_pvals; idx_pvals++)
    {
        if(pvalues[idx_pvals] <= bin_edges[idx_bin_edges]){
            bins[idx_bin_edges]++;
        } else{
            while(pvalues[idx_pvals] > bin_edges[idx_bin_edges]){
                idx_bin_edges++;
                bins[idx_bin_edges] = 0;
            }
        }
    }

}




typedef double (*GoF)(const double *, unsigned long long);

void GoF_pvals(const char* src_file, int sample_size, int repetitions, int GoF_idx, double* resulted_pvals){
    int i, j, num_pvals;
    double *pvalues, *sample;
    GoF GoF_functions[] = {&dieharder_pvalue, &dieharder_pvalue_kuiper, &nist_pvalue, &testu01_pvalue_snpair_ClosePairs, &testu01_pvalue_sknuth_MaxOft,
                           &testu01_pvalue_ksp, &testu01_pvalue_ksm, &testu01_pvalue_ks, &testu01_pvalue_ad, &testu01_pvalue_cm, &testu01_pvalue_wg,
                           &testu01_pvalue_wu};


    sample = (double *) malloc(sample_size*sizeof(double));
    read_pvals(src_file, &pvalues,&num_pvals);


    const double *tst = sample;
    unsigned long long tst_num = (unsigned long long)sample_size;


    GoF gof_func = GoF_functions[GoF_idx];
    for(i = 0; i < repetitions; i++)
    {
        random_sample(pvalues, num_pvals, sample, sample_size);
        resulted_pvals[i] = gof_func(tst, tst_num);
    }

//    printf("Num pvals %llu\n", tst_num);
//    printf("DIEHARDER         %.8f\n", dieharder_pvalue(tst, tst_num));
//    printf("DIEHARDER_KUIPER  %.8f\n", dieharder_pvalue_kuiper(tst, tst_num));
//    printf("NIST              %.8f\n", nist_pvalue(tst, tst_num));
//    printf("TESTU01_ClosePair %.8f\n", testu01_pvalue_snpair_ClosePairs(tst, tst_num));
//    printf("TESTU01_MaxOft    %.8f\n", testu01_pvalue_sknuth_MaxOft(tst, tst_num));
//    printf("TESTU01_KSP       %.8f\n", testu01_pvalue_ksp(tst, tst_num));
//    printf("TESTU01_KSM       %.8f\n", testu01_pvalue_ksm(tst, tst_num));
//    printf("TESTU01_KS        %.8f\n", testu01_pvalue_ks(tst, tst_num));
//    printf("TESTU01_AD        %.8f\n", testu01_pvalue_ad(tst, tst_num));
//    printf("TESTU01_CM        %.8f\n", testu01_pvalue_cm(tst, tst_num));
//    printf("TESTU01_WG        %.8f\n", testu01_pvalue_wg(tst, tst_num));
//    printf("TESTU01_WU        %.8f\n", testu01_pvalue_wu(tst, tst_num));


}
