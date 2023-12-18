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
    int decimals[10] = {1,10,100, 1000, 10000, 100000, 1000000, 10000000, 100000000};

    //read file to memory
    fp = fopen(filename,"rb");
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
    *pvalues = file_doubles = (double *) malloc(*num_pvalues*sizeof(double));
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

            file_doubles[idx] = tmp/decimals[decimal_counter];
            idx += 1;
            tmp = 0;
        }
    }
}

void random_sample(const double* values, int num_values, double* sample, int sample_size){
    unsigned int idx;
    uint32_t i;

    for (i = 0; i < sample_size; i++) {
        idx = xorshift32() % num_values;                    // random block from
        sample[i] = values[idx];
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
