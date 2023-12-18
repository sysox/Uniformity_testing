//
// Created by user on 17/12/2023.
//
#ifndef UNIFORMITY_TESTING_UTILS_H
#define UNIFORMITY_TESTING_UTILS_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "libpvals/libpvals.h"

void seed_xorshift32(uint32_t seed);
uint32_t xorshift32();

//file must ends with '\n'
void read_pvals(const char* filename, double** pvalues, int *num_pvalues);
void random_sample(const double* values, int num_values, double* sample, int sample_size);
void GoF_pvals(const char* src_file, int sample_size, int repetitions, int GoF_idx, double* resulted_pvals);

int compare_function(const void *a,const void *b);
void hist(int num_pvals, double* pvalues, double* bin_edges, int* bins);

#endif //UNIFORMITY_TESTING_UTILS_H
