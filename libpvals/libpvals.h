#ifndef __LIBPVALS_H__
#define __LIBPVALS_H__

#include "dieharder/dieharder.h"
#include "nist-sts/nist-sts.h"
#include "testu01/testu01.h"

/* Kolmogorov's distribution, from Marsaglia */
#include "k.h"

/* exported utils */
void GoF_pvals(const char* src_file, int sample_size, int repetitions, int GoF_idx, double* resulted_pvals);

#endif