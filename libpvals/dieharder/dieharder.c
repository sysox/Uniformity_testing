#include <stdint.h>

#include "dieharder.h"

#define CONST_CAST(x) (x)(uintptr_t)

double kstest(double *pvalue, int count);
double kstest_kuiper(double *pvalue,int count);

double dieharder_pvalue(const double *P, unsigned long long N)
{
    return kstest(CONST_CAST(double*)(P), (int)N);
}

double dieharder_pvalue_kuiper(const double *P, unsigned long long N)
{
    return kstest_kuiper(CONST_CAST(double*)(P), (int)N);
}
