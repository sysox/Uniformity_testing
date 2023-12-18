#ifndef __TESTU01_H__
#define __TESTU01_H__

double testu01_pvalue_snpair_ClosePairs(const double *P, unsigned long long N);
double testu01_pvalue_sknuth_MaxOft(const double *P, unsigned long long N);

double testu01_pvalue_ksp(const double *P, unsigned long long N);
double testu01_pvalue_ksm(const double *P, unsigned long long N);
double testu01_pvalue_ks(const double *P, unsigned long long N);
double testu01_pvalue_ad(const double *P, unsigned long long N);
double testu01_pvalue_cm(const double *P, unsigned long long N);
double testu01_pvalue_wg(const double *P, unsigned long long N);
double testu01_pvalue_wu(const double *P, unsigned long long N);

#endif
