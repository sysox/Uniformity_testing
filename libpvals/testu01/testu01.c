#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "testu01.h"
#include "gofs.h"
#include "gofw.h"
#include "fbar.h"
#include "tables.h"

typedef double (*gofs_CFUNC) (double [], long int);
typedef double (*fbar_DFUNC) (long int, double);

void debug_array(const double *V, unsigned long long N, const char *s)
{
    printf("%s ----\n", s);
    for (int i = 1; i <= N; i++)
        printf("%f\n", V[i]);
    printf("----\n");
}

/* Copy and shift, TestU01 do not use value at index 0, values are stored 1..N */
static void tables_CopyTabD_1 (const double T1[], double T2[], int n1, int n2)
{
    for (int i = n1; i <= n2; i++)
        T2[i] = T1[i-1];
}

static double generic_pvalue(const double *P, unsigned long long N,
                             gofs_CFUNC gofs_F, fbar_DFUNC fbar_F)
{
    double d, p;
    double *T;

    T = calloc(N+1, sizeof(double));
    if (!T) {
        printf("testu01: out of memory\n");
        return nan("NAN");
    }
    tables_CopyTabD_1(P, T, 1, N);
    tables_QuickSortD(T, 1, N);
    d = gofs_F(T, N);
    p = fbar_F(N, d);

    free(T);
    return p;
}

static double ks_pvalue(const double *P, unsigned long long N,
                             double *pplus, double *pminus)
{
    double p, DP, DM, D;
    double *T;

    T = calloc(N+1, sizeof(double));
    if (!T) {
        printf("testu01: out of memory\n");
        return nan("NAN");
    }
    tables_CopyTabD_1(P, T, 1, N);
    tables_QuickSortD(T, 1, N);
    gofs_KS(T, N, &DP, &DM, &D);

    if (pplus)
        *pplus  = fbar_KSPlus(N, DP);
    if (pminus)
        *pminus = fbar_KSPlus(N, DM);
    p = fbar_KSPlus(N, D);

    free(T);
    return p;
}

/* m-NP */
double testu01_pvalue_snpair_ClosePairs(const double *P, unsigned long long N)
{
    return generic_pvalue(P, N, gofs_AndersonDarling, fbar_AndersonDarling);
}

/* AD A2 */
double testu01_pvalue_sknuth_MaxOft(const double *P, unsigned long long N)
{
    double p;
    double *T;
    gofw_TestArray sVal;

    T = calloc(N+1, sizeof(double));
    if (!T) {
        printf("testu01: out of memory\n");
        return nan("NAN");
    }

    tables_CopyTabD_1(P, T, 1, N);
    tables_QuickSortD(T, 1, N);
    /* optimization - overwriting T[] in-place */
    gofs_ContUnifTransform (T, N, wdist_Unif, NULL, T);
    gofw_Tests0 (T, N, sVal);
    p = fbar_AndersonDarling (N, sVal[gofw_AD]);

    free(T);
    return p;
}

double testu01_pvalue_ksp(const double *P, unsigned long long N)
{
    double p, pminus;
    p = ks_pvalue(P, N, NULL, &pminus);
    return isnan(p) ? p : pminus;
}
double testu01_pvalue_ksm(const double *P, unsigned long long N)
{
    double p, pplus;
    p = ks_pvalue(P, N, &pplus, NULL);
    return isnan(p) ? p : pplus;
}

double testu01_pvalue_ks(const double *P, unsigned long long N)
{
    return ks_pvalue(P, N, NULL, NULL);
}

double testu01_pvalue_ad(const double *P, unsigned long long N)
{
    return generic_pvalue(P, N, gofs_AndersonDarling, fbar_AndersonDarling);
}

double testu01_pvalue_cm(const double *P, unsigned long long N)
{
    return generic_pvalue(P, N, gofs_CramerMises, fbar_CramerMises);
}

double testu01_pvalue_wg(const double *P, unsigned long long N)
{
    return generic_pvalue(P, N, gofs_WatsonG, fbar_WatsonG);
}

double testu01_pvalue_wu(const double *P, unsigned long long N)
{
    return generic_pvalue(P, N, gofs_WatsonU, fbar_WatsonU);
}

/*
      pVal[gofw_KSP] = fbar_KSPlus (N, sVal[gofw_KSP]);
      pVal[gofw_KSM] = fbar_KSPlus (N, sVal[gofw_KSM]);
      pVal[gofw_KS] = fbar_KS1 (N, sVal[gofw_KS]);
      pVal[gofw_AD] = fbar_AndersonDarling (N, sVal[gofw_AD]);
      pVal[gofw_CM] = fbar_CramerMises (N, sVal[gofw_CM]);
      pVal[gofw_WG] = fbar_WatsonG (N, sVal[gofw_WG]);
      pVal[gofw_WU] = fbar_WatsonU (N, sVal[gofw_WU]);

      wdist_Normal, wdist_ChiSquare, FDistProd, wdist_Unif, FDistMeans
*/
