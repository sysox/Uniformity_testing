/*
 The C program to compute Kolmogorov's distribution

             K(n,d) = Prob(D_n < d),         where

      D_n = max(x_1-0/n,x_2-1/n...,x_n-(n-1)/n,1/n-x_1,2/n-x_2,...,n/n-x_n)

    with  x_1<x_2,...<x_n  a purported set of n independent uniform [0,1)
    random variables sorted into increasing order.
    See G. Marsaglia, Wai Wan Tsang and Jingbo Wong, J.Stat.Software.
    https://www.jstatsoft.org/article/view/v008i18
*/
#include "k.h"

    /* prototypes*/
double K(int n,double d) ;                        //Kolmogorov distribution
static void mMultiply(double *A,double *B,double *C,int m) ;      //Matrix product
static void mPower(double *A,int eA,double *V,int *eV,int m,int n); //Matrix power

double D_n(double *pvalues, int n); //KS statistic with expected pdf = [0,1] uniformity
double KS_stat_to_pval(int n,double d, int tails);
double KS_pvalue(double *pvalues, int n, int tails);

double K(int n,double d)
{
   int k,m,i,j,g,eH,eQ;
   double h,s,*H,*Q;
    /* OMIT NEXT TWO LINES IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL*/
s=d*d*n;
if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
   k=(int)(n*d)+1;
   m=2*k-1;
   h=k-n*d;
   H=(double*)malloc((m*m)*sizeof(double));
   Q=(double*)malloc((m*m)*sizeof(double));
       for(i=0;i<m;i++)
         for(j=0;j<m;j++)
           if(i-j+1<0) H[i*m+j]=0;
          else     H[i*m+j]=1;
    for(i=0;i<m;i++)
    {
    H[i*m]-=pow(h,i+1);
    H[(m-1)*m+i]-=pow(h,(m-i));
    }
    H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
    for(i=0;i<m;i++)
    for(j=0;j<m;j++)
    if(i-j+1>0)
        for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
    eH=0;
    mPower(H,eH,Q,&eQ,m,n);
    s=Q[(k-1)*m+k-1];
    for(i=1;i<=n;i++)
    {
    s=s*i/n;
    if(s<1e-140){s*=1e140; eQ-=140;}
    }
    s*=pow(10.,eQ);
    free(H);
    free(Q);
    return s;
}
static void mMultiply(double *A,double *B,double *C,int m)
{
    int i,j,k;
    double s;
    for(i=0;i<m;i++)
        for(j=0;j<m;j++)
        {
            s=0.;
            for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j];
            C[i*m+j]=s;
        }
}
static void mPower(double *A,int eA,double *V,int *eV,int m,int n)
{
    double *B;int eB,i;
    if(n==1)
    {
        for(i=0;i<m*m;i++) V[i]=A[i];
        *eV=eA;
        return;
    }
    mPower(A,eA,V,eV,m,n/2);
    B=(double*)malloc((m*m)*sizeof(double));
    mMultiply(V,V,B,m);
    eB=2*(*eV);
    if(n%2==0) {for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
    else   {mMultiply(A,B,V,m);*eV=eA+eB;}
    if(V[(m/2)*m+(m/2)]>1e140)
    {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140; *eV+=140;}
    free(B);
}


//New stuff
static int compare_function(const void *a,const void *b) {
    double *x = (double *) a;
    double *y = (double *) b;
    if ((double)*x < (double)*y) return -1;
    else{
        if ((double)*x > (double)*y) return 1;
        else return 0;
    }
}
double D_n(double *pvalues, int n){
    int i;
    double d_minus, d_plus, d_max, step;

    qsort (pvalues, n, sizeof(double ), compare_function);

    d_minus = d_plus = 0.0;
    step = 1.0 / (n);

    for (i = 0; i < n; i++) {
        d_minus = fmax(d_minus, fabs(pvalues[i] - i*step));
        d_plus = fmax(d_plus, fabs(pvalues[i] - (i+1)*step));
    }
    return fmax(d_minus, d_plus);
}


double KS_stat_to_pval(int n,double d, int tails){
    double cdf;

    cdf = K(n,d);
    switch(tails){
        case -1:
            return cdf;
            break;
        case 1:
            return 1-cdf;
            break;
        case 2:
            return 2*fmin(cdf, 1-cdf);
            break;
    }
}
double KS_tails(double *pvalues, int n, int tails){
    double d, cdf;
    d = D_n(pvalues, n);
    return KS_stat_to_pval(n, d, tails);
}

double KS_both(double *pvalues, int n){
    return KS_tails(pvalues, n, 2);
}

double KS_left(double *pvalues, int n){
    return KS_tails(pvalues, n, -1);
}

double KS_right(double *pvalues, int n){
    return KS_tails(pvalues, n, 1);
}
