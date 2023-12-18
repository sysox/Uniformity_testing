// extracted from assess.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "../include/decls.h"
//#include "../include/cephes.h"  
//#include "../include/utilities.h"

#define ALPHA 0.01 /* SIGNIFICANCE LEVEL */

#define TEST_FREQUENCY          1
#define TEST_BLOCK_FREQUENCY    2
#define TEST_CUSUM              3
#define TEST_RUNS               4
#define TEST_LONGEST_RUN        5
#define TEST_RANK               6
#define TEST_FFT                7
#define TEST_NONPERIODIC        8
#define TEST_OVERLAPPING        9
#define TEST_UNIVERSAL         10
#define TEST_APEN              11
#define TEST_RND_EXCURSION     12
#define TEST_RND_EXCURSION_VAR 13
#define TEST_SERIAL            14
#define TEST_LINEARCOMPLEXITY  15

#include "nist-sts-cephes.h"
#include "nist-sts.h"

const char *testNames[16] = {
    " ",
    "Frequency",
    "BlockFrequency",
    "CumulativeSums",
    "Runs",
    "LongestRun",
    "Rank",
    "FFT",
    "NonOverlappingTemplate",
    "OverlappingTemplate",
    "Universal",
    "ApproximateEntropy",
    "RandomExcursions",
    "RandomExcursionsVariant",
    "Serial",
    "LinearComplexity"
};

static int
cmp(const double *a, const double *b)
{
	if ( *a < *b )
		return -1;
	if ( *a == *b )
		return 0;
	return 1;
}

int nist_compute(const double *P, int numOfBitStreams, int test, double *pval)
{
	int		j, pos, count, passCount, sampleSize, expCount, proportion_threshold_min, proportion_threshold_max;
	int		freqPerBin[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double	*A, *T, chi2, uniformity, p_hat;
	// FIXME: why they lose precission here, it should be double?
	float	c;

	if ( (A = (double *)calloc(numOfBitStreams, sizeof(double))) == NULL ) {
		printf("Final Analysis Report aborted due to insufficient workspace\n");
		return 0;
	}

	/* Compute Metric 1: Proportion of Passing Sequences */

	count = 0;
	sampleSize = numOfBitStreams;

	if ( (test == TEST_RND_EXCURSION) || (test == TEST_RND_EXCURSION_VAR) ) { /* Special Case: Random Excursion Tests */
		if ( (T = (double *)calloc(numOfBitStreams, sizeof(double))) == NULL ) {
			printf("Final Analysis Report aborted due to insufficient workspace\n");
			return 0;
		}
		for ( j=0; j<sampleSize; j++ ) {
			//fscanf(fp, "%f", &c);
			c = P[j];
			if ( c > 0.000000 )
				T[count++] = c;
		}
		
		if ( (A = (double *)calloc(count, sizeof(double))) == NULL ) {
			printf("Final Analysis Report aborted due to insufficient workspace\n");
			return 0;
		}
		
		for ( j=0; j<count; j++ )
			A[j] = T[j];
		
		sampleSize = count;
		count = 0;
		for ( j=0; j<sampleSize; j++ )
			if ( A[j] < ALPHA )
				count++;
		free(T);
	}
	else {
		if ( (A = (double *)calloc(sampleSize, sizeof(double))) == NULL ) {
			printf("Final Analysis Report aborted due to insufficient workspace\n");
			return 0;
		}
		for ( j=0; j<sampleSize; j++ ) {
			//fscanf(fp, "%f", &c);
			c = P[j];
			if ( c < ALPHA )
				count++;
			A[j] = c;
		}
	}

	if ( sampleSize == 0 )
		passCount = 0;
	else
		passCount = sampleSize - count;
	
	p_hat = 1.0 - ALPHA;
	proportion_threshold_max = (p_hat + 3.0 * sqrt((p_hat*ALPHA)/sampleSize)) * sampleSize;
	proportion_threshold_min = (p_hat - 3.0 * sqrt((p_hat*ALPHA)/sampleSize)) * sampleSize;

	/* Compute Metric 2: Histogram */
	
	qsort((void *)A, sampleSize, sizeof(double), (void *)cmp);
	for ( j=0; j<sampleSize; j++ ) {
		pos = (int)floor(A[j]*10);
		if ( pos == 10 )
			pos--;
		freqPerBin[pos]++;
	}
	chi2 = 0.0;
	expCount = sampleSize/10;
	if ( expCount == 0 )
		uniformity = 0.0;
	else {
		for ( j=0; j<10; j++ )
			chi2 += pow(freqPerBin[j]-expCount, 2)/expCount;
		uniformity = cephes_igamc(9.0/2.0, chi2/2.0);
	}

#if 0 /* Disable NIST STS default stat output */
	FILE	*summary = stdout;
	for ( j=0; j<10; j++ )			/* DISPLAY RESULTS */
		fprintf(summary, "%3d ", freqPerBin[j]);
	
	if ( expCount == 0 )
		fprintf(summary, "    ----    ");
	else if ( uniformity < 0.0001 )
		fprintf(summary, " %8.6f * ", uniformity);
	else
		fprintf(summary, " %8.6f   ", uniformity);
	
	if ( sampleSize == 0 )
		fprintf(summary, " ------     %s\n", testNames[test]);
	//	else if ( proportion < 0.96 )
	else if ( (passCount < proportion_threshold_min) || (passCount > proportion_threshold_max))
		fprintf(summary, "%4d/%-4d *  %s\n", passCount, sampleSize, testNames[test]);
	else
		fprintf(summary, "%4d/%-4d    %s\n", passCount, sampleSize, testNames[test]);
#endif
	//fclose(fp);
	free(A);

	if (pval)
		*pval = uniformity;

	/* Silence warning for not yet unused variables */
#define UNUSED(x) (void)(x)
	UNUSED(proportion_threshold_max);
	UNUSED(proportion_threshold_min);
	UNUSED(passCount);

	return sampleSize;
}

double nist_pvalue(const double *P, unsigned long long N)
{
	double p  = 0;
	nist_compute(P, (int)N, 0, &p);
	return p;
}
