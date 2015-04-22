#include "mpuLikelihoods.h"
#include "mpuFile.h"
#include "Constant.h"
#include "Error.h"

mpuLikelihood::mpuLikelihood(int count, mpuFile * mpuPointers, bool contam)
{
  n = count;
  mpu = mpuPointers;
  hasContam = contam;
  
  sexes = new char [n];

  for (int i = 0; i < n; i++)
    sexes[i] = SEX_FEMALE;
  
  chromosomeType = CT_AUTOSOME;
}

mpuLikelihood::~mpuLikelihood()
{
  if (sexes != NULL)
    delete [] sexes;
}

void mpuLikelihood::SetAlleles(int al1, int al2)
{
  allele1 = al1;
  allele2 = al2;

  geno11 = GenotypeIndex(allele1, allele1);
  geno12 = GenotypeIndex(allele1, allele2);
  geno22 = GenotypeIndex(allele2, allele2);
}

double mpuLikelihood::Evaluate(double freq)
{
  double prior11 = freq * freq;
  double prior12 = freq * (1.0 - freq) * 2.0;
  double prior22 = (1.0 - freq) * (1.0 - freq);
  
  double prior1 = freq;
  double prior2 = 1.0 - freq;
  
  double likelihood = 0.0;
  double* lks;

   switch (chromosomeType)
      {
      case CT_MITO :
	prior11 = prior1;
	prior12 = 0.0;
	prior22 = prior2;
      case CT_AUTOSOME :
	for (int i = 0; i < n; i++) {
	  lks = mpu[i].LKs;
	  likelihood += log(prior11 * lks[geno11] +
			    prior12 * lks[geno12] +
			    prior22 * lks[geno22] +
			    1e-30);
	}
	break;
      case CT_CHRY :
	for (int i = 0; i < n; i++) {
	  lks = mpu[i].LKs;
	  if (sexes[i] == SEX_MALE)
	    likelihood += log(prior1 * lks[geno11] +
			     prior2 * lks[geno22] +
			     1e-30);
       }
       break;
      case CT_CHRX :
	for (int i = 0; i < n; i++) {
	  lks = mpu[i].LKs;
	  if (sexes[i] == SEX_MALE)
	    likelihood += log(prior1 * lks[geno11] +
			      prior2 * lks[geno22] +
			      1e-30);
	  else
	    likelihood += log(prior11 * lks[geno11] +
				prior12 * lks[geno12] +
				prior22 * lks[geno22] +
                                 1e-30);
	}
      }
   return likelihood;
}

void mpuLikelihood::GetPriors(double * priors, double freq, int i)
{
  if (sexes[i] == SEX_MALE)
    GetMalePriors(priors, freq);
  else
    GetFemalePriors(priors, freq);
}

void mpuLikelihood::GetMalePriors(double * priors, double freq)
{
  for (int i = 0; i < 10; i++)
    priors[i] = 0.0;
  
  switch (chromosomeType)
    {
    case CT_AUTOSOME :
      priors[geno11] = freq * freq;
      priors[geno12] = 2 * (1. - freq) * freq;
      priors[geno22] = (1. - freq) * (1. - freq);
      break;
    case CT_CHRY :
      priors[geno11] = freq;        /* would be zero for females */
      priors[geno12] = 0.0;
      priors[geno22] = 1. - freq;   /* would be zero for females */
      break;
    case CT_CHRX :
      priors[geno11] = freq;        /* would be freq * freq for females */
      priors[geno12] = 0.;          /* would be 2 * (1. - freq) * freq for females */
      priors[geno22] = 1.  - freq;  /* would be (1. - freq) * (1. - freq) for females */
      break;
    case CT_MITO :
      priors[geno11] = freq;
      priors[geno12] = 0;
      priors[geno22] = 1. - freq;
      break;
    }
}

void mpuLikelihood::GetFemalePriors(double * priors, double freq)
{
  for (int i = 0; i < 10; i++)
    priors[i] = 0.0;
  
  switch (chromosomeType)
    {
    case CT_AUTOSOME :
      priors[geno11] = freq * freq;
      priors[geno12] = 2 * (1. - freq) * freq;
      priors[geno22] = (1. - freq) * (1. - freq);
      break;
    case CT_CHRY :
      priors[geno11] = 0.0;            /* would be freq for males */
      priors[geno12] = 0.0;
      priors[geno22] = 0.0;            /* would be 1. - freq for males */
      break;
    case CT_CHRX :
      priors[geno11] = freq * freq;               /* would be freq for males */
      priors[geno12] = 2 * (1. - freq) * freq;    /* would be 0 for males */
      priors[geno22] = (1. - freq) * (1. - freq); /* would be 1. - freq for males */
      break;
    case CT_MITO :
      priors[geno11] = freq;
      priors[geno12] = 0;
      priors[geno22] = 1. - freq;
      break;
    }
}

double mpuLikelihood::OptimizeFrequency()
{
  for(int i=0; i < n; ++i) { 
    mpu[i].computeMixGL(0); 
  }

  a = 0.00001; fa = f(a);
  b = 0.4; fb = f(b);
  c = 0.99999; fc = f(c);

  if ( fc == fb ) { 
    min = 0.99999999;
    fmin = fb; 
    return min; 
  }

  Brent(0.00001);

  if ( hasContam ) {
    //double org = min;
    //notice("%s:%d\t%lf\t%d", mpu[0].chrom.c_str(), mpu[0].pos, min, 0);
    double prev;
    int iter = 0;
    do {
      prev = min;
      for(int i=0; i < n; ++i) { 
	mpu[i].computeMixGL(1.-min); 
      }      

      a = 0.00001; fa = f(a);
      b = min; fb = f(min);
      c = 0.99999; fc = f(c);

      Brent(0.00001);
      ++iter;
    } while ( ( iter < 20 ) && ( fabs(min - prev) > 0.0001 ) );
    //notice("%s:%d\t%lf\t%lf\t%d", mpu[0].chrom.c_str(), mpu[0].pos, org, min, iter);
  }

  return min;
}
