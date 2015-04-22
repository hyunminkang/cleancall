#ifndef __MPU_SET_LIKELIHOODS_H__
#define __MPU_SET_LIKELIHOODS_H__

#include "mpuSet.h"
#include "MathGold.h"
#include "fPed.h"

#include <math.h>

// Chromosome types
#define CT_AUTOSOME     0
#define CT_CHRX         1
#define CT_CHRY         2
#define CT_MITO         3

class mpuSetLikelihood : public ScalarMinimizer {
 public:
  mpuSet     * ms;
  std::vector<int>      pairs;
  std::vector<uint8_t>  sexes;
  std::vector<double>   fmixs;
  int     chromosomeType;
  bool    hasPair;
  bool    hasMix;

  // likelihoods under no contamination - independent of AF, fMix
  double * nGLs;
  double * nLKs;
  // likelihoods for paired sample contamination - dep. on fMix, indep. of AF
  double * pGLs;
  double * pLKs;
  // likelihoods for mixture - dep. on fMix, AF
  double * mGLs;
  double * mLKs; 

  double * mGPs;
  double hweAF, hwdAF1, hwdAF2, brentAF, FIC, SLRT, ABL;

  mpuSetLikelihood(mpuSet* pms, const char* pedf = NULL, const char* datf = NULL, const char* fmixcol = NULL, const char* paircol = NULL);
  virtual ~mpuSetLikelihood();
  
  void SetAlleles(int al1, int al2);
  
  virtual double f(double freq) { return -Evaluate(freq); }
  
  virtual double Evaluate(double freq);
  
  void GetPriors(double * priors, double freq, int i);
  void GetMalePriors(double * priors, double freq);
  void GetFemalePriors(double * priors, double freq);
  
  double OptimizeFrequency();
  void OptimizeFrequencyEM();
  
  static int GenotypeIndex(int base1, int base2) {
    // 00 10 11 20 21 22 30 31 32 33
    // i*(i+1)/2 + j
    return base1 < base2 ? base2*(base2+1)/2+base1 : base1*(base1+1)/2+base2;
      //return base1 < base2 ? (base1) *(9 - base1) / 2 + (base2 - base1) :
      //(base2) *(9 - base2) / 2 + (base1 - base2);
  }
  
  inline double probBaseGivenGeno(int geno, uint8_t nref, uint8_t nalt, uint8_t qual, uint8_t base) {
    if ((base != nref) && (base != nalt)) 
      return phredConv.phred2Err[qual]/3.0;
    else if (geno==1) 
      return 1/2.0 - phredConv.phred2Err[qual]/3.0;
    else if ((geno==0 && base==nref) || (geno==2 && base==nalt)) 
      return phredConv.phred2Mat[qual];
    else
      return phredConv.phred2Err[qual]/3.0;
  }

  static void int2str(int x, std::string& s) {
    char buf[255];
    sprintf(buf,"%d",x);
    s = buf;
  }

  static void makeGeno(int g1, int g2, std::string& s);

  void ReportSNP(wFile& baseCalls, int al1 = 0, int al2 = 1);
  void ReportGenotypes(int al1, int al2, std::string& info, std::string &genotypes);
  std::pair<int,int> GetBestRatio(const double priors[], int self, int pair = -1);

  void computePureGLs();
  void computeIBDGLs();
  void computeMixGLs(double af = 0);

  void getInfoGL(std::string &info, int al1 = 0, int al2 = 1);
  
 protected:
  int allele1, allele2;
  int geno11, geno12, geno22;
};
#endif
