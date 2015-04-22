#include "mpuSetLikelihoods.h"
#include "mpuFile.h"
#include "Constant.h"
#include "Error.h"

#include <map>
#include <string>

mpuSetLikelihood::mpuSetLikelihood(mpuSet* pms, const char* pedf, const char* datf, const char* fmixcol, const char* paircol) : ms(pms), chromosomeType(CT_AUTOSOME), hweAF(0), hwdAF1(0), hwdAF2(0), brentAF(0) {
  int i;

  int n = pms->nInds;
  // create an ID map to make sure that every ID is covered
  std::map<std::string,int> idmap;
  for(i=0; i < n; ++i) {
    if ( idmap.find(pms->inds[i]) != idmap.end() ) {
      error("Duplicated ID %s in the mpuSet input",pms->inds[i].c_str());
    }
    idmap[pms->inds[i]] = i;
  }

  // check whether PED matches with the 
  pairs.resize(n,-1);
  sexes.resize(n,0);
  fmixs.resize(n,0);
  if ( pedf != NULL ) {
    fPed ped(pedf, datf);
    std::map<std::string,int> cidmap = idmap;
    for(i=0; i < n; ++i) {
      const std::string& id = ped.getPheno("IND_ID",i);
      if ( idmap.find(id) == idmap.end() ) {
	error("Individual ID %s exists in PED file, but not in the PILEUPs",id.c_str());
      }
      int vid = idmap[id];
      idmap.erase(id);
      if ( fmixcol != NULL ) {
	const std::string& s = ped.getPheno(fmixcol,i);
	if ( ( s == "." ) || ( s == "NA" ) ) fmixs[vid] = 0;
	else fmixs[vid] = atof(s.c_str());
      }
      if ( paircol != NULL ) {
	const std::string& s = ped.getPheno(paircol,i);
	if ( ( s == "." ) || ( s == "NA" ) ) pairs[vid] = -1;
	else {
	  if ( cidmap.find(s) == cidmap.end() ) {
	    warning("Paired sample %s does not exist. Ignoring..",s.c_str());
	    pairs[vid] = -1;
	  }
	  else {
	    pairs[vid] = cidmap[s];
	  }
	}
      }
    }
    // make sure that all pair INFO is consistent
    for(i=0; i < n; ++i) {
      if ( pairs[i] >= 0 ) {
	if ( pairs[pairs[i]] != i ) {
	  error("The pair information in PED file is not consistent between %s and %s",pms->inds[i].c_str(),pms->inds[pairs[i]].c_str());
	}
      }
    }
  }

  nGLs = new double[n*10];
  nLKs = new double[n*10];
  pGLs = new double[n*9];
  pLKs = new double[n*9];
  mGLs = new double[n*3];
  mLKs = new double[n*3];
  mGPs = new double[n*3];
}

mpuSetLikelihood::~mpuSetLikelihood()
{
  delete [] nGLs;
  delete [] nLKs;
  delete [] pGLs;
  delete [] pLKs;
  delete [] mGLs;
  delete [] mLKs;
}

void mpuSetLikelihood::SetAlleles(int al1, int al2)
{
  allele1 = al1;
  allele2 = al2;

  geno11 = GenotypeIndex(allele1, allele1);
  geno12 = GenotypeIndex(allele1, allele2);
  geno22 = GenotypeIndex(allele2, allele2);
}

double mpuSetLikelihood::Evaluate(double freq)
{
  double prior11 = freq * freq;
  double prior12 = freq * (1.0 - freq) * 2.0;
  double prior22 = (1.0 - freq) * (1.0 - freq);

  //notice("%lg %lg %lf %lg",freq,prior11,prior12,prior22);
  //if ( hasMix ) computeMixGLs(1.-freq);
  
  double prior1 = freq;
  double prior2 = 1.0 - freq;
  
  double likelihood = 0.0;
  double* lks;
  int n = ms->nInds;
  double log10 = log(10.);

   switch (chromosomeType)
   {
      case CT_MITO :
	prior11 = prior1;
	prior12 = 0.0;
	prior22 = prior2;
      case CT_AUTOSOME : // now only evaluate autosomes
	for (int i = 0; i < n; i++) {
	  if ( freq == 1 ) {
	    likelihood += nGLs[i*10+geno11]*log10;
	  }
	  else if ( fmixs[i] == 0 ) {
	    lks = &nLKs[i*10];
	    likelihood += log(prior11 * lks[geno11] +
			      prior12 * lks[geno12] +
			      prior22 * lks[geno22] +
			      1e-30);
	  }
	  else if ( pairs[i] < 0 ) {
	    lks = &mLKs[i*3];
	    likelihood += log(prior11 * lks[geno11] +
			      prior12 * lks[geno12] +
			      prior22 * lks[geno22] +
			      1e-30);
	  }
	  else if ( pairs[i] < i ) {
	    double* lks1 = &pLKs[i*9];
	    double* lks2 = &pLKs[pairs[i]*9];
	    likelihood += log(prior11 * prior11 * lks1[0] * lks2[0] +
			      prior11 * prior12 * lks1[1] * lks2[3] + 
			      prior11 * prior22 * lks1[2] * lks2[6] +
			      prior12 * prior11 * lks1[3] * lks2[1] +
			      prior12 * prior12 * lks1[4] * lks2[4] +
			      prior12 * prior22 * lks1[5] * lks2[7] +
			      prior22 * prior11 * lks1[6] * lks2[2] +
			      prior22 * prior12 * lks1[7] * lks2[5] +
			      prior22 * prior22 * lks1[8] * lks2[8] +
			      1e-30);
	  }
	}
	break;
      case CT_CHRY :
	for (int i = 0; i < n; i++) {
	  lks = &nLKs[i*10];
	  if (sexes[i] == SEX_MALE)
	    likelihood += log(prior1 * lks[geno11] +
			     prior2 * lks[geno22] +
			     1e-30);
       }
       break;
      case CT_CHRX :
	for (int i = 0; i < n; i++) {
	  lks = &nLKs[i*10];
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

void mpuSetLikelihood::GetPriors(double * priors, double freq, int i)
{
  if (sexes[i] == SEX_MALE)
    GetMalePriors(priors, freq);
  else
    GetFemalePriors(priors, freq);
}

void mpuSetLikelihood::GetMalePriors(double * priors, double freq)
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

void mpuSetLikelihood::GetFemalePriors(double * priors, double freq)
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

void mpuSetLikelihood::computePureGLs() {
  int i,j,k,l,m,o,p,q,b;
  int n = ms->nInds;

  //notice("3 %d",ms->nInds);

  // initialize GLs and LKs

  std::fill(nGLs, nGLs+n*10, 0);
  std::fill(nLKs, nLKs+n*10, 1);

  // calculate per-sample likelihood
  for(i=0; i < n; ++i) {
    o = i*(ms->maxDP+1);
    p = i*10;
    for(j=0; j < (int)ms->nBases[i]; ++j) {
      b = ms->bases[o+j];
      q = ms->bQs[o+j];
      for(k=0,m=0; k < 4; ++k) {
	for(l=0; l <= k; ++l, ++m) {
	  //for(k=0,m=0; k < 4; ++k) {
	  //for(l=k; l < 4; ++l, ++m) {
	  if ( k == b ) {
	    if ( l == b ) {
	      nGLs[p+m] += phredConv.phred2LogMat[q];
	    }
	    else {
	      nGLs[p+m] += phredConv.phred2HalfLogMat3[q];
	    }
	  }
	  else if ( l == b ) {
	    nGLs[p+m] += phredConv.phred2HalfLogMat3[q];
	  }
	  else {
	    nGLs[p+m] += (-0.1*q-phredConv.log3);
	  }
	}
      }
    }
    double maxGL = nGLs[p];
    for(j=p+1; j < p+10; ++j) {
      if ( maxGL < nGLs[j] ) maxGL = nGLs[j];
    }
    for(j=p; j < p+10; ++j) {
      nGLs[j] -= maxGL;
      if ( nGLs[j] < MINGL ) nGLs[j] = MINGL;
      nLKs[j] = pow(10,nGLs[j]);
    }
  }
}

// IBD GLs are independent of AF, but dependent on fmix
void mpuSetLikelihood::computeIBDGLs() {
  int i,j,o,p,q,b;
  int n = ms->nInds;

  // initialize GLs and LKs
  std::fill(pGLs, pGLs+n*9, 0);
  std::fill(pLKs, pLKs+n*9, 1);

  for(i=0; i < n; ++i) {
    o = i*(ms->maxDP+1);
    p = i*9;
    if ( fmixs[i] == 0 ) {
      for(j=0; j < (int)ms->nBases[i]; ++j) {
	b = ms->bases[o+j];
	q = ms->bQs[o+j];
	if ( b == allele1 ) {
	  pGLs[p+0] += phredConv.phred2LogMat[q];
	  pGLs[p+3] += phredConv.phred2HalfLogMat3[q];
	  pGLs[p+6] += (-0.1*q-phredConv.log3);
	}
	else if ( b == allele2 ) {
	  pGLs[p+0] += (-0.1*q-phredConv.log3);
	  pGLs[p+3] += phredConv.phred2HalfLogMat3[q];
	  pGLs[p+6] += phredConv.phred2LogMat[q];
	}
	else {
	  pGLs[p+0] += (-0.1*q-phredConv.log3);
	  pGLs[p+3] += (-0.1*q-phredConv.log3);
	  pGLs[p+6] += (-0.1*q-phredConv.log3);
	}
      }
      pGLs[p+2] = pGLs[p+1] = pGLs[p+0];
      pGLs[p+5] = pGLs[p+4] = pGLs[p+3];
      pGLs[p+8] = pGLs[p+7] = pGLs[p+6];
      
      double maxGL = pGLs[p+0];
      if ( maxGL < pGLs[p+3] ) maxGL = pGLs[p+3];
      if ( maxGL < pGLs[p+6] ) maxGL = pGLs[p+6];
      
      for(j=0; j < 9; ++j) { 
	pGLs[p+j] -= maxGL;
	if ( pGLs[p+j] < MINGL ) pGLs[p+j] = MINGL;
	pLKs[p+j] = pow(10,pGLs[p+j]);
      }
    }
    else {
      int geno, cgeno;
      double pg[3];

      for(j=0; j < (int)ms->nBases[i]; ++j) {
	b = ms->bases[o+j];
	q = ms->bQs[o+j];
	pg[0] = probBaseGivenGeno(0,allele1,allele2,q,b);
	pg[1] = probBaseGivenGeno(1,allele1,allele2,q,b);
	pg[2] = probBaseGivenGeno(2,allele1,allele2,q,b);
	for(geno=0; geno<3; ++geno) {
	  for(cgeno=0; cgeno<3; ++cgeno) {
	    pGLs[p+geno*3+cgeno] += log((1.-fmixs[i]) * pg[geno] + (fmixs[i]) * pg[cgeno]);
	  }
	}
      }
      double maxGL = pGLs[p+0];
      for(j=1; j < 9; ++j) { 
	if ( maxGL < pGLs[p+j] ) maxGL = pGLs[p+j];
      }
      for(j=0; j < 9; ++j) { 
	pGLs[p+j] -= maxGL;
	if ( pGLs[p+j] < MINGL ) pGLs[p+j] = MINGL;
	pLKs[p+j] = pow(10,pGLs[p+j]);
      }
    }
  }
}


// IBD GLs are independent of AF, but dependent on fmix
void mpuSetLikelihood::computeMixGLs(double af) {
  int i,m,p;
  int n = ms->nInds;

  // initialize GLs and LKs
  // initialize GLs and LKs
  std::fill(mGLs, mGLs+n*3, 0);
  std::fill(mLKs, mLKs+n*3, 1);

  double frqs[3];
  frqs[0] = (1-af)*(1-af);
  frqs[1] = 2*af*(1-af);
  frqs[2] = af*af;

  for(i=0; i < n; ++i) {
    if ( ( fmixs[i] == 0 ) || ( af == 0 ) ) {
      //o = i*(ms->maxDP+1);
      m = i*3;
      p = i*9;
      mGLs[m+0] = pGLs[p+0];
      mGLs[m+1] = pGLs[p+3];
      mGLs[m+2] = pGLs[p+6];
      mLKs[m+0] = pLKs[p+0];
      mLKs[m+1] = pLKs[p+3];
      mLKs[m+2] = pLKs[p+6];
    }
    else {
      //o = i*(ms->maxDP+1);
      m = i*3;
      p = i*9;
      mGLs[m+0] = log( exp(pGLs[p+0])*frqs[0] + exp(pGLs[p+1])*frqs[1] + exp(pGLs[p+2])*frqs[2] );
      mGLs[m+1] = log( exp(pGLs[p+3])*frqs[0] + exp(pGLs[p+4])*frqs[1] + exp(pGLs[p+5])*frqs[2] );
      mGLs[m+2] = log( exp(pGLs[p+6])*frqs[0] + exp(pGLs[p+7])*frqs[1] + exp(pGLs[p+8])*frqs[2] );
      
      double maxGL = mGLs[m+0];
      if ( maxGL < mGLs[m+1] ) maxGL = mGLs[m+1];
      if ( maxGL < mGLs[m+2] ) maxGL = mGLs[m+2];
      
      mGLs[m+0] -= maxGL;
      mGLs[m+1] -= maxGL;
      mGLs[m+2] -= maxGL;
      
      if ( mGLs[m+0] < MINGL ) mGLs[m+0] = MINGL;
      if ( mGLs[m+1] < MINGL ) mGLs[m+1] = MINGL;
      if ( mGLs[m+1] < MINGL ) mGLs[m+2] = MINGL;
      
      mLKs[m+0] = pow(10,mGLs[m+0]);
      mLKs[m+1] = pow(10,mGLs[m+1]);
      mLKs[m+2] = pow(10,mGLs[m+2]);
    }
  }
}

double mpuSetLikelihood::OptimizeFrequency()
{
  //notice("OptimizeFrequency called. hasmix = %d",hasMix);

  computePureGLs();   //notice("ComputePureGLs ended");
  computeIBDGLs();    //notice("ComputeIBDGLs ended");
  computeMixGLs(0);   //notice("ComputeMixGLs ended");

  /*
  notice("%lf %lf %lf",nGLs[0],nGLs[1],nGLs[2]);
  notice("%lf %lf %lf",pGLs[0],mGLs[3],mGLs[6]);
  notice("%lg %lg %lg",nLKs[0],nLKs[1],nLKs[2]);
  notice("%lg %lg %lg",pLKs[0],mLKs[3],mLKs[6]);
  error("%lf %lf %lf",mGLs[0],mGLs[1],mGLs[2]);
  */
  
  a = 0.00001; fa = f(a);
  b = 0.4; fb = f(b);
  c = 0.99999; fc = f(c);
  
  Brent(0.00001);

  //notice("Evaluation ended");

  if ( hasMix ) {
    double prev;
    int iter = 0;
    do {
      prev = min;
      computeMixGLs(1.-min);

      a = 0.00001; fa = f(a);
      b = min; fb = f(min);
      c = 0.99999; fc = f(c);

      Brent(0.00001);
      ++iter;
    } while ( ( iter < 20 ) && ( fabs(min - prev) > 0.0001 ) );
    //notice("%s:%d\t%lf\t%lf\t%d", mpu[0].chrom.c_str(), mpu[0].pos, org, min, iter);
  }

  //notice("%s:%d\tmin = %lg",ms->chrom.c_str(),ms->pos,min);
  brentAF = 1.-min;

  return min;
}

void mpuSetLikelihood::ReportSNP(wFile& baseCalls, int al1, int al2) {
  SetAlleles(al1, al2);
  char alleles[] = { 'A', 'C', 'G', 'T', 0 };

  // #CHROM\tPOS\tID
  //notice("%s\t%d\t%s", ms->chrom.c_str(), ms->pos, ms->id.c_str());
  baseCalls.printf("%s\t%d\t%s", ms->chrom.c_str(), ms->pos, ms->id.c_str());
  
  // REF
  if ( ( !ms->alt.empty() ) && ( ms->alt != "." ) ) {
    baseCalls.printf("\t%s\t%s",ms->ref.c_str(),ms->alt.c_str());
  }
  else {
    baseCalls.printf("\t%c\t%c",alleles[al1], alleles[al2]);
  }

  OptimizeFrequency();
  OptimizeFrequencyEM();

  double lRef = Evaluate(1);        // Pr(Data|Ref)
  double lVar = Evaluate(min);      // Pr(Data|Var)
  // Pr(Var|Data) = Pr(Data|Var)*Pr(Var)/[Pr(Data|Var)*Pr(Var)+Pr(Data|Ref)*Pr(Ref)]
  double qual = -10 * log10(0.999/(exp(lVar-lRef)*0.001+0.999));

  //double qual = (lVar-lRef)/log(10.)*10;
  if ( qual < 0 ) qual = 0;
  if ( qual > 255 ) qual = 255;

  baseCalls.printf("\t%.2lf",qual);
  baseCalls.printf("\t%s",qual < 3 ? "LowQual" : "PASS");
  baseCalls.printf("\t");
  
  std::string info, genotypes;
  ReportGenotypes(al1, al2, info, genotypes);
  
  baseCalls.printf("%s\t",info.c_str());
  baseCalls.printf("GT:GD:GQ:GL:LM");
  baseCalls.printf("%s\n",genotypes.c_str());
}

void mpuSetLikelihood::makeGeno(int g1, int g2, std::string& s) {
  char buf[255];
  sprintf(buf,"%d/%d",g1,g2);
  s = buf;
}

void mpuSetLikelihood::ReportGenotypes(int al1, int al2, std::string& info, std::string &genotypes) {
  int n = ms->nInds;

  info.clear();
  genotypes.clear();

  double priors[2][10];
  char buf[255];

  GetFemalePriors(priors[0], min);
  GetMalePriors(priors[1], min);

  //error("%lg %lg %lg %lg",min,priors[0][0],priors[0][1],priors[0][2]);

  int geno11 = GenotypeIndex(al1, al1);
  int geno12 = GenotypeIndex(al1, al2);
  int geno22 = GenotypeIndex(al2, al2);
  
  int label1 = 0; //al1 == refAllele ? 0 : 1;
  int label2 = 1; //al2 == refAllele ? 0 : al1 == al2 ? label1 : label1 + 1;
  
  //int genoRR = 0, genoR1 = 0, genoR2 = 0;
  
  if (label2 == 2) {
    error("Currently multi-allelic variants are not supported");    
  }

  std::string label11[2], label12[2], label22[2];

  if (chromosomeType == CT_CHRY)
    label11[0] = label12[0] = label22[0] = ".";
  else if (chromosomeType == CT_MITO) {
    int2str(label1, label11[0]);
    label12[0] = ".";
    int2str(label2, label22[0]);
  }
  else { /* CT_AUTO, CT_CHRX */
    makeGeno(label1, label1, label11[0]);
    makeGeno(label1, label2, label12[0]);
    makeGeno(label2, label2, label22[0]);
  }
  
  if (chromosomeType != CT_AUTOSOME) {
    int2str(label1, label11[1]);
    label12[1] = ".";
    int2str(label2, label22[1]);
  }
  else {
    makeGeno(label1, label1, label11[1]);
    makeGeno(label1, label2, label12[1]);
    makeGeno(label2, label2, label22[1]);
  }
  
  int ns         = 0;
  int dp         = 0;
  int ac[4]      = {0, 0, 0, 0};
  //double Acount  = 0.5;
  //double ABcount = 1.0;
  
  // iterate over n individuals
  
  for (int i = 0; i < n; i++) {
    int sex = sexes[i] == SEX_MALE ? 1 : 0;

    std::pair<int,int> r = GetBestRatio(priors[sex], i, pairs[i]);
    int best = r.first;
    int quality = r.second;
      
    std::string & label = best == geno11 ? label11[sex] : best == geno12 ? label12[sex] : label22[sex];
    bool    nocall = label[0] == '.';
    int     depth  = ms->nBases[i];
    dp += depth;
      
    sprintf(buf,"\t%s:%d:%d",label.c_str(), depth, nocall ? 0 : quality);
    genotypes += buf;
      
    if (label[0] != '.') {
      if ( depth > 0 ) {
	ns++;
	ac[label[0] - '0']++;
      }
      
      if (label.size() > 1 && label[2] != '.') {
	if ( depth > 0 ) ac[label[2] - '0']++;
      }
      
      if (label1 == 0 && label2 == 0)
	continue;

      if (label2 >= 2) 
	error("label2=%s, Currently multi-allelic variant is not supported",label2);

      if ( fmixs[i] == 0 ) {
	double * GLs = &nGLs[i*10];
	sprintf(buf,":%.2lf,%.2lf,%.2lf:.",GLs[geno11],GLs[geno12],GLs[geno22]);
      }
      else if ( pairs[i] < 0 ) {
	double * GLs = &mGLs[i*3];
	sprintf(buf,":%.2lf,%.2lf,%.2lf:.",GLs[0],GLs[1],GLs[2]);
      }
      else {
	double* GLs = &mGLs[i*3];
	double* LMs = &pGLs[i*9];
	sprintf(buf,":%.2lf,%.2lf,%.2lf:%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf",GLs[0],GLs[1],GLs[2],LMs[0],LMs[1],LMs[2],LMs[3],LMs[4],LMs[5],LMs[6],LMs[7],LMs[8]);
      }
      genotypes += buf;
    }
  }

  //double AB = Acount / (ABcount + 1e-30);
  
  sprintf(buf,"DP=%d;NS=%d", dp, ns); info += buf;
  sprintf(buf,";AN=%d", ac[0] + ac[1] + ac[2] + ac[3]); info += buf;
  if (label1 == 0 && label2 == 0)
    return;
  
  //if (label2 < 2) {
  sprintf(buf,";AC=%d;AF=%.6lf", ac[1], 1. - min); info += buf;
  getInfoGL(info, al1,al2);
  //}
  //else {
  //sprintf(buf,";AC=%d,%d;AF=%.6lf;%.6lf", ac[1], ac[2],lk.min, 1. - min); info += buf;
  //}
  //sprintf(buf,";AB=%.4lf", AB); info += buf;
}

std::pair<int,int> mpuSetLikelihood::GetBestRatio(const double priors[], int self, int pair) {
  //error("%d %d",self,pair);
  if ( pair < 0 ) {
    double* lks = &mLKs[self*3];
    double p0 = lks[0] * priors[geno11];
    double p1 = lks[1] * priors[geno12];
    double p2 = lks[2] * priors[geno22];
    int best = (p0 < p1) ? (p1 < p2 ? geno22 : geno12) : (p0 < p2 ? geno22 : geno11);
    double sum = p0+p1+p2;
    double e = ((best == geno11) ? (p1+p2) : (best == geno12 ? (p0+p2) : (p0+p1)))/sum;
    //notice("%d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf",self,mGLs[self*3],mGLs[self*3+1],mGLs[self*3+2],nGLs[self*10],nGLs[self*10+1],nGLs[self*10+2]);
    //if ( best == 1 ) abort();
    if (e < 3.16e-26 ) 
      return std::pair<int,int>(best,255);
    else 
      return std::pair<int,int>(best,int (-log10(e) * 10 + 0.5));
  }
  else {  // paired sample exists
    double* lks1 = &pLKs[self*9];
    double* lks2 = &pLKs[pair*9];
    double lks[9], mlks[3];

    double sum = 0.0;
    double bestlk = 0;
    int best = 0;
    
    for( int i=0; i < 3; ++i) {
      for(int j=0; j < 3; ++j) {
	lks[3*i+j] = lks1[i*3+j] * lks2[j*3+i] * priors[i] * priors[j];
	sum += lks[i*3+j];
      }
      mlks[i] = lks[3*i+0]+lks[3*i+1]+lks[3*i+2];
      if ( bestlk < mlks[i] ) {
	bestlk = mlks[i];
	best = i;
      }
    }
    
    double e = (sum > 0) ? (1.0 - bestlk/sum) : 1.0;
    //int qual = (error < 0.0000000001) ? 100 : int(-log10(error) * 10 + 0.5);
    int qual = (e < 3.16e-26) ? 255 : int(-log10(e) * 10 + 0.5);
    return std::pair<int,int>(best,qual);
  }
}

void mpuSetLikelihood::OptimizeFrequencyEM() {
  // fit AF with/without HWE
  double p = .1 + rand()/(RAND_MAX+1.)*.8; // start from random AF [.1-.9]
  double f0,f1,f2, fsum, sum, sum0, sum1, sum2;
  bool convE = false, convD = false;
  
  double q = 1.-p;
  double p0 = q * q;
  double p1 = 2. * p * q;
  double p2 = p * p;
  int i;

  int geno11 = 0, geno12 = 1, geno22 = 2;
  double eps = 1e-6;
  int n = ms->nInds;
  int rounds;
  int maxiter = 20;
  double depth, scale, minimum, delta, nref = 0;
  double Acount = 0.05, ABcount = 0.1;
  double llkHWE = 0, llkHWD = 0, obsHET = 0;

  if ( brentAF == 0 ) {
    computePureGLs();   //notice("ComputePureGLs ended");
    computeIBDGLs();    //notice("ComputeIBDGLs ended");
    computeMixGLs(0);   //notice("ComputeMixGLs ended");
  }
  else {
    computeMixGLs(0);
  }
  
  for(rounds = 0; rounds < maxiter; ++rounds) {
    if ( ( hasMix ) && ( rounds > 0 ) ) {
      computeMixGLs(hweAF);
    }

    sum = sum0 = sum1 = sum2 = 0;
    for(i=0; i < n; ++i) {
      const double * lks = &mLKs[i*3];
      if ( !convE ) {
	f0 = q * q * lks[geno11];
	f1 = 2. * p * q * lks[geno12];
	f2 = p * p * lks[geno22];
	fsum = f0+f1+f2; 
	sum += f1/fsum; 
	sum += (2. * f2/fsum);
      }
      
      if ( !convD ) {
	f0 = p0 * lks[geno11];
	f1 = p1 * lks[geno12];
	f2 = p2 * lks[geno22];
	fsum = f0+f1+f2;
	sum0 += f0/fsum; 
	sum1 += f1/fsum; 
	sum2 += f2/fsum; 
      }
    }
    
    if ( !convE ) {
      p = sum / (2*n);
      if ( fabs(p + q - 1.) < eps ) convE = true;
      q = 1.-p;
    }
    if ( !convD ) {
      if ( ( fabs(p1 - sum1/n) < eps ) && ( fabs(p2 - sum2/n) < eps ) ) convD = true;
      p0 = sum0 / (n);
      p1 = sum1 / (n);
      p2 = sum2 / (n);
    }
    
    if ( convE && convD ) break;
  }

  for(i=0; i < n; ++i) {
    const double* lks = &mLKs[i*3];
    const double* gls = &mGLs[i*3];
    depth = ms->nBases[i];
    
    f0 = q * q * lks[geno11];
    f1 = 2. * p * q * lks[geno12];
    f2 = p * p * lks[geno22];
    fsum = f0+f1+f2+1e-30; // sum_g Pr(g)Pr(Data|g)
    llkHWE += log(fsum);

    mGPs[i*3+0] = f0/fsum;
    mGPs[i*3+1] = f1/fsum;
    mGPs[i*3+2] = f2/fsum;

    if ( ( mGPs[i*3+1] > .1 ) && ( depth > 0 ) ) {
      scale = -10.0*(gls[geno22] + gls[geno11] - 2 * gls[geno12]) + 6 * depth;
      minimum = abs(10*(gls[geno22] - gls[geno11]));
      
      if (scale < 4) scale = 4;
      if (scale < minimum) scale = minimum;
      
      delta = (-10*(gls[geno22] - gls[geno11])) / (scale + 1e-30);
      nref = 0.5 * depth * (1.0 + delta);
      
      Acount += (mGPs[i*3+1] * nref);
      ABcount += (mGPs[i*3+1] * depth);
    }

    f0 = p0 * lks[geno11];
    f1 = p1 * lks[geno12];
    f2 = p2 * lks[geno22];
    fsum = f0+f1+f2+1e-30; // sum_g Pr(g)Pr(Data|g)
    llkHWD += log(fsum);
    obsHET += f1/fsum;
  }
  hweAF = p;
  hwdAF1 = p1;
  hwdAF2 = p2;

  FIC = 1.-obsHET/(2.*p*q*n);
  SLRT = (FIC > 0 ? 2 : -2 )*(llkHWD-llkHWE);
  ABL = Acount/ABcount;

  //notice("p=%lf\tq=%lf\tAcount=%lf\tABcount=%lf\tABL=%lf",p,q,Acount,ABcount,ABL);
}

/*
void mpuSetLikelihood::getInfoGL(std::string &info) {
  // calculate statistics
  int sa[4] = {1,1,1,1};
  int ta[4] = {1,1,1,1};
  int q1 = 60, q2 = 1800, c1 = 100, c2 = 5000, a1 = 1, qa = 30, ca = 50, nr = 0, na = 0, ne = 0, m1 = 0, m2 = 0, ma1 = 0, ma2 = 0, lm = 0;

  int vsa[4] = {1,1,1,1};
  int vta[4] = {1,1,1,1};

  // among variant-carrying individuals
  // count statistics.
  
  int b, bq, s, c, m, o;
  for(i=0; i < n; ++i) {
    o = i*(ms->maxDP+1);
    for(j=0; j < (int)ms->nBases[i]; ++j) {
      b = ms->bases[o+j];
      bq = ms->bQs[o+j];
      m = ms->mQs[o+j];
      c = ms->cycles[o+j];
      s = ms->strands[o+j];

      ++m1;
      m2 += (m*m);
      if ( m <= 13 ) { ++lm; };
      if ( b <= 1 ) {
	++sa[s*2+b];
	++ta[(c<10 || c>70)*2+b];
	q1 += bq;
	q2 += (bq*bq);
	c1 += c;
	c2 += (c*c);
	a1 += b;
	qa += bq*b;
	ca += c*b;
	if ( b == 1 ) {
	  ++ma1;
	  ma2 += (m*m);
	}
	if ( bq >= 13 ) {
	  if ( b == 0 ) ++nr;
	  else ++na;
	}
      }
      else {
	if ( bq >= 13 ) {
	  ++ne;
	}
      }
    }
  }

  int sasum = sa[0]+sa[1]+sa[2]+sa[3];
  double STR = (double)(sa[0]*sa[3]-sa[1]*sa[2])/sqrt((double)(sa[0]+sa[1])*(sa[2]+sa[3])*(sa[0]+sa[2])*(sa[1]+sa[3]));
  double STZ = STR * sqrt((double)sasum);
  double BQR = (double)(qa - q1*a1/(double)sasum)/sqrt((double)(q2-q1*q1/(double)sasum)*(a1-a1*a1/(double)sasum));
  double BQZ = BQR * sqrt((double)sasum);
  double CBR = (double)(ca - c1*a1/(double)sasum)/sqrt((double)(c2-c1*c1/(double)sasum)*(a1-a1*a1/(double)sasum));
  double CBZ = CBR * sqrt((double)sasum);
  double TBR = (double)(ta[0]*ta[3]-ta[1]*ta[2])/sqrt((double)(ta[0]+ta[1])*(ta[2]+ta[3])*(ta[0]+ta[2])*(ta[1]+ta[3]));
  double TBZ = TBR * sqrt((double)sasum);
  double FOB = (double)(ne+.1)/(double)(na+.1);
  int MQ = (int)sqrt(m2/(m1+1e-6));
  int AMQ = (int)sqrt((ma2)/(ma1+1e-6));
  double LMQ = lm/(m1+1e-6);

  char buf[65536];
  sprintf(buf,";HWEAF=%.6lf;HWDAF=%.6lf,%.6lf;FIC=%.3lf;SLRT=%.3lf;ABL=%.4lf;STR=%.3lf;STZ=%.3lf;BQR=%.3lf;BQZ=%.3lf;CBR=%.3lf;CBZ=%.3lf;TBR=%.3lf;TBZ=%.3lf;FOB=%.3lf;MMQ=%d;AMQ=%d;LMQ=%.3lf",pHWE,pHWD1,pHWD2,inbreedingCoeff,signedLRT,Acount/ABcount,STR,STZ,BQR,BQZ,CBR,CBZ,TBR,TBZ,FOB,MQ,AMQ,LMQ);
  info += buf;
}
*/

void mpuSetLikelihood::getInfoGL(std::string &info, int al1, int al2) {
  int n = ms->nInds;
  // calculate statistics
  int sa[4] = {1,1,1,1};
  int ta[4] = {1,1,1,1};
  int q1 = 60, q2 = 1800, c1 = 100, c2 = 5000, a1 = 1, qa = 30, ca = 50, nr = 0, na = 0, ne = 0, m1 = 0, m2 = 0, ma1 = 0, ma2 = 0, lm = 0;

  // among variant-carrying individuals
  // count statisticss

  // What statistics do we need?
  // Per individual
  // # of REF/ALT/OTH bases (for allele balance)
  double sgpb[9] = {0,0,0,0,0,0,0,0,0};
  int ns[3] = {0,0,0};
  int nMaxAlt = 0, nMaxRef = 0, nRefAlt = 0;
  
  int i, j, b, bq, s, c, m, o; // base, BQ, strand, cycle, MQ, offset
  for(i=0; i < n; ++i) {
    double gpb[9] = {0,0,0,0,0,0,0,0,0};
    o = i*(ms->maxDP+1);
    if ( mGPs[i*3+0] > 0.5 ) ++ns[0];
    else if ( mGPs[i*3+1] > 0.5 ) ++ns[1];
    else if ( mGPs[i*3+2] > 0.5 ) ++ns[2];
    for(j=0; j < (int)ms->nBases[i]; ++j) {
      b = ms->bases[o+j];
      bq = ms->bQs[o+j];
      m = ms->mQs[o+j];
      c = ms->cycles[o+j];
      s = ms->strands[o+j];
      ++m1;        // # of MQs
      m2 += (m*m); // sumsqMQ

      if ( b == al1 ) b = 0;
      else if ( b == al2 ) b = 1;
      else b = 3;

      if ( m <= 13 ) { ++lm; };
      if ( b <= 1 ) {
	++sa[s*2+b];
	++ta[(c<10 || c>70)*2+b];
	q1 += bq;
	q2 += (bq*bq);
	c1 += c;
	c2 += (c*c);
	a1 += b;
	qa += bq*b;
	ca += c*b;
	if ( b == 1 ) {
	  ++ma1;
	  ma2 += (m*m);
	}
	if ( mGPs[i*3+0] > 0.5 ) ++gpb[3*0+b];
	else if ( mGPs[i*3+1] > 0.5 ) ++gpb[3*1+b];
	else if ( mGPs[i*3+2] > 0.5 ) ++gpb[3*2+b];
	//gpb[3*0+b] += mGPs[i*3+0];
	//gpb[3*1+b] += mGPs[i*3+1];
	//gpb[3*2+b] += mGPs[i*3+2];
	if ( bq >= 13 ) {
	  if ( b == 0 ) ++nr;
	  else ++na;
	}
      }
      else {
	if ( mGPs[i*3+0] > 0.5 ) ++gpb[3*0+2];
	else if ( mGPs[i*3+1] > 0.5 ) ++gpb[3*1+2];
	else if ( mGPs[i*3+2] > 0.5 ) ++gpb[3*2+2];
	//gpb[3*0+2] += mGPs[i*3+0];
	//gpb[3*1+2] += mGPs[i*3+1];
	//gpb[3*2+2] += mGPs[i*3+2];
	if ( bq >= 13 ) {
	  ++ne;
	}
      }
    }

    if ( nMaxRef < gpb[0]+gpb[3]+gpb[6] ) nMaxRef = gpb[0]+gpb[3]+gpb[6];
    else if ( nMaxAlt < gpb[1]+gpb[4]+gpb[7] ) {
      nMaxAlt = gpb[1]+gpb[4]+gpb[7];
      nRefAlt = gpb[0]+gpb[3]+gpb[6];
    }
    for(j=0; j < 9; ++j) { sgpb[j] += gpb[j]; }
  }

  int sasum = sa[0]+sa[1]+sa[2]+sa[3];
  double STR = (double)(sa[0]*sa[3]-sa[1]*sa[2])/sqrt((double)(sa[0]+sa[1])*(sa[2]+sa[3])*(sa[0]+sa[2])*(sa[1]+sa[3]));
  double STZ = STR * sqrt((double)sasum);
  double BQR = (double)(qa - q1*a1/(double)sasum)/sqrt((double)(q2-q1*q1/(double)sasum)*(a1-a1*a1/(double)sasum));
  double BQZ = BQR * sqrt((double)sasum);
  double CBR = (double)(ca - c1*a1/(double)sasum)/sqrt((double)(c2-c1*c1/(double)sasum)*(a1-a1*a1/(double)sasum));
  double CBZ = CBR * sqrt((double)sasum);
  double TBR = (double)(ta[0]*ta[3]-ta[1]*ta[2])/sqrt((double)(ta[0]+ta[1])*(ta[2]+ta[3])*(ta[0]+ta[2])*(ta[1]+ta[3]));
  double TBZ = TBR * sqrt((double)sasum);
  double FOB = (double)(ne+.1)/(double)(na+.1);
  int MQ = (int)sqrt(m2/(m1+1e-6));
  int AMQ = (int)sqrt((ma2)/(ma1+1e-6));
  double LMQ = lm/(m1+1e-6);
  double ABE = (sgpb[3*1+0]+0.5)/(sgpb[3*1+0]+sgpb[3*1+1]+1.);
  //double NRO = (gpb[5]+gpb[8]+1)*(gpb[0]+gpb[1]+gpb[2]+1)/(gpb[3]+gpb[4]+gpb[5]+gpb[6]+gpb[7]+gpb[8]+1)/(gpb[2]+1);
  double NRO = (sgpb[5]/(ns[1]+.01) + sgpb[8]/(ns[2]+.01) + 1)/(sgpb[2]/(ns[0]+.01) + 1);

  char buf[65536];
  sprintf(buf,";HWEAF=%.5lf;HWDAF=%.5lf,%.5lf;FIC=%.3lf;SLRT=%.3lf;ABL=%.4lf;STR=%.3lf;STZ=%.3lf;BQR=%.3lf;BQZ=%.3lf;CBR=%.3lf;CBZ=%.3lf;TBR=%.3lf;TBZ=%.3lf;FOB=%.3lf;MMQ=%d;AMQ=%d;LMQ=%.3lf;ABE=%.3lf;NRO=%.3lf;MAXALT=%d;MAXREF=%d;ALTFRAC=%.3lf", //;GPB=%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf"
	  hweAF, hwdAF1, hwdAF2, FIC,SLRT,ABL,STR,STZ,BQR,BQZ,CBR,CBZ,TBR,TBZ,FOB,MQ,AMQ,LMQ,ABE,NRO,nMaxAlt,nMaxRef,(double)nMaxAlt/(double)(nMaxAlt+nRefAlt+1e-10));//,gpb[0]/ns[0],gpb[1]/ns[0],gpb[2]/ns[0],gpb[3]/ns[1],gpb[4]/ns[1],gpb[5]/ns[1],gpb[6]/ns[2],gpb[7]/ns[2],gpb[8]/ns[2]);
  info += buf;
}
