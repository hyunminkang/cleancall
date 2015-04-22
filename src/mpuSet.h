#ifndef __MPU_SET_H__
#define __MPU_SET_H__

// MPUSET is supposed to contain
// #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT IND1 IND2 IND3

#include "pFile.h"
#include "fVcf.h"
#include "PhredHelper.h"
#include "BaseAsciiMap.h"

#include <vector>
#include <cstdlib>
#include <cmath>

#define ENDPOS 1000000000
#define MINGL -30
#define MINLK 1e-30
#define DEFAULT_MAX_DP 255

// when mpuSet is being used, 
// we assume that chromosomal order is the same between samples
class mpuSet {
 public:
  // information shared across individuals
  pFile tf;          // file to read
  int nInds;         // number of individuals
  std::string chrom; // chromosome name
  int pos;           // 1-based position
  std::string id;
  std::string ref; 
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;
  int maxDP;         // max depth per individual
  char* line;        // line
  bool nobgzf;
  bool hasHeader;
  bool hasMQ;
  bool hasCY;

  // individual level info
  std::vector<uint32_t> nBases;
  //std::vector<uint32_t> cumBases;
  std::vector<std::string> inds;

  // base level info
  std::vector<uint8_t>  bases;
  std::vector<uint8_t>  bQs;
  std::vector<uint8_t>  mQs;
  std::vector<uint8_t>  cycles;
  std::vector<bool>     strands;

  mpuSet() : maxDP(DEFAULT_MAX_DP), nobgzf(false), hasHeader(false) {}

  mpuSet(int _ninds, int _maxDP = DEFAULT_MAX_DP, bool _nobgzf = false) : nobgzf(_nobgzf), hasHeader(false) {
    setParams(_ninds, _maxDP);
  }

  mpuSet(const char* filename, const char* region = NULL, bool printHeader = false) : maxDP(DEFAULT_MAX_DP), nobgzf(false), hasHeader(false) {
    load(filename, region, printHeader);
  }

  static int nchr(const std::string& s) {
    int n = atoi(s.c_str());
    if ( n == 0 ) {
      if ( s == "X" ) return 23;
      else if ( s == "Y" ) return 24;
      else if ( s == "XY" ) return 25;
      else if ( s == "MT" ) return 26;
      else error("Cannot recognize chromosome %s",s.c_str());
    }
    else {
      return n;
    }
    return 0;
  }

  int compare(const mpuSet& m) {
    if ( isEOF() ) {
      if ( m.isEOF() ) return 0;
      else return 1;
    }
    else if ( m.isEOF() ) return -1;
    else if ( chrom == m.chrom ) return (pos - m.pos);
    else return ( nchr(chrom) - nchr(m.chrom) );
  }

  void setParams(int _ninds, int _maxDP = DEFAULT_MAX_DP) {
    nInds = _ninds;
    maxDP = _maxDP;
    nBases.resize(nInds);

    char buf[255];
    inds.resize(nInds);
    for(int i=0; i < nInds; ++i) {
      if ( inds[i].empty() ) {
	sprintf(buf,"ID%d",i+1);
	inds[i] = buf;
      }
    }

    bases.resize((maxDP+1) * nInds,0);
    bQs.resize((maxDP+1) * nInds,0);
    mQs.resize((maxDP+1) * nInds,255);
    cycles.resize((maxDP+1) * nInds,255);
    strands.resize((maxDP+1) * nInds,false);
  }

  int parseInds(char* line, int startIdx = 9) {
    char* pch = line;
    char* nch = NULL;
    int i, j;

    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	inds.push_back(id);
	++j;
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    return j;
  }

  bool load(const char* filename, const char* region = NULL, bool printHeader = false) {
    tf.load(filename, region, printHeader, nobgzf);

    while ( (line = (char*)tf.getLine()) != NULL ) {
      if ( line[0] == '#' ) {
	hasHeader = true;
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line);
	}
	else if ( line[1] == '#' ) {
	  // meta lines - ignore
	}
      }
      else {
	break;
      }
    }

    if ( !hasHeader ) nInds = 1; 

    setParams(nInds, maxDP);

    return parseMarker();
  }

  bool parseMarker() {
    if ( isEOF() ) return false;

    hasMQ = hasCY = false;

    if ( hasHeader ) { // if header exists, parse as VCF-like format
      char *pch, *nch, *p;
      int i,j,k,o,nbase,onbase;
      int nref = 0;
      std::string s;
      int startIdx = 9;
      bool hasAlt = false;
    
      if ( line != NULL ) {
	pch = line;
	for(i=0, j=0; pch != NULL; ++i) {
	  nch = strchr(pch,'\t');
	  if ( i < startIdx ) {
	    if ( nch == NULL ) s.assign(pch);
	    else s.assign(pch, nch - pch);
	    switch(i) {
	    case 0:
	      chrom = s; break;
	    case 1:
	      pos = atoi(pch); 
	      break;
	    case 2:
	      id = s; break;
	    case 3:
	      ref = s; 
	      nref = BaseAsciiMap::base2int[(int)*pch];
	      break;
	    case 4:
	      alt = s; 
	      //notice("%s:%d\t%s\t%s",chrom.c_str(),pos,ref.c_str(),alt.c_str());
	      if ( ( !alt.empty() ) && ( alt != "." ) ) hasAlt = true;
	      break;
	    case 5:
	      qual = s; break;
	    case 6:
	      filter = s; break;
	    case 7:
	      info = s; break;
	    }
	    pch = ( nch == NULL ) ? NULL : nch+1;
	  }
	  else {
	    // parse nBase first
	    onbase = atoi(pch);
	    if ( onbase > maxDP ) nbase = maxDP;
	    else nbase = onbase;
	    nBases[i-startIdx] = nbase;
	    //cumBases[i-startIdx] = (i == startIdx) ? nbase : ( cumBases[i-startIdx] + nbase );

	    if ( nbase == 0 ) {
	      p = strchr(pch, '\t');
	      pch = (p == NULL) ? NULL : p+1;
	      continue;
	    }

	    // parse bases
	    p = strchr(pch, ':') + 1;
	    j = 0;
	    o = (i-startIdx)*(maxDP+1);
	    if ( hasAlt ) {
	      for(k=0; j < nbase; ++k) {
		switch(p[k]) {
		case 'A':
		  bases[o+j] = 1;
		  strands[o+j] = true;
		  ++j;
		  break;
		case 'a':
		  bases[o+j] = 1;
		  strands[o+j] = false;
		  ++j;
		  break;
		case 'R':
		  bases[o+j] = 0;
		  strands[o+j] = true;
		  ++j;
		  break;
		case 'r':
		  bases[o+j] = 0;
		  strands[o+j] = false;
		  ++j;
		  break;
		case '*': case 'n': case 'N':
		  bases[o+j] = 4;
		  strands[o+j] = false;
		  ++j;
		  break;
		default:
		  error("Cannot recognize %c in the pileups at %s:%d i=%d, j= %d, k=%d, nbase=%d\n%s\n%s\n",p[k],chrom.c_str(),pos,i,j,k,nbase,pch,p+k);
		}
	      }
	    }
	    else {
	      for(k=0; j < nbase; ++k) {
		switch(p[k]) {
		case '.': 
		  bases[o+j] = nref;
		  strands[o+j] = true;
		  ++j;
		  break;
		case ',':
		  bases[o+j] = nref;
		  strands[o+j] = false;
		  ++j;
		  break;
		case 'A':
		  bases[o+j] = 0;
		  strands[o+j] = true;
		  ++j;
		  break;
		case 'a':
		  bases[o+j] = 0;
		  strands[o+j] = false;
		  ++j;
		  break;
		case 'C':
		  bases[o+j] = 1;
		  strands[o+j] = true;
		  ++j;
		  break;
		case 'c':
		  bases[o+j] = 1;
		  strands[o+j] = false;
		  ++j;
		  break;
		case 'G':
		  bases[o+j] = 2;
		  strands[o+j] = true;
		  ++j;
		  break;
		case 'g':
		  bases[o+j] = 2;
		  strands[o+j] = false;
		  ++j;
		  break;
		case 'T':
		  bases[o+j] = 3;
		  strands[o+j] = true;
		  ++j;
		  break;
		case 't':
		  bases[o+j] = 3;
		  strands[o+j] = false;
		  ++j;
		  break;
		case '*': case 'n': case 'N':
		  bases[o+j] = 4;
		  strands[o+j] = false;
		  ++j;
		  break;
		default:
		  error("Cannot recognize %c in the pileups at i=%d, j= %d, k=%d, nbase=%d\n%s\n%s\n",p[k],i,j,k,pch,p+k);
		}
	      }
	    }

	    //p += k;
	    p += onbase;

	    while( ( *p != ':' ) && ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	    if ( ( *p == '\t' ) || ( *p == '\0' ) ) {
	      error("Cannot find BQ field");
	    }

	    ++p;
	    for(j=0; j < onbase; ++j) {
	      if ( j > 0 ) {
		if ( *p != ',' ) return parseError("Comma",p);
		++p;
	      }
	      if ( ( *p < '0' ) || ( *p > '9' ) ) return parseError("Number",p);
	      if ( j < nbase ) bQs[j+o] = strtol(p, &p, 10);
	      else strtol(p, &p, 10);
	    }
	    if ( *p == ':' ) hasMQ = true;
	    else if ( ( *p != '\t' ) && ( *p != '\0' ) ) {
	      pch = ( *p == '\0' ) ? NULL : p+1;
	      continue;
	    }
	    else return parseError("Colon or whitespace",p);

	    ++p;
	    for(j=0; j < onbase; ++j) {
	      if ( j > 0 ) {
		if ( *p != ',' ) return parseError("Comma",p);
		++p;
	      }
	      if ( ( *p < '0' ) || ( *p > '9' ) ) return parseError("Number",p);
	      if ( j < nbase ) mQs[j+o] = strtol(p, &p, 10);
	      else strtol(p, &p, 10);
	    }
	    if ( *p == ':' ) hasCY = true;
	    else if ( ( *p != '\t' ) && ( *p != '\0' ) ) {
	      pch = ( *p == '\0' ) ? NULL : p+1;
	      continue;
	    }
	    else return parseError("Colon or whitespace",p);

	    ++p;
	    for(j=0; j < onbase; ++j) {
	      if ( j > 0 ) {
		if ( *p != ',' ) return parseError("Comma",p);
		++p;
	      }
	      if ( ( *p < '0' ) || ( *p > '9' ) ) return parseError("Number",p);
	      if ( j < nbase ) cycles[j+o] = strtol(p, &p, 10);
	      else strtol(p, &p, 10);
	    }
	    if ( ( *p != '\t' ) && ( *p != '\0' ) ) return parseError("Whitespace",p);
	    else { pch = ( *p == '\0' ) ? NULL : p+1; }

	    // parse bQs;	  
	    //for(j=0; j < nbase; ++j) {
	    //bQ[j+o] = strtol(p, &p, 10);
	    //if ( ( *p != ':' ) && ( *p != '\t' ) && ( *p != '\0' ) ) break;
	    //++p;
	      /*
	      if ( j > 0 ) {
		while( (*p != ',') && (*p != '\0') ) ++p;
	      }
	      ++p;
	      bQs[j+o] = atoi(p);*/
	      //}
	    //while( ( *p != ':' ) && ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	    //if ( ( *p == '\t' ) || ( *p == '\0' ) ) {
	    //  pch = ( *p == '\0' ) ? NULL : p+1;
	    //  continue;
	    //}

	    /*
	    hasMQ = true;
	    // parse mQs;
	    for(j=0; j < nbase; ++j) {
	      if ( j > 0 ) {
		//while(*p != ',') ++p;
		while( (*p != ',') && (*p != '\0') ) ++p;
	      }
	      ++p;
	      mQs[j+o] = atoi(p);
	    }
	    while( ( *p != ':' ) && ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	    if ( ( *p == '\t' ) || ( *p == '\0' ) ) {
	      pch = ( *p == '\0' ) ? NULL : p+1;
	      continue;
	    }

	    hasCY = true;
	    // parse cycles;
	    for(j=0; j < nbase; ++j) {
	      if ( j > 0 ) {
		//while(*p != ',') ++p;
		while( (*p != ',') && (*p != '\0') ) ++p;
	      }
	      ++p;
	      cycles[j+o] = atoi(p);
	    }
	    while( ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	    pch = ( *p == '\0' ) ? NULL : p+1;
	    */
	  }
	}
	if ( i != startIdx + nInds ) {
	  char buf[255];
	  sprintf(buf,"%d columns",startIdx+nInds);
	  return parseError(buf,NULL);
	}
      }
    }
    else {  // this is a single-sample pileup
      char *pch, *nch, *p;
      int j,k,l,o,nbase,nref;

      if ( line != NULL ) {
	pch = line;

	nch = strchr(pch,'\t');
	chrom.assign( pch, nch-pch ); // CHROM
	pch = nch+1;

	nch = strchr(pch,'\t');
	pos = atoi(pch); // POS
	pch = nch+1;

	nch = strchr(pch,'\t');	
	ref.assign( pch, nch-pch ); // CHROM
	nref = BaseAsciiMap::base2int[(int)*pch];
	pch = nch+1;

	nbase = atoi(pch);
	if ( nbase > maxDP ) nbase = maxDP;
	nBases[0] = nbase;

	//notice("%s %d %c %d",chrom.c_str(), pos, ref[0], nBases[0]);

	if ( nbase == 0 ) return true;
	// parse bases
	nch = strchr(pch, '\t');
	p = pch = nch+1;
	j = 0;
	o = 0;
	for(k=0; j < nbase; ++k) {
	  switch(p[k]) {
	  case '$': // do nothing
	    break;
	  case '^': // ignore the next character
	    ++k;
	    break;
	  case '-': case '+':
	    l = atoi(&p[k+1]);
	    k += (l+(l<10 ? 1 : (l < 100 ? 2 : 3)));
	    break;
	  case '.': 
	    bases[o+j] = nref;
	    strands[o+j] = true;
	    ++j;
	    break;
	  case ',':
	    bases[o+j] = nref;
	    strands[o+j] = false;
	    ++j;
	    break;
	  case 'A':
	    bases[o+j] = 0;
	    strands[o+j] = true;
	    ++j;
	    break;
	  case 'a':
	    bases[o+j] = 0;
	    strands[o+j] = false;
	    ++j;
	    break;
	  case 'C':
	    bases[o+j] = 1;
	    strands[o+j] = true;
	    ++j;
	    break;
	  case 'c':
	    bases[o+j] = 1;
	    strands[o+j] = false;
	    ++j;
	    break;
	  case 'G':
	    bases[o+j] = 2;
	    strands[o+j] = true;
	    ++j;
	    break;
	  case 'g':
	    bases[o+j] = 2;
	    strands[o+j] = false;
	    ++j;
	    break;
	  case 'T':
	    bases[o+j] = 3;
	    strands[o+j] = true;
	    ++j;
	    break;
	  case 't':
	    bases[o+j] = 3;
	    strands[o+j] = false;
	    ++j;
	    break;
	  case 'N':
	    bases[o+j] = 4;
	    strands[o+j] = true;
	    ++j;
	  case '*': case 'n':
	    bases[o+j] = 4;
	    strands[o+j] = false;
	    ++j;
	    break;
	  default:
	    error("Cannot recognize %c in the pileups at i=%d, j= %d, nbase=%d",p[k],j,k);
	  }
	}

	p += k;
	while( ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	if ( *p == '\0' )  {
	  error("Cannot find BQ field");
	}

	++p;
	for(j=0; j < nbase; ++j, ++p) {
	  bQs[j] = *p-33;
	}
	while( ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	if ( *p == '\0' ) return true;

	hasMQ = true;

	++p;
	for(j=0; j < nbase; ++j, ++p) {
	  mQs[j] = (*p)-33;
	}
	while( ( *p != '\t' ) && ( *p != '\0' ) ) ++p;

	if ( *p == '\0' ) return true;

	//error("CY: %s",p);
	hasCY = true;

	// parse cycles;
	for(j=0; j < nbase; ++j) {
	  if ( j > 0 ) {
	    if ( *p != ',' ) return parseError("Comma",p);
	    ++p;
	  }
	  cycles[j] = strtol(p,&p,10);
	  /*
	  if ( j > 0 ) {
	    while(*p != ',') ++p;
	  }
	  ++p;
	  cycles[j] = atoi(p);
	  */
	}
      }
      //notice("%d\t%s\t%d\t%c\t%d",bp,chrom.c_str(),pos,ref,nbase);
    }
    return true;
  }

  bool parseError(const char* expected, const char* observed) {
    warning("%s expected but observed %c(%d) while parsing a line at %s:%d. Stop reading...",expected,(observed == NULL) ? 0 : *observed,(observed == NULL) ? -1 : *observed,chrom.c_str(),pos);
    line = NULL;
    pos = ENDPOS;
    for(int i=0; i < nInds; ++i) {
      nBases[i] = 0;
    }
    return false;
  }

  bool next() {
    if ( (line = (char*)tf.getLine()) != NULL ) {
      return parseMarker();
    }
    else {
      pos = ENDPOS;
      for(int i=0; i < nInds; ++i) {
	nBases[i] = 0;
      }
      return false;
    }
  }

  bool isEOF() const {
    return ( pos == ENDPOS );
  }

  bool advanceTo(const char* chr, int bp) {
    if ( ( chrom.compare(chr) != 0 ) && ( tf.getType() == 2 ) ) {
      char buf[256];
      sprintf(buf,"%s:%d",chr,bp);
      tf.updateRegion(buf);
      chrom = chr;
      pos = 0;
    }

    while ( ( bp > pos ) && ( line != NULL ) ) {
      next();
      char* p = strchr(line,'\t')+1;
      pos = strtol(p,&p,10); //atoi(p);
    }

    if ( line == NULL ) { 
      pos = ENDPOS; 
      for(int i=0; i < nInds; ++i) nBases[i] = 0;
      return false; 
    }
    else if ( pos > bp ) { return false; }
    else { return parseMarker(); }
  }

  void writeHeader(wFile& wf) {
    wf.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(int i=0; i < nInds; ++i) {
      wf.printf("\t%s",inds[i].c_str());
    }
    wf.printf("\n");
  }

  static void writeTabStr(wFile& wf, std::string& s) {
    if ( s.empty() ) { wf.printf("\t."); }
    else { wf.printf("\t%s",s.c_str()); }
  }

  void writeMarker(wFile& wf, bool tosingle = false) {
    if ( tosingle ) {
      int i,j,k,o,nb,nref;

      wf.printf("%s\t%d",chrom.c_str(),pos);
      writeTabStr(wf, ref);
      nref = BaseAsciiMap::base2int[(int)ref[0]];
      nb = 0;
      for(i=0; i < nInds; ++i) nb += nBases[i];
      wf.printf("\t%d",nb);

      wf.printf("\t");
      for(i=0; i < nInds; ++i) {
	o = i*(maxDP+1);
	for(j=0; j < (int)nBases[i]; ++j) {
	  if ( bases[o+j] == nref ) {
	    if ( strands[o+j] ) wf.printf(".");
	    else wf.printf(",");
	  }
	  else {
	    wf.printf("%c",strands[o+j] ? BaseAsciiMap::int2base[(int)bases[o+j]] : BaseAsciiMap::int2basel[(int)bases[o+j]]);
	  }
	}
      }

      wf.printf("\t");
      for(i=0; i < nInds; ++i) {
	o = i*(maxDP+1);
	for(j=0; j < (int)nBases[i]; ++j) {
	  wf.printf("%c",(char)(33+bQs[o+j]));
	}
      }

      if ( hasMQ ) {
	wf.printf("\t");
	for(i=0; i < nInds; ++i) {
	  o = i*(maxDP+1);
	  for(j=0; j < (int)nBases[i]; ++j) {
	    wf.printf("%c",(char)(33+mQs[o+j]));
	  }
	}
      }

      if ( hasCY ) {
	wf.printf("\t");
	for(i=0, k=0; i < nInds; ++i) {
	  o = i*(maxDP+1);
	  for(j=0; j < (int)nBases[i]; ++j, ++k) {
	    if ( k > 0 ) wf.printf(",");
	    wf.printf("%d",cycles[o+j]);
	  }
	}
      }

      wf.printf("\n");
    }
    else {
      int i,j,o;

      wf.printf("%s\t%d",chrom.c_str(),pos);
      writeTabStr(wf, id);
      writeTabStr(wf, ref);
      writeTabStr(wf, alt);
      writeTabStr(wf, qual);
      writeTabStr(wf, filter);
      writeTabStr(wf, info);
      wf.printf("\tN:RD:BQ");
      if ( hasMQ ) wf.printf(":MQ");
      if ( hasCY ) wf.printf(":CY");

      for(i=0; i < nInds; ++i) {
	if ( nBases[i] > 0 ) {
	  wf.printf("\t%d",nBases[i]);
	  o = i*(maxDP+1);
	  wf.printf(":");
	  for(j=0; j < (int)nBases[i]; ++j) {
	    wf.printf("%c",strands[o+j] ? BaseAsciiMap::int2base[(int)bases[o+j]] : BaseAsciiMap::int2basel[(int)bases[o+j]]);
	  }
	  wf.printf(":");
	  for(j=0; j < (int)nBases[i]; ++j) {
	    if ( j > 0 ) wf.printf(",");
	    wf.printf("%d",bQs[o+j]);
	  }
	  if ( hasMQ ) {
	    wf.printf(":");
	    for(j=0; j < (int)nBases[i]; ++j) {
	      if ( j > 0 ) wf.printf(",");
	      wf.printf("%d",mQs[o+j]);
	    }
	  }
	  if ( hasCY ) {
	    wf.printf(":");
	    for(j=0; j < (int)nBases[i]; ++j) {
	      if ( j > 0 ) wf.printf(",");
	      wf.printf("%d",cycles[o+j]);
	    }
	  }
	}
	else {
	  wf.printf("\t.");
	}
      }
      wf.printf("\n");
    }
  }
};
#endif
