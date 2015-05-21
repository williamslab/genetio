// Library for I/O of genetic data
// Author: Amy Williams <alw289 @ cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "marker.h"


////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<Marker *> Marker::_allMarkers(600000);
dynarray<int> Marker::_omitMarkers;
bool Marker::_readOnlyOneChrom = false;
int Marker::_firstStoredMarkerIdx = -1; // by default -1 => not applicable
int Marker::_numMarkersInFile = 0; // gets updated as we read the file
dynarray<char *> Marker::_chromNames;
dynarray<int> Marker::_firstMarkerNum;
dynarray<int> Marker::_lastMarkerNum;
dynarray<int> Marker::_firstHapChunk;
dynarray<int> Marker::_lastHapChunk;
int Marker::_numHapChunks = 0;
dynarray<int>    Marker::_hapWindowEnds;
dynarray<float> Marker::_hapWindowMapCenter;

// Read Reich lab format .snp file
void Marker::readSNPFile(const char *snpFile, const char *onlyChr, int startPos,
			 int endPos) {
  FILE *in = fopen(snpFile, "r");
  if (!in) {
    fprintf(stderr, "\n\nERROR: Couldn't open SNP file %s\n", snpFile);
    perror(snpFile);
    exit(2);
  }

  readMarkers(in, onlyChr, /*type=*/ 1, startPos, endPos);
  fclose(in);
}

// Read PLINK format .map file
void Marker::readMapFile(const char *mapFile, const char *onlyChr, int startPos,
			 int endPos) {
  FILE *in = fopen(mapFile, "r");
  if (!in) {
    fprintf(stderr, "\n\nERROR: Couldn't open map file %s\n", mapFile);
    perror(mapFile);
    exit(2);
  }

  readMarkers(in, onlyChr, /*type=*/ 2, startPos, endPos);
  fclose(in);
}

// Read PLINK format .bim file
void Marker::readBIMFile(const char *bimFile, const char *onlyChr, int startPos,
			 int endPos) {
  FILE *in = fopen(bimFile, "r");
  if (!in) {
    fprintf(stderr, "\n\nERROR: Couldn't open BIM file %s\n", bimFile);
    perror(bimFile);
    exit(2);
  }

  readMarkers(in, onlyChr, /*type=*/ 3, startPos, endPos);
  fclose(in);
}

// Read markers from VCF file (which is required to be gzipped)
void Marker::readVCFFile(htsFile *vcfIn, tbx_t *index, hts_itr_t *itr,
			 int startPos, int endPos) {
  std::string markerName;
  std::string alleles;
  std::string tmpStr;
  int chromIdx = -1;
  Marker *prevMarker = NULL;
  int physPos;

  int numMarkersCurChrom = 0;

  // Start parsing the VCF using HTSlib; gives one line at a time:
  // Go through all the lines in the query region
  while (tbx_itr_next(vcfIn, index, itr, &vcfIn->line) >= 0) {
    // string for the current line is in vcfIn->line.s
    // parse the parts of the string having to do with the marker:

    // read chromosome/contig name:
    int c;
    tmpStr.clear();
    for(c = 0; vcfIn->line.s[c] != '\t'; c++) // to end of the field
      tmpStr += vcfIn->line.s[c];

    // Note: because we're using the HTS iterator and have given it a region to
    // iterate over (in PersonIO::readVCF()), we don't need to check whether the
    // chromosome is relevant here

    // same chromosome as the previous marker?
    if (chromIdx < 0 || strcmp(tmpStr.c_str(), _chromNames[chromIdx]) != 0) {
      // new chromosome
      if (chromIdx >= 0 && _readOnlyOneChrom) {
	printf("\n\n");
	printf("ERROR: markers present from multiple chromosomes/contigs.\n");
	printf("Please specify a chromosome to process with --chr\n");
	exit(1);
      }

      int length = c+1; // + 1 due to '\0'
      char *newChrom = new char[length];
      strcpy(newChrom, tmpStr.c_str());
      _chromNames.append(newChrom);
      chromIdx = _chromNames.length() - 1;

      _firstMarkerNum.append(0);
      _lastMarkerNum.append(-1);
      _firstHapChunk.append(0);
      _lastHapChunk.append(-1);
    }

    // read physical position:
    int s = c + 1; // shift to the start of this field
    tmpStr.clear();
    for(c = 0; vcfIn->line.s[s+c] != '\t'; c++) // to end of the field
      tmpStr += vcfIn->line.s[s+c];
    s += c + 1; // update shift value to the start of the next field
    physPos = atoi(tmpStr.c_str());

    // Again, due to iterators, we are guaranteed to be in range of startPos and
    // endPos
    assert(physPos >= startPos && physPos <= endPos);

    // read marker name:
    markerName.clear();
    for(c = 0; vcfIn->line.s[s+c] != '\t'; c++) // to end of the field
      markerName += vcfIn->line.s[s+c];
    s += c + 1; // update shift

    // read reference allele:
    alleles.clear();
    for(c = 0; vcfIn->line.s[s+c] != '\t'; c++)
      alleles += vcfIn->line.s[s+c];
    alleles += ' ';
    s += c + 1; // update shift

    // read alternate alleles:
    int numAlleles = 1;
    for(c = 0; vcfIn->line.s[s+c] != '\t'; c++) {
      if (vcfIn->line.s[s+c] == ',')
	alleles += ' ';
      else
	alleles += vcfIn->line.s[s+c];
      numAlleles++;
    }
    s += c + 1;

    // TODO: check if the Person type P can handle this number of alleles

    ///////////////////////////////////////////////////////////////////////////
    // Have data read, now make the Marker and do bookkeeping as needed

    // Note: for VCFs, since we're using iterators, we don't need
    // _firstStartMarkerIdx?

    // Note: don't have genetic map distance here, so set from physical:
    float morganDistToPrev = 1.0f;
    if (prevMarker != NULL && chromIdx == prevMarker->_chromIdx) {
      morganDistToPrev = (physPos - prevMarker->getPhysPos()) / (100 * 1000000);
    }

    // TODO: need the ability to read genetic map from other files, require it
    // for VCFs

    Marker *m = new Marker(markerName.c_str(), chromIdx, /*mapPos=*/ 0.0f,
			   morganDistToPrev, physPos, alleles.c_str(),
			   numAlleles);

    if (prevMarker != NULL) {
      int prevChromIdx = prevMarker->_chromIdx;
      if (chromIdx == prevChromIdx) {
	if (physPos == prevMarker->_physPos) {
	  fprintf(stderr,
		  "WARNING: marker %s has same position as previous marker\n",
		  markerName.c_str());
	}
      }
      else { // Have a valid prev chrom?  Update marker/chunk counts
	// Note: numMarkersCurChrom actually applies to the previous chrom
	updateInfoPrevChrom(prevChromIdx, numMarkersCurChrom);

	numMarkersCurChrom = 0; // reset

	int curMarkerNum = _allMarkers.length();
	// about to append the first marker for this chrom to this index:
	_firstMarkerNum[chromIdx] = curMarkerNum;
      }
    }
    
    _allMarkers.append(m);
    numMarkersCurChrom++;
    prevMarker = m;
  }

  // Set starting chunk for final chromosome:
  if (prevMarker != NULL) {
    int prevChromIdx = prevMarker->_chromIdx;
    updateInfoPrevChrom(prevChromIdx, numMarkersCurChrom);
  }
}

// Prints a .phsnp/.snp (PackedAncestryMap/Eigenstrat format) file
void Marker::printSNPFile(FILE *out) {
  int numMarkers = _allMarkers.length();
  for(int m = 0; m < numMarkers; m++) {
    Marker *cur = _allMarkers[m];

    if (cur->getNumAlleles() > 2)
      // these .snp files only support biallelics: skip this marker
      continue;

    // print SNP id
    int padSpaces = 20 - strlen(cur->_name);
    for(int i = 0; i < padSpaces; i++)
      fprintf(out, " ");
    fprintf(out, "%s", cur->_name);

    // print chromosome
    padSpaces = 4 - strlen(cur->getChromName());
    for(int i = 0; i < padSpaces; i++)
      fprintf(out, " ");
    fprintf(out, "%s", cur->getChromName());

    // print genetic position
    fprintf(out, "        %1.12f", cur->_mapPos);

    // print physical position
    fprintf(out, "       ");
    for(int cap = 100000000; cap > 0 && cur->_physPos < cap; cap /= 10)
      fprintf(out, " ");
    fprintf(out, "%d", cur->_physPos);

    // print alleles
    fprintf(out, " %s\n", cur->_alleles);
  }
}

// Prints the first 5 columns of an IMPUTE2 format .haps file (i.e., SNP
// information)
void Marker::printImpute2Prefix(FILE *out, int markerNum) {
  Marker *cur = _allMarkers[markerNum];

  if (cur->getNumAlleles() > 2)
    // IMPUTE2 format only supports biallelic SNPs: skip
    return;

  // IMPUTE2 format uses the opposite numerical encoding to Eigenstrat and
  // packed Ancestrymap formats, so we flip the allele order here:
  fprintf(out, "%s %s %d %s", cur->getChromName(), cur->getName(),
	  cur->getPhysPos(), cur->_alleles);
}

// Print, in gzipped format, the first 5 columns of an IMPUTE2 format .haps
// file (i.e., SNP information)
void Marker::printGzImpute2Prefix(gzFile out, int markerNum) {
  Marker *cur = _allMarkers[markerNum];

  if (cur->getNumAlleles() > 2)
    // IMPUTE2 format only supports biallelic SNPs: skip
    return;

  // IMPUTE2 format uses the opposite numerical encoding to Eigenstrat and
  // packed Ancestrymap formats, so we flip the allele order here:
  gzprintf(out, "%s %s %d %s", cur->getChromName(), cur->getName(),
	   cur->getPhysPos(), cur->_alleles);
}

// TODO: comment, including on return value
bool skipWhitespace(char *curBuf, int &bind, size_t nread, const int BUF_SIZE) {
  for ( ; (curBuf[bind] == ' ' ||curBuf[bind] == '\t') && bind < nread; bind++);
  if (bind == nread) {
    if (nread < BUF_SIZE) // done reading file
      return -1; // EOF (or potentially error) reached
    else { // more to read
      // TODO!
//      nread = fread(curBuf, sizeof(char), BUF_SIZE, in);
      return 0;
    }
  }
  else if (curBuf[bind] == EOF) // done reading file? end loop
    return -1;

  return 1;
}

// Helper function for setting appropriate null character points
// TODO : naming of function
// Switch "buf_ind " to "bind"
int readDoubleBuffer(FILE *in, char *field, char *curBuf, char *nextBuf, int BUF_SIZE, int buf_ind, size_t nread){
    // First skip leading whitespace...
    int status = skipWhitespace(curBuf, buf_ind, nread, BUF_SIZE);
    if (status < 0) return status; // We have reached EOF
    int mstart = buf_ind;
    field = &curBuf[mstart];
    for ( ; !isspace(curBuf[buf_ind]) && buf_ind < nread; buf_ind++)
    ;	 
    if (buf_ind == nread) {
      // reached end of curBuf and markerName is incomplete; copy into
      // <nextBuf> and then read more into that buffer
      int numCpy = nread - mstart;
      strncpy(nextBuf, field, nread - mstart);
      fread(&nextBuf[numCpy], sizeof(char), BUF_SIZE - numCpy, in);
      field = &nextBuf[0];

      char *tmpBuf = curBuf;
      curBuf = nextBuf;
      nextBuf = tmpBuf;
      // now get the end of the field
      buf_ind = numCpy;
      for ( ; !isspace(curBuf[buf_ind]) && buf_ind < nread; buf_ind++);
    }

		// TODO: print error message if this doesn't hold?
    assert(buf_ind < nread); 

    // null terminate field by inserting '\0' in curBuf:
    curBuf[buf_ind] = '\0';  

    // Return the current buffer index so that we can pick up where we left off
    return buf_ind;
}


// Read marker/genetic map definition file of the following formats:
// If type == 1, reads Reich lab format .snp file
// If type == 2, reads PLINK format .map file
// If type == 3, reads PLINK format .bim file
void Marker::readMarkers(FILE *in, const char *onlyChr, int type, int startPos,
			 int endPos) {
  const int BUF_SIZE = 2048;
  char buf1[BUF_SIZE], buf2[BUF_SIZE];
  char *curBuf, *nextBuf;
  size_t nread; // number of chars read into <curBuf>
  int bind = 0; // current buffer index in <curBuf> (during parsing below)
  // TODO: we could realistically only have one "field pointer right?"
  // aab227 : set all to null so no warnings
  char *markerNameX = NULL;
  char *chromNameX = NULL;
  char *tmpStrX = NULL;
  std::string markerName;
  std::string chromName;
  std::string tmpStr;
  int chromIdx = -1;
  Marker *prevMarker = NULL;
  float mapPos;
  int physPos;
  char alleles[4] = "Z Z"; // initially

  int numMarkersCurChrom = 0;

  // set genetic positions from physical?  Yes if all genetic positions are 0
  int setGenetFromPhys = -1;

  curBuf = buf1;
  nextBuf = buf2;
  nread = fread(curBuf, sizeof(char), BUF_SIZE, in);

  while (1) {
    // Note: I assume the map positions are in Morgans per the spec of both
    // the Reich lab SNP file format and the spec of the PLINK .map file format
    if (type == 1) {

      // get the marker name
      bind = readDoubleBuffer(in, markerNameX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      // Actually setting the markername from the char pointer
      markerName.assign(markerNameX);

      // get the chromosome name
      bind = readDoubleBuffer(in, chromNameX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      // setting the chomosome name
      chromName.assign(chromNameX);

      // get the genetic map position
      bind = readDoubleBuffer(in, tmpStrX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      // Reading in the float for map posittion
      tmpStr.assign(tmpStrX);
      mapPos = atof(tmpStr.c_str());

      // Get the physical position of the variant
      bind = readDoubleBuffer(in, tmpStrX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      tmpStr.assign(tmpStrX);
      physPos = atoi(tmpStr.c_str());

      // // // readToken(in, markerName); // read marker name
      // // if (fgetc(in) == EOF) break; // done reading file 

      // readToken(in, chromName);  // read chromosome name
      // readToken(in, tmpStr);     // read map position, then convert to float
      // mapPos = atof(tmpStr.c_str());
      // readToken(in, tmpStr);     // read physical position, then convert to int
      // physPos = atoi(tmpStr.c_str());

      // this loop ends when the allele gets read:
      while(isspace(alleles[0] = fgetc(in))); // read first allele
      char c = fgetc(in);
      if (c != ' ' && c != '\t') {
				fprintf(stderr, "ERROR: alleles expected to be single characters\n");
				fprintf(stderr, "At marker %s\n", markerName.c_str());
				exit(1);
      }
      // Note: index 1 is a space
      while(isspace(alleles[2] = fgetc(in))); // read second allele

      // read to the end of the line
      while (isspace(c = fgetc(in)) && c != '\n' && c != EOF);

      if (c == EOF) break; // done reading file
      else if (c != '\n') {
				fprintf(stderr, "ERROR: extra characters on line for marker %s\n",
					markerName.c_str());
				exit(1);
      }
    }
    else if (type == 2 || type == 3) {
      // TODO : remove dependence on this single character...
      char c; 
      // c = readToken(in, chromName); // read chromosome
      // if (c == EOF) // done reading file

      // read in the chromosome name
      bind = readDoubleBuffer(in, chromNameX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      chromName.assign(chromNameX);

      // read in the marker name
      bind = readDoubleBuffer(in, markerNameX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      markerName.assign(markerNameX);

      // readToken(in, markerName); // read marker name
      // readToken(in, tmpStr);     // read map position, then convert to float
      // mapPos = atof(tmpStr.c_str());

      // read in genetic map position
      bind = readDoubleBuffer(in, tmpStrX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      tmpStr.assign(tmpStrX);
      mapPos = atof(tmpStr.c_str());

      // readToken(in, tmpStr);     // read physical position, then convert to int
      // physPos = atoi(tmpStr.c_str());

      // read in physical position
      bind = readDoubleBuffer(in, tmpStrX, curBuf, nextBuf, BUF_SIZE, bind, nread);
      if (bind < 0) break;
      tmpStr.assign(tmpStrX);
      physPos = atof(tmpStr.c_str());


      if (type == 3) { // for .bim files, must read alleles
	// this loop ends when the allele gets read:
	while(isspace(alleles[0] = fgetc(in))); // read first allele
	c = fgetc(in);
	if (c != ' ' && c != '\t') {
	  fprintf(stderr, "ERROR: alleles expected to be single characters\n");
	  fprintf(stderr, "At marker %s\n", markerName.c_str());
	  exit(1);
	}
	// Note: index 1 is a space
	while(isspace(alleles[2] = fgetc(in))); // read second allele
      }

      // read to the end of the line
      while (isspace(c = fgetc(in)) && c != '\n' && c != EOF);

      if (c == EOF)
	break; // done reading file
      else if (c != '\n') {
	fprintf(stderr, "ERROR: extra characters on line for marker %s\n",
		markerName.c_str());
	exit(1);
      }
    }
    else {
      fprintf(stderr, "ERROR: unknown marker file type %d!\n", type);
      exit(2);
    }

    _numMarkersInFile++;

    if (onlyChr != NULL && strcmp(chromName.c_str(), onlyChr) != 0) {
      // only keeping markers on <onlyChr> -- skip
      continue;
    }

    // Is the physical position missing?  If so, omit this marker
    if (physPos == 0) {
      fprintf(stderr, "Omitting marker %s: missing physical position\n",
	      markerName.c_str());
      if (_allMarkers.length() > 0) {
	// only track omit markers when we have some markers prior to it that
	// we *are* storing (i.e., when _allMarkers.length() > 0); those
	// previous to the markers that we are storing will be omitted using
	// the mechanism that skips markers we don't store
	int curIndex = _allMarkers.length();
	_omitMarkers.append(curIndex);
      }
      continue;
    }

    if (physPos < startPos || physPos > endPos) {
      // marker outside of range to be inspected
      continue;
    }

    // same chromosome as the previous marker?
    if (chromIdx < 0 || strcmp(chromName.c_str(), _chromNames[chromIdx]) != 0) {
      // new chromosome
      if (chromIdx >= 0 && _readOnlyOneChrom) {
	printf("\n\n");
	printf("ERROR: markers present from multiple chromosomes/contigs.\n");
	printf("Please specify a chromosome to process with --chr\n");
	exit(1);
      }

      char *newChrom = new char[ chromName.length() + 1 ]; // + 1 due to '\0'
      strcpy(newChrom, chromName.c_str());
      _chromNames.append(newChrom);
      chromIdx = _chromNames.length() - 1;

      _firstMarkerNum.append(0);
      _lastMarkerNum.append(-1);
      _firstHapChunk.append(0);
      _lastHapChunk.append(-1);
    }

    if (_allMarkers.length() == 0) {
      _firstStoredMarkerIdx = _numMarkersInFile - 1;
    }

    float morganDistToPrev = 1.0f;
    if (prevMarker != NULL && chromIdx == prevMarker->_chromIdx) {
      // on second marker? update setGenetFromPhys
      if (_allMarkers.length() == 1) {
	if (prevMarker->getMapPos() == 0.0f && mapPos == 0.0f) {
	  setGenetFromPhys = 1;
	  printf("WARNING: Setting genetic position from physical\n");
	}
	else
	  setGenetFromPhys = 0;
      }
      if (setGenetFromPhys) {
	if (mapPos != 0.0f) {
	  fprintf(stderr, "ERROR: some markers have non-zero genetic position and cannot be used\n");
	  fprintf(stderr, "To force them to be dropped, set their physical position to 0\n");
	  exit(1);
	}
	// 1 cM per Mb
	morganDistToPrev = (physPos - prevMarker->getPhysPos()) /
								(100 * 1000000);
      }
      else {
	morganDistToPrev = mapPos - prevMarker->getMapPos();
      }
    }
    // Note: if alleles == NULL, the numAlleles argument will be ignored
    Marker *m = new Marker(markerName.c_str(), chromIdx, mapPos,
			   morganDistToPrev, physPos,
			   (type == 2) ? NULL : alleles, /*numAlleles=*/ 2);

    if (prevMarker != NULL) {
      int prevChromIdx = prevMarker->_chromIdx;
      if (chromIdx == prevChromIdx) {
	if (mapPos < prevMarker->_mapPos || physPos < prevMarker->_physPos) {
	  fprintf(stderr,
		  "ERROR: marker %s has position before previous marker!\n",
		  markerName.c_str());
	  exit(1);
	}
	else if (physPos == prevMarker->_physPos) {
	  fprintf(stderr,
		  "WARNING: marker %s has same position as previous marker\n",
		  markerName.c_str());
	}

	// For hotspot-based windows selection
//	float genetDistToPrev = mapPos - prevMarker->_mapPos;
//	int distToPrev = physPos - prevMarker->_physPos;
//	float MorgPerMb = genetDistToPrev / ((float) distToPrev / 1000000);
//	if (MorgPerMb >= .1) { // >= 10 cM/Mb?
//	  // hotspot!
//	  int curMarkerNum = _allMarkers.length();
//	  // note: window *won't* include the current marker
//	  setNumMarkersInWindow(windowStartIdx, curMarkerNum - windowStartIdx);
//	  float curWindowDist = mapPos - windowStartMorg;
//	  fprintf(stderr, "%lf %d\n", curWindowDist,
//		  curMarkerNum - windowStartIdx);
//	  windowStartIdx = curMarkerNum;
//	  windowStartMorg = mapPos;
//	}
      }
      else { // Have a valid prev chrom?  Update marker/chunk counts
	// Note: numMarkersCurChrom actually applies to the previous chrom
	updateInfoPrevChrom(prevChromIdx, numMarkersCurChrom);

	numMarkersCurChrom = 0; // reset

	int curMarkerNum = _allMarkers.length();
	// about to append the first marker for this chrom to this index:
	_firstMarkerNum[chromIdx] = curMarkerNum;
      }
    }

    _allMarkers.append(m);
    numMarkersCurChrom++;
    prevMarker = m;
  }

  // Set starting chunk for final chromosome:
  if (prevMarker != NULL) {
    int prevChromIdx = prevMarker->_chromIdx;
    updateInfoPrevChrom(prevChromIdx, numMarkersCurChrom);
  }
}


// Updates the size and location of windows so that the first window starts
// at <initOffset> and each subsequent has size <windowNumMarkers>.
void Marker::updateWindows(int initOffset, int windowNumMarkers) {
  _hapWindowEnds.clear();
  _hapWindowMapCenter.clear();
  if (initOffset > 0)
    setNumMarkersInWindow(0, /*numMarkers=*/ initOffset);
  // NOTE: this is *not* chromosome-aware, but that's OK since we only phase one
  // chromosome at a time.
  int totalMarkers = _allMarkers.length();
  for(int m = initOffset; m < totalMarkers; m += windowNumMarkers) {
    if (m + windowNumMarkers > totalMarkers)
      setNumMarkersInWindow(m, totalMarkers - m);
    else
      setNumMarkersInWindow(m, windowNumMarkers);
  }
}

// This is no longer used:
// Updates the size and location of windows based on genetic distance, requiring
// each window to be at least windowLengthMorgans long in genetic distance and
// to have a minimum of minNumMarkers.
void Marker::updateWindowsMap(int initOffset, float windowLengthMorgans,
			      int minNumMarkers) {
  _hapWindowEnds.clear();
  _hapWindowMapCenter.clear();
  if (initOffset > 0)
    setNumMarkersInWindow(0, /*numMarkers=*/ initOffset);

  int windowStartIdx = 0;
  float windowStartMorg = 0.0f;
  const Marker *prevMarker = NULL;

  int totalMarkers = _allMarkers.length();
  for(int m = initOffset; m < totalMarkers; m++) {
    const Marker *curMarker = Marker::getMarker(m);
    float curMapPos = curMarker->getMapPos();

    if (m == initOffset) {
      windowStartMorg = curMapPos;
    }
    else {
      int prevChromIdx = prevMarker->_chromIdx;
      if (curMarker->_chromIdx != prevChromIdx) {
	assert(curMarker->_chromIdx > prevChromIdx);

	// note: window *won't* include the current marker
	setNumMarkersInWindow(windowStartIdx, m - windowStartIdx);
	windowStartIdx = m;
	windowStartMorg = curMapPos;
      }
    }

    if (curMapPos - windowStartMorg > windowLengthMorgans &&
					  m - windowStartIdx >= minNumMarkers) {
      int numMarkers = m - windowStartIdx;

      if ((uint32_t) numMarkers > 2 * BITS_PER_CHUNK) {
	// find the ideal number of windows to divided this into
	int numWins = 2;
	for( ; ; numWins++) {
	  if ((uint32_t) numMarkers / numWins <= 2 * BITS_PER_CHUNK &&
	      (uint32_t) numMarkers - (numWins - 1) * (numMarkers / numWins) <=
							     2 * BITS_PER_CHUNK)
	    break;  // find ideal number of windows (numWins)
	}
	int markersPerWin = numMarkers / numWins;
	for(int i = 0; i < numWins - 1; i++) {
	  setNumMarkersInWindow(windowStartIdx, markersPerWin);
	  windowStartIdx += markersPerWin;
	}
	// last window has different number:
	int markersLastWin = numMarkers - (numWins - 1) * markersPerWin;
	setNumMarkersInWindow(windowStartIdx, markersLastWin);
	windowStartIdx += markersLastWin;

	assert(windowStartIdx == m);
	windowStartMorg = curMapPos;
      }
      else {
	// note: window *won't* include the current marker
	setNumMarkersInWindow(windowStartIdx, numMarkers);
//	float curWindowDist = prevMarker->_mapPos - windowStartMorg;
//	fprintf(stderr, "%lf %d\n", curWindowDist, numMarkers);
	windowStartIdx = m;
	windowStartMorg = curMapPos;
      }
    }

    prevMarker = curMarker;
  }

  // make window for last few markers:
  int curMarkerNum = _allMarkers.length();
  setNumMarkersInWindow(windowStartIdx, curMarkerNum - windowStartIdx);
}

// Sets the last marker number for <prevChromIdx> along with the first and last
// chunk numbers
void Marker::updateInfoPrevChrom(int prevChromIdx, int numMarkersPrevChrom) {
  _lastMarkerNum[prevChromIdx] = _allMarkers.length() - 1;
  // Note: _numHapChunks is presently 1 more than the final chunk index for the
  // chromosome before <prevChromIdx>, i.e., exactly the first index we need
  // here:
  _firstHapChunk[prevChromIdx] = _numHapChunks;
  // Update the total number of haplotype chunks:
  _numHapChunks += getNumHapChunksFor(numMarkersPrevChrom);
  _lastHapChunk[prevChromIdx] = _numHapChunks - 1;
}

// Stores the number of markers within the 0.25cM window: corrects for marker
// density/sparsity in a region and is a somewhat hacky way of correcting for
// LD (each 0.25cM block is treated as one marker).
void Marker::setNumMarkersInWindow(int startMarkerNum, int numMarkers) {
  int endMarker = startMarkerNum + numMarkers - 1;
  for(int i = startMarkerNum; i <= endMarker; i++) {
    _allMarkers[i]->_numSNPsWindow = numMarkers;
  }
  // add end point for this window:
  _hapWindowEnds.append(endMarker);
  float mapCenter = (Marker::getMarker(startMarkerNum)->getMapPos() +
			      Marker::getMarker(endMarker)->getMapPos()) / 2.0f;
  _hapWindowMapCenter.append(mapCenter);
}

// Returns the number of haplotype chunks required to store <numMarkers>
int Marker::getNumHapChunksFor(int numMarkers) {
  int numChunks = numMarkers / BITS_PER_CHUNK;
  if (numMarkers % BITS_PER_CHUNK > 0)
    numChunks++;

  return numChunks;
}

// Returns the number of markers not divisible by the num of bits in a chunk
int Marker::getChunkModMarkers(int numMarkers) {
  return numMarkers % BITS_PER_CHUNK;
}

int Marker::getFirstMarkerNumForChunk(int chromIdx, int chunkNum) {
  int numChunksIntoChrom = chunkNum - getFirstHapChunk(chromIdx);
  return getFirstMarkerNum(chromIdx) + (BITS_PER_CHUNK * numChunksIntoChrom);
}

// Returns the total physical length of the markers that were input
//uint32_t Marker::getTotalPhysLength(bool analyzeChrX) {
//  uint32_t total = 0;
//  // Note: we only include the autosomes and optionally the X chromosome in this
//  // we don't include Y, PAR, or MT
//  int lastChr = CHR_LAST_AUTOSOME;
//  if (analyzeChrX)
//    lastChr = CHR_X;
//  for (int chrom = 1; chrom <= lastChr; chrom++) {
//    if (Marker::getNumChromMarkers(chrom) == 0)
//      continue;
//
//    const Marker *firstMarker =
//			  Marker::getMarker( Marker::getFirstMarkerNum(chrom) );
//    const Marker *lastMarker =
//			  Marker::getMarker( Marker::getLastMarkerNum(chrom) );
//    total += lastMarker->getPhysPos() - firstMarker->getPhysPos();
//  }
//  return total;
//}

// Returns the total genetic length (in Morgans) of the markers that were input
//float Marker::getTotalGenetLength(bool analyzeChrX) {
//  float total = 0.0f;
//  // Note: we only include the autosomes and optionally the X chromosome in this
//  // we don't include Y, PAR, or MT
//  int lastChr = CHR_LAST_AUTOSOME;
//  if (analyzeChrX)
//    lastChr = CHR_X;
//  for (int chrom = 1; chrom <= lastChr; chrom++) {
//    if (Marker::getNumChromMarkers(chrom) == 0)
//      continue;
//
//    const Marker *firstMarker =
//			  Marker::getMarker( Marker::getFirstMarkerNum(chrom) );
//    const Marker *lastMarker =
//			  Marker::getMarker( Marker::getLastMarkerNum(chrom) );
//    total += lastMarker->getMapPos() - firstMarker->getMapPos();
//  }
//  return total;
//}

Marker::Marker(const char *markerName, int chromIdx, float mapPos,
	       float morganDistToPrev, int physPos, const char *alleles,
	       int numAlleles) {
  _name = new char[strlen(markerName)+1];
  strcpy(_name, markerName);
  _chromIdx = chromIdx;
  _physPos = physPos;
  _mapPos = mapPos;

  if (alleles != NULL) {
    _alleles = new char[ strlen(alleles) + 1 ]; // +1 for '\0'
    strcpy(_alleles, alleles);
    _numAlleles = numAlleles;
  }
  else {
    _numAlleles = 0;
  }
}

// Compute allele frequency stats
void Marker::setAlleleFreq(int alleleCount, int totalGenoWithData, 
					bool nonStandardGeno) {
	if (nonStandardGeno) {
	   _logAlleleFreq = _logVarAlleleFreq = 1; // illegal log value
	  return;
  }
  // calling allele 1 the variant allele, though it need not be:
  float variantFrequency = (float) alleleCount / (2 * totalGenoWithData);
  float referenceFrequency = 1 - variantFrequency;
  if (totalGenoWithData != 0)
    // if we have any data, we should have a valid frequency; otherwise it's
    // nan
    assert(0.0 <= variantFrequency && variantFrequency <= 1.0);

  _logAlleleFreq = log( referenceFrequency );
  _logVarAlleleFreq = log( variantFrequency );
}