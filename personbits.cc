// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "marker.h"
#include "personbits.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
int PersonBits::_numDuos = 0;
int PersonBits::_numTrioKids = 0;
dynarray<PersonBits *> PersonBits::_allIndivs;
Hashtable<char *, PersonBits *> PersonBits::_idToPerson(2003, stringHash,
							stringcmp);

PersonBits::PersonBits(char *id, char sex, int popIndex,
		       short familyIdLength) :
		       SuperPerson(id, sex, popIndex, familyIdLength) {
  if (!_ignore) {
    // Will update _trioDuoType and _tdData later as we read in relationships
    _trioDuoType = UNRELATED;
    _tdData = NULL;

    int numHapChunks = Marker::getNumHapChunks();
    assert(numHapChunks > 0);
    _homozyLoci = new chunk[numHapChunks];
    _knownHap = new chunk[numHapChunks];
    _missingLoci = new chunk[numHapChunks];
//    _expDiffBit = 0;

    // initialize these arrays:
    for(int i = 0; i < numHapChunks; i++) {
      _homozyLoci[i] = _knownHap[i] = _missingLoci[i] = 0;
    }

    _resolvedHaplotype[0] = _resolvedHaplotype[1] = NULL;

    _sampledHaplotypes = NULL;

    if (_idToPerson.lookup(_id)) {
      fprintf(stderr, "\nERROR: multiple individuals with id %s!\n", _id);
      exit(3);
    }

    _idToPerson.add(_id, this);
  }
}

PersonBits::~PersonBits() {
  if (!_ignore) {
    delete [] _homozyLoci;
    delete [] _knownHap;
    delete [] _missingLoci;
    assert(_tdData == NULL || getTrioDuoType() == TRIO_CHILD);
//    if (_tdData != NULL) {
//      delete [] _tdData->_tdKnownLoci;
//      if (_tdData->_childIsHet != NULL)
//	delete [] _tdData->_childIsHet;
//      delete _tdData;
//    }
    assert(_resolvedHaplotype[0] == NULL);
    assert(_sampledHaplotypes == NULL);
  }
}

// Deletes genotype data for <this>.  Used for trio children after their data
// were already used to infer the known phase of their parents
void PersonBits::empty() {
  delete [] _homozyLoci;
  delete [] _knownHap;
  delete [] _missingLoci;
  _homozyLoci = _knownHap = _missingLoci = NULL;
}

// Given the parents of <this>, sets the relationships, and infers the trio or
// duo phase.
// See comment above inferTrioDuoHaplotypes() regarding <numMendelError> and
// <numMendelCounted>
void PersonBits::setParents(char *familyid, PersonBits *parents[2],
			    int numParents, bool &warningPrinted, FILE *log,
			    int *numMendelError, int *numMendelCounted) {
  if (numParents == 2) {
    setTrioDuoType(TRIO_CHILD);
    _numTrioKids++;
  }
  else {
    setTrioDuoType(DUO_CHILD);
    _numDuos++;
  }

  for(int p = 0; p < 2; p++) {
    if (parents[p] == NULL)
      continue;

    // ensure this parent isn't part of another trio/duo relationship:
    if (parents[p]->getTrioDuoType() != UNRELATED) {
      if (!warningPrinted) {
	printf("\n");
	if (log != NULL) fprintf(log, "\n");
	warningPrinted = true;
      }
      fprintf(stderr, "ERROR: parent %s is a member of another trio or duo\n",
	      parents[p]->getId());
      exit(1);
    }

    if (numParents == 1)
      parents[p]->setTrioDuoType(PARENT_0);
    else {
      if (p == 0)
	parents[p]->setTrioDuoType(PARENT_0);
      else
	parents[p]->setTrioDuoType(PARENT_1);
    }
  }

  if (numParents == 1 && parents[0] == NULL) {
    parents[0] = parents[1];
    parents[1] = NULL;
    assert(parents[0] != NULL);
  }

  inferTrioDuoHaplotypes(this, parents, numMendelError, numMendelCounted);

  if (numMendelError != NULL || numMendelCounted != NULL) {
    assert(numMendelError != NULL || numMendelCounted != NULL);

    // Don't want errors for the purpose of counting non-Mendelian errors when
    // other relationships for the same person (e.g., parents of several
    // children, etc.) come along
    this->resetTrioDuoType();
    for(int p = 0; p < 2; p++) {
      if (parents[p] != NULL)
	parents[p]->resetTrioDuoType();
    }
  }
}

// Infers trio and duo haplotypes for <child> and its parent(s):
// Note that no duos reach here at present since we warn about and ignore
// them in the code that parses the .fam file
// If non-NULL, <numMendelError> is set to be the count of the number of
// Mendelian errors, and <numMendelCounted> is the number of relationships with
// non-missing data that were examined for Mendelian errors. These values are
// assumed to point to arrays with size of the number of markers; counts are
// cumulative across calls to inferTrioDuoHaplotypes(), so this method does
// not reset these values, but only increments them.
void PersonBits::inferTrioDuoHaplotypes(PersonBits *child,
					PersonBits *parents[2],
					int *numMendelError,
					int *numMendelCounted) {
  int numMarkers = Marker::getNumMarkers();

  int chromIdx = Marker::getMarker(/*marker=*/ 0)->getChromIdx();

  bool isTrio = parents[1] != NULL;
  if (isTrio) {
    parents[0]->allocTrioDuoData(isTrio, Marker::getNumHapChunks(),
				 parents[1]);
    parents[1]->_tdData = parents[0]->_tdData;
    // Note: cheap hack -- having this point to the parent, a PersonBits object,
    // rather than a TrioDuoData object.  Doing this so that we don't need an
    // extra field in all individuals
    child->_tdData = (TrioDuoData *) parents[0];
  }
  else { // a duo
    parents[0]->allocTrioDuoData(isTrio, Marker::getNumHapChunks(),
				 child);
    child->_tdData = parents[0]->_tdData;
    assert(child->getTrioDuoType() == DUO_CHILD);
  }

  // because we haven't changed known loci at all due to trios and duos,
  // the known loci are the homozygous loci.
  chunk childHomozyLoci = child->getHomozyLoci(/*curHapChunk=*/ 0);
  chunk childHomozyGenos = child->getKnownHaplotype(/*curHapChunk=*/ 0);
  chunk childMissingLoci = child->getMissingLoci(/*curHapChunk=*/ 0);
  chunk parentsHomozyLoci[2], parentsHomozyGenos[2], parentsMissingLoci[2];
  for(int p = 0; p < 2; p++) {
    if (parents[p] == NULL) {
      // all missing data for any NULL parent:
      parentsHomozyLoci[p] = 0;
      parentsHomozyGenos[p] = 0;
      parentsMissingLoci[p] = ALL_CHUNK_BITS_SET;
      continue;
    }
    parentsHomozyLoci[p] = parents[p]->getHomozyLoci(/*curHapChunk=*/ 0);
    parentsHomozyGenos[p] = parents[p]->getKnownHaplotype(/*curHapChunk=*/ 0);
    parentsMissingLoci[p] = parents[p]->getMissingLoci(/*curHapChunk=*/ 0);
  }

  if (isTrio)
    parents[0]->_tdData->_childIsHet[0] =
      ~(child->getHomozyLoci(/*curHapChunk=*/ 0) |
				    child->getMissingLoci(/*curHapChunk=*/ 0));

  // curHapChunk: which haplotype chunk are we on? (BITS_PER_CHUNK bit chunks)
  // curChunkIdx: which bit/locus within the chunk are we on?
  uint32_t curHapChunk = 0, curChunkIdx = 0;
  for(int m = 0; m < numMarkers; m++, curChunkIdx++) {
    if (Marker::getLastMarkerNum(chromIdx) == m - 1) {
      chromIdx = Marker::getMarker(m)->getChromIdx();

      if (curChunkIdx < BITS_PER_CHUNK) {
	// clear out the final bits in _childIsHet that aren't defined because
	// the chromosome ends:
	// The following will clear bits (curChunkIdx..BITS_PER_CHUNK-1); since
	// curChunkIdx is one more than the last bit of the previous chrom, this
	// is exactly what we want.
	int numClearBits = BITS_PER_CHUNK - curChunkIdx;
	chunk clearBits = setLastNumBits(numClearBits);
	if (isTrio)
	  parents[0]->_tdData->_childIsHet[curHapChunk] &= ~clearBits;
      }

      // Now on next chromosome; update chunk indices
      if (curChunkIdx != 0) { // markers from prev chrom on current chunk?
	curHapChunk++; // markers for current chrom are on next chunk number
	curChunkIdx = 0;
      }
    }
    else if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;

      childHomozyLoci = child->getHomozyLoci(curHapChunk);
      childHomozyGenos = child->getKnownHaplotype(curHapChunk);
      childMissingLoci = child->getMissingLoci(curHapChunk);
      for(int p = 0; p < 2; p++) {
	if (parents[p] == NULL)
	  continue;
	parentsHomozyLoci[p] = parents[p]->getHomozyLoci(curHapChunk);
	parentsHomozyGenos[p] = parents[p]->getKnownHaplotype(curHapChunk);
	parentsMissingLoci[p] = parents[p]->getMissingLoci(curHapChunk);
      }

      if (isTrio)
	parents[0]->_tdData->_childIsHet[curHapChunk] =
	  ~(child->getHomozyLoci(curHapChunk) |
					    child->getMissingLoci(curHapChunk));
    }

    int childIsHomozy = (childHomozyLoci >> curChunkIdx) & 1;
    int childHomozyAllele = (childHomozyGenos >> curChunkIdx) & 1;
    int childMissing = (childMissingLoci >> curChunkIdx) & 1;
    int parentsHomozy[2], parentsHomozyAllele[2], parentsMissing[2];
    for(int p = 0; p < 2; p++) {
      parentsHomozy[p] = (parentsHomozyLoci[p] >> curChunkIdx) & 1;
      parentsHomozyAllele[p] = (parentsHomozyGenos[p] >> curChunkIdx) & 1;
      parentsMissing[p] = (parentsMissingLoci[p] >> curChunkIdx) & 1;
    }


    //////////////////////////////////////////////////////////////////////////
    // check for Mendelian errors

    if (!childMissing && (!parentsMissing[0] || !parentsMissing[1]))
      // for rate calculation: should only count observations where child and
      // at least one parent is non-missing
      numMendelCounted[m]++;

    // Can only find Mendelian errors if child is non-missing and at least one
    // parent is homozygous:
    if (!childMissing && (parentsHomozy[0] || parentsHomozy[1])) {
      if (parentsHomozy[0] && parentsHomozy[1]) {
	int childGeno = parentsHomozyAllele[0] + parentsHomozyAllele[1];
	if (childIsHomozy) {
	  if (childGeno == 1) {
	    // Mendelian error!
	    setTrioDuoMissing(child, parents, curHapChunk, curChunkIdx);
	    numMendelError[m]++;
	    continue;
	  }
	  else if (childGeno / 2 != childHomozyAllele) {
	    // Mendelian error!
	    setTrioDuoMissing(child, parents, curHapChunk, curChunkIdx);
	    numMendelError[m]++;
	    continue;
	  }
	}
	else if (childGeno != 1) {
	  // Mendelian error!
	  setTrioDuoMissing(child, parents, curHapChunk, curChunkIdx);
	  numMendelError[m]++;
	  continue;
	}
      }
      else {
	// only one homozygous parent; can only get Mendelian error if child is
	// homozygous

	if (childIsHomozy) {

	  int parHomozyAllele;
	  if (parentsHomozy[0]) {
	    parHomozyAllele = parentsHomozyAllele[0];
	  }
	  else {
	    assert(parentsHomozy[1]);
	    parHomozyAllele = parentsHomozyAllele[1];
	  }

	  if (parHomozyAllele != childHomozyAllele) {
	    // Mendelian error!
	    setTrioDuoMissing(child, parents, curHapChunk, curChunkIdx);
	    numMendelError[m]++;
	    continue;
	  }
	}
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // No Mendelian error: infer phase

    if (parentsHomozy[0] && parentsHomozy[1]) { // both parents homozygous
      continue; // phase of both parents is known
    }

    // one or both parents' phase unknown

    if (childIsHomozy) {
      // child homozygous: don't need parents to infer:
      for(int p = 0; p < 2; p++) {
	if (parents[p] == NULL)
	  continue;
	// assign phase
	parents[p]->_knownHap[curHapChunk] |=
				    ((chunk) childHomozyAllele) << curChunkIdx;
      }
      parents[0]->_tdData->_knownLoci[curHapChunk] |= 1ul << curChunkIdx;
      continue;
    }

    // child is heterozygous or missing; try to infer phase:

    if ((!parentsHomozy[0] || parentsMissing[0]) &&
	(!parentsHomozy[1] || parentsMissing[1])) {
      // both parents either heterozygous or missing => ambiguous: skip
      continue;
    }

    if (!isTrio) {
      // Duo!
      // At this point, we know the child is heterozygous or missing, but at
      // least one parent is homozygous.  Since we only have data for one
      // parent, that parent must be homozygous, and that resolves the
      // transmitted haplotype for the child
      assert(parentsHomozy[0]);
      child->_knownHap[curHapChunk] |=
				((chunk) parentsHomozyAllele[0]) << curChunkIdx;
      child->_tdData->_knownLoci[curHapChunk] |= 1ul << curChunkIdx;
      continue;
    }

    // Only trios from here down!

    if (childMissing) {
      // Note: could impute the child's genotype if both parents are
      // homozygous, but we aren't phasing the (trio) child anyway, so can just
      // move on.  The parents don't need to be phased at homozygous sites.
      continue;
    }

    // child is heterozygous with one parent homozygous

    int homozyPar = 0;
    for(; homozyPar < 2; homozyPar++) {
      if (parentsHomozy[homozyPar])
	break;
    }
    assert(homozyPar < 2);

    int otherParent = homozyPar ^ 1;
    assert(!parentsHomozy[otherParent]); // ensure other parent not homozygous

    // child is heterozygous, so its genotype is 1 and it received whatever
    // allele from the heterozygous parent that the homozygous parent didn't
    // transmit:
    int otherParTransAllele = 1 - parentsHomozyAllele[homozyPar];
    parents[otherParent]->_knownHap[curHapChunk] |=
				  ((chunk) otherParTransAllele) << curChunkIdx;
    parents[otherParent]->_tdData->_knownLoci[curHapChunk] |= 1ul <<curChunkIdx;
  }

  if (curChunkIdx < BITS_PER_CHUNK) { // identical to code above
    // clear out the final bits in _childIsHet that aren't defined because
    // the chromosome ends:
    // The following will clear bits (curChunkIdx..BITS_PER_CHUNK-1); since
    // curChunkIdx is one more than the last bit of the previous chrom, this
    // is exactly what we want.
    int numClearBits = BITS_PER_CHUNK - curChunkIdx;
    chunk clearBits = setLastNumBits(numClearBits);
    if (isTrio)
      parents[0]->_tdData->_childIsHet[curHapChunk] &= ~clearBits;
  }
}

// For sites in trios and duos that are Mendelian errors, sets the locus to
// missing in all individuals
void PersonBits::setTrioDuoMissing(PersonBits *child, PersonBits *parents[2],
				   int chunkNum, int chunkIdx) {
  child->setMissing(chunkNum, chunkIdx);
  parents[0]->setMissing(chunkNum, chunkIdx);
  if (parents[1] != NULL)
    parents[1]->setMissing(chunkNum, chunkIdx);
}

// To be run prior to removeIgnoreIndivs() in PersonIO::readData():
// After parsing and performing trio-based phasing in
// setParents(), trio children are no longer needed and
// their memory can be reclaimed. Method either sets these children's
// population to Ignore (so the Person object will be removed), or clears
// their genotype data so they take up minimal space. The latter occurs
// if we're going to later print their inferred haplotypes
void PersonBits::cleanUpPostParse(bool keepTrioKids) {
  if (!keepTrioKids) {
    // no kids after this and removeIgnoreIndivs() run
    _numTrioKids = 0;
  }

  int length = _allIndivs.length();
  for(int p = 0; p < length; p++) {
    PersonBits *cur = _allIndivs[p];
    if (cur->_ignore)
      continue;
    if (cur->getTrioDuoType() == TRIO_CHILD) {
      if (keepTrioKids)
	// need not retain genotype data for trio child, even though we are
	// keeping the child's entry for later printing
	cur->empty();
      else
	cur->setIgnore();
    }
  }
}


// Prints parent's haplotypes for trios/duos.  For each parent, the first
// haplotype printed is the transmitted haplotype and the second is
// untransmitted.  At present, this method only gets called before population
// phasing takes place, so the triple het sites are untrustworthy and are
// set to missing in inferTrioDuoHaplotypes().
void PersonBits::printTrioEigenstratPhased(FILE *out) {
  int numMarkers = Marker::getNumMarkers();
  int numIndivs = _allIndivs.length();

  int chromIdx = Marker::getMarker(/*marker=*/ 0)->getChromIdx();

  // curHapChunk: which haplotype chunk are we on? (BITS_PER_CHUNK bit chunks)
  // curChunkIdx: which bit/locus within the chunk are we on?
  for(int m = 0, curHapChunk = 0, curChunkIdx = 0; m < numMarkers;
							  m++, curChunkIdx++) {
    if (Marker::getLastMarkerNum(chromIdx) == m - 1) {
      // Now on next chromosome; update chunk indices
      if (curChunkIdx != 0) { // markers from prev chrom on current chunk?
	curHapChunk++; // markers for current chrom are on next chunk number
	curChunkIdx = 0;
      }
      chromIdx = Marker::getMarker(m)->getChromIdx();

    }
    if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;
    }

    for(int i = 0; i < numIndivs; i++) {
      PersonBits *thePerson = _allIndivs[i];

      // only print phase for the parents of trios/duos:
      if (thePerson->getTrioDuoType() == UNRELATED)
	continue;

      chunk haplotype;
      chunk knownLoci;

      // minimally the person's haplotype will contain the alleles at
      // homozygous sites (this matters at sites where the child is
      // missing data -- there _trioDuoHaplotype->defined won't be set to 1):
      haplotype = thePerson->_knownHap[curHapChunk];
      knownLoci = thePerson->_homozyLoci[curHapChunk] |
		  (thePerson->_tdData->_knownLoci[curHapChunk] &
				    (~thePerson->getMissingLoci(curHapChunk)));

      if (((~knownLoci) >> curChunkIdx) & 1) {
	// site is missing/ambiguous:
	fprintf(out, "99");
	continue;
      }

      bool subtract = false;
      for(int h = 0; h < 2; h++) {
	if (subtract) {
	  haplotype ^= (~thePerson->_homozyLoci[curHapChunk]);
	}

	int hapAllele = (haplotype >> curChunkIdx) & 1;

	fprintf(out, "%d", hapAllele);

	subtract = true;
      }
    }
    fprintf(out, "\n");

  }
}

// Allocates space for sampled haplotypes and initializes them using random
// values for heterozygous sites
void PersonBits::initRandSampledHaps() {
  int numSamples = _allIndivs.length();
  int numHapChunks = Marker::getNumHapChunks();
  for(int id = 0; id < numSamples; id++) {
    PersonBits *thePerson = _allIndivs[id];

    // To ensure that the transmitted haplotypes match between the parent
    // and duo child (necessary to avoid a case in which the randomly
    // initialized haplotypes are infeasible according to the relationship
    // and genotypes), and the case where the transmitted haplotypes of two
    // trio parents do not match the trio child, we initialize parent-offspring
    // duos and the two trio parents together.
    PersonBits *otherPerson = NULL;
    TDTYPE otherType = UNRELATED;

    if (thePerson->getTrioDuoType() == PARENT_0) {
      otherPerson = thePerson->_tdData->_otherParentOrChild;
      otherType = otherPerson->getTrioDuoType();
    }
    else if (thePerson->getTrioDuoType() != UNRELATED)
      // skip DUO_CHILD and PARENT_1 indivs, they are initialized with PARENT_0
      continue;


    thePerson->_sampledHaplotypes = new chunk*[numHapChunks];
    if (otherPerson != NULL)
      otherPerson->_sampledHaplotypes = new chunk*[numHapChunks];

    for(int curHapChunk = 0; curHapChunk < numHapChunks; curHapChunk++) {
      thePerson->_sampledHaplotypes[curHapChunk] =
					  new chunk[NUM_HAPLOTYPES_TO_SAMPLE*2];
      if (otherPerson != NULL)
	otherPerson->_sampledHaplotypes[curHapChunk] =
					  new chunk[NUM_HAPLOTYPES_TO_SAMPLE*2];

      for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
	thePerson->_sampledHaplotypes[curHapChunk][s*2+0] =
					      thePerson->_knownHap[curHapChunk];
	chunk knownSites = thePerson->getKnownLoci(curHapChunk);
	chunk missingSites = thePerson->getMissingLoci(curHapChunk);
	for(uint32_t curChunkIdx = 0; curChunkIdx < BITS_PER_CHUNK;
								curChunkIdx++) {
	  int isKnown = (knownSites >> curChunkIdx) & 1;
	  int isMiss = (missingSites >> curChunkIdx) & 1;

	  int randGeno = rand() / (RAND_MAX / 2 + 1); // +1 so result < 2
	  // only randomize if the haplotype value is not known and the
	  // genotype is not missing:
	  //
	  // Note: we won't bother randomizing missing data sites since the
	  // first window size has only 4 markers, so we will seed all possible
	  // values for missing data sites (since 4 < MAX_NUM_MISSING_FULL_SEED)
	  chunk genoVal = (1 - isKnown) & (1 - isMiss) & randGeno;

	  thePerson->_sampledHaplotypes[curHapChunk][s*2+0]
						    |= (genoVal << curChunkIdx);
	}

	// invert for haplotype 1:
	chunk homozySites = thePerson->getHomozyLoci(curHapChunk);
	chunk hetSites = ~(homozySites | missingSites);
	thePerson->_sampledHaplotypes[curHapChunk][s*2+1] =
		  thePerson->_sampledHaplotypes[curHapChunk][s*2+0] ^ hetSites;

	if (otherPerson != NULL) {
	  // Set haplotype 0 (different for duo children and other trio parents:
	  if (otherType == DUO_CHILD) {
	    // haplotype 0 identical to parent:
	    otherPerson->_sampledHaplotypes[curHapChunk][s*2+0] =
			      thePerson->_sampledHaplotypes[curHapChunk][s*2+0];
	  }
	  else {
	    assert(otherType == PARENT_1);

	    chunk thePersonTrans =
			      thePerson->_sampledHaplotypes[curHapChunk][s*2+0];
	    chunk trioChildHet = thePerson->getTrioChildHet(curHapChunk);
	    otherPerson->_sampledHaplotypes[curHapChunk][s*2+0] =
		      otherPerson->_knownHap[curHapChunk] |
					    ((~thePersonTrans) & trioChildHet);
	    chunk knownSites = otherPerson->getKnownLoci(curHapChunk) |
					  trioChildHet;
	    chunk missingSites = thePerson->getMissingLoci(curHapChunk);
	    for(uint32_t curChunkIdx = 0; curChunkIdx < BITS_PER_CHUNK;
								curChunkIdx++) {
	      int isKnown = (knownSites >> curChunkIdx) & 1;
	      int isMiss = (missingSites >> curChunkIdx) & 1;

	      int randGeno = rand() / (RAND_MAX / 2 + 1); // +1 so result < 2
	      chunk genoVal = (1 - isKnown) & (1 - isMiss) & randGeno;

	      otherPerson->_sampledHaplotypes[curHapChunk][s*2+0]
						    |= (genoVal << curChunkIdx);
	    }
	  }

	  // haplotype 1 inverted at heterozygous sites:
	  chunk otherHomozySites = otherPerson->getHomozyLoci(curHapChunk);
	  chunk otherMissingSites = otherPerson->getMissingLoci(curHapChunk);
	  chunk otherHetSites = ~(otherHomozySites | otherMissingSites);
	  otherPerson->_sampledHaplotypes[curHapChunk][s*2+1] =
	    otherPerson->_sampledHaplotypes[curHapChunk][s*2+0] ^ otherHetSites;
	}
      }
    }

  }
}

void PersonBits::setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
			     int chrMarkerIdx, int geno[2]) {
  if (_ignore)
    return; // no need/nowhere to store genotypes
  if (geno[0] < 0) {
    // missing data
    _missingLoci[hapChunkNum] += 1ul << chunkIdx;
    return;
  }
  if ((geno[0] == 0 && geno[1] == 1) || (geno[0] == 1 && geno[1] == 0))
    return; // heterozygote => do nothing

  assert(geno[0] <= 1 && geno[1] <= 1);

  // have homozygote: need to set the corresponding bit:
  chunk homozyVal = 1ul << chunkIdx;
  // want lowest order bit of genotype (0 or 1) to be set at chunkIdx in the
  // chunk:
  chunk genoVal = ((unsigned long) geno[0]) << chunkIdx;

  _homozyLoci[hapChunkNum] += homozyVal;
  _knownHap[hapChunkNum] += genoVal;
}

// Sets <this>'s genotype to missing at <hapChunkNum>, <chunkIdx>
void PersonBits::setMissing(int hapChunkNum, int chunkIdx) {
  chunk bitToSetMiss = 1ul << chunkIdx;
  // make sure the locus is not set to be homozygous and that the homozygous
  // genotype value is 0:
  _homozyLoci[hapChunkNum] &= ~bitToSetMiss;
  _knownHap[hapChunkNum] &= ~bitToSetMiss;
  // set missing:
  _missingLoci[hapChunkNum] |= bitToSetMiss;
}

// Sets all heterozygous sites to missing.  Used for males on the X chromosome
// where heterozygous sites are errors.
// If non-NULL, sets *numHets and *numCalls to be the number of heterozygous
// sites and the number of called (non-missing) sites for <this>.  Those
// counts are *before* setting heterozygous sites to missing.
// Code assumes that the entire dataset is the X chromosome
void PersonBits::setXHetToMissing(int *numHets, int *numCalls) {
//  assert(Marker::getMarker(0)->getChrom() == CHR_X);
//  assert(Marker::getNumChromMarkers(CHR_X) == Marker::getNumMarkers());
  assert(getSex() == 'M');

  if (numHets != NULL) {
    assert(numCalls != NULL);
    *numHets = 0;
    *numCalls = 0;
  }

  int numMissing = 0;

  for(int curChunk = 0; curChunk < Marker::getNumHapChunks(); curChunk++) {
    chunk homozyBits = getHomozyLoci(curChunk);
    chunk missingBits = getMissingLoci(curChunk);

    numMissing += countBitsSet(missingBits);
    chunk hets = (~homozyBits) & (~missingBits);

    if (curChunk+1 == Marker::getNumHapChunks()) {
      // Don't count as heterozygous those sites that appear beyond 
      int numValid = Marker::getNumMarkers() % BITS_PER_CHUNK;
      if (numValid > 0)
	hets &= setBitsToIdx( numValid - 1);
    }
    if (numHets != NULL)
      *numHets += countBitsSet(hets);

    assert((_missingLoci[curChunk] & hets) == 0ul);
    _missingLoci[curChunk] |= hets;
  }

  if (numCalls != NULL) {
    *numCalls = Marker::getNumMarkers() - numMissing;
  }
}

// prints the portion of the haplotype for <this> that is defined by its
// homozygous loci, with ? for heterozygous loci
void PersonBits::printChunkHap(FILE *out, int chunkNum) {
  printHap(out, _knownHap[chunkNum], _homozyLoci[chunkNum]);
}

// Sets all sampled haplotypes to 0 so that we can call orSampledHaplotypes()
// to set newly sampled haplotypes
void PersonBits::clearSampledHaplotypes() {
  int numHapChunks = Marker::getNumHapChunks();
  for(int c = 0; c < numHapChunks; c++) {
    for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
      for(int h = 0; h < 2; h++) {
	_sampledHaplotypes[c][s*2+h] = 0;
      }
    }
  }
}


void PersonBits::orSampledHaplotype(int sampNum, int homolog, int chunkNum,
			chunk haplotype) {
  _sampledHaplotypes[chunkNum][sampNum*2+homolog] |= haplotype;
}

// Allocates space in which to store the final (Viterbi decoded) haplotypes for
// this individual
void PersonBits::initFinalHaplotype() {
  // don't need the genotype data anymore, so steal their memory:
  _resolvedHaplotype[0] = _knownHap;
  _resolvedHaplotype[1] = _homozyLoci;
  int numHapChunks = Marker::getNumHapChunks();
  for(int h = 0; h < 2; h++) {
    for(int c = 0; c < numHapChunks; c++) {
      _resolvedHaplotype[h][c] = 0;
    }
  }
  _knownHap = _homozyLoci = NULL;
}

void PersonBits::orFinalHaplotype(int homolog, int chunkNum, chunk haplotype) {
  _resolvedHaplotype[homolog][chunkNum] |= haplotype;
}
