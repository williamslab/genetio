// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <time.h>
#include "nuclearfamily.h"
#include "util.h"

NuclearFamily::fam_ht NuclearFamily::_families;

// Given two <parents>, looks up or creates a corresponding NuclearFamily and
// adds <child> to their list of children
void NuclearFamily::addChild(PersonBulk *parents[2], PersonBulk *child) {
  par_pair_real parentsKey(parents[0], parents[1]);

  fam_ht_iter it = _families.find( &parentsKey );
  if (it == _families.end()) {
    // No family yet for these parents; create and insert it
    NuclearFamily *newFamily = new NuclearFamily();
    newFamily->_parents = new par_pair_real(parents[0], parents[1]);
    newFamily->_children.append(child); // add first child
    par_pair parentsKey = newFamily->_parents;
    _families[ parentsKey ] = newFamily;
  }
  else {
    // Have a family: append child to list of children
    // <*it> is a pair with <it->first> being the key and <it->second> the value
    // so <it->second> is the NuclearFamily
    it->second->_children.append(child);
  }
}

// Prints the haplotypes for <this> to <out> in human readable text format
void NuclearFamily::printHapTxt(FILE *out, int chrIdx) {
  int numChildren = _children.length();

  /////////////////////////////////////////////////////////////////////////
  // Print header
  int famIdLen = _parents->first->getFamilyIdLength();
  if (famIdLen == 0) {
    fprintf(out, "%s %s", _parents->first->getId(), _parents->second->getId());
    for(int c = 0; c < numChildren; c++)
      fprintf(out, " %s", _children[c]->getId());
    fprintf(out, "\n");
  }
  else {
    fprintf(out, "(%.*s) %s %s", famIdLen, _parents->first->getId(),
	    &_parents->first->getId()[famIdLen+1],
	    &_parents->second->getId()[famIdLen+1]);
    for(int c = 0; c < numChildren; c++)
      fprintf(out, " %s", &_children[c]->getId()[famIdLen+1]);
    fprintf(out, "\n");
  }

  /////////////////////////////////////////////////////////////////////////
  // Print haplotypes
  // TODO: better name?
  const char ambigMissType[4] = { '/', '_', '?', '!' };
  const char markerType[3] = { '0', '1', 'B' };
  const char ivLabel[4] = { 'A', 'B', 'a', 'b' };

  const char *chrName = Marker::getChromName(chrIdx);

  int prevM_OK = -1;
  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);
  for(int m = firstMarker; m <= lastMarker; m++) {
    const char *theAlleles = Marker::getMarker(m)->getAlleleStr();
    // For use with the printGeno() method: allows printing an unknown genotype
    // for untransmitted alleles at positions where we attempt to impute the
    // parents' genotype: alleles[1] is '0' or unknown.
    char alleles[3] = { theAlleles[0], '0', theAlleles[2] };

    fprintf(out, "%-2s %-5d ", chrName, m - firstMarker);

    PhaseStatus status = _phase[m].status;
    uint8_t parentData, hetParent, untrans;
    uint64_t childrenData, iv, ambigMiss, ivFlippable;
    uint8_t imputeParGeno = G_MISS; // by default no imputation
    uint8_t swapHet = 0; // by default no swapping hets
    uint8_t parAlleleIdx[2][2];
    uint8_t homGenotypes[2] = { G_HOM0, G_HOM1 };
    bool bothParHomozy = status == PHASE_AMBIG && _phase[m].ambigBothParHomozy;
    switch(status) {
      case PHASE_UNINFORM:
	// can impute at uninformative markers, not the others, using the
	// <homParentGeno> field:
	imputeParGeno = _phase[m].homParentGeno;
	// only will need/know to swap heterozygous markers for uninformative
	// markers
	swapHet = _phase[m].uninfHetSwap;
      case PHASE_AMBIG:
      case PHASE_ERROR:
      case PHASE_ERR_RECOMB:
	/////////////////////////////////////////////////////////////////////
	// Not phased / trivially phased cases:

	// print the parent's genotypes
	parentData = _phase[m].parentData;
	untrans = _phase[m].untransParHap;
	for(int p = 0; p < 2; p++) {
	  uint8_t thisParGeno = parentData & 3;
	  uint8_t thisUntrans = untrans & 3;
	  // <isMissing> is 1 iff thisParGeno == 1:
	  uint8_t isMissing = (thisParGeno & 1) & ~(thisParGeno >> 1);
	  // Will print either the parent genotype if it is not missing or
	  // the imputed genotype if it is missing (note that the default is
	  // for <imputeParGeno> is missing as well):
	  imputeParGeno = (1 - bothParHomozy) * imputeParGeno +
			      bothParHomozy * homGenotypes[p];
	  uint8_t genoToPrint = (1 - isMissing) * thisParGeno +
					      isMissing * imputeParGeno;
	  printGeno(out, alleles, genoToPrint, ambigMissType[ isMissing ],
		    isMissing * thisUntrans);
	  parentData >>= 2;
	  untrans >>= 2;
	  if (p == 0)
	    fprintf(out, "   ");
	}
	switch(status) {
	  case PHASE_UNINFORM:
	    fprintf(out, " |");
	    break;
	  case PHASE_AMBIG:
	    fprintf(out, " ?");
	    break;
	  case PHASE_ERROR:
	    fprintf(out, " E");
	    break;
	  case PHASE_ERR_RECOMB:
	    fprintf(out, " R");
	    break;
	  case PHASE_OK: // can't get here but put case in to prevent warning
	  case NUM_PHASE_STATUS:
	    break;
	}

	// print the children's genotypes
	childrenData = _phase[m].iv;
	for(int c = 0; c < numChildren; c++) {
	  fprintf(out, " ");
	  printGeno(out, alleles, childrenData & 3, /*sep=*/ '/',
		    /*untrans=NA=*/ 0, swapHet);
	  fprintf(out, "   ");
	  childrenData >>= 2;
	}
	fprintf(out, "\n");
	break;

      case PHASE_OK:
	/////////////////////////////////////////////////////////////////////
	// Standard phased marker

	// determine which alleles each parent has on each haplotype;
	getParAlleles(m, parAlleleIdx);

	hetParent = _phase[m].hetParent;

	// print parent's haplotypes
	uint8_t parIsMiss[2];
	parIsMiss[0] = _phase[m].parentData & 3;
	parIsMiss[1] = (_phase[m].parentData >> 2) & 3;
	fprintf(out, "%c%c%c %c %c%c%c ",
		alleles[ parAlleleIdx[0][0] ], ambigMissType[ parIsMiss[0] ],
		alleles[ parAlleleIdx[0][1] ],
		markerType[hetParent],
		alleles[ parAlleleIdx[1][0] ], ambigMissType[ parIsMiss[1] ],
		alleles[ parAlleleIdx[1][1] ]);

	// print phase status
	if (_phase[m].ambigParPhase | _phase[m].ambigParHet |
						      _phase[m].arbitraryPar) {
	  if (_phase[m].arbitraryPar)
	    fprintf(out, "P");

	  // Will set this value to be only the bits relative to
	  // <_phase[m].hetParent> if this first assignment is non-zero:
	  uint8_t thisParHetAmbigPhase = _phase[m].ambigParPhase;
	  if (thisParHetAmbigPhase) {
	    // update to the bits relevant to this hetParent
	    thisParHetAmbigPhase >>= 2 * hetParent;
	    uint8_t isPI = hetParent >> 1;
	    thisParHetAmbigPhase &= isPI * 15 + (1-isPI) * 3;
	  }
	  if (thisParHetAmbigPhase) {
	    fprintf(out, "S");
	    if (hetParent == 2) {
	      char swapTypes[4] = { '?', '0', '1', 'B' };
	      uint8_t diff;
	      switch(thisParHetAmbigPhase) {
		// only a single bit in the ambiguous value cases:
		case 1:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 2:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 4:
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 8:
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		// two bits in the ambiguous value cases:
		case 3:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 5:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 6:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 9:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 10:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 12:
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		// three bits in the ambiguous value cases:
		case 7:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 11:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 13:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 14:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		default:
		  // 0:  impossible to get to this branch
		  // 15: includes the current value which shouldn't happen
		  fprintf(out, "ERROR: ambigParPhase status %d\n",
			  thisParHetAmbigPhase);
		  break;
	      }
	    }
	  }
	  if (_phase[m].ambigParHet) {
	    switch (_phase[m].ambigParHet) {
	      case 1:
		fprintf(out, "H0");
		break;
	      case 2:
		fprintf(out, "H1");
		break;
	      case 3:
		fprintf(out, "H01");
		break;
	      case 4:
		fprintf(out, "H2");
		break;
	      case 5:
		fprintf(out, "H02");
		break;
	      case 6:
		fprintf(out, "H12");
		break;

	      default:
		fprintf(out, "ERROR: ambigParHet status %d\n",
			_phase[m].ambigParHet);
		break;
	    }
	  }
	}
	else
	  fprintf(out, "|");

	// now print children's haplotypes
	iv = _phase[m].iv;
	ivFlippable = _phase[m].ivFlippable;
	ambigMiss = _phase[m].ambigMiss;
	for(int c = 0; c < numChildren; c++) {
	  uint8_t curIV = iv & 3;
	  uint8_t curIVflip = ivFlippable & 3;
	  uint8_t curAmbigMiss = ambigMiss & 3;
	  int ivs[2] = { curIV & 1, curIV >> 1 };
	  // IV letter is '?' if either parent is flippable
	  // mechanically this happens by adding 2 to the index to <ivLabel>.
	  // This occurs when different states yield minimum recombinant
	  // haplotypes and differ in their inheritance vector values
	  // TODO: better name for <ivLetters> and <ivLabel>
	  // TODO: add to the IV output
	  char curIVchars[2] = { ivLabel[ ivs[0] + ((curIVflip & 1) << 1) ],
				 ivLabel[ ivs[1] + (curIVflip & 2) ] };
	  fprintf(out, " %c%c%c %c%c", alleles[ parAlleleIdx[0][ ivs[0] ] ],
		  ambigMissType[curAmbigMiss],
		  alleles[ parAlleleIdx[1][ ivs[1] ] ],
		  curIVchars[0], curIVchars[1]);
	  iv >>= 2;
	  ivFlippable >>= 2;
	  ambigMiss >>= 2;
	}
	// TODO: remove at some point
	if (prevM_OK >= 0) {
	  uint8_t count = popcount(_phase[m].iv ^ _phase[prevM_OK].iv);
	  if (_phase[m].numRecombs != count) {
	    fprintf(out, " %2d!!!! ", count);
	  }
	}
	fprintf(out, " %2d\n", _phase[m].numRecombs);

	// TODO: remove
	prevM_OK = m;
	break;

      default:
	fprintf(stderr, "ERROR: marker status %d\n", _phase[m].status);
	exit(1);
	break;
    }
  }
}

// Prints the parents' and optionally children's haplotypes for <this> nuclear
// family to <out> in JSON format. Assumes that a containing object (hash) is
// the current context in the JSON file.
void NuclearFamily::printHapJson(FILE *out, bool withChildren) {
  fprintf(out, "\"%s-%s\":", _parents->first->getId(),
			     _parents->second->getId());
  // The order of the JSON format output is the transpose of other ways we print
  // to save redundant work, we store strings representing the haplotypes of the
  // parents and any children and then printing is simple.

  // Because HAPI2 is limited to 32 children, we'll allocate an array of fixed
  // (maximum) size. Easy to change if HAPI changes.
  char *hapStrs[64 + 4];

  int numHaps = 2 * 2;
  int numChildren = _children.length();
  int numMarkers = Marker::getNumMarkers();
  if (withChildren) {
    numHaps += 2 * numChildren;
    assert(numChildren <= 32);
  }

  for(int h = 0; h < numHaps; h++)
    hapStrs[h] = new char[numMarkers];

  for(int m = 0; m < numMarkers; m++) {
    const char *theAlleles = Marker::getMarker(m)->getAlleleStr();
    char alleles[3] = { theAlleles[0], '0', theAlleles[2] };

    PhaseStatus status = _phase[m].status;
    uint8_t imputeParGeno = G_MISS; // by default no imputation
    uint8_t swapHet = 0; // by default no swapping hets
    uint8_t parAlleleIdx[2][2];
    uint8_t homGenotypes[2] = { G_HOM0, G_HOM1 };
    bool bothParHomozy = status == PHASE_AMBIG && _phase[m].ambigBothParHomozy;
    switch(status) {
      case PHASE_UNINFORM:
	// can impute at uninformative markers, not the others, using the
	// <homParentGeno> field:
	imputeParGeno = _phase[m].homParentGeno;
	// only will need/know to swap heterozygous markers for uninformative
	// markers
	swapHet = _phase[m].uninfHetSwap;
      case PHASE_AMBIG:
      case PHASE_ERROR:
      case PHASE_ERR_RECOMB:
	///////////////////////////////////////////////////////////////////
	// Not phased / trivially phased cases:

	{
	  for(int par = 0; par < 2; par++) {
	    uint8_t thisParGeno = (_phase[m].parentData >> (par * 2)) & 3;
	    uint8_t thisUntrans = (_phase[m].untransParHap >> (par * 2)) & 3;
	    // <isMissing> is 1 iff thisParGeno == 1:
	    uint8_t isMissing = (thisParGeno & 1) & ~(thisParGeno >> 1);
	    // Will print either the parent genotype if it is not missing or
	    // the imputed genotype if it is missing (note that the default is
	    // for <imputeParGeno> is missing as well):
	    imputeParGeno = (1 - bothParHomozy) * imputeParGeno +
			      bothParHomozy * homGenotypes[par];
	    uint8_t genoToPrint = (1 - isMissing) * thisParGeno +
					      isMissing * imputeParGeno;
	    for(int hap = 0; hap < 2; hap++) {
	      uint8_t curUntrans = (thisUntrans >> hap) & 1;
	      // if untransmitted, then missing:
	      uint8_t curGenoToPrint = (1 - curUntrans) * genoToPrint +
							  curUntrans * G_MISS;
	      switch (curGenoToPrint) {
		case G_HOM0:
		  hapStrs[par * 2 + hap][m] = alleles[0];
		  break;
		case G_MISS:
		  hapStrs[par * 2 + hap][m] = '0';
		  break;
		case G_HET:
		  hapStrs[par * 2 + hap][m] = alleles[ hap * 2 ];
		  break;
		case G_HOM1:
		  hapStrs[par * 2 + hap][m] = alleles[2];
		  break;
	      }
	    } // haplotype
	  } // parent

	  if (!withChildren)
	    // not printing children: done
	    break;

	  // also printing children
	  uint64_t childrenData = _phase[m].iv;
	  for(int c = 0; c < numChildren; c++) {
	    uint8_t childGenotype = childrenData & 3;

	    uint8_t hap0;
	    switch (childGenotype) {
	      case G_HOM0:
		hapStrs[4 + c*2 + 0][m] = hapStrs[4 + c*2 + 1][m] = alleles[0];
		break;
	      case G_MISS:
		hapStrs[4 + c*2 + 0][m] = hapStrs[4 + c*2 + 1][m] = '0';
		break;
	      case G_HET:
		// if <swapHet> == 0, hap0 == 0. Otherwise, hap0 == 1 (swapped)
		hap0 = swapHet;
		hapStrs[4 + c*2 + hap0][m] = alleles[0];
		hapStrs[4 + c*2 + (1 - hap0)][m] = alleles[2];
		break;
	      case G_HOM1:
		hapStrs[4 + c*2 + 0][m] = hapStrs[4 + c*2 + 1][m] = alleles[2];
		break;
	    }

	    childrenData >>= 2;
	  }
	}
	break;

      case PHASE_OK:
	///////////////////////////////////////////////////////////////////
	// Standard phased marker

	getParAlleles(m, parAlleleIdx);

	for(int par = 0; par < 2; par++)
	  for(int har = 0; har < 2; har++)
	    hapStrs[par * 2 + har][m] = alleles[ parAlleleIdx[par][har] ];

	if (!withChildren)
	  // not printing children: done
	  break;

	// children's haplotypes
	for(int c = 0; c < numChildren; c++) {
	  uint8_t curIV = (_phase[m].iv >> (2*c)) & 3;

	  int ivs[2] = { curIV & 1, curIV >> 1 };
	  hapStrs[4 + c*2 + 0][m] = alleles[ parAlleleIdx[0][ ivs[0] ] ];
	  hapStrs[4 + c*2 + 1][m] = alleles[ parAlleleIdx[1][ ivs[1] ] ];

	}
	break;

      default:
	fprintf(stderr, "ERROR: marker status %d\n", _phase[m].status);
	exit(1);
	break;
    }
  }

  // print the parent haplotypes across all chromosomes
  //    will have an array of length 2 indexing the haplotypes for each parent
  fprintf(out, "{\"parhaps\":[");
  for(int par = 0; par < 2; par++) {
    if (par == 1)
      fprintf(out, ","); // separate the parent values in "parhaps" by commas

    // Place the two haplotypes each parent possesses into another array 
    fprintf(out, "[");

    for(int hap = 0; hap < 2; hap++) {
      if (hap == 1)
	fprintf(out, ","); // separate the two haplotypes by commas

      // Haplotype array containing <numMarkers> elements
      fprintf(out, "[");
      for(int m = 0; m < numMarkers; m++) {
	if (m > 0)
	  fprintf(out, ","); // separate haplotype alleles by commas

	// Allele at marker <m>:
	fprintf(out, "\"%c\"", hapStrs[par * 2 + hap][m]);
      }
      fprintf(out, "]"); // end haplotype array
    }
    fprintf(out, "]"); // the pair of haplotypes for <par>
  }
  fprintf(out, "]"); // "parhaps" value


  if (withChildren) {
    // print the children's haplotypes across all chromosomes
    fprintf(out, ",\"childhaps\":{");
    for(int c = 0; c < numChildren; c++) {
      if (c > 0)
	fprintf(out, ","); // separate the child values in "childhaps" by commas

      // top level is a dictionary indexed by child id; value of each entry is
      // an array with two values for the two haplotypes
      fprintf(out, "\"%s\":[", _children[c]->getId());

      for(int hap = 0; hap < 2; hap++) {
	if (hap == 1)
	  fprintf(out, ","); // separate the two haplotypes by commas

	// Haplotype array containing <numMarkers> elements
	fprintf(out, "[");
	for(int m = 0; m < numMarkers; m++) {
	  if (m > 0)
	    fprintf(out, ","); // separate haplotype alleles by commas

	  // Allele at marker <m>:
	  fprintf(out, "\"%c\"", hapStrs[4 + c * 2 + hap][m]);
	}
	fprintf(out, "]"); // end haplotype array
      }

      fprintf(out, "]"); // the pair of haplotypes for child <c>
    }
    fprintf(out, "}"); // end "childhaps" value
  }

  // Code array containing <numMarkers> elements
  // -> Gives codes that indicate phase state and any ambiguities
  fprintf(out, ",\"codes\":[");
  for(int m = 0; m < numMarkers; m++) {
    if (m > 0)
      fprintf(out, ","); // separate marker codes by commas

    PhaseStatus status = _phase[m].status;
    switch(status) {
      case PHASE_UNINFORM:
	fprintf(out, "null");
	break;
      case PHASE_AMBIG:
	fprintf(out, "[\"?\"]");
	break;
      case PHASE_ERROR:
	fprintf(out, "[\"E\"]");
	break;
      case PHASE_ERR_RECOMB:
	fprintf(out, "[\"R\"]");
	break;
      case PHASE_OK:
	if (!(_phase[m].arbitraryPar || _phase[m].ambigParPhase ||
	      _phase[m].ambigParHet)) {
	  // well phased:
	  fprintf(out, "null");
	  break;
	}
	{
	  fprintf(out, "[");
	  bool somethingPrinted = false;
	  if (_phase[m].arbitraryPar) {
	    fprintf(out, "\"P\"");
	    somethingPrinted = true;
	  }

	  // Will set this value to be only the bits relative to
	  // <_phase[m].hetParent> if this first assignment is non-zero:
	  uint8_t thisParHetAmbigPhase = _phase[m].ambigParPhase;
	  if (thisParHetAmbigPhase) {
	    // update to the bits relevant to this hetParent
	    thisParHetAmbigPhase >>= 2 * _phase[m].hetParent;
	    uint8_t isPI = _phase[m].hetParent >> 1;
	    thisParHetAmbigPhase &= isPI * 15 + (1-isPI) * 3;
	  }
	  if (thisParHetAmbigPhase) {
	    if (somethingPrinted)
	      fprintf(out, ",");

	    if (_phase[m].hetParent != 2) {
	      fprintf(out, "\"S\"");
	    }
	    else { // both parent's heterozygous, print details of swap types
	      char swapTypes[4] = { '?', '0', '1', 'B' };
	      uint8_t diff;
	      switch(thisParHetAmbigPhase) {
		// only a single bit in the ambiguous value cases:
		case 1:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 2:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 4:
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 8:
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		// two bits in the ambiguous value cases:
		case 3:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 5:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 6:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 9:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 10:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 12:
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		// three bits in the ambiguous value cases:
		case 7:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 11:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 13:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		case 14:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\",", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "\"S%c\"", swapTypes[diff]);
		  break;
		default:
		  // 0:  impossible to get to this branch
		  // 15: includes the current value which shouldn't happen
		  fprintf(out, "ERROR: ambigParPhase status %d\n",
			  thisParHetAmbigPhase);
		  break;
	      }
	    }
	    somethingPrinted = true;
	  }
	  if (_phase[m].ambigParHet) {
	    if (somethingPrinted)
	      fprintf(out, ",");

	    switch (_phase[m].ambigParHet) {
	      case 1:
		fprintf(out, "\"H0\"");
		break;
	      case 2:
		fprintf(out, "\"H1\"");
		break;
	      case 3:
		fprintf(out, "\"H01\"");
		break;
	      case 4:
		fprintf(out, "\"H2\"");
		break;
	      case 5:
		fprintf(out, "\"H02\"");
		break;
	      case 6:
		fprintf(out, "\"H12\"");
		break;

	      default:
		fprintf(out, "ERROR: ambigParHet status %d\n",
			_phase[m].ambigParHet);
		break;
	    }
	  }

	  fprintf(out, "]");
	}
	break;
      default:
	break;
    }
  }
  // end the "codes" value and the "<parent0_id>-<parent1_id>" value
  fprintf(out, "]}");

  for(int h = 0; h < numHaps; h++)
    delete [] hapStrs[h];
}

// Determines the alleles the two parents carry at <marker>. The caller has a
// character array with 3 entries: index 0 is allele 0, index 1 is the missing
// allele code (likely '0'), and index 2 is allele 1
void NuclearFamily::getParAlleles(int marker, uint8_t parAlleleIdx[2][2]) {
  uint8_t hetParent = _phase[marker].hetParent;
  uint8_t defaultPhase = _phase[marker].parentPhase;

  if (hetParent == 2) {
    // Both parents heterozygous
    int defaultParPhase[2] = { defaultPhase & 1, defaultPhase >> 1 };
    int idx0[2];

    // assuming no ambiguity, these are the parent allele indexes:
    for(int p = 0; p < 2; p++) {
      idx0[p] = 0 * (1 - defaultParPhase[p]) + 2 * defaultParPhase[p];
      parAlleleIdx[p][0] = idx0[p];
      parAlleleIdx[p][1] = 2 - idx0[p];
    }

    // any ambig par het means this parent could be homozygous:
    if (_phase[marker].ambigParHet) {
      for(int p = 0; p < 2; p++) {
	// If the ambiguous heterozygosity can include only <p>:
	uint8_t phaseBitsIfHetP = (_phase[marker].ambigParPhase >> (2*p)) & 3;
	if (phaseBitsIfHetP) {
	  // First: effects on <p>
	  if (phaseBitsIfHetP == 3) {
	    // Both phase types valid, so neither allele is certain:
	    parAlleleIdx[p][0] = parAlleleIdx[p][1] = 1; // 1 is missing idx
	  }
	  else {
	    // get the bit that's assigned in <phaseBitsIfHetP>;
	    // In this case it's trivial to do:
	    uint8_t phaseIfHetPar = phaseBitsIfHetP - 1;
	    if (phaseIfHetPar != defaultParPhase[p])
	      // opposite phase assignment when heterozygous only for the
	      // current parent, so both alleles are ambiguous:
	      parAlleleIdx[p][0] = parAlleleIdx[p][1] = 1; // 1 is missing idx
	  }

	  // Second: this phase possibility implies that <1-p> can be homozygous
	  // must set the allele that doesn't match the homozygous genotype
	  // to missing
	  int otherP = 1 - p;
	  for(int h = 0; h < 2; h++)
	    if (parAlleleIdx[otherP][h] != (_phase[marker].homParentGeno /3) *2)
	      parAlleleIdx[otherP][h] = 1;
	}
      }
    }
  }
  else {
    // One parent homozygous; we should know which genotype:
    assert(_phase[marker].homParentGeno != G_MISS);

    // assuming no ambiguity, these are the parent allele indexes:
    int idx0 = 0 * (1 - defaultPhase) + 2 * defaultPhase;
    parAlleleIdx[hetParent][0] = idx0;
    parAlleleIdx[hetParent][1] = 2 - idx0;
    parAlleleIdx[1 - hetParent][0] = parAlleleIdx[1 - hetParent][1] =
					(_phase[marker].homParentGeno / 3) * 2;

    // ambig par het? if so, haplotype allele index could be uncertain
    if (_phase[marker].ambigParHet) {
      // If the ambiguous heterozygosity is with only parent (1 - <hetParent>)
      // being heterozygous (H0 or H1):
      uint8_t otherP = 1 - hetParent;
      uint8_t phaseBitsIfHetOtherP =
			    (_phase[marker].ambigParPhase >> (2 * otherP)) & 3;
      if (phaseBitsIfHetOtherP) {
	// there exists a possibility with <hetParent> homozygous;
	// must set the allele that doesn't match the homozygous genotype
	// to missing
	for(int h = 0; h < 2; h++)
	  if (parAlleleIdx[hetParent][h] != (_phase[marker].homParentGeno /3)*2)
	    parAlleleIdx[hetParent][h] = 1; // 1 is missing idx

	// in this case, <otherP> is heterozygous; set uncertain alleles to
	// missing:
	if (phaseBitsIfHetOtherP == 3) {
	  // Both phase types valid, so neither allele is certain:
	  parAlleleIdx[otherP][0] = parAlleleIdx[otherP][1] = 1;
	}
	else {
	  // get the bit that's assigned in <phaseBitsIfHetOtherP>;
	  // In this case it's trivial to do:
	  uint8_t phaseIfOtherPHet = phaseBitsIfHetOtherP - 1;
	  int otherPHet_idx0 = 0 * (1 - phaseIfOtherPHet) + 2 *phaseIfOtherPHet;
	    // check both haplotypes in <otherP>, and set missing sites that
	    // differ
	    if (parAlleleIdx[otherP][0] != otherPHet_idx0)
	      parAlleleIdx[otherP][0] = 1; // 1 is missing idx
	    // TODO: else?
	    if (parAlleleIdx[otherP][1] != 2 - otherPHet_idx0)
	      parAlleleIdx[otherP][1] = 1; // 1 is missing idx
	}
      }

      // If the ambiguous heterozygosity is with both parents being
      // heterozygous (H2):
      uint8_t phaseBitsIfHet2 = _phase[marker].ambigParPhase >> (2 * 2);
      if (phaseBitsIfHet2) {
	uint8_t phaseSwapType = 0; // one of the possible swap types
	uint8_t phaseSwapTypeAssigned = 0; // has <phaseSwapType> been assigned?
	// is there more than one valid swap types? bits in <ambigSwap> indicate
	// which bits vary, with bit 0 corresponding to phase for parent 0, and
	// bit 1 for parent 1.
	uint8_t ambigSwap = 0;
	for(int curSwap = 0; curSwap < 4; curSwap++) {
	  uint8_t curIsValid = (phaseBitsIfHet2 >> curSwap) & 1;
	  phaseSwapType |= (1 - phaseSwapTypeAssigned) * curIsValid * curSwap;
	  phaseSwapTypeAssigned |= curIsValid;
	  ambigSwap |= curIsValid * (phaseSwapType ^ curSwap);
	}

	if (ambigSwap == 3) {
	  // all alleles completely uncertain:
	  parAlleleIdx[0][0] = parAlleleIdx[0][1] =
	    parAlleleIdx[1][0] = parAlleleIdx[1][1] = 1; // 1 is missing idx
	}
	else {
	  if (ambigSwap) {
	    // uncertain phase for one parent, set sites missing for that parent
	    uint8_t ambigPar = ambigSwap - 1;
	    parAlleleIdx[ambigPar][0] = parAlleleIdx[ambigPar][1] = 1;

	    uint8_t unambigPar = 1 - ambigPar;
	    uint8_t bothHetPhase_unambigPar = (phaseSwapType >> unambigPar) & 1;
	    int bothHet_unambig_idx0 = 0 * (1 - bothHetPhase_unambigPar) +
						    2 * bothHetPhase_unambigPar;
	    // check both haplotypes in <unambigPar>, and set missing sites that
	    // differ
	    if (parAlleleIdx[unambigPar][0] != bothHet_unambig_idx0)
	      parAlleleIdx[unambigPar][0] = 1; // 1 is missing idx
	    // TODO: else?
	    if (parAlleleIdx[unambigPar][1] != 2 - bothHet_unambig_idx0)
	      parAlleleIdx[unambigPar][1] = 1; // 1 is missing idx
	  }
	  else {
	    // only one phase type; determine whether the <hetParent> phase
	    // remains the same as the default
	    uint8_t bothHetPhase_defaultHetPar =(phaseSwapType >> hetParent) &1;
	    if (bothHetPhase_defaultHetPar != defaultPhase) 
	      // opposite phase assignment when heterozygous for both parents,
	      // so both alleles are ambiguous:
	      parAlleleIdx[hetParent][0] = parAlleleIdx[hetParent][1] = 1;

	    uint8_t homPar = 1 - hetParent;
	    uint8_t bothHetPhase_homPar = (phaseSwapType >> homPar) & 1;
	    int bothHet_hom_idx0 = 0 * (1 - bothHetPhase_homPar) +
						  2 * bothHetPhase_homPar;
	    // check both haplotypes in <homPar>, and set missing sites that
	    // differ
	    if (parAlleleIdx[homPar][0] != bothHet_hom_idx0)
	      parAlleleIdx[homPar][0] = 1; // 1 is missing idx
	    // TODO: else?
	    if (parAlleleIdx[homPar][1] != 2 - bothHet_hom_idx0)
	      parAlleleIdx[homPar][1] = 1; // 1 is missing idx
	  }
	}
      }
    } // ambiguous parent heterozygosity
  } // one parent heterozygous case
}

// Prints the haplotypes for <this> nuclear family to <out> in VCF format
void NuclearFamily::printPhasedVCF(FILE *out, const char *programName) {
  // print meta-information
  fprintf(out, "##fileformat=VCFv4.3\n");
  time_t time_raw_format;
  struct tm * ptr_time;
  time ( &time_raw_format );
  ptr_time = localtime ( &time_raw_format );
  fprintf(out, "##fileDate=%d%02d%02d\n", ptr_time->tm_year + 1900,
	  ptr_time->tm_mon + 1, ptr_time->tm_mday);
  fprintf(out, "##source=\"%s\"\n", programName);

  fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

  // print initial header
  fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  // print ids for dad, mom, and then children
  fprintf(out, "\t%s\t%s", _parents->first->getId(), _parents->second->getId());
  int numChildren = _children.length();
  for(int c = 0; c < numChildren; c++) {
    fprintf(out, "\t%s", _children[c]->getId());
  }
  fprintf(out, "\n");

  // print the haplotypes across all chromosomes
  for(int chrIdx = 0; chrIdx < Marker::getNumChroms(); chrIdx++) {
    const char *chrName = Marker::getChromName(chrIdx);

    int firstMarker = Marker::getFirstMarkerNum(chrIdx);
    int lastMarker = Marker::getLastMarkerNum(chrIdx);
    for(int m = firstMarker; m <= lastMarker; m++) {
      const Marker *theMarker = Marker::getMarker(m);
      const char *theAlleles = theMarker->getAlleleStr();

      fprintf(out, "%s\t%d\t%s\t%c\t%c\t.\tPASS\t.\tGT",
	      chrName, theMarker->getPhysPos(), theMarker->getName(),
	      theAlleles[0], theAlleles[2]);

      // For use with the printGeno() method: allows printing an unknown
      // genotype for untransmitted alleles at positions where we attempt to
      // impute the parents' genotype: alleles[1] is '.' or missing.
      // Note that in VCF format, we want to print an index to either ref or
      // alt alleles, so we use '0' and '1' for these values.
      char alleles[3] = { '0', '.', '1' };

      PhaseStatus status = _phase[m].status;
      switch(status) {
	case PHASE_AMBIG:
	case PHASE_ERROR:
	case PHASE_ERR_RECOMB:
	  // not phased -- code as missing
	  // print dad and mom's haplotypes and then children's
	  fprintf(out, "\t.|.\t.|.");
	  for(int c = 0; c < numChildren; c++) {
	    fprintf(out, "\t.|.");
	  }
	  break;

	case PHASE_UNINFORM:
	  {
	    // can impute at uninformative markers, not the others, using the
	    // <homParentGeno> field:
	    uint8_t imputeParGeno = _phase[m].homParentGeno;

	    // print the parent's genotypes (which are homozygous so phased)
	    uint8_t parentData = _phase[m].parentData;
	    uint8_t untrans = _phase[m].untransParHap;
	    for(int p = 0; p < 2; p++) {
	      uint8_t thisParGeno = parentData & 3;
	      uint8_t thisUntrans = untrans & 3;
	      // <isMissing> is 1 iff thisParGeno == 1:
	      uint8_t isMissing = (thisParGeno & 1) & ~(thisParGeno >> 1);
	      // VCF format seems not to support printing genotypes with one
	      // allele missing such as .|1, so we'll only impute when we can
	      // impute both alleles. That corresponds to the case where
	      // both were transmitted
	      // The following is only 1 if <thisUntrans> == 0, i.e., if both
	      // alleles were transmitted. Otherwise it's 0
	      uint8_t fullyTrans = (3 ^ thisUntrans) / 3;
	      uint8_t doImpute = isMissing * fullyTrans;
	      uint8_t genoToPrint = (1 - doImpute) * thisParGeno +
						      doImpute * imputeParGeno;
	      fprintf(out, "\t");
	      printGeno(out, alleles, genoToPrint, /*sep=*/ '|',
			isMissing * thisUntrans);
	      parentData >>= 2;
	      untrans >>= 2;
	    }

	    // print the children's genotypes
	    uint64_t childrenData = _phase[m].iv;
	    for(int c = 0; c < numChildren; c++) {
	      fprintf(out, "\t");
	      printGeno(out, alleles, childrenData & 3, /*sep=*/ '|',
			/*untrans=NA=*/ 0, /*swapHet=*/ _phase[m].uninfHetSwap);
	      childrenData >>= 2;
	    }
	  }
	  break;

	case PHASE_OK:
	  {
	    ///////////////////////////////////////////////////////////////////
	    // Standard phased marker

	    // Note about ambiguities:
	    // See large comment in printPhasedPed()

	    // first determine which alleles each parent has on each haplotype;
	    // Note that the <alleles> string has alleles at index 0 and 2
	    uint8_t hetParent = _phase[m].hetParent;
	    uint8_t parentPhase = _phase[m].parentPhase;
	    int parAlleleInds[2][2]; // Allele indexes: 0 for ref, 1 for alt
	    if (hetParent == 0 || hetParent == 1) {
	      int idx0 = 0 * (1 - parentPhase) + parentPhase;
	      parAlleleInds[hetParent][0] = idx0;
	      parAlleleInds[hetParent][1] = 1 - idx0;
	      assert(_phase[m].homParentGeno != G_MISS);
	      parAlleleInds[1 - hetParent][0] =
		parAlleleInds[1 - hetParent][1] = (_phase[m].homParentGeno / 3);
	    }
	    else {
	      int idx0[2] = { 0 * (1 - (parentPhase & 1)) + (parentPhase & 1),
			      0 * (1 - (parentPhase >>1)) + (parentPhase >>1) };
	      for(int p = 0; p < 2; p++) {
		parAlleleInds[p][0] = idx0[p];
		parAlleleInds[p][1] = 1 - idx0[p];
	      }
	    }

	    if (_phase[m].ambigParHet) {
	      // unclear which parent is heterozygous
	      // will show parents as missing
	      fprintf(out, "\t.|.\t.|.");
	    }
	    else {
	      // print parent's haplotypes
	      fprintf(out, "\t%d|%d\t%d|%d",
		      parAlleleInds[0][0], parAlleleInds[0][1],
		      parAlleleInds[1][0], parAlleleInds[1][1]);
	    }

	    // now print children's haplotypes
	    uint64_t iv = _phase[m].iv;
	    uint64_t ambigMiss = _phase[m].ambigMiss;
	    for(int c = 0; c < numChildren; c++) {
	      uint8_t curIV = iv & 3;
	      uint8_t curAmbigMiss = ambigMiss & 3;
	      int ivs[2] = { curIV & 1, curIV >> 1 };

	      if (_phase[m].ambigParHet || curAmbigMiss >= 2) {
		// Either which parent is heterozygous is unknown or
		// child's phase is ambiguous: print missing
		fprintf(out, "\t.|.");
	      }
	      else {
		fprintf(out, "\t%d|%d", parAlleleInds[0][ ivs[0] ],
		    parAlleleInds[1][ ivs[1] ]);
	      }
	      iv >>= 2;
	      ambigMiss >>= 2;
	    }
	  }
	  break;

	default:
	  fprintf(stderr, "ERROR: marker status %d\n", _phase[m].status);
	  exit(1);
	  break;
      }
      fprintf(out, "\n");
    }
  }
}

// Prints the haplotypes for <this> nuclear family to <out> in PLINK ped format
void NuclearFamily::printPhasedPed(FILE *out) {
  int famIdLen = _parents->first->getFamilyIdLength();

  // print the dad:
  fprintf(out, "%s\t%s\t0\t0\t%c\t0",
	  (famIdLen == 0) ? "0" : _parents->first->getId(),
	  &_parents->first->getId()[famIdLen+1], '1');
  printOnePedHap(out, /*p=*/0, /*c=*/-1);

  // print the mom:
  fprintf(out, "%s\t%s\t0\t0\t%c\t0",
	  (famIdLen == 0) ? "0" : _parents->second->getId(),
	  &_parents->second->getId()[famIdLen+1], '2');
  printOnePedHap(out, /*p=*/1, /*c=*/-1);

  // print the children:
  int numChildren = _children.length();
  for(int c = 0; c < numChildren; c++) {
    fprintf(out, "%s\t%s\t%s\t%s\t%c\t0",
	    (famIdLen == 0) ? "0" : _children[c]->getId(),
	    &_children[c]->getId()[famIdLen+1],
	    &_parents->first->getId()[famIdLen+1],
	    &_parents->second->getId()[famIdLen+1],
	    (_children[c]->getSex() == 'M') ? '1' : '2');
    printOnePedHap(out, /*p=*/-1, /*c=*/c);
  }
}

// Print the genotypes for a given parent <p> (0 or 1) or a given child <c> in
// PLINK ped format, terminating with an endline character.
void NuclearFamily::printOnePedHap(FILE *out, int p, int c) {
  assert(p < 0 || c < 0); // can only process one person, a parent or child

  int numMarkers = Marker::getNumMarkers();
  for(int m = 0; m < numMarkers; m++) {
    const char *theAlleles = Marker::getMarker(m)->getAlleleStr();
    // For use with the printGeno() method: allows printing an unknown genotype
    // for untransmitted alleles at positions where we attempt to impute the
    // parents' genotype: alleles[1] is '0' or unknown.
    char alleles[3] = { theAlleles[0], '0', theAlleles[2] };

    PhaseStatus status = _phase[m].status;
    switch(status) {
      case PHASE_AMBIG:
      case PHASE_ERROR:
      case PHASE_ERR_RECOMB:
	// not phased -- code as missing
	fprintf(out, "\t0\t0");
	break;

      case PHASE_UNINFORM:
	{
	  // can impute at uninformative markers, not the others, using the
	  // <homParentGeno> field:
	  uint8_t imputeParGeno = _phase[m].homParentGeno;

	  if (p >= 0) {
	    // print the parent's genotypes (which are homozygous so phased)
	    uint8_t parentData = _phase[m].parentData;
	    uint8_t untrans = _phase[m].untransParHap;

	    // must shift by 2 bits to get data for p == 1:
	    uint8_t thisParGeno = (parentData >> (2*p)) & 3;
	    uint8_t thisUntrans = (untrans >> (2*p)) & 3;
	    // <isMissing> is 1 iff thisParGeno == 1:
	    uint8_t isMissing = (thisParGeno & 1) & ~(thisParGeno >> 1);
	    // Will print either the parent genotype if it is not missing or
	    // the imputed genotype if it is missing (note that the default is
	    // for <imputeParGeno> is missing as well):
	    uint8_t genoToPrint = (1 - isMissing) * thisParGeno +
					      isMissing * imputeParGeno;
	    fprintf(out, "\t");
	    printGeno(out, alleles, genoToPrint, /*sep=*/ '\t',
		      isMissing * thisUntrans);
	  }
	  else {
	    // print child's haplotype
	    // must shift by 2*c bits to get this child's data:
	    uint8_t childGeno = (_phase[m].iv >> (2*c)) & 3;
	    fprintf(out, "\t");
	    printGeno(out, alleles, childGeno, /*sep=*/ '\t',
		      /*untrans=NA=*/ 0, /*swapHet=*/ _phase[m].uninfHetSwap);
	  }
	}
	break;

      case PHASE_OK:
	{
	  /////////////////////////////////////////////////////////////////////
	  // Standard phased marker

	  // Note about ambiguities:
	  // - arbitraryPar: which parent is which was chosen arbitarily. This
	  //   can only occur when neither parent has data. We'll just print the
	  //   phase as is, and note that if this happens internally, it can
	  //   lead to a switch error in which the parent's haplotypes get mixed
	  //   up with each other.
	  // - ambigParPhase: the parent's phase is potentially swappable. We'll
	  //   print the arbitrarily chosen phase here without warning
	  // - ambigParHet: which parent is heterozygous is ambiguous. We won't
	  //   print any haplotypes (print missing data instead).
	  // - ivFlippable: happens when different states yield minimum
	  //   recombinant haplotypes and differ in their inheritance vector
	  //   values. We'll stick with the convention adopted throughout these
	  //   ambiguities of choosing a phase arbitrarily among the
	  //   possibilities, and will therefore print the values corresponding
	  //   to the state stored here.
	  // - ambigMiss: samples that are ambiguous (or both ambiguous and
	  //   missing) do not have their haplotypes printed (print missing data
	  //   instead).

	  // first determine which alleles each parent has on each haplotype;
	  // Note that the <alleles> string has alleles at index 0 and 2
	  uint8_t hetParent = _phase[m].hetParent;
	  uint8_t parentPhase = _phase[m].parentPhase;
	  char parAlleles[2][2];
	  if (hetParent == 0 || hetParent == 1) {
	    int idx0 = 0 * (1 - parentPhase) + 2 * parentPhase;
	    parAlleles[hetParent][0] = alleles[idx0];
	    parAlleles[hetParent][1] = alleles[ 2 - idx0 ];
	    assert(_phase[m].homParentGeno != G_MISS);
	    parAlleles[1 - hetParent][0] = parAlleles[1 - hetParent][1] =
				   alleles[ (_phase[m].homParentGeno / 3) * 2 ];
	  }
	  else {
	    int idx0[2] = { 0 * (1 - (parentPhase & 1)) + 2 *(parentPhase & 1),
			  0 * (1 - (parentPhase >> 1)) + 2*(parentPhase >> 1) };
	    for(int p = 0; p < 2; p++) {
	      parAlleles[p][0] = alleles[ idx0[p] ];
	      parAlleles[p][1] = alleles[ 2 - idx0[p] ];
	    }
	  }

	  if (p >= 0) {
	    if (_phase[m].ambigParHet) {
	      // unclear which parent is heterozygous
	      // will show parents as missing
	      fprintf(out, "\t0\t0");
	    }
	    else {
	      // print parent's haplotypes
	      fprintf(out, "\t%c\t%c",
		      parAlleles[p][0], parAlleles[p][1]);
	    }
	  }
	  else {
	    // print child's haplotype
	    // must shift by 2*c bits to get this child's data:
	    uint8_t curIV = (_phase[m].iv >> (2*c)) & 3;
	    uint8_t curAmbigMiss = (_phase[m].ambigMiss >> (2*c)) & 3;
	    int ivs[2] = { curIV & 1, curIV >> 1 };

	    if (_phase[m].ambigParHet || curAmbigMiss >= 2) {
	      // Either which parent is heterozygous is unknown or
	      // child's phase is ambiguous: print missing
	      fprintf(out, "\t0\t0");
	    }
	    else {
	      fprintf(out, "\t%c\t%c", parAlleles[0][ ivs[0] ],
		      parAlleles[1][ ivs[1] ]);
	    }
	  }
	}
	break;

      default:
	fprintf(stderr, "ERROR: marker status %d\n", _phase[m].status);
	exit(1);
	break;
    }
  }
  fprintf(out, "\n");
}


void NuclearFamily::printIvCSV(FILE *out, int chrIdx) {
  int numChildren = _children.length();

  /////////////////////////////////////////////////////////////////////////
  // Print header
  int famIdLen = _parents->first->getFamilyIdLength();
  if (famIdLen == 0) {
    fprintf(out, ",%s,", _parents->first->getId());
    for(int c = 0; c < numChildren; c++)
      fprintf(out, "P-%s,", _children[c]->getId());
    fprintf(out, ",%s,", _parents->second->getId());
    for(int c = 0; c < numChildren; c++)
      fprintf(out, "M-%s,", _children[c]->getId());
  }
  else {
    // Print family id in first column, omit from all other samples
    fprintf(out, "%.*s,%s,", famIdLen, _parents->first->getId(),
	    &_parents->first->getId()[famIdLen+1]);
    for(int c = 0; c < numChildren; c++)
      fprintf(out, "P-%s,", &_children[c]->getId()[famIdLen+1]);
    fprintf(out, "%s,", &_parents->second->getId()[famIdLen+1]);
    for(int c = 0; c < numChildren; c++)
      fprintf(out, "M-%s,", &_children[c]->getId()[famIdLen+1]);
  }
  // Print column header for recombination count and the family id
  fprintf(out, "N_recomb,AmbigPar\n");

  /////////////////////////////////////////////////////////////////////////
  // Print inheritance vectors
  const char *chrName = Marker::getChromName(chrIdx);

  bool sawPrevIV = false;
  uint64_t prevIV = 0;
  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);
  for(int m = firstMarker; m <= lastMarker; m++) {
    fprintf(out, "%s,%s,", chrName, Marker::getMarker(m)->getName());

    PhaseStatus status = _phase[m].status;
    uint64_t iv = _phase[m].iv;
    uint64_t ivFlippable = _phase[m].ivFlippable;
    uint64_t ambigMiss = _phase[m].ambigMiss;
    uint8_t hetParent = _phase[m].hetParent;
    uint8_t numRecombs = _phase[m].numRecombs;
    char status_indicator[NUM_PHASE_STATUS] = { '|', '|', '?', 'E', 'R' };
    // Row indicates 0/1 whether the marker is informative
    // Column indicates which haplotype was transmitted with '?' used for
    // flippable <iv> values where ambiguities among states lead to varying
    // haplotype transmissions and leave the true transmitted haplotype
    // uncertain.
    char ivLabel[2][4] = { { 'a', 'b', '?', '?' },
			   { 'A', 'B', '?', '?' } };
    switch(status) {
      case PHASE_AMBIG:
      case PHASE_ERROR:
      case PHASE_ERR_RECOMB:
	for (int p = 0; p < 2; p++) {
	  for(int c = 0; c < numChildren; c++) {
	    fprintf(out, "%c,", status_indicator[status]);
	  }
	  if (p == 0)
	    fprintf(out, ",");
	}
	fprintf(out, "0,\n");
	break;
      case PHASE_UNINFORM:
	if (!sawPrevIV) {
	  for (int p = 0; p < 2; p++) {
	    for(int c = 0; c < numChildren; c++) {
	      fprintf(out, "0,");
	    }
	    if (p == 0)
	      fprintf(out, ",");
	  }
	  fprintf(out, "0,\n");
	  break;
	}
	// else: fall through to PHASE_OK case to print <prevIV>
	iv = prevIV;
	numRecombs = 0;
      case PHASE_OK:
	for (int p = 0; p < 2; p++) {
	  // If the marker is informative for p, print either 'A' or 'B'
	  // otherwise print lowercase 'a' or 'b'
	  int informative = 0;
	  if (status == PHASE_OK && (hetParent == 2 || hetParent == p))
	    informative = 1;

	  for(int c = 0; c < numChildren; c++) {
	    uint8_t curIV = (iv >> (2 * c + p)) & 1;
	    uint8_t curIVflip = (ivFlippable >> (2 * c + p)) & 1;
	    uint8_t curAmbigMiss = (ambigMiss >> (2 * c)) & 3;
	    uint8_t missing = curAmbigMiss & 1;
	    uint8_t ambig = curAmbigMiss >> 1;

	    // Note: 'a' + 1 == 'b' and 'A' + 1 == 'B' for ASCII chars
	    char curIVchar = ivLabel[informative][ curIV + 2 * curIVflip ];
	    // Will print '0' if the child is missing data
	    char printChar = missing * '0' + (1 - missing) * curIVchar;
	    // Will append a '?' if the child's IV is ambiguous
	    fprintf(out, "%c%.*s,", printChar, ambig, "?");
	  }
	  if (p == 0)
	    fprintf(out, ",");
	}
	fprintf(out, "%d,", numRecombs);
	// print any ambiguous status values
	if (status == PHASE_OK && //<ambigParHet> field only defined in PHASE_OK
	    (_phase[m].ambigParPhase | _phase[m].ambigParHet |
						      _phase[m].arbitraryPar)) {
	  if (_phase[m].arbitraryPar) {
	    fprintf(out, "ParArbitrary");
	  }

	  // Will set this value to be only the bits relative to
	  // <_phase[m].hetParent> if this first assignment is non-zero:
	  uint8_t thisParHetAmbigPhase = _phase[m].ambigParPhase;
	  if (thisParHetAmbigPhase) {
	    // update to the bits relevant to this hetParent
	    thisParHetAmbigPhase >>= 2 * _phase[m].hetParent;
	    uint8_t isPI = _phase[m].hetParent >> 1;
	    thisParHetAmbigPhase &= isPI * 15 + (1-isPI) * 3;
	  }
	  if (thisParHetAmbigPhase) {
	    if (_phase[m].hetParent != 2) {
	      fprintf(out, "Swap");
	    }
	    else { // both parent's heterozygous, print details of swap types
	      fprintf(out, "S");
	      char swapTypes[4] = { '?', '0', '1', 'B' };
	      uint8_t diff;
	      switch(thisParHetAmbigPhase) {
		// only a single bit in the ambiguous value cases:
		case 1:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 2:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 4:
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 8:
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		// two bits in the ambiguous value cases:
		case 3:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 5:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 6:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 9:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 10:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 12:
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		// three bits in the ambiguous value cases:
		case 7:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 11:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 13:
		  diff = 0 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		case 14:
		  diff = 1 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 2 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  diff = 3 ^ _phase[m].parentPhase;
		  fprintf(out, "%c", swapTypes[diff]);
		  break;
		default:
		  // 0:  impossible to get to this branch
		  // 15: includes the current value which shouldn't happen
		  fprintf(out, "ERROR: ambigParPhase status %d\n",
			  thisParHetAmbigPhase);
		  break;
	      }
	    }
	    if (_phase[m].ambigParHet)
	      fprintf(out, "_");
	  }
	  if (_phase[m].ambigParHet) {
	    switch (_phase[m].ambigParHet) {
	      case 1:
		fprintf(out, "HetPar0");
		break;
	      case 2:
		fprintf(out, "HetPar1");
		break;
	      case 3:
		fprintf(out, "HetPar0_1");
		break;
	      case 4:
		fprintf(out, "HetPar2");
		break;
	      case 5:
		fprintf(out, "HetPar0_2");
		break;
	      case 6:
		fprintf(out, "HetPar1_2");
		break;

	      default:
		fprintf(out, "ERROR: ambigParHet status %d\n",
			_phase[m].ambigParHet);
		break;
	    }
	  }
	}
	fprintf(out, "\n");

	prevIV = iv;
	sawPrevIV = true;
	break;

      case NUM_PHASE_STATUS:
      default:
	fprintf(out, "ERROR: marker status %d\n", _phase[m].status);
	break;
    }
  }
}
