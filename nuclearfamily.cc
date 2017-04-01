// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

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

// Prints the haplotypes for <this> to <out>
void NuclearFamily::printHaplotypes(FILE *out) {
  int numChildren = _children.length();

  /////////////////////////////////////////////////////////////////////////
  // Print header
  fprintf(out, "%s %s", _parents->first->getId(), _parents->second->getId());
  for(int c = 0; c < numChildren; c++) {
    fprintf(out, " %s", _children[c]->getId());
  }
  fprintf(out, "\n");

  /////////////////////////////////////////////////////////////////////////
  // Print haplotypes
  // TODO: better name?
  const char ambigMissType[4] = { '/', '_', '?', '!' };
  const char markerType[3] = { '0', '1', 'B' };

  for(int c = 0; c < Marker::getNumChroms(); c++) {
    const char *chrName = Marker::getChromName(c);

    int prevM_OK = -1;
    int lastMarker = Marker::getLastMarkerNum(c);
    for(int m = Marker::getFirstMarkerNum(c); m <= lastMarker; m++) {
      const char *alleles = Marker::getMarker(m)->getAlleleStr();

      fprintf(out, "%2s %-5d ", chrName, m);

      PhaseStatus status = _phase[m].status;
      // TODO: can move declarations down if we encapsulate cases in { }
      uint8_t parentData, hetParent, parentPhase;
      uint64_t childrenData, iv, ambigMiss;
      char parAlleles[2][2];
      switch(status) {
	case PHASE_UNINFORM:
	case PHASE_AMBIG:
	case PHASE_ERROR:
	case PHASE_ERR_RECOMB:
	  /////////////////////////////////////////////////////////////////////
	  // Not phased / trivially phased cases:

	  // print the parent's genotypes
	  parentData = _phase[m].parentData;
	  printGeno(out, alleles, parentData & 3);
	  fprintf(out, "   ");
	  printGeno(out, alleles, parentData >> 2);
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
	      break;
	  }

	  // print the children's genotypes
	  childrenData = _phase[m].iv;
	  for(int c = 0; c < numChildren; c++) {
	    fprintf(out, " ");
	    printGeno(out, alleles, childrenData & 3);
	    fprintf(out, "   ");
	    childrenData >>= 2;
	  }
	  fprintf(out, "\n");
	  break;

	case PHASE_OK:
	  /////////////////////////////////////////////////////////////////////
	  // Standard phased marker

	  // first determine which alleles each parent has on each haplotype;
	  // Note that the <alleles> string has alleles at index 0 and 2
	  hetParent = _phase[m].hetParent;
	  parentPhase = _phase[m].parentPhase;
	  if (hetParent == 0 || hetParent == 1) {
	    int ind0 = 0 * (1 - parentPhase) + 2 * parentPhase;
	    parAlleles[hetParent][0] = alleles[ind0];
	    parAlleles[hetParent][1] = alleles[ 2 - ind0 ];
	    // TODO: eventually require/assert that this value not be G_MISS
	    if (_phase[m].homParentGeno == G_MISS)
	      parAlleles[1 - hetParent][0] = parAlleles[1 - hetParent][1] = '0';
	    else
	      parAlleles[1 - hetParent][0] = parAlleles[1 - hetParent][1] =
				   alleles[ (_phase[m].homParentGeno / 3) * 2 ];
	  }
	  else {
	    assert(hetParent == 2);
	    int ind0[2] = { 0 * (1 - (parentPhase & 1)) + 2 *(parentPhase & 1),
			  0 * (1 - (parentPhase >> 1)) + 2*(parentPhase >> 1) };
	    for(int p = 0; p < 2; p++) {
	      parAlleles[p][0] = alleles[ ind0[p] ];
	      parAlleles[p][1] = alleles[ 2 - ind0[p] ];
	    }
	  }

	  // print parent's haplotypes
	  // TODO: indicate when parent's genotypes were missing
	  fprintf(out, "%c/%c %c %c/%c |", parAlleles[0][0], parAlleles[0][1],
		  markerType[hetParent], parAlleles[1][0], parAlleles[1][1]);

	  // now print children's haplotypes
	  iv = _phase[m].iv;
	  ambigMiss = _phase[m].ambigMiss;
	  for(int c = 0; c < numChildren; c++) {
	    uint8_t curIV = iv & 3;
	    uint8_t curAmbigMiss = ambigMiss & 3;
	    int ivs[2] = { curIV & 1, curIV >> 1 };
	    fprintf(out, " %c%c%c %c%c", parAlleles[0][ ivs[0] ],
		    ambigMissType[curAmbigMiss], parAlleles[1][ ivs[1] ],
		    (char) ('A' + ivs[0]), (char) ('A' + ivs[1]));
	    iv >>= 2;
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
	  fprintf(out, "ERROR: marker status %d\n", _phase[m].status);
	  break;
      }
    }
  }
}
