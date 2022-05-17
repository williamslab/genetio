// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <unordered_map>
#include "personbulk.h"
#include "dynarray.h"

#ifndef NUCLEARFAMILY_H
#define NUCLEARFAMILY_H

// See comments for the State structure in phaser.h from the hapi2 repository
// for more details on the meaning of these values
enum PhaseStatus {
  PHASE_OK,
  PHASE_UNINFORM,
  PHASE_AMBIG,
  PHASE_ERROR,
  PHASE_ERR_RECOMB,
  PHASE_X_SPECIAL,
  NUM_PHASE_STATUS
};

struct PhaseVals {
  // Inheritance vector for PHASE_OK status and contains the children's
  // genotype data in PLINK bed format (2 bits per child) otherwise:
  uint64_t iv;
  uint64_t ambigMiss;
  // For ambiguous states, which <iv> values differ among these states at the
  // current marker
  uint64_t ivFlippable;
  // Stores parent data for non-PHASE_OK. For PHASE_OK, indicates if missing
  uint8_t parentData;    // Can fit in 4 bits
  uint8_t hetParent;     // Can fit in 2 bits
  // For PHASE_UNINFORM, when one or both parents are missing data and we can
  // impute the missing parent(s)'s homozygous genotypes, the following is that
  // genotype:
  uint8_t homParentGeno; // Can fit in 2 bits
  // For PHASE_UNINFORM, for the imputation procedure, we determine which
  // parent haplotypes were untransmitted; the data are ambiguous as to the
  // genotype of any such haplotypes, and so we don't impute those alleles
  uint8_t untransParHap; // Can fit in 4 bits
  // For PHASE_UNINFORM, both parents are homozygous, so when a child is
  // heterozygous, it is possible to use one (or both) of the parent's
  // genotypes to deduce what the phase of heterozygous children is.
  // If <uninfHetSwap> == 0, the default print order for heterozygous genotypes
  // gives the correct phase, otherwise this should be swapped.
  uint8_t uninfHetSwap;  // Can fit in 1 bit
  // For PHASE_AMBIG, are the parents homozygous? Based on the IVs at nearby
  // informative markers
  bool ambigBothParHomozy;
  uint8_t parentPhase;   // Can fit in 2 bits
  PhaseStatus status;    // Can fit in 3 bits
  // TODO: recalculate this when needed?
  uint8_t numRecombs;    // Since previous marker
  uint8_t ambigParHet;   // Can fit in 3 bit
  // Each bit indicates presence or absence of different equal probability
  // <parentPhase> values; represents the phase possibilities for all three
  // possible parent heterozygosity statuses. See Phaser::ambigParPhase in the
  // HAPI2 repository
  uint8_t ambigParPhase; // Requires 8 bits
  uint8_t arbitraryPar;  // Can fit in 1 bit
  // TODO; either compress the above so they fit in 3 words or remove numRecombs
};

class NuclearFamily {
  public:
    //////////////////////////////////////////////////////////////////
    // type definitions
    //////////////////////////////////////////////////////////////////

    typedef typename std::pair<PersonBulk*,PersonBulk*>* par_pair; //parent pair
    typedef typename std::pair<PersonBulk*,PersonBulk*> par_pair_real;
    struct EqParPair {
      bool operator()(const par_pair k1, const par_pair k2) const {
	return k1 == k2 || (k1 && k2 && k1->first == k2->first &&
						      k1->second == k2->second);
      }
    };
    struct HashParPair {
      size_t operator()(NuclearFamily::par_pair const key) const {
	// make a better hash function?
	return std::hash<PersonBulk*>{}(key->first) +
	       std::hash<PersonBulk*>{}(key->second);
      }
    };
    typedef std::unordered_map<par_pair, NuclearFamily *,
			       HashParPair, EqParPair> fam_ht;
    typedef fam_ht::const_iterator fam_ht_iter;


    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void addChild(PersonBulk *parents[2], PersonBulk *child);

    static fam_ht_iter familyIter() { return _families.begin(); }
    static fam_ht_iter familyIterEnd() { return _families.end(); }
    static fam_ht::size_type numFamilies() { return _families.size(); }

    // not needed -- only delete when program done: OS will manage
    static void cleanUp() {
      for(fam_ht_iter it = _families.begin(); it != _families.end(); it++) {
        delete it->first;
        delete it->second;
      }
    }

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    // not needed -- only delete when program done: OS will manage
    ~NuclearFamily() {
      if (_phase)  delete [] _phase;
    }

    int numChildren() { return _children.length(); }
    void initPhase() { _phase = new PhaseVals[ Marker::getNumMarkers() ]; }
    void deletePhase() { delete [] _phase; _phase = NULL; }

    // As below, <missing> should have the lower order bit for each child that
    // is missing set to 1.
    void setStatus(int marker, PhaseStatus status, uint8_t parentData,
		   uint64_t childrenData, uint64_t missing) {
      // should only use this method to set status besides PHASE_OK and
      // PHASE_UNINFORM
      assert(status != PHASE_OK && status != PHASE_UNINFORM);
      _phase[marker].iv = childrenData;
      _phase[marker].ambigMiss = missing;
      _phase[marker].parentData = parentData;
      _phase[marker].status = status;
    }

    PhaseStatus getStatus(int marker) {
      return _phase[marker].status;
    }

    // For setting a marker's phase status to be uninformative; we require a
    // value for <homParentGeno> which is used for imputing the genotypes of
    // the parents when one or both are missing data
    void setUninform(int marker, uint8_t parentData, uint64_t childrenData,
		     uint64_t missing, uint8_t homParentGeno,
		     uint8_t uninfHetSwap) {
      _phase[marker].iv = childrenData;
      _phase[marker].ambigMiss = missing;
      _phase[marker].parentData = parentData;
      _phase[marker].homParentGeno = homParentGeno;
      _phase[marker].uninfHetSwap = uninfHetSwap;
      _phase[marker].status = PHASE_UNINFORM;
    }

    // For setting a marker's phase status to be the special X type. See
    // HAPI2's Phaser::getMarkerTypeX() method for more details
    void setSpecialX(int marker, uint8_t parentData, uint64_t childrenData,
		     uint64_t missing, uint8_t homParentGeno,
		     uint8_t uninfHetSwap) {
      _phase[marker].iv = childrenData;
      _phase[marker].ambigMiss = missing;
      _phase[marker].parentData = parentData;
      _phase[marker].homParentGeno = homParentGeno;
      _phase[marker].uninfHetSwap = uninfHetSwap;
      _phase[marker].status = PHASE_X_SPECIAL;
    }

    // <ambig> should have the higher order bit (of the two bits allotted each
    // child) set to 1 if the corresponding child has ambiguous phase.
    // <missing> should have the lower order bit set to 1 if the corresponding
    // child is missing data.
    void setPhase(int marker, uint64_t iv, uint64_t ambig, uint64_t missing,
		  uint64_t ivFlippable, uint8_t parMissing, uint8_t hetParent,
		  uint8_t homParentGeno, uint8_t parentPhase,
		  uint8_t numRecombs, uint8_t ambigParHet,
		  uint8_t ambigParPhase, uint8_t arbitraryPar) {
      _phase[marker].iv = iv;
      _phase[marker].ambigMiss = ambig | missing;
      _phase[marker].ivFlippable = ivFlippable;
      _phase[marker].parentData = parMissing;
      _phase[marker].hetParent = hetParent;
      _phase[marker].homParentGeno = homParentGeno;
      _phase[marker].parentPhase = parentPhase;
      _phase[marker].status = PHASE_OK;
      _phase[marker].numRecombs = numRecombs;
      _phase[marker].ambigParHet = ambigParHet;
      _phase[marker].ambigParPhase = ambigParPhase;
      _phase[marker].arbitraryPar = arbitraryPar;
    }

    // Sets the <untransParHap> field in <_phase[marker]>. This indicates which
    // parent haplotypes were transmitted to the children at the corresponding
    // marker. This used to help impute the parent(s)'s genotypes at
    // uninformative markers.
    // Also sets <ambigBothParHomozy>.
    void setUntransPar(int marker, uint8_t untransParHap,
		       bool ambigBothParHomozy) {
      _phase[marker].untransParHap = untransParHap;
      _phase[marker].ambigBothParHomozy = ambigBothParHomozy;
    }

    // Same as above but only assigns untransParHap
    void setUntransPar(int marker, uint8_t untransParHap) {
      _phase[marker].untransParHap = untransParHap;
    }

    // Support for assigning one haplotype transmission regions (associated
    // values sometimes need to be modified after their initial assignment)
    void updateOHTVals(int marker, uint64_t ivFlippable, uint8_t ambigParHet,
		       uint8_t ambigParPhase) {
      _phase[marker].ivFlippable |= ivFlippable;
      _phase[marker].ambigParHet |= ambigParHet;
      _phase[marker].ambigParPhase |= ambigParPhase;
    }

    const PhaseVals &getPhase(int marker) {
      return _phase[marker];
    }

    void printHapTxt(FILE *out, int chrIdx);
    void printHapJson(FILE *out, bool withChildren);
    void printPhasedPed(FILE *out);
    void printPhasedVCF(FILE *out, const char *progamName);
    void printIvCSV(FILE *out, int chrIdx);

    //////////////////////////////////////////////////////////////////
    // public variables
    //////////////////////////////////////////////////////////////////

    // TODO: make private?
    par_pair _parents;
    dynarray<PersonBulk *> _children;

  private:
    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // Hash table indexed by a par_pair (pair of parent PersonBulk* objects)
    // with the value being a dynarray containin the PesonBulk* objects for
    // the children of the couple
    static fam_ht _families;

    //////////////////////////////////////////////////////////////////
    // private methods
    //////////////////////////////////////////////////////////////////

    void getParAlleles(int marker, uint8_t parAlleleIdx[2][2], bool onChrX);
    void printOnePedHap(FILE *out, int p, int c);

    // alleles is expected to contain in elements 0 and 2 the two alleles for
    // this marker (single characters) and in element 1, a missing data
    // character (such as '0'). The latter enables for indicating that the
    // allele is unknown when the parent did not transmit a potentially imputed
    // allele.
    void printGeno(FILE *out, const char *alleles, uint8_t genotype,
		   char sep, uint8_t untrans, uint8_t swapHet = 0) {
      switch(genotype) {
	case G_HOM0:
	  fprintf(out, "%c%c%c", alleles[untrans & 1], sep,
		  alleles[(untrans & 2) >> 1]);
	  break;
	case G_MISS:
	  fprintf(out, "%c%c%c", alleles[1], sep, alleles[1]);
	  break;
	case G_HET:
	  // if <swapHet> == 0:
//	  fprintf(out, "%c%c%c", alleles[untrans & 1], sep,
//		  alleles[2 - ((untrans & 2) >> 1)]);
	  {
	    int idxs[2] = { (untrans & 1) + (1 - untrans) * swapHet * 2,
			  (untrans & 1) + (1 - untrans) * (1 - swapHet) * 2 };
	    fprintf(out, "%c%c%c", alleles[ idxs[0] ], sep, alleles[ idxs[1] ]);
	  }
	  break;
	case G_HOM1:
	  fprintf(out, "%c%c%c", alleles[2 - (untrans & 1)], sep,
		  alleles[2 - ((untrans & 2) >> 1)]);
	  break;
      }
    }

    //////////////////////////////////////////////////////////////////
    // private variables
    //////////////////////////////////////////////////////////////////

    // Stores phase for the family in the same format as used to phase it
    PhaseVals *_phase;
};

#endif // NUCLEARFAMILY_H
