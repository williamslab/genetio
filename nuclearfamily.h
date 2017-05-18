// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

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
  uint8_t homParentGeno; // Can fit in 2 bits
  uint8_t parentPhase;   // Can fit in 2 bits
  PhaseStatus status;    // Can fit in 2 bits
  // TODO: recalculate this when needed?
  uint8_t numRecombs;    // Since previous marker
  uint8_t ambigParHet;   // Can fit in 1 bit
  // Each bit indicates presence or absence of different equal probability
  // <parentPhase> values
  uint8_t ambigParPhase; // Can fit in 4 bits
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
    struct eqParPair {
      bool operator()(const par_pair k1, const par_pair k2) const {
	return k1 == k2 || (k1 && k2 && k1->first == k2->first &&
						      k1->second == k2->second);
      }
    };
    struct hashParPair {
      size_t operator()(NuclearFamily::par_pair const key) const {
	// make a better hash function?
	return std::tr1::hash<PersonBulk*>{}(key->first) +
	       std::tr1::hash<PersonBulk*>{}(key->second);
      }
    };
    typedef typename google::dense_hash_map<par_pair, NuclearFamily *,
					    hashParPair, eqParPair> fam_ht;
    typedef typename fam_ht::const_iterator fam_ht_iter;


    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void init() {
      _families.set_empty_key(NULL);
    }

    static void addChild(PersonBulk *parents[2], PersonBulk *child);

    static fam_ht_iter familyIter() { return _families.begin(); }
    static fam_ht_iter familyIterEnd() { return _families.end(); }
    static fam_ht::size_type numFamilies() { return _families.size(); }

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    int numChildren() { return _children.length(); }
    void initFam() { _phase = new PhaseVals[ Marker::getNumMarkers() ]; }

    // As below, <missing> should have the lower order bit for each child that
    // is missing set to 1.
    void setStatus(int marker, PhaseStatus status, uint8_t parentData,
		   uint64_t childrenData, uint64_t missing) {
      // should only use this method to set bad status
      assert(status != PHASE_OK);
      _phase[marker].iv = childrenData;
      _phase[marker].ambigMiss = missing;
      _phase[marker].parentData = parentData;
      _phase[marker].status = status;
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

    const PhaseVals &getPhase(int marker) {
      return _phase[marker];
    }

    void printHapTxt(FILE *out, int chrIdx);
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

    void printGeno(FILE *out, const char *alleles, uint8_t genotype) {
      switch(genotype) {
	case G_HOM0:
	  fprintf(out, "%c/%c", alleles[0], alleles[0]);
	  break;
	case G_MISS:
	  fprintf(out, "0/0");
	  break;
	case G_HET:
	  fprintf(out, "%c/%c", alleles[0], alleles[2]);
	  break;
	case G_HOM1:
	  fprintf(out, "%c/%c", alleles[2], alleles[2]);
	  break;
      }
    }

    //////////////////////////////////////////////////////////////////
    // public variables
    //////////////////////////////////////////////////////////////////

    // Stores phase for the family in the same format as used to phase it
    PhaseVals *_phase;
};

#endif // NUCLEARFAMILY_H
