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
  PHASE_AMBIG,
  PHASE_ERROR,
};

struct PhaseVals {
  uint64_t iv;
  uint64_t ambig;
  uint8_t hetParent : 2;
  uint8_t parentPhase : 2;
  PhaseStatus status : 2;
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
	return k1 == k2 || (k1->first == k2->first && k1->second == k2->second);
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

    void setStatus(int marker, PhaseStatus status) {
      // should only use this method to set bad status
      assert(status != PHASE_OK);
      _phase[marker].status = status;
    }

    void setPhase(int marker, uint64_t iv, uint64_t ambig, uint8_t hetParent,
		  uint8_t parentPhase) {
      _phase[marker].iv = iv;
      _phase[marker].ambig = ambig;
      _phase[marker].hetParent = hetParent;
      _phase[marker].parentPhase = parentPhase;
      _phase[marker].status = PHASE_OK;
    }

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
    // public variables
    //////////////////////////////////////////////////////////////////

    // Stores phase for the family in the same format as used to phase it
    PhaseVals *_phase;
};

#endif // NUCLEARFAMILY_H
