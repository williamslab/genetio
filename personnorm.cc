// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include "personnorm.h"
#include "marker.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<PersonNorm *> PersonNorm::_allIndivs;
Hashtable<char *, PersonNorm *> PersonNorm::_idToPerson(2003, stringHash,
							stringcmp);
PersonNorm::fam_ht PersonNorm::_families;

PersonNorm::PersonNorm(char *id, char sex,int popIndex,
		       short familyIdLength, bool allocData) :
		     SuperPerson(id, sex, popIndex, familyIdLength) {
  if (!_ignore) {
    if (allocData) {
      int numChroms = Marker::getNumChroms();
      _geno = new Genotype *[numChroms];
      for(int c = 0; c < numChroms; c++) {
	int numChrMarkers = Marker::getNumChromMarkers(c);
	assert(numChrMarkers > 0);
	_geno[c] = new Genotype[numChrMarkers];
      }
    }

    if (_idToPerson.lookup(_id)) {
      fprintf(stderr, "\nERROR: multiple individuals with id %s!\n", _id);
      exit(3);
    }

    _idToPerson.add(_id, this);
  }
}

PersonNorm::~PersonNorm() {
  if (!_ignore) {
    int numChroms = Marker::getNumChroms();
    for(int c = 0; c < numChroms; c++) {
      delete [] _geno[c];
    }
  }
}

void PersonNorm::setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
			     int chromMarkerIdx, int geno[2]) {
  if (_ignore)
    return; // no need/nowhere to store genotypes

  for(int h = 0; h < 2; h++)
    _geno[chromIdx][chromMarkerIdx][h] = geno[h];
}

// Given the parents of <this>, uses a hashtable to store the family
// relationships (for later phasing).
void PersonNorm::setParents(char *familyid, PersonNorm *parents[2],
			    int numParents, bool &warningPrinted, FILE *log,
			    int *numMendelError, int *numMendelCounted) {
  // This method does not identify non-Mendelian errors, so
  // ensure that the caller isn't requesting these counts
  assert(numMendelError == NULL && numMendelCounted == NULL);

  if (parents[0] == NULL || parents[1] == NULL)
    // Only process families where the ids of both parents is known; otherwise
    // can't identify siblings
    return;

  par_pair_real parentsKey(parents[0], parents[1]);

  fam_ht_iter it = _families.find( &parentsKey );
  if (it == _families.end()) {
    // No family yet for these parents; create and insert it
    dynarray<PersonNorm *> *newChildren = new dynarray<PersonNorm *>;
    newChildren->append(this); // add first child
    par_pair newParentsKey = new par_pair_real(parents[0], parents[1]);
    _families[ newParentsKey ] = newChildren;
  }
  else {
    // Have a family: append child to list of children
    // <*it> is a pair with <it->first> being the key and <it->second> the value
    it->second->append(this);
  }
}

// TODO: comment
void PersonNorm::setXHetToMissing(FILE *log, int *numHets, int *numCalls) {
  // TODO!
}
