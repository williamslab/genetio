// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include "personbulk.h"
#include "marker.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<PersonBulk *> PersonBulk::_allIndivs;
Hashtable<char *, PersonBulk *> PersonBulk::_idToPerson(2003, stringHash,
							stringcmp);
PersonBulk::fam_ht PersonBulk::_families;
uint8_t *PersonBulk::_data;
int      PersonBulk::_bytesPerMarker;

PersonBulk::PersonBulk(char *id, char sex, int popIndex, uint32_t sampNum,
		       short familyIdLength) :
		     SuperPerson(id, sex, popIndex, familyIdLength) {
  if (!_ignore) {
    if (_idToPerson.lookup(_id)) {
      fprintf(stderr, "\nERROR: multiple individuals with id %s!\n", _id);
      exit(3);
    }

    _sampNum = sampNum;

    _idToPerson.add(_id, this);
  }
}

PersonBulk::~PersonBulk() {
}

// Given the parents of <this>, uses a hashtable to store the family
// relationships (for later phasing).
void PersonBulk::setParents(char *familyid, PersonBulk *parents[2],
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
    dynarray<PersonBulk *> *newChildren = new dynarray<PersonBulk *>;
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
void PersonBulk::setXHetToMissing(FILE *log, int *numHets, int *numCalls) {
  // TODO!
}
