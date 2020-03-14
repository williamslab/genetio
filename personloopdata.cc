// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include "personloopdata.h"
#include "marker.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<PersonLoopData *> PersonLoopData::_allIndivs;
std::unordered_map<const char *, PersonLoopData *, HashString, EqualString>
						    PersonLoopData::_idToPerson;

PersonLoopData::PersonLoopData(char *id, char sex, int popIndex,
				   uint32_t sampNum, short familyIdLength) :
		     SuperPerson(id, sex, popIndex, familyIdLength) {
  if (!_ignore) {
    if (_idToPerson.find(_id) != _idToPerson.end()) {
      fprintf(stderr, "\nERROR: multiple individuals with id %s!\n", _id);
      exit(3);
    }

    _sampNum = sampNum;

    _idToPerson[_id] = this;
  }
}

PersonLoopData::~PersonLoopData() {
}

// Parent ids not used, so basically does nothing
void PersonLoopData::setParents(char *familyid, PersonLoopData *parents[2],
				  int numParents, bool &warningPrinted,
				  FILE *log, int *numMendelError,
				  int *numMendelCounted) {
  // This method does not identify non-Mendelian errors, so
  // ensure that the caller isn't requesting these counts
  assert(numMendelError == NULL && numMendelCounted == NULL);
}

// TODO: comment
void PersonLoopData::setXHetToMissing(FILE *log, int *numHets,
					int *numCalls) {
  // TODO!
}
