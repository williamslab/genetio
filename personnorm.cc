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


PersonNorm::PersonNorm(char *id, char gender,int popIndex,
		       short familyIdLength) : 
		     SuperPerson(id, gender, popIndex, familyIdLength) {
  if (!_ignore) {
    for(int c = 1; c <= LAST_CHROM; c++) {
      int numChrMarkers = Marker::getNumChromMarkers(c);
      if (numChrMarkers > 0) {
	_geno[c] = new Genotype[numChrMarkers];
      }
      else
	_geno[c] = NULL;
    }

    _parents[0] = _parents[1] = NULL;

    if (_idToPerson.lookup(_id)) {
      fprintf(stderr, "\nERROR: multiple individuals with id %s!\n", _id);
      exit(3);
    }

    _idToPerson.add(_id, this);
  }
}

PersonNorm::~PersonNorm() {
  if (!_ignore) {
    for(int c = 1; c <= LAST_CHROM; c++) {
      if (_geno[c] != NULL)
	delete [] _geno[c];
    }
  }
}

void PersonNorm::setGenotype(int hapChunkNum, int chunkIdx, int chrom,
			     int chromMarkerIdx, int geno[2]) {
  if (_ignore)
    return; // no need/nowhere to store genotypes

  for(int h = 0; h < 2; h++)
    _geno[chrom][chromMarkerIdx][h] = geno[h];
}

// Given the parents of <this>, sets pointers to those parents in
// <this>
void PersonNorm::setParents(char *familyid, PersonNorm *parents[2],
			    int numParents, bool &warningPrinted, FILE *log,
			    int *numMendelError, int *numMendelCounted) {
  // This method does not identify non-Mendelian errors, so
  // ensure that the caller isn't requesting these counts
  assert(numMendelError == NULL && numMendelCounted == NULL);

  _parents[0] = parents[0];
  _parents[1] = parents[1];
}

// TODO: comment
void PersonNorm::setXHetToMissing(int *numHets, int *numCalls) {
  // TODO!
}
