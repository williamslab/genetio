// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include "personhapbits.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<PersonHapBits *> PersonHapBits::_allIndivs;
std::unordered_map<const char *, PersonHapBits *, HashString, EqualString>
						    PersonHapBits::_idToPerson;

PersonHapBits::PersonHapBits(char *id, char sex, int popIndex, uint32_t sampNum,
			     short familyIdLength) :
		       SuperPerson(id, sex, popIndex, familyIdLength) {
  if (!_ignore) {
    int numHapChunks = Marker::getNumHapChunks();
    assert(numHapChunks > 0);
    for(int h = 0; h < 2; h++) {
      _hap[h] = new chunk[numHapChunks];

      // init:
      for(int i = 0; i < numHapChunks; i++) {
	_hap[h][i] = 0;
      }
    }

    if (_idToPerson.find(_id) != _idToPerson.end()) {
      fprintf(stderr, "\nERROR: multiple individuals with id %s!\n", _id);
      exit(3);
    }

    _idToPerson[_id] = this;
  }
}

PersonHapBits::~PersonHapBits() {
  if (!_ignore) {
    for(int h = 0; h < 2; h++)
      delete [] _hap[h];
  }
}

void PersonHapBits::setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
				int chromMarkerIdx, int geno[2]) {
  if (_ignore)
    return; // no need/nowhere to store genotypes
  assert(geno[0] >= 0); // should not have any missing data
  assert(geno[0] <= 1 && geno[1] <= 1); // biallelic variants

  for(int h = 0; h < 2; h++) {
    if (geno[h] == 0)
      continue; // already init'd to 0, so done

    // set corresponding bit to 1:
    _hap[h][hapChunkNum] += 1ul << chunkIdx;
  }
}

// This class does not support family relationships so attempt to set
// parents prints error and exits
void PersonHapBits::setParents(char *familyid, PersonHapBits *parents[2],
			       int numParents, bool &warningPrinted, FILE *log,
			       int *numMendelError, int *numMendelCounted) {
  FILE *outs[2] = { stdout, log };
  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    if (out == NULL)
      continue;
    fprintf(out, "\n\nERROR: Class that stores haplotype data does not currently support families\n");
  }
  exit(5);
}

// Checks that there are no heterozygous sites on the X chromosome for
// <this> (which is required to be male).
// Gives error and quits if there are heterozygous sites: this class requires
// complete data so setting a value to missing is not allowed.
// Assumes that we only have data for the X chromosome.
void PersonHapBits::setXHetToMissing(FILE *log, int *numHets, int *numCalls) {
  assert(getSex() == 'M');
  FILE *outs[2] = { stdout, log };

  for(int curChunk = 0; curChunk < Marker::getNumHapChunks(); curChunk++) {
    if (_hap[0][curChunk] != _hap[1][curChunk]) {
      for (int o = 0; o < 2; o++) {
	FILE *out = outs[o];
	if (out == NULL)
	  continue;
	fprintf(out, "\n\nERROR: Male %s has heterozygous variants on X\n",
		getId());
	fprintf(out, "       Set these variants to missing and re-phase to correct\n");
	exit(2);
      }
    }
  }
}
