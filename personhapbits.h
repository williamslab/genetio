// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdint.h>
#include <unordered_map>
#include "marker.h"
#include "superperson.h"
#include "personio.h"
#include "util.h"

#ifndef PERSONHAPBITS_H
#define PERSONHAPBITS_H

class PersonHapBits : public SuperPerson {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static PersonHapBits * lookupId(char *id) {
      auto entry = _idToPerson.find(id);
      if (entry == _idToPerson.end())
	return NULL;
      else
	return entry->second;
    }

    static void getBulkContainers(uint8_t **&bulk_data, int *&bytesPerMarker) {
      // This class does not support bulk data storage
      bulk_data = NULL;
      bytesPerMarker = NULL;
    }

    friend class PersonIO<PersonHapBits>;

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    int getGenotype(int chunkNum, int chunkIdx, int chromIdx,
		    int chromMarkerIdx)
		    { return getGenotype(chunkNum, chunkIdx); }
    int getGenotype(int chunkNum, int chunkIdx) {
      int alleleSum = 0;
      for(int h = 0; h < 2; h++) {
	alleleSum += (_hap[h][chunkNum] >> chunkIdx) & 1;
      }

      return alleleSum;
    }
    int getHapAllele(int homolog, int chunkNum, int chunkIdx, int chromIdx,
		     int chromMarkerIdx) {
      return (_hap[homolog][chunkNum] >> chunkIdx) & 1;
    }

    bool isPhased() { return true; }
    bool isUnrelated() { return true; }

    //////////////////////////////////////////////////////////////////
    // public static variables
    //////////////////////////////////////////////////////////////////

    static dynarray<PersonHapBits *> _allIndivs;
    // Hash from PersonHapBits ids to PersonHapBits *
    static std::unordered_map<const char *, PersonHapBits *,
			      HashString, EqualString> _idToPerson;

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    // Nothing to clean up after parsing all the persons; the one
    // parameter to this function applies to a different Person class
    static void cleanUpPostParse(bool keepTrioKids) { }

    //////////////////////////////////////////////////////////////////
    // private methods
    //////////////////////////////////////////////////////////////////

    PersonHapBits(char *id, char sex, int popIndex, uint32_t sampNum,
		  short familyIdLength = 0);
    ~PersonHapBits();

    void setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
		     int chromMarkerIdx, int geno[2]);
    void setParents(char *familyid, PersonHapBits *parents[2],
		    int numParents, bool &warningPrinted, FILE *log,
		    int *numMendelError, int *numMendelCounted);

    // Used for males on the X chromosome: heterozygous sites are erroneous
    void setXHetToMissing(FILE *log, int *numHets = NULL, int *numCalls = NULL);

    //////////////////////////////////////////////////////////////////
    // private variables
    //////////////////////////////////////////////////////////////////

    // An array of bit vectors where each bit corresponds to a marker and the
    // two alleles are encoded as either 0/1
    chunk *_hap[2];
};

#endif // PERSONHAPBITS_H
