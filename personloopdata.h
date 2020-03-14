// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <unordered_map>
#include "marker.h"
#include "superperson.h"
#include "personio.h"
#include "dynarray.h"
#include "util.h"

#ifndef PERSONLOOPDATA_H
#define PERSONLOOPDATA_H

class PersonLoopData : public SuperPerson {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static PersonLoopData * lookupId(char *id) {
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

    friend class PersonIO<PersonLoopData>;

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    // Note: these next several methods are currently only used by PersonIO:
    // isUnrelated() really doesn't apply for PersonLoopData, which is
    // intended to be used by IBIS
    int getGenotype(int chunkNum, int chunkIdx, int chromIdx,
		    int chromMarkerIdx) {
      fprintf(stderr, "ERROR: PersonLoopData objects do not store genotype data; cannot get genotype values\n");
      exit(9);
    }
    int getHapAllele(int homolog, int chunkNum, int chunkIdx,
		     int chromIdx, int chromMarkerIdx) {
      fprintf(stderr, "ERROR: PersonLoopData objects do not store genotype data; cannot get haplotype values\n");
      exit(9);
    }
    bool isUnrelated() { return false; }
    bool isPhased() { return false; }

    //////////////////////////////////////////////////////////////////
    // public static variables
    //////////////////////////////////////////////////////////////////

    static dynarray<PersonLoopData *> _allIndivs;
    // Hash from PersonLoopData ids to PersonLoopData *
    static std::unordered_map<const char *, PersonLoopData *,
			      HashString, EqualString> _idToPerson;

  private:
    //////////////////////////////////////////////////////////////////
    // private methods
    //////////////////////////////////////////////////////////////////

    void setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
		     int chromMarkerIdx, int geno[2]) {
      fprintf(stderr, "ERROR: PersonLoopData objects do not store genotype data; cannot set genotype values\n");
      exit(9);
    }

    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    // Nothing to clean up after parsing all the persons; the one
    // parameter to this function applies to a different Person class
    static void cleanUpPostParse(bool dontcare) { }

    //////////////////////////////////////////////////////////////////
    // private methods
    //////////////////////////////////////////////////////////////////

    PersonLoopData(char *id, char sex, int popIndex, uint32_t sampNum,
		     short familyIdLength = 0);
    ~PersonLoopData();

    void setParents(char *familyid, PersonLoopData *parents[2],
		    int numParents, bool &warningPrinted, FILE *log,
		    int *numMendelError, int *numMendelCounted);
    // Used for males on the X chromosome: heterozygous sites are erroneous
    void setXHetToMissing(FILE *log, int *numHets = NULL, int *numCalls = NULL);

    //////////////////////////////////////////////////////////////////
    // private variables
    //////////////////////////////////////////////////////////////////

    // Stores the sample number of <this> (this is the 0-indexed row number of
    // the individual in the .fam/.ind file).
    // This is used to index into the marker data: starting from the first
    // genotype for a given marker, the data for <this> is stored 2*_sampNum
    // bits from that location.
    uint32_t _sampNum;
};

#endif // PERSONLOOPDATA_H
