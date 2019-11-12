// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <sparsehash/dense_hash_map>
#include <tr1/hashtable.h>
#include "marker.h"
#include "superperson.h"
#include "personio.h"
#include "dynarray.h"

#ifndef PERSONBULK_H
#define PERSONBULK_H

// For reference, the following are the 2-bit values of all the genotypes
// in PLINK format data:
// 0 - homozygous for allele 0
// 1 - missing
// 2 - heterozygous
// 3 - homozygous for allele 1
enum Geno {
  G_HOM0 = 0,
  G_MISS = 1,
  G_HET  = 2,
  G_HOM1 = 3
};

class PersonBulk : public SuperPerson {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void init();

    static PersonBulk * lookupId(char *id) { return _idToPerson.lookup(id); }

    static void getBulkContainers(uint8_t **&bulk_data, int *&bytesPerMarker) {
      bulk_data = &_data;
      bytesPerMarker = &_bytesPerMarker;
    }

    friend class PersonIO<PersonBulk>;

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    uint8_t getBitGeno(int marker) {
      // Samples that don't have data (i.e., parents who had a PersonBulk entry
      // created for them though no genotype data exists) have _sampNum == -1:
      if (_sampNum == (uint32_t) -1)
	return G_MISS;

      uint64_t bitNum = _sampNum * 2;
      int byteShift = bitNum / 8;
      uint32_t index = (uint32_t) marker * _bytesPerMarker + byteShift;
      // extract two bits at the relevant position in this byte:
      return (_data[index] >> (bitNum % 8)) & 3;
    }

    int getGenotype(int chunkNum, int chunkIdx, int chromIdx,
		    int chromMarkerIdx) {
      // This function is used for outputing Eigenstrat format data. It is
      // trivial to implement (especially if runtime performance isn't a
      // concern) but likely won't be used.
      fprintf(stderr, "ERROR: support for printing Eigenstrat output from PersonBulk stored data not yet implemented\n");
      exit(9);
    }
    int getHapAllele(int homolog, int chunkNum, int chunkIdx, int chromIdx,
		     int chromMarkerIdx) {
      // TODO: need to implement this once phasing code is completed
      assert(false);
    }

    // Note: these next two methods are currently only used by PersonIO:
    // isUnrelated() really doesn't apply for PersonBulk, which is
    // intended to be used for family-based phasing; samples are not
    // from a population
    bool isUnrelated() { return false; }
    bool isPhased() { return false; }

    bool hasData() { return _sampNum != (uint32_t) -1; }

    //////////////////////////////////////////////////////////////////
    // public static variables
    //////////////////////////////////////////////////////////////////

    static dynarray<PersonBulk *> _allIndivs;
    // Hash from PersonBulk ids to PersonBulk *
    static Hashtable<char *, PersonBulk *> _idToPerson;

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    // Nothing to clean up after parsing all the persons; the one
    // parameter to this function applies to a different Person class
    static void cleanUpPostParse(bool dontcare) { }

    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // Stores all the genotype data in PLINK bed format: two bits per genotype
    // in SNP-major mode. The PersonBulk objects store the location of each
    // genotype.
    static uint8_t *_data;

    // How many bytes for each marker are contained in <data>? Need this to
    // get data for a specific marker: data[ marker * byptesPerMarker ] points
    // to the beginning of this data in memory
    static int _bytesPerMarker;

    //////////////////////////////////////////////////////////////////
    // private methods
    //////////////////////////////////////////////////////////////////

    PersonBulk(char *id, char sex, int popIndex, uint32_t sampNum,
	       short familyIdLength = 0);
    ~PersonBulk();

    void setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
		     int chromMarkerIdx, int geno[2]) {
      fprintf(stderr, "ERROR: PersonBulk objects store genotype data in bulk; cannot set specific genotype values\n");
      exit(9);
    }
    void setParents(char *familyid, PersonBulk *parents[2],
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

#endif // PERSONBULK_H
