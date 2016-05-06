// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include "dynarray.h"

#ifndef SUPERPERSON_H
#define SUPERPERSON_H

#define MAX_PERSON_ID	126

// Simple super class for Person* objects
class SuperPerson {
  public:
    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    SuperPerson(char *id, char sex, int popIndex, short familyIdLength);
    virtual ~SuperPerson();

    void setIgnore() {
      _ignore = true;
      if (_sex == 'M')
	_numMales--;
      else if (_sex == 'U')
	_numSexUnknown--;
    }

    const char *getId() { return _id; }
    const short getFamilyIdLength() { return _familyIdLength; }
    const char *getPopLabel() { return _popLabels[_popIndex];}
    char getSex() { return _sex; }
    int getPopIndex() { return _popIndex; }
    virtual int getGenotype(int chunkNum, int chunkIdx, int chromIdx,
			    int chromMarkerIdx) = 0;
    virtual int getHapAllele(int homolog, int chunkNum, int chunkIdx,
			     int chromIdx, int chromMarkerIdx) = 0;

    bool isIgnore() { return _ignore; }
    virtual bool isPhased() = 0;
    virtual bool isUnrelated() = 0;

    //////////////////////////////////////////////////////////////////
    // public static variables
    //////////////////////////////////////////////////////////////////

    static dynarray<const char *> _popLabels;
    static int _maxPersonIdLength;
    // need these two numbers for X chromosome phasing/checking
    static int _numMales, _numSexUnknown;

  protected:
    //////////////////////////////////////////////////////////////////
    // protected methods
    //////////////////////////////////////////////////////////////////

    virtual void setGenotype(int hapChunkNum, int chunkIdx, int chromIdx,
			     int chromMarkerIdx, int geno[2]) = 0;
    virtual void setXHetToMissing(FILE *log, int *numHets = NULL,
				  int *numCalls = NULL) = 0;

    //////////////////////////////////////////////////////////////////
    // protected variables
    //////////////////////////////////////////////////////////////////

    // String id for person
    char _id[MAX_PERSON_ID];

    // Sex for person -- 'M', 'F', or 'U'
    char _sex;

    // True if the individual should be ignored.  To be ignored,
    // PersonIO<PersonBits>::removeIgnoreIndivs() should be called on the list
    // of indivs
    bool _ignore;

    // Length of the family portion of the id (the IMPUTE2 file format separates
    // the family and individual ids, so we need to be able to print the part of
    // the string corresponding to the family id separate from the individual)
    short _familyIdLength;

    // The index number of the population for this PersonBits
    short _popIndex;
};

#endif // SUPERPERSON_H
