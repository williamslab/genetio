// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "dynarray.h"
#include "hashtable.h"

#ifndef PERSONIO_H
#define PERSONIO_H

template <class P>
class PersonIO {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void readData(const char *genoFile, const char *markerFile,
			 const char *indFile, int onlyChr,
			 int startPos, int endPos, int analyzeChrX,
			 int noFamilyId, int printTrioKids,
			 bool printGenetLength = false, FILE *log = NULL);

    static void readIndivs(FILE *in);
    static bool readPedOrFamFile(FILE *in, bool omitFamilyId,
				 bool knowIsFam = false);
    static void parsePedGenotypes(FILE *in, P *thePerson);
    static void findRelationships(FILE *in, FILE *log, bool omitFamilyId,
				  int *numMendelError = NULL,
				  int *numMendelCounted = NULL);
    static void removeIgnoreIndivs();
    static void printEigenstratGeno(FILE *out);
    static void printEigenstratPhased(FILE *out, int numSamples = -1);
    static void printGzEigenstratPhased(gzFile out);
    static void printPhasedIndFile(FILE *out, bool trioDuoOnly = false);
    static void printImpute2Haps(FILE *out);
    static void printGzImpute2Haps(gzFile out);
    static void printImpute2SampleFile(FILE *out, bool trioDuoOnly = false);
    static void parsePackedAncestryMapFormat(FILE *in);
    static void parseEigenstratFormat(FILE *in);
    static void parsePlinkBedFormat(FILE *in);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void parsePackedGenotypes(FILE *in, int recordLen, char *buf,
				     int numIndivs, int type);
};

#endif // PERSONIO_H
