// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
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
			 const char *indFile, const char *onlyChr,
			 int startPos, int endPos, const char *XchrName,
			 int noFamilyId, bool vcfInput,
			 int printTrioKids = false, FILE *log = NULL,
			 bool phased = false, int **numMendelError = NULL,
			 int **numMendelCounted = NULL);
    static void readVCF(const char *vcfFile, const char *onlyChr, int startPos,
			int endPos, const char *XcharName, FILE *log = NULL);

    static void readIndivs(FILE *in, FILE *log, bool phased);
    static bool readPedOrFamFile(FILE *in, bool omitFamilyId,
				 bool knowIsFam = false);
    static void makePersonsFromIds(char **ids, uint32_t numIds);
    static void parsePedGenotypes(FILE *in, P *thePerson);
    static void findRelationships(FILE *in, FILE *log, bool omitFamilyId,
				  int *numMendelError = NULL,
				  int *numMendelCounted = NULL);
    static void removeIgnoreIndivs();

    static void printEigenstratGeno(FILE *out);
    static void printEigenstratPhased(FILE *out, int numSamples = -1);
    static void printGzEigenstratPhased(gzFile out);
    static void printPed(FILE *out);
    static void printPhasedIndFile(FILE *out, bool trioDuoOnly = false);
    static void printImpute2Haps(FILE *out);
    static void printGzImpute2Haps(gzFile out);
    static void printImpute2SampleFile(FILE *out, bool trioDuoOnly = false);

    static void parsePackedAncestryMapFormat(FILE *in);
    static void parseEigenstratFormat(FILE *in, bool phased);
    static void parsePlinkBedFormat(FILE *in);
    static void parseVCFGenotypes(htsFile *vcfIn, tbx_t *index, hts_itr_t *itr,
				  const char *vcfFile, FILE *outs[2]);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void parsePackedGenotypes(FILE *in, int recordLen, char *buf,
				     int numIndivs, int type);
};

#endif // PERSONIO_H
