// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "dynarray.h"
#include "hashtable.h"

#ifdef VCF
#include <htslib/hts.h>
#include <htslib/tbx.h>
#endif

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
			 bool printTrioKids = false, FILE *log = NULL,
			 bool phased = false, int **numMendelError = NULL,
			 int **numMendelCounted = NULL,
			 bool allowEmptyParents = false, bool bulkData = false,
			 bool loopData = false, bool useParents = true);
    static void readData(const char *genoFile, const char *markerFile,
			 const char *indFile, const char *onlyChr,
			 int startPos, int endPos, const char *XchrName,
			 int noFamilyId, FILE *log, bool allowEmptyParents,
			 bool bulkData, bool loopData = false,
			 bool useParents = true);

#ifdef VCF
    static void readVCF(const char *vcfFile, const char *onlyChr, int startPos,
			int endPos, const char *XcharName, FILE *log = NULL);
#endif

    static int  readGenoRow(uint8_t * &data, int bytesPerMarker);
    static void closeGeno();

    static void printEigenstratGeno(FILE *out);
    static void printEigenstratPhased(FILE *out, int numSamples = -1);
    static void printGzEigenstratPhased(gzFile out);
    static void printPed(FILE *out);
    static void printPhasedIndFile(FILE *out, bool trioDuoOnly = false);
    static void printImpute2Haps(FILE *out);
    static void printGzImpute2Haps(gzFile out);
    static void printImpute2SampleFile(FILE *out, bool trioDuoOnly = false);

    // not needed -- only delete when program done: OS will manage
//    static void cleanUp() {
//      int len = P::_allIndivs.length();
//      for(int p = 0; p < len; p++) {
//	delete P::_allIndivs[p];
//      }
//    }

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static int  getGenoFileType(FILE *genoIn, bool phased, FILE *outs[2]);

    static void readIndivs(FILE *in, FILE *log, bool phased);
    static bool readPedOrFamFile(FILE *in, bool omitFamilyId,
				 bool knowIsFam = false);
    static void makePersonsFromIds(char **ids, uint32_t numIds);
    static void parsePedGenotypes(FILE *in, P *thePerson);
    static void findRelationships(FILE *in, FILE *log, bool omitFamilyId,
				  int *numMendelError, int *numMendelCounted,
				  bool createMissingParents);
    static void removeIgnoreIndivs();

    static void parsePackedAncestryMapFormat(FILE *in);
    static void parseEigenstratFormat(FILE *in, bool phased);
    static void parsePlinkBedFormat(FILE *in, FILE *outs[2]);
    static void readPlinkBedBulk(FILE *in, FILE *outs[2]);
    static void checkPlinkHeader(FILE *in, FILE *outs[2]);
#ifdef VCF
    static void parseVCFGenotypes(htsFile *vcfIn, tbx_t *index, hts_itr_t *itr,
				  const char *vcfFile, FILE *outs[2]);
#endif
    static void parsePackedGenotypes(FILE *in, int recordLen, char *buf,
				     int numIndivs, int type);

    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // for use with readGenoRow()
    static FILE *_loopGenoIn;
    static int   _curLoopMarker;
    static int   _curOmitIdx;
};

#endif // PERSONIO_H
