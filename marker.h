// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdint.h>
#include <zlib.h>
#include <string>
#include "hapi-ur-util.h"
#include "dynarray.h"

#ifdef VCF
#include <htslib/hts.h>
#include <htslib/tbx.h>
#endif

#ifndef MARKER_H
#define MARKER_H


// Determines the size of the _alleles character string.  Want this to be
// somewhat sizeable so that we can store indel character strings.  A value of
// 41 below makes sizeof(Marker) == 72; must be two-word aligned, so any larger
// makes this 80 bytes.
#define MAX_ALLELES	   41

class Marker {
  public:
//    ~Marker() { // not needed -- only delete when program done: OS will manage
//      if (_numAlleles)
//	delete [] _alleles;
//      delete [] _name;
//    }

    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void setReadOnlyOneChrom() { _readOnlyOneChrom = true; }
    static bool getReadOnlyOneChrom() { return _readOnlyOneChrom; }
    static void readSNPFile(const char *snpFile, const char *onlyChr,
			    int startPos, int endPos);
    static void readMapFile(const char *mapFile, const char *onlyChr,
			    int startPos, int endPos);
    static void readBIMFile(const char *bimFile, const char *onlyChr,
			    int startPos, int endPos);
#ifdef VCF
    static void readVCFFile(htsFile *vcfIn, tbx_t *index, hts_itr_t *itr,
			    int startPos, int endPos);
#endif
    static void printSNPFile(FILE *out);
    static void printMapFile(FILE *out);
    static void printImpute2Prefix(FILE *out, int markerNum);
    static void printGzImpute2Prefix(gzFile out, int markerNum);

    static int getNumChroms()  { return _chromNames.length(); }
    static int getNumMarkers() { return _allMarkers.length(); }
    static int getNumMarkersInFile() { return _numMarkersInFile; }
    static int getFirstStoredMarkerFileIdx() { return _firstStoredMarkerIdx; }
    static int getFirstMarkerNum(int chromIdx)
					  { return _firstMarkerNum[chromIdx]; }
    static int getLastMarkerNum(int chromIdx)
					  { return _lastMarkerNum[chromIdx]; }
//    static const Marker * getFirstMarker(int chrom)
//			    { return getMarker( getFirstMarkerNum(chrom) ); }
    static int getNumChromMarkers(int chromIdx)
	    { return _lastMarkerNum[chromIdx] - _firstMarkerNum[chromIdx] + 1; }

    static int getNumHapChunks() { return _numHapChunks; }
    static int getFirstHapChunk(int chromIdx) {return _firstHapChunk[chromIdx];}
    static int getLastHapChunk(int chromIdx)  {return _lastHapChunk[chromIdx];}
    static int getNumHapChunksFor(int numMarkers);
    // Returns the number of markers not divisible by the num of bits in a chunk
    static int getChunkModMarkers(int numMarkers);

    static int getFirstMarkerNumForChunk(int chromIdx, int chunkNum);

    static const char * getChromName(int chrIdx) {return _chromNames[ chrIdx ];}
    static const Marker * getMarker(int num) { return _allMarkers[num]; }

    // cheating in order to set allele frequencies and/or add alleles for .ped:
    static Marker * getMarkerNonConst(int num) { return _allMarkers[num]; }

    static void convertMapTocM() {
      int numMarkers = getNumMarkers();
      for(int m = 0; m < numMarkers; m++) {
	_allMarkers[m]->_mapPos *= 100;
      }
    }

    static int    getNumWindows() { return _hapWindowEnds.length(); }
    static int    getWindowEndMarker(int wind) { return _hapWindowEnds[wind]; }
    static float  getWindowMapCenter(int wind)
					  { return _hapWindowMapCenter[wind]; }
    static void   updateWindows(int initOffset, int windowNumMarkers);
    static void   updateWindowsMap(int initOffset, float windowLengthMorgans,
				   int minNumMarkers);
//    static uint32_t getTotalPhysLength(bool analyzeChrX);
//    static float    getTotalGenetLength(bool analyzeChrX);
    
    static const dynarray<int> & getMarkersToOmit() { return _omitMarkers; }

    // not needed -- only delete when program done: OS will manage
//    static void cleanUp() {
//      int len = _chromNames.length();
//      for(int c = 0; c < len; c++) {
//	delete [] _chromNames[c];
//      }
//
//      len = _allMarkers.length();
//      for(int m = 0; m < len; m++) {
//	delete _allMarkers[m];
//      }
//    }

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    const char * getName() const  { return _name; }
    int   getChromIdx() const     { return _chromIdx; }
    const char * getChromName() const { return _chromNames[ _chromIdx ]; }
    float getMapPos() const       { return _mapPos; }
    int   getPhysPos() const      { return _physPos; }
    const char * getAlleleStr() const  { return _alleles; }
    short getNumAlleles() const   { return _numAlleles; }
    // Note: these values do not get initialized for .ped files for which SNPs
    // are read in an awkward order for this calculation to be done while
    // reading
    float getLogAlleleFreq() const      { return _logAlleleFreq; }
    float getLogVarAlleleFreq() const   { return _logVarAlleleFreq; }
    float getNumMarkersInWindow() const { return _numSNPsWindow; }
    void setAlleleFreq(int alleleCount, int totalGenoWithData,
		       bool nonStandardGenotype = false);


  private:
    Marker(const char *markerName, int chromIdx, float mapPos,
	   float morganDistToPrev, int physPos, const char *alleles,
	   int numAlleles);

    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void readMarkers(FILE *in, const char *onlyChr, int type,
			    int startPos, int endPos);
    static void updateInfoPrevChrom(int prevChromIdx, int numMarkersPrevChrom);
    static void setNumMarkersInWindow(int startMarkerNum, int numMarkers);
    static int  skipWhitespace(char *curBuf, size_t &bind, size_t &nread,
			       const size_t BUF_SIZE);
    static void replaceBuffer(FILE *in, char *&curBuf, char *&nextBuf,
			      size_t &bind, size_t &nread,
			      const size_t BUF_SIZE);
    static int readDoubleBuffer(FILE *in, char *&field, char *&curBuf,
				char *&nextBuf, size_t &bind, size_t &nread,
				const size_t BUF_SIZE);

    // marker name (usually SNP rs id)
    char *_name;

    // physical position:
    int _physPos;

    // genetic position (map distance):
    float _mapPos;

    // Alleles; stored as a single string with a space separating the types.
    char *_alleles;
    uint8_t _numAlleles;

    // chromosome
    // TODO: add check to ensure uint8_t is big enough
    uint8_t _chromIdx;

    // Number of SNPs in this window (changes as HAPI-UR program runs)
    short _numSNPsWindow;

    // Allele frequencies for 0 allele (_logAlleleFreq) and 1 allele
    // (_logVarAlleleFreq)
    // TODO: probably should remove these to save space.  Only need these
    // in rare instances
    float _logAlleleFreq;
    float _logVarAlleleFreq;

    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // List of all the markers read in:
    static dynarray<Marker *> _allMarkers;

    // List of marker indexes for markers that should be omitted.  These markers
    // have either physical or genetic positions of 0 and therefore cannot be
    // placed with confidence.
    static dynarray<int> _omitMarkers;

    // Should we read only one chromosome?  By default, no, but
    // setReadOnlyOneChrom() can change this.
    static bool _readOnlyOneChrom;

    // For when we're only processing one chromosome (i.e., when
    // CmdLineOpts::onlyChr != 0), this value indicates the number of markers
    // that should be skipped when reading the genotype data in order to get
    // to the first marker for the chromosome being processed
    // default value is -1 which means that all markers are read
    static int _firstStoredMarkerIdx;

    // The number of markers that were present in the marker file (for
    // checking the reported values in a packed ancestrymap format file):
    static int _numMarkersInFile;

    // Stores the chromosome/contig names for all that have been read in
    static dynarray<char *> _chromNames;

    // Stores the first marker number on the corresponding chromosome number
    // 1..23, 0 if no markers on chrom
    static dynarray<int> _firstMarkerNum;

    // Stores the last marker number on the corresponding chromosome number
    // 1..23, 0 if no markers on chrom
    static dynarray<int> _lastMarkerNum;

    // Stores the starting haplotype chunk number for each chromosome, 0 if
    // no markers on chrom
    static dynarray<int> _firstHapChunk;

    // Stores the last haplotype chunk number for each chromosome, 0 if no
    // markers on chrom
    static dynarray<int> _lastHapChunk;

    // Number of haplotype chunks to store the data.  We store the haplotypes as
    // ulint type variables, with one locus per bit, so this number is the
    // number of markers divided by the number of bits per ulint type:
    static int _numHapChunks;

    // List of marker numbers for haplotype window breakpoints
    static dynarray<int> _hapWindowEnds;

    // The Morgan distance to the center of this window from the start of the
    // chromosome
    static dynarray<float> _hapWindowMapCenter;
};


#endif // MARKER_H
