// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <limits.h>
#include <stdio.h>
#include <stdint.h>

#ifndef HAPI_UR_UTIL_H
#define HAPI_UR_UTIL_H

// For haplotype data:
typedef  unsigned long int   chunk;

#define BITS_PER_CHUNK          ((unsigned int) (sizeof(chunk) * 8))

// Need to be able to store bit indexes into chunks in a memory-efficient
// manner; typically, BITS_PER_CHUNK will be 64, but surely never more than
// 256; hence 8 bits will be enough.  (We can always change this if computers
// start having native 256 bit values.)
#define BITS_FOR_CHUNK_IDX      8
static_assert( (1 << BITS_FOR_CHUNK_IDX) >= BITS_PER_CHUNK,
	       "Not enough bits to store a chunk index?" );

// All bits set for a chunk:
#define ALL_CHUNK_BITS_SET      ULONG_MAX

// The following is assumed in setBitsToIdx()
static_assert( ((chunk) -1) == ALL_CHUNK_BITS_SET,
	      "(chunk) -1 != every bit set in a chunk?" );

// Returns a chunk with bits 0..<index> set to 1, others 0
inline chunk setBitsToIdx(int index) {
  // The following is the same as (1ul << (index + 1)) - 1 except that
  // when index == 63, it returns ALL_CHUNK_BITS_SET set instead of returning 0.
  // For some reason (1ul << 64) == 1 -- the shift wraps around instead of
  // overflowing.
  return (((1ul << index) - 1) << 1) + 1;
}

// Returns a chunk with bits (LAST_BIT-numBits+1)..LAST_BIT set to 1, others 0
// Note: if numBits == 0, this actually returns all bits set since this is the
// semantics the caller needs.
//inline chunk setLastNumBits(int numBits) {
//  // The following is the same as
//  // setBitsToIdx(numBits) << (BITS_PER_CHUNK - numBits) except that when
//  // numBits == 0, it returns all the bits set instead of returning 1
//  return ((setBitsToIdx(numBits-1) << 1) + 1) << (BITS_PER_CHUNK - numBits);
//}

// Returns a chunk with bits (LAST_BIT-numBits+1)..LAST_BIT set to 1, others 0
// Note: if numBits == 0 or numBits == 64, retuns 0
inline chunk setLastNumBits(int numBits) {
  return (setBitsToIdx(numBits-1) << (BITS_PER_CHUNK - numBits - 1)) << 1;
}


// Returns the bit value (0 or 1) at bit number <bitNum> in <val>
inline int getBit(chunk val, int bitNum) {
  return (val >> bitNum) & 1;
}

uint32_t countBitsSet(chunk value);
uint32_t countBitsSetDense(chunk value); // when value is likely to have many bits
uint32_t getHighestOrderBitIdx(chunk value);
uint32_t getLowestOrderBitIdx(chunk value);
chunk    getHighestOrderBit(chunk value);
chunk    getLowestOrderBit(chunk value);
uint32_t countHighOrderUnsetBits(chunk value);
uint32_t countLowOrderUnsetBits(chunk value);

int  chunkHashFunc(const chunk &key);
bool chunkEqualFunc(const chunk &v1, const chunk &v2);

// prints the portion the haplotype underlying genos & loci, using ? for
// heterozygous loci
inline void printHap(FILE *out, chunk genos, chunk loci) {
  for(uint32_t i = 0; i < BITS_PER_CHUNK; i++) {
    int geno = (genos >> i) & 1;
    int defined = (loci >> i) & 1;
    if (!defined)
      fprintf(out, "?");
    else
      fprintf(out, "%d", geno);
  }
}

inline void swap(chunk &v1, chunk &v2) {
  chunk tmp = v1;
  v1 = v2;
  v2 = tmp;
}

#endif // HAPI_UR_UTIL_H
