// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
std::mt19937 RandGen::v;
std::uniform_real_distribution<> RandGen::dis01(0,1);

int intHashFunc(const int &key) {
  // Take each 2 byte portion of the chunk and combine them:
  int twoByteMask = (1 << 16) - 1;
  int total = 0;

  for(uint16_t i = 0; i < (sizeof(int) / 2); i++) {
    // add in get the ith set of two bytes
    total += (key >> (16 * i)) & twoByteMask;
    total *= 3;
  }

  return total;
}

bool intEqualFunc(const int &v1, const int &v2) {
  return v1 - v2 == 0;
}

int pairIntHashFunc(const PairIdx<int> &key) {
  return 7 * key[0] + 13 * key[1];
}

bool pairIntEqualFunc(const PairIdx<int> &v1, const PairIdx<int> &v2) {
  return v1[0] == v2[0] && v1[1] == v2[1];
}

int stringHash(char * const &key) {
  int sum = 0;
  int strlength = strlen(key);
  for(int i = 0; i < strlength; i++) {
    sum *= 7;
    sum += key[i];
  }
  return sum;
}

bool stringcmp(char * const &s1, char * const &s2) {
    return strcmp(s1, s2) == 0;
}

void RandGen::seed(bool autoSrand, std::mt19937::result_type &randSeed) {
  if (!autoSrand) {
    v.seed(randSeed); // seed with a given value
    return;
  }

  // need to produce our own seed
  randSeed = std::random_device().entropy();

  if (randSeed == 0 && std::random_device().entropy() == 0) {
    // std::random_device is not a real random number generator: fall back on
    // using time to generate a seed:
    timeval tv;
    gettimeofday(&tv, NULL);
    randSeed = tv.tv_sec + tv.tv_usec * 100000;
  }

  v.seed(randSeed); // do the seeding
}
