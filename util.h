// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <random>
#include <bitset>
#include <functional>
#include <cstring>

#ifndef UTIL_H
#define UTIL_H

#define  FILENAME_LEN   2048

template<class A, class B>
struct Pair {
  Pair(A theA, B theB) { a = theA; b = theB; }
  Pair() { }
  A a;
  B b;
};

template<class A>
struct PairIdx {
  PairIdx(A one, A two) { a[0] = one; a[1] = two; }
  PairIdx() { }
  A &operator[] (int i) { return a[i]; }
  const A &operator[] (int i) const { return a[i]; }
  A a[2];
};


// Returns the smaller of <a> and <b>
inline int min(int a, int b) {
  return (a < b) ? a : b;
}

// Returns the larger of <a> and <b>
inline int max(int a, int b) {
  return (a > b) ? a : b;
}

inline double sumLogLikelihood(double a, double b) {
  if (a > b) {
    return a + log1p(exp(b - a));
  }
  else {
    return b + log1p(exp(a - b));
  }
}

inline double minusLogLikelihood(double a, double b) {
  if (a > b) {
    return a + log(1 - exp(b - a));
  }
  else {
    return b + log(1 - exp(a - b));
  }
}


int  intHashFunc(const int &key);
bool intEqualFunc(const int &v1, const int &v2);
int  pairIntHashFunc(const PairIdx<int> &key);
bool pairIntEqualFunc(const PairIdx<int> &v1, const PairIdx<int> &v2);
bool stringcmp(char * const &s1, char * const &s2);

struct LesserString {
  bool operator() (const char *lhs, const char *rhs) const {
    return strcmp(lhs, rhs) < 0;
  }
};

struct HashString {
  size_t operator() (const char *s) const {
    // got the hash from stack overflow and changed macros to constants
    // https://stackoverflow.com/questions/8317508/hash-function-for-a-string
    const int A = 54059; /* a prime */
    const int B = 76963; /* another prime */
    //const int C = 86969; /* yet another prime */
    const int FIRSTH = 37; /* also prime */

    size_t h = FIRSTH;
    while (*s) {
      h = (h * A) ^ (s[0] * B);
      s++;
    }

    return h; // or return h % C;
  }
};

struct EqualString {
  bool operator() (const char *lhs, const char *rhs) const {
    return !strcmp(lhs, rhs);
  }
};

inline int getNumDigits(int val) {
  int ct;
  for(ct = 0; val > 0; ct++) {
    val /= 10;
  }
  return ct;
}


// Google Unit of Least Precision or Unit in Last Place for a possible algorithm
// (available in Java) for determining what the epsilon value should be.
#define EPSILON         1e-9

inline bool doubleEq(double a, double b) {
  return (fabs(a - b) < EPSILON);
}

inline bool doubleGt(double a, double b, bool orEq = false) {
  if (fabs(a - b) < EPSILON)
    return orEq;
  return a > b;
}

inline bool doubleLt(double a, double b, bool orEq = false) {
  if (fabs(a - b) < EPSILON)
    return orEq;
  return a < b;
}

inline void swap(int &v1, int &v2) {
  int tmp = v1;
  v1 = v2;
  v2 = tmp;
}

struct RandGen {
  static void seed(bool autoSrand, std::mt19937::result_type &randSeed);
  static double real01() { return dis01(v); }

  static std::mt19937 v; // variable to generate random numbers
  static std::uniform_real_distribution<> dis01;
};

// Attempts to open <filename> for reading
// If successful, returns a FILE * to the input stream
// On failure, prints an error message to the streams <outs> and exits
inline FILE *openRead(const char *filename, const char *fileDescriptor,
		      FILE *outs[2]) {
  FILE *in = fopen(filename, "r");
  if (!in) {
    for (int o = 0; o < 2; o++) {
      FILE *out = outs[o];
      if (out == NULL)
	continue;
      fprintf(out, "\n\nERROR: Couldn't open %s file %s\n", fileDescriptor,
	      filename);
    }
    perror(filename);
    exit(2);
  }
  
  return in;
}

// Prints the same string to two output streams in <outs>
inline void mult_printf(FILE *outs[2], const char *msg) {
  for (int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    if (out == NULL)
      continue;
    fprintf(out, "%s", msg);
  }
}

inline size_t popcount(uint64_t val) {
    // Currently using std::bitset, but the conversion is probably not free, so
    // optimize? TODO
    std::bitset<64> to_count(val);
    return to_count.count();
}

#endif // UTIL_H
