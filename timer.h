// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#ifndef TIMER_H
#define TIMER_H

#define CLOCK_REALTIME 0 
#define CLOCK_MONOTONIC 0

#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

#ifdef __MACH__
//clock_gettime is not implemented on OSX
int clock_gettime(int /*clk_id*/, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

class Timer {
  public:
    inline Timer();
    inline void reset();
    inline double getElapsedSec();
    inline void printElapsedTime(FILE *out);

  private:
    struct timespec _clk_startTime; // wall clock time using clock_gettime()
//    struct timeval _startTime;      // wall clock time using gettimeofday()
//    clock_t _startTicks;            // cycle time
};

inline Timer::Timer() {
  clock_gettime(CLOCK_REALTIME, &_clk_startTime);
//  gettimeofday(&_startTime, NULL);  // get wall clock time
//  _startTicks = clock();	      // get cycles
}

inline void Timer::reset() {
  clock_gettime(CLOCK_REALTIME, &_clk_startTime);
}

inline double Timer::getElapsedSec() {
  // using clock_gettime() -- purportedly more precise than getimeofday()
  struct timespec curTime;
  int r = clock_gettime(CLOCK_REALTIME, &curTime);
  assert(r == 0);

  double wallDelta = (curTime.tv_sec - _clk_startTime.tv_sec) +
		     // 1 000 000 000 nanoseconds in a second
		     (curTime.tv_nsec - _clk_startTime.tv_nsec) / 1000000000.0;
  return wallDelta;
}

// Prints the total elapsed wall clock and cycle time from the creation of
// this Timer object.
void Timer::printElapsedTime(FILE *out) {
  double wallDelta = getElapsedSec();

  fprintf(out, "Total wall clock time:\t%.3lf seconds.\n", wallDelta);

  // using gettimeofday():
//  struct timeval curTime;
//  gettimeofday(&curTime, NULL);
//
//  double wallDelta = (curTime.tv_sec - _startTime.tv_sec) +
//		     // 1 000 000 microseconds in a second
//		     (curTime.tv_usec - _startTime.tv_usec) / 1000000.0;
//
//  fprintf(out, "Total wall clock time:\t%.3lf seconds.\n", wallDelta);


  // using clock() -- cycle time
//  clock_t curTicks = clock();
//  double cycleDelta = (double) (curTicks - _startTicks) / CLOCKS_PER_SEC;
//
//  fprintf(out, "Total cycle time:\t%.3lf seconds.\n", cycleDelta);
}

#endif // TIMER_H
