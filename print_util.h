// See LICENSE.txt for license details.

#ifndef PRINT_UTIL_H_
#define PRINT_UTIL_H_

#include <cinttypes>
#include <stdio.h>
#include <string>

#include "timer.h"


void PrintTime(const std::string &s, double seconds) {
  printf("%-21s%3.5lf\n", (s + ":").c_str(), seconds);
}

void PrintStep(int step, double seconds, int64_t count=-1) {
  if (count != -1)
    printf("%5d%11lld  %10.5lf\n", step, count, seconds);
  else
    printf("%5d%23.5lf\n", step, seconds);
}

void PrintStep(const std::string &s, double seconds, int64_t count=-1) {
  if (count != -1)
    printf("%5s%11lld  %10.5lf\n", s.c_str(), count, seconds);
  else
    printf("%5s%23.5lf\n", s.c_str(), seconds);
}


#define TIME_PRINT(label, op) {    \
  Timer t_;                       \
  t_.Start();                     \
  (op);                           \
  t_.Stop();                      \
  PrintTime(label, t_.Seconds()); \
}


#endif  // PRINT_UTIL_H_
