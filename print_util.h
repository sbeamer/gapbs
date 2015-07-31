// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef PRINT_UTIL_H_
#define PRINT_UTIL_H_

#include <stdio.h>
#include <cinttypes>
#include <string>

#include "timer.h"


/*
GAP Benchmark Suite
Author: Scott Beamer

Useful output functions to ensure formmating consistency
*/


void PrintTime(const std::string &s, double seconds) {
  printf("%-21s%3.5lf\n", (s + ":").c_str(), seconds);
}

void PrintStep(int step, double seconds, int64_t count = -1) {
  if (count != -1)
    printf("%5d%11" PRId64 " %10.5lf\n", step, count, seconds);
  else
    printf("%5d%23.5lf\n", step, seconds);
}

void PrintStep(const std::string &s, double seconds, int64_t count = -1) {
  if (count != -1)
    printf("%5s%11" PRId64 "  %10.5lf\n", s.c_str(), count, seconds);
  else
    printf("%5s%23.5lf\n", s.c_str(), seconds);
}

// Runs op and prints the time it took to execute labelled by label
#define TIME_PRINT(label, op) {   \
  Timer t_;                       \
  t_.Start();                     \
  (op);                           \
  t_.Stop();                      \
  PrintTime(label, t_.Seconds()); \
}

#endif  // PRINT_UTIL_H_
