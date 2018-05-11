// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef TIMER_H_
#define TIMER_H_


/*
GAP Benchmark Suite
Class:  Timer
Authors: Scott Beamer, Michael Sutton

Simple timer that wraps gettimeofday or uses chrono for Windows
*/

#if defined _WIN32

#include <chrono>

class Timer {
 public:
  Timer() {}

  void Start() {
      t2 = t1 = std::chrono::high_resolution_clock::now();
  }

  void Stop() {
      t2 = std::chrono::high_resolution_clock::now();
  }

  double Seconds() const {
    return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
  }

  double Millisecs() const {
    return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0;
  }

  double Microsecs() const {
    return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1.0;
  }

 private:
  std::chrono::high_resolution_clock::time_point t1, t2;
};

#else

#include <sys/time.h>

class Timer {
 public:
  Timer() {}

  void Start() {
    gettimeofday(&start_time_, NULL);
  }

  void Stop() {
    gettimeofday(&elapsed_time_, NULL);
    elapsed_time_.tv_sec  -= start_time_.tv_sec;
    elapsed_time_.tv_usec -= start_time_.tv_usec;
  }

  double Seconds() const {
    return elapsed_time_.tv_sec + elapsed_time_.tv_usec/1e6;
  }

  double Millisecs() const {
    return 1000*elapsed_time_.tv_sec + elapsed_time_.tv_usec/1000;
  }

  double Microsecs() const {
    return 1e6*elapsed_time_.tv_sec + elapsed_time_.tv_usec;
  }

 private:
  struct timeval start_time_;
  struct timeval elapsed_time_;
};

#endif

// Times op's execution using the timer t
#define TIME_OP(t, op) { t.Start(); (op); t.Stop(); }

#endif  // TIMER_H_
