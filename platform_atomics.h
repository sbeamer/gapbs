// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef PLATFORM_ATOMICS_H_
#define PLATFORM_ATOMICS_H_


/*
GAP Benchmark Suite
File:   Platform Atomics
Author: Scott Beamer

Wrappers for compiler intrinsics for atomic memory operations (AMOs)
 - If not using OpenMP (serial), provides serial fallbacks
*/


#if defined _OPENMP

  #if defined __GNUC__

    // gcc/clang/icc instrinsics

    template<typename T, typename U>
    T fetch_and_add(T &x, U inc) {
      return __sync_fetch_and_add(&x, inc);
    }

    template<typename T, typename U, typename V>
    bool compare_and_swap(T &x, U old_val, V new_val) {
      return __sync_bool_compare_and_swap(&x, old_val, new_val);
    }

  #else   // defined __GNUC__

    #error No atomics available for this compiler but using OpenMP

  #endif  // else defined __GNUC__

#else   // defined _OPENMP

  // serial fallbacks

  template<typename T, typename U>
  T fetch_and_add(T &x, U inc) {
    T orig_val = x;
    x += inc;
    return orig_val;
  }

  template<typename T, typename U, typename V>
  bool compare_and_swap(T &x, U old_val, V new_val) {
    if (x == old_val) {
      x = new_val;
      return true;
    }
    return false;
  }

#endif  // else defined _OPENMP

#endif  // PLATFORM_ATOMICS_H_
