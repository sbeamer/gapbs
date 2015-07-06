// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BITMAP_H_
#define BITMAP_H_

#include <algorithm>
#include <cinttypes>

#include "platform_atomics.h"

/*

Parallel bitmap motivated by need to be able to be able to safely set
bits in parallel (unlike std::vector<bool>).

*/


const uint64_t bits_per_word = 64;
uint64_t word_offset(size_t n) { return n/bits_per_word; }
uint64_t bit_offset(size_t n) { return n & (bits_per_word-1); }

class Bitmap {
 public:
  Bitmap(size_t size) {
    uint64_t num_words = (size + bits_per_word - 1) / bits_per_word;
    start_ = new uint64_t[num_words];
    end_ = start_ + num_words;
  }

  ~Bitmap() {
    delete[] start_;
  }

  void reset() {
    std::fill(start_, end_, 0);
  }

  void set_bit(size_t pos) {
    start_[word_offset(pos)] |= ((uint64_t) 1l << bit_offset(pos));
  }

  void set_bit_atomic(size_t pos) {
    uint64_t old_val, new_val;
    do {
      old_val = start_[word_offset(pos)];
      new_val = old_val | ((uint64_t) 1l << bit_offset(pos));
    } while (!compare_and_swap(start_[word_offset(pos)], old_val, new_val));
  }

  bool get_bit(size_t pos) {
    return (start_[word_offset(pos)] >> bit_offset(pos)) & 1l;
  }

  void or_in(const Bitmap &other) {
    uint64_t *p_local = start_;
    uint64_t *p_other = other.start_;
    while (p_local != end_) {
      *p_local |= *p_other;
      p_local++;
      p_other++;
    }
  }

  void swap(Bitmap &other) {
    std::swap(start_, other.start_);
    std::swap(end_, other.end_);
  }

 private:
  uint64_t *start_;
  uint64_t *end_;
};

#endif  // BITMAP_H_
