// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUCKET_H_
#define BUCKET_H_

#include <algorithm>
#include <cinttypes>
#include <iterator>
#include <vector>


/*
GAP Benchmark Suite
Class:  Bucket
Author: Scott Beamer

Parallel container designed for constant time appends
 - Threads should fill thread-local std::vector first
 - When done with vector, call swap_vector_in(vector)
 - Once started reading with iterator, should not modify or append anymore
 - Like other iterators, comparing iterators for different objects is undefined
 - Internally, built as a vector of vectors but appears contiguous by iterators
*/


template <typename T_>
class Bucket {
 public:
  size_t size() {
    return num_elements_;
  }

  void clear() {
    chunks_.resize(0);
    num_elements_ = 0;
  }

  bool empty() {
    return size() == 0;
  }

  void push_back(T_ to_add) {
    if (chunks_.empty())
      chunks_.emplace_back();
    chunks_.back().push_back(to_add);
    num_elements_++;
  }

  void swap_vector_in(std::vector<T_> &v) {
    if (!v.empty()) {
      #pragma omp critical
      {
        num_elements_ += v.size();
        chunks_.emplace_back();
        chunks_.back().swap(v);
      }
    }
  }

  void swap(Bucket<T_> &other) {
    chunks_.swap(other.chunks_);
    std::swap(num_elements_, other.num_elements_);
  }


  // Doesn't define every operator, but more than enough for OpenMP
  class iterator : public std::iterator<std::random_access_iterator_tag, T_> {
   public:
    iterator(size_t index, size_t offset, std::vector<std::vector<T_>> &chunks)
        : chunk_index_(index), chunk_offset_(offset), chunks_ref_(chunks) {}

    T_ operator*() const {
      return chunks_ref_[chunk_index_][chunk_offset_];
    }

    T_& operator*() {
      return chunks_ref_[chunk_index_][chunk_offset_];
    }

    const iterator &operator++() {
      chunk_offset_++;
      if (chunk_offset_ == chunks_ref_[chunk_index_].size()) {
        chunk_offset_ = 0;
        chunk_index_++;
      }
      return *this;
    }

    iterator operator++(int) {
      iterator copy(*this);
      chunk_offset_++;
      if (chunk_offset_ == chunks_ref_[chunk_index_].size()) {
        chunk_offset_ = 0;
        chunk_index_++;
      }
      return copy;
    }

    iterator & operator +=(int64_t to_add) {
      if (to_add > 0) {
        while ((to_add != 0) && (chunk_index_ < chunks_ref_.size())) {
          chunk_offset_ += to_add;
          if (chunk_offset_ >= chunks_ref_[chunk_index_].size()) {
            to_add = chunk_offset_ - chunks_ref_[chunk_index_].size();
            chunk_offset_ = 0;
            chunk_index_++;
          } else {
            to_add = 0;  // success
          }
        }
      } else {
        while ((to_add != 0) && (chunk_index_ >= 0)) {
          chunk_offset_ += to_add;
          if (chunk_offset_ < 0) {
            chunk_index_--;
            to_add = chunk_offset_;
            chunk_offset_ = chunks_ref_[chunk_index_].size() - 1;
          } else {
            to_add = 0;  // success
          }
        }
      }
      return *this;
    }

    int64_t operator -(const iterator &other) const {
      if (chunk_index_ == other.chunk_index_)
        return chunk_offset_ - other.chunk_offset_;
      size_t op_index = chunk_index_;
      int64_t total = 0;
      if (op_index > other.chunk_index_) {
        total += chunk_offset_;
        while (op_index > other.chunk_index_) {
          op_index--;
          total += chunks_ref_[op_index].size();
        }
        return total - other.chunk_offset_;
      } else {
        total += chunks_ref_[op_index].size() - chunk_offset_;
        op_index++;
        while (op_index < other.chunk_index_) {
          total += chunks_ref_[op_index].size();
          op_index++;
        }
        return total + other.chunk_offset_;
      }
    }

    bool operator <(const iterator &other) const {
      if (chunk_index_ == other.chunk_index_)
        return chunk_offset_ < other.chunk_offset_;
      else
        return chunk_index_ < other.chunk_index_;
    }

    bool operator==(const iterator &other) const {
      return (chunk_index_ == other.chunk_index_) &&
             (chunk_offset_ == other.chunk_offset_) &&
             (chunks_ref_ == other.chunks_ref_);
    }

    bool operator!=(const iterator &other) const {
      return !(operator==(other));
    }

   private:
    size_t chunk_index_;
    size_t chunk_offset_;
    std::vector<std::vector<T_>> &chunks_ref_;
  };


  iterator begin() { return iterator(0, 0, chunks_); }
  iterator end()   { return iterator(chunks_.size(), 0, chunks_); }


 private:
  std::vector<std::vector<T_>> chunks_;
  size_t num_elements_ = 0;
};

#endif  // BUCKET_H_
