// See LICENSE.txt for license details.

#ifndef PVECTOR_H_
#define PVECTOR_H_

#include <algorithm>

#include "platform_atomics.h"

/*

Parallel vector class motivated by need to avoid serial initialization
done by std::vector. When resized, new elements are unitialized.

*/


template <typename T_>
class pvector {
 public:
  typedef T_* iterator;

  pvector() : start_(nullptr), end_size_(nullptr), end_capacity_(nullptr) {}

  pvector(size_t num_elements) {
    start_ = new T_[num_elements];
    end_size_ = start_ + num_elements;
    end_capacity_ = end_size_;
  }

  pvector(size_t num_elements, T_ init_val) : pvector(num_elements) {
    fill(init_val);
  }

  pvector(iterator copy_begin, iterator copy_end)
      : pvector(copy_end - copy_begin) {
    #pragma omp parallel for
    for (size_t i=0; i < capacity(); i++)
      start_[i] = copy_begin[i];
  }

  // don't want this to be copied, too much data to move
  pvector(const pvector &other) = delete;

  // prefer move because too much data to copy
  pvector(pvector &&other)
      : start_(other.start_), end_size_(other.end_size_),
        end_capacity_(other.end_capacity_) {
    other.start_ = nullptr;
    other.end_size_ = nullptr;
    other.end_capacity_ = nullptr;
  }

  // want move assignment
  pvector& operator= (pvector &&other) {
    start_ = other.start_;
    end_size_ = other.end_size_;
    end_capacity_ = other.end_capacity_;
    other.start_ = nullptr;
    other.end_size_ = nullptr;
    other.end_capacity_ = nullptr;
    return *this;
  }

  ~pvector() {
    if (start_ != nullptr)
      delete[] start_;
  }

  // not thread-safe
  void reserve(size_t num_elements) {
    if (num_elements > capacity()) {
      T_ *new_range = new T_[num_elements];
      #pragma omp parallel for
      for (size_t i=0; i < size(); i++)
        new_range[i] = start_[i];
      end_size_ = new_range + size();
      delete[] start_;
      start_ = new_range;
      end_capacity_ = start_ + num_elements;
    }
  }

  bool empty() {
    return end_size_ == start_;
  }

  void clear() {
    end_size_ = start_;
  }

  void resize(size_t num_elements) {
    reserve(num_elements);
    end_size_ = start_ + num_elements;
  }

  T_& operator[](size_t n) {
    return start_[n];
  }

  const T_& operator[](size_t n) const {
    return start_[n];
  }

  void push_back(T_ val) {
    if (size() == capacity()) {
      size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;
      reserve(new_size);
    }
    *end_size_ = val;
    end_size_++;
  }

  void push_back_atomic(T_ val) {
    *fetch_and_add(end_size_, sizeof(T_)) = val;
  }

  void fill(T_ init_val) {
    #pragma omp parallel for
    for (T_* ptr=start_; ptr < end_size_; ptr++)
      *ptr = init_val;
  }

  size_t capacity() {
    return end_capacity_ - start_;
  }

  size_t size() {
    return end_size_ - start_;
  }

  iterator begin() {
    return start_;
  }

  iterator end() {
    return end_size_;
  }

  T_* data() {
    return start_;
  }

  void swap(pvector &other) {
    std::swap(start_, other.start_);
    std::swap(end_size_, other.end_size_);
    std::swap(end_capacity_, other.end_capacity_);
  }


 private:
  T_* start_;
  T_* end_size_;
  T_* end_capacity_;
  static const size_t growth_factor = 2;
};

#endif  // PVECTOR_H_
