// See LICENSE.txt for license details.

#ifndef SLIDING_QUEUE_H_
#define SLIDING_QUEUE_H_

#include <algorithm>

#include "platform_atomics.h"

/*

Double-buffered queue so appends aren't seen until SlideWindow() called.
Use QueueBuffer when used in parallel to avoid false sharing by doing
bulk appends from thread-local storage.

*/

template <typename T>
class QueueBuffer;

template <typename T>
class SlidingQueue {
  T *shared;
  size_t shared_in;
  size_t shared_out_start;
  size_t shared_out_end;
  friend class QueueBuffer<T>;

 public:
  SlidingQueue(size_t shared_size) {
    shared = new T[shared_size];
    Reset();
  }

  ~SlidingQueue() {
    delete[] shared;
  }

  void Push(T to_add) {
    shared[shared_in++] = to_add;
  }

  bool Empty() {
    return shared_out_start == shared_out_end;
  }

  void Reset() {
    shared_out_start = 0;
    shared_out_end = 0;
    shared_in = 0;
  }

  void SlideWindow() {
    shared_out_start = shared_out_end;
    shared_out_end = shared_in;
  }

  typedef T* iterator;

  iterator begin() {
    return shared + shared_out_start;
  }

  iterator end() {
    return shared + shared_out_end;
  }

  size_t size() {
    return end() - begin();
  }
};


template <typename T>
class QueueBuffer {
  size_t in;
  T *local_queue;
  SlidingQueue<T> &sq;
  const size_t local_size;

 public:
  QueueBuffer(SlidingQueue<T> &master, size_t given_size=16384)
      : sq(master), local_size(given_size) {
    in = 0;
    local_queue = new T[local_size];
  }

  ~QueueBuffer() {
    delete[] local_queue;
  }

  void Push(T to_add) {
    if (in == local_size)
      Flush();
    local_queue[in++] = to_add;
  }

  void Flush() {
    T *shared_queue = sq.shared;
    size_t copy_start = fetch_and_add(sq.shared_in, in);
    std::copy(local_queue, local_queue+in, shared_queue+copy_start);
    in = 0;
  }
};

#endif  // SLIDING_QUEUE_H_
