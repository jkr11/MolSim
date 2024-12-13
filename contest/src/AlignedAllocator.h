//
// Created by mcarn on 12/13/24.
//

#pragma once

#include <vector>
#include <memory>
#include <cstdlib>
#include <iostream>

template <typename T>
class AlignedAllocator {
public:
  using value_type = T;

  AlignedAllocator() noexcept = default;

  template <typename U>
  constexpr AlignedAllocator(const AlignedAllocator<U>&) noexcept {}

  T* allocate(std::size_t n) {
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 64, n * sizeof(T)) != 0)
      throw std::bad_alloc();

    return static_cast<T*>(ptr);
  }

  void deallocate(T* p, std::size_t) noexcept {
    free(p);
  }
};

template <typename T, typename U>
bool operator==(const AlignedAllocator<T>&, const AlignedAllocator<U>&) {
  return true;
}

template <typename T, typename U>
bool operator!=(const AlignedAllocator<T>&, const AlignedAllocator<U>&) {
  return false;
}
