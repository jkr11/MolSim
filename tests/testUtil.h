#pragma once
#include "../src/defs/Particle.h"
#include "gtest/gtest.h"

/**
 * @brief checks if vector a is close to vector b to a desired tolerance epsilon
 * @note this is basiclly like numpy.allclose
 *
 * @param a vector a
 * @param b vector b
 * @param error error message on failure
 * @param eps tolerance, default = 0.001
 */
inline void DVEC3_NEAR(const dvec3& a, const dvec3& b, const std::string& error,
                       const double eps = 0.001) {
  for (int i = 0; i < 3; ++i) EXPECT_NEAR(a[i], b[i], eps) << error;
}

/**
 * @brief Check if two vectors are the same
 * @param list1 first list
 * @param list2 second list
 * @tparam T generic object type
 */
template <typename T>
void EXPECT_VECTOR_EQ(const std::vector<T>& list1,
                      const std::vector<T>& list2) {
  EXPECT_EQ(list1.size(), list2.size()) << "Lists are of different sizes";

  for (size_t i = 0; i < list1.size(); ++i) {
    EXPECT_EQ(list1[i], list2[i])
        << "Elements at index " << i << " are not equal";
  }
}