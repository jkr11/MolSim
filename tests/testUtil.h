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

/**
 * @brief compares two ivec3's
 * @param a ivec3 source
 * @param b ivec3 target
 */
inline void EXPECT_IVEC3_EQ(const ivec3& a, const ivec3& b) {
  EXPECT_EQ(a.size(), b.size()) << "Lists are of different sizes";
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_EQ(a[i], b[i]) << "Elements at index " << i << " are not equal";
  }
}

/**
 * @brief checks if a list contains an element
 * @tparam T template type T
 * @param elt query element
 * @param list list possibly containing elt
 */
template <typename T>
void VEC_CONTAINS(const T elt, const std::vector<T>& list) {
  EXPECT_TRUE(std::find(list.begin(), list.end(), elt));
}