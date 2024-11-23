//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once
#include <stdexcept>

#include "utils/SpdWrapper.h"

using dvec3 = std::array<double, 3>;
using ivec3 = std::array<int, 3>;

/**
 * @brief translates a vector from the xml parser to a valid "standard" c++ type
 * @tparam SVec source vector as one of the xsd types, e.g. IVec3Type
 * @tparam TVec target in c++-space, e.g ivec3 (:= std::array<int,3>)
 * @param source the vector passed form the xml parser
 * @param paramName the name of the parameter we are dealing with for throwing
 * an error
 * @return the vector in c++-space types.
 */
template <typename SVec, typename TVec>
TVec unwrapVec(const SVec& source, const std::string& paramName) {
  try {
    return TVec{source.x(), source.y(), source.z()};
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to unwrap vector " + paramName + ": " +
                             e.what());
  }
}

#endif  // SIMULATION_H
