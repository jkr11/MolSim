//
// Created by jkr on 11/16/24.
//

#ifndef XMLREADER_H
#define XMLREADER_H
#pragma once
#include "defs/Particle.h"
#include "io/CLArgumentParser.h"

/**
 * @brief class for reading xml files specified by input.xsd in /input/
 *
 * Currently supports reading cuboids and spheroids simultaneously, metadata is
 * optional as it is currently handled by CLArgumentParser
 */
class XmlReader {
 public:
  XmlReader() = default;

  /**
   * @brief supports reading cuboids and sphereoids as specified in input.xsd
   * @param particles particles in which we store the configuration passed
   * @param filepath path at which the xml file is located
   * @param simulation_parameters passed form MolSim.cpp
   */
  static void read(std::vector<Particle>& particles,
                   const std::string& filepath,
                   Arguments& simulation_parameters);

  static void loadCheckpoint(const std::string& filepath, std::vector<Particle>& particles);
};

/**
 * @brief converts the xsd type to a ::Arguments:: type
 * @param boundary_type is the choice type passed from the xsd schema
 * specification
 * @return the enum type for arguments
 */
template <typename BT>
LinkedCellsConfig::BoundaryType toBoundaryType(const BT& boundary_type);

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
#endif  // XMLREADER_H
