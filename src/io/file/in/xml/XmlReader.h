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

  /**
   * @brief validate that reflective boundaries are always pairs
   * @param boundary the boundary configuration
   */
  static inline void validateBoundaries(
      const LinkedCellsConfig::BoundaryConfig& boundary);

  /**
   * @brief Reads a xml checkpoint file as specified in
   * MolSim/src/io/file/out/checkpoint-schema.xsd
   * @param filepath path pointing to checkpoint // TODO default location
   * @param particles particles to which the checkpoint is read
   */
  static void loadCheckpoint(const std::string& filepath,
                             std::vector<Particle>& particles);

  static void loadCheckpointMembrane(const std::string& filepath,
                                     const ivec3& dimensions,
                                     std::vector<Particle>& particles);
};

/**
 * @brief validates that the number of bins are correct and not smaller than 1
 * and the output time is larger than 0
 */
inline void validateStatisticsInput(const StatisticsConfig& stats) {
  if (stats.x_bins < 1 || stats.y_bins < 1) {
    SpdWrapper::get()->error(
        "Number of x-bins ({}) or y-bins ({}) is too small", stats.x_bins,
        stats.y_bins);
    throw std::runtime_error("Number of x-bins or y-bins is too small");
  }

  if (stats.output_interval < 1) {
    SpdWrapper::get()->error(
        "Output interval for the statistics is too small ({})",
        stats.output_interval);
    throw std::runtime_error("Output interval for the statistics is too small");
  }
}

/**
 * @brief ensures the input files are valid (not empty, extension, path)
 * @param path the path pointing to the file
 * @param extension the desired file extension
 * @param type the type of file you are checking (Input, checkpoint, etc)
 */
inline void validatePath(const std::filesystem::path& path,
                          const std::string& extension,
                          const std::string& type) {
  if (!exists(path)) {
    throw std::runtime_error(type + " file not found: " + path.string());
  }
  if (path.extension() != extension) {
    throw std::invalid_argument(
        type + " file extension is not supported: " + path.string());
  }
  if (std::filesystem::is_empty(path)) {
    throw std::runtime_error(type + "file is empty " + path.string());
  }
}

/**
 * @brief converts the xsd type to a ::Arguments:: type
 * @param boundary_type is the choice type passed from the xsd schema
 * specification
 * @return the enum type for arguments
 */
template <typename Bt>
LinkedCellsConfig::BoundaryType toBoundaryType(const Bt& boundary_type);

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
