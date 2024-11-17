//
// Created by jkr on 11/16/24.
//

#ifndef XMLREADER_H
#define XMLREADER_H
#pragma once
#include "io/file/in/FileReader.h"
#include "defs/Simulation.h"

/**
 * @brief class for reading xml files specified by input.xsd in /input/
 *
 * Currently supports reading cuboids and spheroids simultaneously, metadata is
 * optional as it is currently handled by CLArgumentParser
 */
class XmlReader final : public FileReader {
  private:
  Simulation simulation_parameters;
 public:
  explicit XmlReader(const Simulation &_simulation_parameters) {
    simulation_parameters = _simulation_parameters;
  };

  /**
   * @brief supports reading cuboids and sphereoids as specified in input.xsd
   * @param particles particles in which we store the configuration passed
   * @param filepath path at which the xml file is located
   */
  void read(std::vector<Particle> &particles, const std::string &filepath) override;

  /**
   * @brief passes the struct to other classes
   * @return struct with Simulation parameters
   */
  [[nodiscard]] Simulation pass() const;
};
#endif  // XMLREADER_H
