//
// Created by jkr on 11/16/24.
//

#ifndef XMLREADER_H
#define XMLREADER_H
#pragma once
#include "defs/Simulation.h"
#include "io/file/in/FileReader.h"

/**
 * @brief class for reading xml files specified by input.xsd in /input/
 *
 * Currently supports reading cuboids and spheroids simultaneously, metadata is
 * optional as it is currently handled by CLArgumentParser
 */
class XmlReader final : public FileReader {
 private:
  Simulation simulation_parameters;
  double delta_t{};
  double t_end{};
  double cutoff_radius{};
  ivec3 domain{};

 public:
  explicit XmlReader(const Simulation &_simulation_parameters) {
    simulation_parameters = _simulation_parameters;
  };

  /**
   * @brief supports reading cuboids and sphereoids as specified in input.xsd
   * @param particles particles in which we store the configuration passed
   * @param filepath path at which the xml file is located
   */
  void read(std::vector<Particle> &particles,
            const std::string &filepath) override;

  /**
   * @brief passes the struct to other classes
   * @return struct with Simulation parameters
   */
  [[nodiscard]] std::tuple<double, double, double, ivec3> pass() const;
};
#endif  // XMLREADER_H
