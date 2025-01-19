//
// Created by maximilian on 19.01.25.
//
#pragma once

#include "defs/containers/LinkedCellsContainer.h"
#include "io/file/out/CSVWriter.h"
class Statistics {
 private:
  std::string density_profile_output_location;
  std::string velocity_profile_output_location;
  int x_bins;
  int y_bins;
  CSVWriter density_profile_writer;
  CSVWriter velocity_profile_writer;
  LinkedCellsContainer& container;
  double x_bin_size;
  double y_bin_size;

 public:
  Statistics(const int x_bins, const int y_bins,
             LinkedCellsContainer& container,
             const std::string& density_profile_output_location,
             const std::string& velocity_profile_output_location)
      : density_profile_output_location(density_profile_output_location),
        velocity_profile_output_location(velocity_profile_output_location),
        x_bins(x_bins),
        y_bins(y_bins),
        density_profile_writer(density_profile_output_location),
        velocity_profile_writer(velocity_profile_output_location),
        container(container) {
    const ivec3 domain = container.getDomain();
    x_bin_size = 1.0 * domain[0] / x_bins;
    y_bin_size = 1.0 * domain[1] / y_bins;
  }

  /**
   * @brief calculates the density and velocity profile for the bins.
   * The bins can also be along the y axis at the same time for a checkerboard
   * pattern. The bins are from left to right (x-axis), and for the lines from
   * down to up (y-axis)
   */
  void writeStatistics() {
    std::vector<std::vector<Particle&>> bins(x_bins * y_bins);
    for (const auto p : container.getParticles()) {
      dvec3 position = p->getX();
      dvec3 velocity = p->getV();

      int x_bin = static_cast<int>(position[0] / x_bin_size);
      int y_bin = static_cast<int>(position[1] / y_bin_size);

    }
  }
};