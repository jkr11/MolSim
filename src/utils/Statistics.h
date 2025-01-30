//
// Created by maximilian on 19.01.25.
//
#pragma once


#include <defs/containers/ParticleContainer.h>

#include "io/file/out/CSVWriter.h"
#include "utils/ArrayUtils.h"

class Statistics {
 private:
  int x_bins_;
  int y_bins_;
  CSVWriter density_profile_writer_;
  CSVWriter velocity_profile_writer_;
  ParticleContainer& container_;
  double x_bin_size_;
  double y_bin_size_;
  double bin_volume_;

 public:
  Statistics(const int x_bins, const int y_bins, ParticleContainer& container,
             const std::string& density_profile_output_location,
             const std::string& velocity_profile_output_location)
      : x_bins_(x_bins),
        y_bins_(y_bins),
        density_profile_writer_(density_profile_output_location),
        velocity_profile_writer_(velocity_profile_output_location),
        container_(container) {
    const ivec3 domain = container.getDomain();
    x_bin_size_ = 1.0 * domain[0] / x_bins;
    y_bin_size_ = 1.0 * domain[1] / y_bins;
    bin_volume_ = x_bin_size_ * y_bin_size_ * domain[2];
  }

  /**
   * @brief calculates the density and velocity profile for the bins.
   * The bins can also be along the y axis at the same time for a checkerboard
   * pattern. The bins are from left to right (x-axis), and for the lines from
   * down to up (y-axis)
   */
  void writeStatistics(const double time) {
    std::vector<std::vector<Particle*>> bins(x_bins_ * y_bins_);
    for (const auto p : container_.getParticles()) {
      if (p->getType() < 0)
        continue;  // ignore walls, because their velocity is not strictly set
                   // to 0 but rather ignored in other calculations

      dvec3 position = p->getX();
      const int x_bin = static_cast<int>(position[0] / x_bin_size_);
      const int y_bin = static_cast<int>(position[1] / y_bin_size_);

      bins[x_bin + y_bin * x_bins_].push_back(p);
    }

    std::vector<std::string> density_data;
    std::vector<std::string> velocity_data;

    for (const auto& bin : bins) {
      const size_t num_particles = bin.size();
      dvec3 sum_velocity = {0.0, 0.0, 0.0};

      for (const auto p : bin) {
        sum_velocity = sum_velocity + p->getV();
      }

      density_data.push_back(
          std::to_string(static_cast<double>(num_particles) / bin_volume_));

      if (num_particles == 0) {
        velocity_data.emplace_back("0.0 0.0 0.0");
      } else {
        velocity_data.push_back(
            std::to_string(sum_velocity[0] /
                           static_cast<double>(num_particles)) +
            " " +
            std::to_string(sum_velocity[1] /
                           static_cast<double>(num_particles)) +
            " " +
            std::to_string(sum_velocity[2] /
                           static_cast<double>(num_particles)));
      }
    }

    density_profile_writer_.writeLine(time, density_data);
    velocity_profile_writer_.writeLine(time, velocity_data);
  }

  void closeFiles() {
    density_profile_writer_.closeFile();
    velocity_profile_writer_.closeFile();

    SpdWrapper::get()->info("Statistics were completed");
  }
};