//
// Created by jkr on 12/26/24.
//

#ifndef INDEXFORCE_H
#define INDEXFORCE_H

#include <vector>

#include "defs/containers/ParticleContainer.h"

class IndexForce {
 private:
  std::vector<int> indeces{};
  double time{};
  dvec3 force_values{};

 public:
  explicit IndexForce(const std::vector<int>& indeces, const double time,
                      const dvec3& force_values)
      : indeces(indeces), time(time), force_values(force_values) {}

  void applyForce(std::vector<Particle>& particles);
};

#endif  // INDEXFORCE_H
