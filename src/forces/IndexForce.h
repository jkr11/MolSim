//
// Created by jkr on 12/26/24.
//

#ifndef INDEXFORCE_H
#define INDEXFORCE_H

#include <vector>

#include "defs/containers/ParticleContainer.h"

class IndexForce {
  std::vector<int> indeces{};
  double time{};
  dvec3 force_values{};

 public:
  explicit IndexForce(const std::vector<int>& ids, const double time,
                      const dvec3& force_values)
      : indeces(ids), time(time), force_values(force_values) {}
  IndexForce() = default;

  dvec3 applyForce(Particle& p, double sim_time) const;

  [[nodiscard]] std::vector<int> getIndeces() const { return indeces; }
};

#endif  // INDEXFORCE_H
