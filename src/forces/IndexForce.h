//
// Created by jkr on 12/26/24.
//

#ifndef INDEXFORCE_H
#define INDEXFORCE_H

#include <vector>

#include "defs/containers/ParticleContainer.h"

class IndexForce {
  std::vector<int> indeces_{};
  double time_{};
  dvec3 force_values_{};

 public:
  explicit IndexForce(const std::vector<int>& ids, const double time,
                      const dvec3& force_values)
      : indeces_(ids), time_(time), force_values_(force_values) {}
  IndexForce() = default;

  dvec3 applyForce(Particle& p, double sim_time) const;

  [[nodiscard]] std::vector<int> getIndeces() const { return indeces_; }
};

#endif  // INDEXFORCE_H
