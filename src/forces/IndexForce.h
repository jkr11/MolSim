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
  ivec3 dimensions{};

 public:
  explicit IndexForce(const std::vector<ivec3>& indeces, const double time,
                      const dvec3& force_values, const ivec3& dimensions)
      : indeces(getIndeces(indeces)),
        time(time),
        force_values(force_values),
        dimensions(dimensions) {}

  void applyForce(std::vector<Particle>& particles);

  [[nodiscard]] std::vector<int> getIndeces(
      const std::vector<ivec3>& two_d_indeces) const {
    std::vector<int> indeces;
    for (auto [fst, snd, thd] : two_d_indeces) {
      indeces.push_back((fst * dimensions[0] * dimensions[2]) +
                        (snd * dimensions[2]) + thd);
    }
    return indeces;
  }
};

#endif  // INDEXFORCE_H
