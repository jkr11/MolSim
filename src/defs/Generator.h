//
// Created by jkr on 10/28/24.
//
#pragma once
#ifndef GENERATOR_H
#define GENERATOR_H
#include "Particle.h"

#include <array>
#include <vector>

class ParticleGenerator {
private:
  std::array<int, 3> dimensions{};
  dvec3 origin{};
  double h{};
  double m{};
  dvec3 velocity{};
  double mvbm = 1.0;

public:
  ParticleGenerator();
  ParticleGenerator(const ParticleGenerator& other);
  ParticleGenerator(
        std::array<double, 3> origin,
        std::array<int, 3> dimensions,
        double h,
        double m,
        std::array<double, 3> initialVelocity,
        double temperature
    ) : dimensions(dimensions), origin(origin), h(h), m(m), velocity({0,0,0}), mvbm(1.0) {}
  std::vector<Particle> generate_cuboid() {
    std::vector<Particle> particles;
    for(int i = 0; i < dimensions[0]; i++) {
      for(int j = 0; j < dimensions[1]; j++) {
        for(int k = 0; k < dimensions[2]; k++) {
          std::array<int,3> pos = {
            origin[0] + h * i,
            origin[1] + h * j,
            origin[2] + h * k
          };
          particles.emplace_back(pos, velocity, m);

  }
};

#endif // GENERATOR_H
