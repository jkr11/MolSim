#pragma once
#include <cstddef>
#include <vector>
#include <functional>

#include "DirectSumContainer.h"
#include "../Particle.h"

  DirectSumContainer::DirectSumContainer() : ParticleContainer() {
    this->particles = {};
  }

  DirectSumContainer::DirectSumContainer(const std::vector<Particle>& particles) : ParticleContainer() {
    for (const auto &particle : particles) {
      this->particles.push_back(particle);
    }
  }

  //DirectSumContainer::~DirectSumContainer();

  void DirectSumContainer::addParticle(const Particle &p) {
    particles.push_back(p);
  }

  void DirectSumContainer::removeParticle(const Particle &p) {
    particles.erase(std::remove_if(particles.begin(), particles.end(),
                                    [&p](const Particle& q) { return &p == &q; }),
                    particles.end());
  }

  std::vector<Particle> DirectSumContainer::getParticles() const {
    return particles; //TODO: is this a copy because it should be
  }

  [[nodiscard]] std::size_t DirectSumContainer::size() const {
    return particles.size();
  }

  void DirectSumContainer::singleIterator(const std::function<void(Particle&)>& f) {
    for (auto &p : particles) {
      f(p);
    }
  }

  void DirectSumContainer::pairIterator(const std::function<void(Particle&, Particle&)>& f)  {
    // note that the upper tri-diag matrix is iterated over
    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = i + 1; j < particles.size(); ++j) {
        f(particles[i], particles[j]);
      }
    }
  }