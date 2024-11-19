//
// Created by mcarn on 10/18/24.
//

#include "defs/containers/DirectSumContainer.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"

DirectSumContainer::DirectSumContainer() : ParticleContainer() {
  this->particles = {};
}

DirectSumContainer::DirectSumContainer(const std::vector<Particle>& particles)
    : ParticleContainer() {
  for (const auto& particle : particles) {
    this->particles.push_back(particle);
  }
}

// DirectSumContainer::~DirectSumContainer();

void DirectSumContainer::addParticle(const Particle& p) {
  particles.push_back(p);
}

void DirectSumContainer::addParticles(const std::vector<Particle>& particles) {
  for (const auto& p : particles) {
    addParticle(p);
  }
}

void DirectSumContainer::removeParticle(const Particle& p) {
  particles.erase(std::remove_if(particles.begin(), particles.end(),
                                 [&p](const Particle& q) { return p == q; }),
                  particles.end());
}

std::vector<Particle*> DirectSumContainer::getParticles() {
  std::vector<Particle*> refs;
  refs.reserve(particles.size());
  for (auto &p : particles) {
    refs.push_back(&p);
  }

  return refs;
}

[[nodiscard]] std::size_t DirectSumContainer::size() const {
  return particles.size();
}

void DirectSumContainer::singleIterator(
    const std::function<void(Particle&)>& f) {
  for (auto& p : particles) {
    f(p);
  }
}

void DirectSumContainer::pairIterator(
    const std::function<void(Particle&, Particle&)>& f) {
  // note that the upper tri-diag matrix is iterated over
  for (size_t i = 0; i < particles.size(); ++i) {
    for (size_t j = i + 1; j < particles.size(); ++j) {
      f(particles[i], particles[j]);
    }
  }
}

void DirectSumContainer::imposeInvariant(){};