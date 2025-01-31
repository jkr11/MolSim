//
// Created by mcarn on 10/18/24.
//

#include "defs/containers/DirectSumContainer.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include "debug/debug_print.h"
#include "defs/Particle.h"
#include "utils/ArrayUtils.h"

DirectSumContainer::DirectSumContainer() : ParticleContainer() {
  DEBUG_PRINT("DirectSumContainer::DirectSumContainer()");

  this->particles_ = {};
}

DirectSumContainer::DirectSumContainer(const std::vector<Particle>& particles)
    : ParticleContainer() {
  DEBUG_PRINT("explicit DirectSumContainer::DirectSumContainer()");
  for (const auto& particle : particles) {
    this->particles_.push_back(particle);
  }
}

void DirectSumContainer::addParticle(const Particle& p) {
  particles_.push_back(p);
}

void DirectSumContainer::addParticles(const std::vector<Particle>& particles) {
  for (const Particle& p : particles) {
    addParticle(p);
  }
}

void DirectSumContainer::removeParticle(const Particle& p) {
  particles_.erase(std::remove_if(particles_.begin(), particles_.end(),
                                  [&p](const Particle& q) { return p == q; }),
                   particles_.end());
}

std::vector<Particle*> DirectSumContainer::getParticles() {
  std::vector<Particle*> refs;
  refs.reserve(particles_.size());
  for (auto& p : particles_) {
    refs.push_back(&p);
  }

  return refs;
}

std::vector<Particle> DirectSumContainer::getParticlesObjects() {
  std::vector<Particle> refs;
  refs.reserve(particles_.size());
  for (auto& p : particles_) {
    refs.push_back(p);
  }
  return refs;
}

[[nodiscard]] std::size_t DirectSumContainer::size() const {
  return particles_.size();
}

void DirectSumContainer::singleIterator(
    const std::function<void(Particle&)>& f) {
  for (auto& p : particles_) {
    f(p);
  }
}

void DirectSumContainer::pairIterator(
    const std::function<void(Particle&, Particle&)>& f) {
  // note that the upper tri-diag matrix is iterated over
  for (size_t i = 0; i < particles_.size(); ++i) {
    for (size_t j = i + 1; j < particles_.size(); ++j) {
      f(particles_[i], particles_[j]);
    }
  }
}

double DirectSumContainer::getKineticEnergy() {
  double e_kin = 0.0;
  singleIterator([&e_kin](const Particle& p) {
    e_kin += 0.5 * p.getM() * ArrayUtils::L2InnerProduct(p.getV());
  });
  return e_kin;
}

void DirectSumContainer::imposeInvariant() {
  SPDLOG_TRACE("DirectSumContainer::imposeInvariant()");
}

ivec3 DirectSumContainer::getDomain() {
  SPDLOG_TRACE("DirectSumContainer::getDomain()");
  return {-1, -1, -1};
}

void DirectSumContainer::computeInteractiveForcesC18(
    const std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces) {
  SpdWrapper::get()->warn(
      "DirectSumContainer does not parallelize and does not use C18");
  for (size_t i = 0; i < particles_.size(); ++i) {
    dvec3 force_accum = {0.0, 0.0, 0.0};
    for (size_t j = i + 1; j < particles_.size(); ++j) {
      dvec3 f12 = {0.0, 0.0, 0.0};
      for (auto& force : interactive_forces) {
        dvec3 ftmp = force->directionalForce(particles_[i], particles_[j]);
        f12 = f12 + ftmp;
      }
      force_accum = force_accum + f12;
      particles_[j].subF(f12);
    }
    particles_[i].addF(force_accum);
  }
}
