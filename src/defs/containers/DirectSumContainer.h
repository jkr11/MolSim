#pragma once
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "defs/containers/ParticleContainer.h"

class DirectSumContainer : public ParticleContainer {
 private:
  std::vector<Particle> particles;

 public:
  DirectSumContainer();

  explicit DirectSumContainer(const std::vector<Particle>& particles);

  ~DirectSumContainer() override = default;

  void addParticle(const Particle& p) override;

  void addParticles(const std::vector<Particle>& particles) override;

  void removeParticle(const Particle& p) override;

  [[nodiscard]] std::vector<Particle> getParticles() const override;

  [[nodiscard]] std::size_t size() const override;

  void singleIterator(const std::function<void(Particle&)>& f) override;

  void pairIterator(
      const std::function<void(Particle&, Particle&)>& f) override;

  void imposeInvariant() override;
};