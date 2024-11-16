#pragma once
#include <cstddef>
#include <vector>
#include <functional>

#include "../containers/ParticleContainer.h"
#include "../Particle.h"

class DirectSumContainer : public ParticleContainer {
private:
  std::vector<Particle> particles;

public:
  DirectSumContainer();

  explicit DirectSumContainer(const std::vector<Particle>& particles);

  ~DirectSumContainer() = default;

  void addParticle(const Particle &p) override;

  void removeParticle(const Particle &p) override;

  std::vector<Particle> getParticles() const override;

  [[nodiscard]] std::size_t size() const override;

  void singleIterator(const std::function<void(Particle&)>& f) override;

  void pairIterator(const std::function<void(Particle&, Particle&)>& f) override;
};