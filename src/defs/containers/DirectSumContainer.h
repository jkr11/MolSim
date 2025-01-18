#pragma once
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "defs/containers/ParticleContainer.h"

/**
 * Direct sum container class, standard n^2 / 2 newton 3 scheme is applied here
 */
class DirectSumContainer final : public ParticleContainer {
 private:
  std::vector<Particle> particles;

 public:
  /**
   * @brief Constructs a empty DirectSumContainer
   */
  DirectSumContainer();

  /**
   * @brief Constructs a DirectSumContainer with the particles
   * @param particles Particles to be added
   */
  explicit DirectSumContainer(const std::vector<Particle>& particles);

  /**
   * @brief Destructor
   */
  ~DirectSumContainer() override = default;

  /**
   * @brief Add a particle to the container
   * @param p Particle to be added
   */
  void addParticle(const Particle& p) override;

  /**
   * @brief Add a vector of particles to the container
   * @param particles  Particles to be added
   */
  void addParticles(const std::vector<Particle>& particles) override;

  /**
   * @brief Remove a particle from the container
   * @param p Particle to be removed
   */
  void removeParticle(const Particle& p) override;

  /**
   * @brief Get a vector of all references to particles in the container
   * @return Vector of references to particles in the container
   */
  [[nodiscard]] std::vector<Particle*> getParticles() override;

  /**
   * @brief Get a vector of all particles in the container
   * @returns Vector of all particles
   */
  [[nodiscard]] std::vector<Particle> getParticlesObjects() override;

  /**
   * @brief Get the count of particles in the container
   * @return Count of particles in the container
   */
  [[nodiscard]] std::size_t size() const override;

  /**
   * @brief Impose the invariant, that the particles are spatially sorted into
   * the correct vectors.
   * @note as this is a directSum container, this does nothing
   */
  void imposeInvariant() override;

  /**
   * @brief Single iterator over all particles in the container
   * @param f Function to be applied
   */
  void singleIterator(const std::function<void(Particle&)>& f) override;

  /**
   * @brief Pair iterator over all distinct particle pairs in the container with
   * distance <= cutoff
   * @param f Function to be applied
   */
  void pairIterator(
      const std::function<void(Particle&, Particle&)>& f) override;

  double getKineticEnergy() override;

  size_t getParticleCount() override { return particles.size(); }

  size_t getSpecialParticleCount() override { return 0;};
};