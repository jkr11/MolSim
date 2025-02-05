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
  std::vector<Particle> particles_;

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
  void addParticle(const Particle& p);

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

  /**
   * @brief does not use parallelization nor c18 coloring, done because of
   * pattern and inheritance
   * @param interactive_forces
   */
  [[deprecated]] void computeInteractiveForcesC18(
      const std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces)
      override;

  /**
   * does not use parallelization, done because of inheritance and pattern
   * @param interactive_forces
   */
  [[deprecated]] void computeInteractiveForcesForceBuffer(
      const std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces)
      override {
    computeInteractiveForcesC18(interactive_forces);
  }

  /**
   * @brief calculates the singular forces on all particles
   * @param singular_forces singular forces to be iterated over
   */
  void computeSingularForces(const std::vector<std::unique_ptr<SingularForce>>&
                                 singular_forces) override {}

  /**
   * Calculates the kinetic energy of the system 1/2 sum m_i abs(v_i^2)
   * @return kinetic energy of the system
   */
  double getKineticEnergy() override;

  /**
   * @brief Get particle count
   * @return particle in the container
   * @note this is technically unnecessary as its just size
   */
  size_t getParticleCount() override { return particles_.size(); }

  /**
   * @brief Get fixed particle count
   * @return number of fixed particles
   */
  size_t getSpecialParticleCount() override { return 0; };

  /**
   * @brief Get simulation domain
   * @return domain of the simulation
   */
  ivec3 getDomain() override;
};