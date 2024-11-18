#pragma once
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "defs/containers/ParticleContainer.h"

typedef std::array<int, 3> ivec3;

class LinkedCellsContainer final : public ParticleContainer {
 private:
  /**
   * x is left - right
   * y is up - down
   * z is back - front
   *
   * position = x * (cellsY * cellsZ) + y * (cellsZ) + z
   */
  std::vector<std::vector<Particle>> cells;

  // number of cells for domain + 2 (halo)
  std::array<int, 3> cellCount{};
  // cell dimensions
  ivec3 cellDim{};
  // cutoff distance
  double cutoff{};

 public:
  LinkedCellsContainer() = default;

  /**
   * @brief Constructs a ParticleContainer which is implemented using
   * linkedCells
   * @param domain Domain of the container
   * @param cutoff maximum distance between valid particle pairs
   */
  explicit LinkedCellsContainer(const ivec3& domain, double cutoff);

  /**
   * @brief Destructor
   */
  ~LinkedCellsContainer() override = default;

  /**
   * @brief Add a particle to the container
   * @param p Particle to be added
   */
  void addParticle(const Particle& p) override;

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
  [[nodiscard]] std::vector<Particle> getParticles() const override;

  /**
   * @brief Get the count of particles in the container
   * @return Count of particles in the container
   */
  [[nodiscard]] std::size_t size() const override;

  /**
   * @brief Impose the invarient, that the particles are spatially sorted into
   * the correct vectors
   */
  void imposeInvariant();

  /**
   * @brief Single iterator over all particles in the container
   * @param f Function to be applied
   * @note Does not impose the invariant automatically!
   */
  void singleIterator(const std::function<void(Particle&)>& f) override;

  /**
   * @brief Pair iterator over all distinct particle pairs in the container with
   * distance <= cutoff
   * @param f Function to be applied
   * @note Does not impose the invariant automatically!
   */
  void pairIterator(
      const std::function<void(Particle&, Particle&)>& f) override;

  /**
   * @brief Single iterator over all particles in the boundary of the container
   * @param f Function to be applied
   * @note Does not impose the invariant automatically!
   */
  void boundaryIterator(const std::function<void(Particle&)>& f);

  /**
   * @brief Single iterator over all particles in the halo of the container
   * @param f Function to be applied
   * @note Does not impose the invariant automatically!
   */
  void haloIterator(const std::function<void(Particle&)>& f);

  /**
   * @brief Get the amount of cells in each dimension
   * @return ivec3 of cells in each dimension
   */
  [[nodiscard]] std::array<int, 3> getCellCount() const { return cellCount; }

  /**
   * @brief Get the dimensions of a all cells in the container
   * @return dvec3 of the dimensions of all cells
   */
  [[nodiscard]] ivec3 getCellDim() const { return cellDim; }

  /**
   * @brief Gets the cell index of a position
   * @param position Position in space
   * @return Associated cell index
   */
  [[nodiscard]] inline std::size_t dvec3ToCellIndex(
      const dvec3& position) const;

  /**
   * @brief Gets the cell index of the specified cell coordinate
   * @param position Cell coordinate in 3 dimensions
   * @return Associated cell index
   */
  [[nodiscard]] inline std::size_t cellCoordToIndex(
      std::array<int, 3> position) const;

  /**
   * @brief Gets the cell coordinate from the cell index
   * @param cellIndex Index of the cell
   * @return Cell coodinate in 3 dimensions
   */
  [[nodiscard]] inline std::array<int, 3> cellIndexToCoord(
      std::size_t cellIndex) const;

  /**
   * @brief Checks if a cell coordinate exists in the container
   * @param coordinate Cell coordinate to be checked
   * @return If cell coordinate is inside the container
   */
  [[nodiscard]] inline bool isValidCellCoordinate(ivec3 coordinate) const;

  /**
   * @brief Checks if a cell coordinate is in the halo of the container
   * @param cellCoord Cell coordinate to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] inline bool isHaloVec(ivec3 cellCoord) const;

  /**
   * @brief Checks if a cell index is in the halo of the container
   * @param cellIndex cell index to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] inline bool isHaloCell(std::size_t cellIndex) const;

  /**
   * @brief Checks if a cell coordinate is in the boundary of the container
   * @param cellCoord Cell coordinate to be checked
   * @return cell is part of the boundary
   */
  [[nodiscard]] inline bool isBoundaryVec(ivec3 cellCoord) const;

  /**
   * @brief Checks if a cell index is in the boundary of the container
   * @param cellIndex cell index to be checked
   * @return If cell is part of the boundary
   */
  [[nodiscard]] inline bool isBoundaryCell(std::size_t cellIndex) const;

  /**
   * @brief Debug method to get direct access to the cells vector
   * @return Reference to the cell vector
   */
  std::vector<std::vector<Particle>>& getCells() { return cells; }
};