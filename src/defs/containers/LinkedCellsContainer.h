#pragma once
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "defs/Simulation.h"
#include "defs/containers/ParticleContainer.h"

/**
 * @brief a particle container with linked cells
 */
class LinkedCellsContainer final : public ParticleContainer {
 private:
  /**
   * @brief
   * x is left - right
   * y is up - down
   * z is back - front
   *
   * position = x * (cellsY * cellsZ) + y * (cellsZ) + z
   */
  std::vector<std::vector<Particle>> cells;

  /**
   * @brief
   * stores the indexes of all halo_cells for faster iteration in the
   * corresponding direction vector
   */
  std::array<std::vector<std::size_t>, 6> halo_direction_cells;

  /**
   * @brief
   * stores the indexes of all halo_cells for faster iteration in the
   * corresponding direction vector
   */
  std::array<std::vector<std::size_t>, 6> boundary_direction_cells;

  /**
   * @brief
   * a more processing friendly storage of LinkedCellsConfig::BoundaryConfig
   */
  std::array<LinkedCellsConfig::BoundaryType, 6> boundaries;

  /**
   * @brief
   * number of cells for domain + 2 (halo)
   */
  ivec3 cell_count{};

  /**
   * @brief
   * cell dimensions
   */
  ivec3 cell_dim{};

  /**
   * @brief
   * cutoff radius
   */
  double cutoff{};

  /**
   * @brief
   * the boundary config of each direction of the simulation
   */
  LinkedCellsConfig::BoundaryConfig boundary_config{};

 public:
  /**
   * Empty constructor
   * TODO: why does this even exist?
   */
  LinkedCellsContainer() = default;

  /**
   * @brief Constructs a ParticleContainer which is implemented using
   * linkedCells
   * @param linked_cells_config configuration struct of the simulation
   */
  explicit LinkedCellsContainer(const LinkedCellsConfig& linked_cells_config);

  /**
   * @brief Destructor
   */
  ~LinkedCellsContainer() override = default;

  /**
   * @brief Add a particle to the container
   * @param p Particle to be added
   * @note Does not impose the invariant automatically!
   */
  void addParticle(const Particle& p) override;

  /**
   * @brief Add a vector of particles to the container
   * @param particles Particles to be added
   * @note Does not impose the invariant automatically!
   */
  void addParticles(const std::vector<Particle>& particles) override;

  /**
   * @brief Remove a particle from the container
   * @param p Particle to be removed
   * @note Does not impose the invariant automatically!
   */
  void removeParticle(const Particle& p) override;

  /**
   * @brief Get a vector of all references to particles in the container
   * @return Vector of references to particles in the container
   */
  [[nodiscard]] std::vector<Particle*> getParticles() override;

  /**
   * @brief Get the count of particles in the container
   * @return Count of particles in the container
   */
  [[nodiscard]] std::size_t size() const override;

  /**
   * @brief Impose the invariant, that the particles are spatially sorted into
   * the correct vectors
   */
  void imposeInvariant() override;

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
  [[nodiscard]] std::array<int, 3> getCellCount() const { return cell_count; }

  /**
   * @brief Get the dimensions of a all cells in the container
   * @return dvec3 of the dimensions of all cells
   */
  [[nodiscard]] ivec3 getCellDim() const { return cell_dim; }

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
  [[nodiscard]] inline ivec3 cellIndexToCoord(std::size_t cellIndex) const;

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
  [[nodiscard]] inline bool isHalo(ivec3 cellCoord) const;

  /**
   * @brief Checks if a cell index is in the halo of the container
   * @param cellIndex cell index to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] inline bool isHalo(std::size_t cellIndex) const;

  /**
   * @brief Checks if a cell coordinate is in the boundary of the container
   * @param cellCoord Cell coordinate to be checked
   * @return cell is part of the boundary
   */
  [[nodiscard]] inline bool isBoundary(ivec3 cellCoord) const;

  /**
   * @brief Checks if a cell index is in the boundary of the container
   * @param cellIndex cell index to be checked
   * @return If cell is part of the boundary
   */
  [[nodiscard]] inline bool isBoundary(std::size_t cellIndex) const;

  /**
   * @brief calculates all directions of the halo cell
   *  empty: no halo cell
   *  0: west
   *  1: east
   *  2: down
   *  3: up
   *  4: south
   *  5: north
   * @param cellIndex cell index to be checked
   * @return the direction of the halo cell
   */
  [[nodiscard]] std::vector<std::size_t> halo_direction(
      std::size_t cellIndex) const;

  /**
   * @brief calculates all directions of the boundary cell
   *
   *  empty: no boundary cell
   *  0: west
   *  1: east
   *  2: down
   *  3: up
   *  4: south
   *  5: north
   * @param cellIndex cell index to be checked
   * @return the directions of the boundary cell
   */
  [[nodiscard]] std::vector<std::size_t> boundary_direction(
      std::size_t cellIndex) const;

  /**
   * @brief Debug method to get direct access to the cells vector
   * @return Reference to the cell vector
   */
  std::vector<std::vector<Particle>>& getCells() { return cells; }
};