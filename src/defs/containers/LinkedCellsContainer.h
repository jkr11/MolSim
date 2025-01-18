#pragma once
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "defs/Simulation.h"
#include "defs/containers/ParticleContainer.h"

/**
 * @brief a particle container with linked cells
 * @image html benchmark/graph.png
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
   * @brief current number of particles
   */
  size_t particle_count{};

  /**
   * @brief number of particles, that are immovable
   */
  size_t special_particle_count{};

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
  std::array<LinkedCellsConfig::BoundaryType, 6> boundaries{};

  /**
   * @brief
   * number of cells per direction for domain + 2 halo cells
   */
  ivec3 cell_count{};

  /**
   * @brief
   * cell dimensions
   */
  dvec3 cell_dim{};

  /**
   * @brief
   * cutoff radius
   */
  double cutoff{};

  /**
   * @brief
   * the domain of the container
   */
  ivec3 domain{};

  /**
   * @brief
   * the boundary config of each direction of the simulation
   */
  LinkedCellsConfig::BoundaryConfig boundary_config{};

  /**
   * @brief apply reflective boundary condition to a dimension
   * @param dimension the problematic dimension
   */
  inline void apply_reflective_boundary(size_t dimension);

  /**
   *@brief index offsets orthogonal to a cell for each dimension, optimized for
   *2D simulations
   *
   */
  std::array<std::array<ivec3, 9>, 3> index_offsets = {{
      // x
      {{{1, 1, 0},
        {1, 0, 0},
        {1, -1, 0},
        // optional for 3D
        {1, 1, 1},
        {1, 0, 1},
        {1, -1, 1},
        {1, 1, -1},
        {1, 0, -1},
        {1, -1, -1}}},
      // y
      {{{1, 1, 0},
        {0, 1, 0},
        {-1, 1, 0},
        // optional for 3D
        {1, 1, 1},
        {0, 1, 1},
        {-1, 1, 1},
        {1, 1, -1},
        {0, 1, -1},
        {-1, 1, -1}}},
      // z, order irrelevant for 2D
      {{{1, 1, 1},
        {1, 0, 1},
        {1, -1, 1},
        {0, 1, 1},
        {0, 0, 1},
        {0, -1, 1},
        {-1, 1, 1},
        {-1, 0, 1},
        {-1, -1, 1}}},
  }};

 public:
  /**
   * 6th root of 2
   */
  static const double sigma_factor;

  /**
   * Empty constructor
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
   * @brief the exact number of current particles, updated accordingly
   * @return the current count of particles left in the simulation
   */
  [[nodiscard]] size_t getParticleCount() const { return particle_count; }

  /**
   * @brief the exact number of current special particles, updated accordingly
   * @return the current count of special particles left in the simulation
   */
  [[nodiscard]] size_t getSpecialParticleCount() const {
    return special_particle_count;
  }

  /**
   * @brief Get a vector of all particles in the container
   * @return Vector of all particles
   */
  [[nodiscard]] std::vector<Particle> getParticlesObjects() override;

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
  [[nodiscard]] dvec3 getCellDim() const { return cell_dim; }

  /**
   * @brief Gets the cell index of a position
   * @param position Position in space
   * @return Associated cell index
   */
  [[nodiscard]] inline std::size_t dvec3ToCellIndex(
      const dvec3& position) const;

  /**
   * @brief API for testing, because gtest does not like inline
   * @param position Position in space
   * @return Associated cell index
   */
  [[nodiscard]] std::size_t dvec3ToCellIndex_testing(
      const dvec3& position) const;

  /**
   * @brief Gets the cell index of the specified cell coordinate
   * @param position Cell coordinate in 3 dimensions
   * @return Associated cell index
   */
  [[nodiscard]] inline std::size_t cellCoordToIndex(ivec3 position) const;

  /**
   * @brief API for testing
   * @param position Cell coordinate in 3 dimensions
   * @return Associated cell index
   */
  [[nodiscard]] std::size_t cellCoordToIndex_testing(ivec3 position) const;

  /**
   * @brief Gets the cell coordinate from the cell index
   * @param cellIndex Index of the cell
   * @return Cell coodinate in 3 dimensions
   */
  [[nodiscard]] inline ivec3 cellIndexToCoord(std::size_t cellIndex) const;

  /**
   * @brief API for testing
   * @param cellIndex Index of the cell
   * @return Cell coodinate in 3 dimensions
   */
  [[nodiscard]] ivec3 cellIndexToCoord_testing(std::size_t cellIndex) const;

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
   * @brief API for testing
   * @param cellIndex cell index to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] bool isHalo_testing(std::size_t cellIndex) const;

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
   * @brief API for testing, because gtest does not like inline
   * @param cellIndex cell index to be checked
   * @return If cell is part of the boundary
   */
  [[nodiscard]] bool isBoundary_testing(std::size_t cellIndex) const;

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
   * @param f function to check if it is a special cell
   * @param lowerMagicNumber lower bound index
   * @param upperMagicNumber upper bound index
   * @return the directions of the boundary cell
   */
  [[nodiscard]] std::vector<std::size_t> special_cell_direction(
      std::size_t cellIndex, const std::function<bool(std::size_t)>& f,
      int lowerMagicNumber, int upperMagicNumber) const;

  /**
   * @brief Debug method to get direct access to the cells vector
   * @return Reference to the cell vector
   */
  std::vector<std::vector<Particle>>& getCells() { return cells; }

  /**
   * @brief warp negative cell index to maximum cell coordinate to enable
   * multiple periodic boundaries in corners. For now this does only work in 2D
   * @param cell_coordinate the cell to be checked
   * @param raw_dimension the dimension axis looked at
   * @return bool: whether it is a valid cell to be checked, ivec3: real cell,
   * dvec3: offset to be applied
   */
  [[nodiscard]] inline std::tuple<bool, ivec3, dvec3> reflective_warp_around(
      ivec3 cell_coordinate, std::size_t raw_dimension) const;

  /**
   * API for testing
   * @param cell_coordinate the cell to be checked
   * @param raw_dimension the dimension axis looked at
   * @return bool: whether it is a valid cell to be checked, ivec3: real cell,
   * dvec3: offset to be applied
   */
  [[nodiscard]] std::tuple<bool, ivec3, dvec3> reflective_warp_around_testing(
      ivec3 cell_coordinate, std::size_t raw_dimension) const;

  double getKineticEnergy() override;

  std::size_t getParticleCount() override { return particle_count; }
  std::size_t getSpecialParticleCount() override {
    return special_particle_count;
  }
};

/**
 * @brief directions for better readability; implicitly cast
 */
enum Directions { xlow, xhigh, ylow, yhigh, zlow, zhigh };