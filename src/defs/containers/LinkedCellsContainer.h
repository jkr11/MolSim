#pragma once
#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "defs/Simulation.h"
#include "defs/containers/ParticleContainer.h"
#include "forces/IndexForce.h"
#include "forces/InteractiveForce.h"
#include "debug/debug_print.h"

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

  std::vector<Particle> particles_;

  std::vector<std::vector<Particle*>> cells_;

  // std::vector<std::vector<Particle*>> cell_orders_;

  std::vector<std::vector<std::size_t>> c18_colours_;

  std::vector<ivec3> c_18_schema_ = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0},
    {-1, 0, 1},   {-1, 1, -1}, {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1},
    {0, -1, 0},   {0, -1, 1},  {0, 0, -1},  {0, 0, 0},   {0, 0, 1},
    {0, 1, -1},   {0, 1, 0},   {0, 1, 1}};

  /**
   * @brief current number of particles
   */
  size_t particle_count_{};

  /**
   * @brief number of particles, that are immovable
   */
  size_t special_particle_count_{};

  /**
   * @brief enables the neighbour calculation for membranes
   */
  bool is_membrane{};

  /**
   * @brief
   * stores the indexes of all halo_cells for faster iteration in the
   * corresponding direction vector
   */
  std::array<std::vector<std::size_t>, 6> halo_direction_cells_;

  /**
   * @brief
   * stores the indexes of all halo_cells for faster iteration in the
   * corresponding direction vector
   */
  std::array<std::vector<std::size_t>, 6> boundary_direction_cells_;

  /**
   * @brief
   * a more processing friendly storage of LinkedCellsConfig::BoundaryConfig
   */
  std::array<LinkedCellsConfig::BoundaryType, 6> boundaries_{};

  /**
   * @brief
   * number of cells per direction for domain + 2 halo cells
   */
  ivec3 cell_count_{};

  /**
   * @brief
   * cell dimensions
   */
  dvec3 cell_dim_{};

  /**
   * @brief
   * cutoff radius
   */
  double cutoff_{};

  /**
   * @brief
   * the domain of the container
   */
  ivec3 domain_{};

  /**
   * @brief
   * the boundary config of each direction of the simulation
   */
  LinkedCellsConfig::BoundaryConfig boundary_config_{};

  /**
   * @brief apply reflective boundary condition to a dimension
   * @param dimension the problematic dimension
   */
  inline void applyReflectiveBoundary(size_t dimension);

  IndexForce index_force{};

  /**
   *@brief index offsets orthogonal to a cell for each dimension, optimized for
   *2D simulations
   *
   */
  std::array<std::array<ivec3, 9>, 3> index_offsets_ = {{
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

  /**
 * @brief Add a particle to the container
 * @param p Particle to be added
 * @note Does not impose the invariant automatically!
 */
  void addParticle(Particle& p);

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
  [[nodiscard]] size_t getParticleCount() const { return particle_count_; }

  /**
   * @brief the exact number of current special particles, updated accordingly
   * @return the current count of special particles left in the simulation
   */
  [[nodiscard]] size_t getSpecialParticleCount() const {
    return special_particle_count_;
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

  void setIndexForce(const IndexForce& index_force);

  /**
   * applies the periodic boundary conditions to the given dimension
   * @param dimension the dimension that the periodic boundary should be applied
   * to
   */
  inline void applyPeriodicBoundary(size_t dimension);

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
   * @brief Compute interactive forces
   */
  void computeInteractiveForces(
      const std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces)
      override;

  /**
   * @brief Compute singular forces
   */
  void computeSingularForces(const std::vector<std::unique_ptr<SingularForce>>&
                                 singular_forces) override;

  /**
   * @brief Get the amount of cells in each dimension
   * @return ivec3 of cells in each dimension
   */
  [[nodiscard]] std::array<int, 3> getCellCount() const { return cell_count_; }

  /**
   * @brief Get the dimensions of a all cells in the container
   * @return dvec3 of the dimensions of all cells
   */
  [[nodiscard]] dvec3 getCellDim() const { return cell_dim_; }

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
  [[nodiscard]] std::size_t cellCoordToIndexTesting(ivec3 position) const;

  /**
   * @brief Gets the cell coordinate from the cell index
   * @param cell_index Index of the cell
   * @return Cell coodinate in 3 dimensions
   */
  [[nodiscard]] inline ivec3 cellIndexToCoord(std::size_t cell_index) const;

  /**
   * @brief API for testing
   * @param cell_index Index of the cell
   * @return Cell coodinate in 3 dimensions
   */
  [[nodiscard]] ivec3 cellIndexToCoordTesting(std::size_t cell_index) const;

  /**
   * @brief Checks if a cell coordinate exists in the container
   * @param coordinate Cell coordinate to be checked
   * @return If cell coordinate is inside the container
   */
  [[nodiscard]] inline bool isValidCellCoordinate(ivec3 coordinate) const;

  /**
   * @brief Checks if a cell coordinate is in the halo of the container
   * @param cell_coord Cell coordinate to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] inline bool isHalo(ivec3 cell_coord) const;

  /**
   * @brief Checks if a cell index is in the halo of the container
   * @param cell_index cell index to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] inline bool isHalo(std::size_t cell_index) const;

  /**
   * @brief API for testing
   * @param cell_index cell index to be checked
   * @return If cell is part of the halo
   */
  [[nodiscard]] bool isHaloTesting(std::size_t cell_index) const;

  /**
   * @brief Checks if a cell coordinate is in the boundary of the container
   * @param cell_coord Cell coordinate to be checked
   * @return cell is part of the boundary
   */
  [[nodiscard]] inline bool isBoundary(ivec3 cell_coord) const;

  /**
   * @brief Checks if a cell index is in the boundary of the container
   * @param cell_index cell index to be checked
   * @return If cell is part of the boundary
   */
  [[nodiscard]] inline bool isBoundary(std::size_t cell_index) const;

  /**
   * @brief API for testing, because gtest does not like inline
   * @param cell_index cell index to be checked
   * @return If cell is part of the boundary
   */
  [[nodiscard]] bool isBoundaryTesting(std::size_t cell_index) const;

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
   * @param cell_index cell index to be checked
   * @param f function to check if it is a special cell
   * @param lower_magic_number lower bound index
   * @param upper_magic_number upper bound index
   * @return the directions of the boundary cell
   */
  [[nodiscard]] std::vector<std::size_t> specialCellDirection(
      std::size_t cell_index, const std::function<bool(std::size_t)>& f,
      int lower_magic_number, int upper_magic_number) const;

  /**
   * @brief Debug method to get direct access to the cells vector
   * @return Reference to the cell vector
   */
  std::vector<std::vector<Particle*>>& getCells() { return cells_; }

  /**
   * @brief warp negative cell index to maximum cell coordinate to enable
   * multiple periodic boundaries in corners. For now this does only work in 2D
   * @param cell_coordinate the cell to be checked
   * @param raw_dimension the dimension axis looked at
   * @return bool: whether it is a valid cell to be checked, ivec3: real cell,
   * dvec3: offset to be applied
   */
  [[nodiscard]] inline std::tuple<bool, ivec3, dvec3> reflectiveWarpAround(
      ivec3 cell_coordinate, std::size_t raw_dimension) const;

  /**
   * API for testing
   * @param cell_coordinate the cell to be checked
   * @param raw_dimension the dimension axis looked at
   * @return bool: whether it is a valid cell to be checked, ivec3: real cell,
   * dvec3: offset to be applied
   */
  [[nodiscard]] std::tuple<bool, ivec3, dvec3> reflectiveWarpAroundTesting(
      ivec3 cell_coordinate, std::size_t raw_dimension) const;

  /**
   * @brief calculates the kinetic energy in the container
   * @return kinetic energy of all particles
   */
  double getKineticEnergy() override;

  /**
   * @brief returns particle count of the container
   * @return particle count
   */
  std::size_t getParticleCount() override { return particle_count_; }

  /**
   * @brief returns the count of immovable particles
   * @return count of immovable particles
   */
  std::size_t getSpecialParticleCount() override {
    return special_particle_count_;
  }

  /**
   * @brief returns the domain of the container
   * @return the domain of the container
   */
  ivec3 getDomain() override { return domain_; }

  /**
   * @brief if true, this corner should not be evaluated because it was already
   * done
   * @param cell_coordinate the coordinate of the cell
   * @param raw_dimension the evaluated dimension
   * @return if cell should be ignored
   */
  [[nodiscard]] inline bool isDoubleCorner(ivec3 cell_coordinate,
                                           std::size_t raw_dimension) const;

  /**
   * @brief since the neighbour references are invalid after adding particles,
   * set them again
   */
  void setNeighbourReferences();

  void initializeC18Schema() {
    INFO_FMT("Cells_size : {}", cells_.size());

    for (auto start_offset : c_18_schema_) {
      // std::vector<std::vector<Particle*>*> iterators;
      std::vector<std::size_t> iterators;
      for (int cx = start_offset[0]; cx <= cell_count_[0]; cx += 2) {
        for (int cy = start_offset[1]; cy <= cell_count_[1]; cy += 3) {
          for (int cz = start_offset[2]; cz <= cell_count_[2]; cz += 3) {
            // INFO_FMT("Cell index {} {} {}", cx, cy, cz);
            if (!isValidCellCoordinate({cx, cy, cz})) {
              continue;
            }
            auto cell_index = cellCoordToIndex({cx, cy, cz});
            iterators.push_back(cell_index);
            // iterators.push_back(&cells_.at(cellCoordToIndex({cx, cy, cz})));
          }
        }
      }
      c18_colours_.push_back(iterators);
    }
  }
};


/**
 * @brief directions for better readability; implicitly cast
 */
enum Directions { xlow, xhigh, ylow, yhigh, zlow, zhigh };