//
// created by mcarn on 11/15/24
//
#include "LinkedCellsContainer.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "debug/debug_print.h"
#include "defs/Particle.h"
#include "forces/LennardJones.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

const double LinkedCellsContainer::sigma_factor = std::pow(2.0, 1.0 / 6.0);

LinkedCellsContainer::LinkedCellsContainer(
    const LinkedCellsConfig &linked_cells_config) {
  domain = linked_cells_config.domain;
  particle_count = 0;
  special_particle_count = 0;

  DEBUG_PRINT("LinkedCellsContainer instantiated");
  SpdWrapper::get()->info("domain size: ({}, {}, {})", domain[0], domain[1],
                          domain[2]);

  cells = {};
  this->cutoff = linked_cells_config.cutoff_radius;

  cell_count = {std::max(static_cast<int>(std::floor(domain[0] / cutoff)), 1),
                std::max(static_cast<int>(std::floor(domain[1] / cutoff)), 1),
                std::max(static_cast<int>(std::floor(domain[2] / cutoff)), 1)};

  cell_dim = {static_cast<double>(domain[0]) / cell_count[0],
              static_cast<double>(domain[1]) / cell_count[1],
              static_cast<double>(domain[2]) / cell_count[2]};

  // safety check that minimum cell count is satisfied so the boundaries work as
  // expected:
  // could be ok, but then reflective calculation has to be disabled, since all
  // possible pairs are already iterated over
  for (std::size_t i = 0; i < 2; i++) {
    if (cell_count[i] < 3) {
      SpdWrapper::get()->error(
          "Cell count is too small if reflective boundaries are used! If this "
          "is not a testing instance, please exit the simulation");
    }
  }

  // add 2 for halo
  cell_count = {cell_count[0] + 2, cell_count[1] + 2, cell_count[2] + 2};

  cells.resize(cell_count[0] * cell_count[1] * cell_count[2]);

  this->boundary_config = linked_cells_config.boundary_config;

  DEBUG_PRINT_FMT("Num Cells: {}", cells.size());

  halo_direction_cells = {};
  boundary_direction_cells = {};
  // precalculate special cells
  for (std::size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
    auto halo_directions = specialCellDirection(
        cell_index, [this](const std::size_t index) { return isHalo(index); },
        -1, 2);
    auto boundary_directions = specialCellDirection(
        cell_index,
        [this](const std::size_t index) { return isBoundary(index); }, 0, 3);

    if (!halo_directions.empty()) {
      for (const unsigned long halo_direction : halo_directions) {
        halo_direction_cells[halo_direction].push_back(cell_index);
      }
    }
    if (!boundary_directions.empty()) {
      for (const unsigned long boundary_direction : boundary_directions) {
        boundary_direction_cells[boundary_direction].push_back(cell_index);
      }
    }
  }

  this->boundaries = {
      linked_cells_config.boundary_config.x_low,
      linked_cells_config.boundary_config.x_high,
      linked_cells_config.boundary_config.y_low,
      linked_cells_config.boundary_config.y_high,
      linked_cells_config.boundary_config.z_low,
      linked_cells_config.boundary_config.z_high,
  };
  SpdWrapper::get()->info("cell dim: {}, {}, {}; cell count: {}, {}, {}",
                          cell_dim[0], cell_dim[1], cell_dim[2], cell_count[0],
                          cell_count[1], cell_count[2]);
}

void LinkedCellsContainer::addParticle(const Particle &p) {
  const std::size_t index = dvec3ToCellIndex(p.getX());
  if (!isValidCellCoordinate(cellIndexToCoord(index))) {
    SpdWrapper::get()->error("Tried to add particle out of bounds");
    exit(1);
  }
  cells[index].emplace_back(p);

  this->particle_count++;

  if (p.getType() < 0) {
    this->special_particle_count++;
  }

  DEBUG_PRINT_FMT("Added particle with coords ({}, {}, {}) into cell index: {}",
                  p.getX()[0], p.getX()[1], p.getX()[2], index)
}

void LinkedCellsContainer::addParticles(
    const std::vector<Particle> &particles) {
  for (const Particle &p : particles) {
    addParticle(p);
  }
}

void LinkedCellsContainer::removeParticle(const Particle &p) {
  const std::size_t index = dvec3ToCellIndex(p.getX());
  std::vector<Particle> &particles = cells[index];

  particles.erase(std::remove_if(particles.begin(), particles.end(),
                                 [&p](const Particle &q) { return p == q; }),
                  particles.end());

  this->particle_count--;
  if (p.getType() < 0) {
    this->special_particle_count--;
  }

  DEBUG_PRINT_FMT(
      "Removed particle with coords ({}, {}, {}) from cell index: {}",
      p.getX()[0], p.getX()[1], p.getX()[2], index)
}

std::vector<Particle *> LinkedCellsContainer::getParticles() {
  std::vector<Particle *> res;
  singleIterator([&res](Particle &p) { res.push_back(&p); });

  return res;
}

std::vector<Particle> LinkedCellsContainer::getParticlesObjects() {
  std::vector<Particle> res;
  singleIterator([&res](const Particle &p) { res.push_back(p); });
  return res;
}

[[nodiscard]] std::size_t LinkedCellsContainer::size() const {
  std::size_t count = 0;
  for (auto &c : cells) {
    count += c.size();
  }
  return count;
}

void LinkedCellsContainer::imposeInvariant() {
  // register in corresponding cell
  for (std::size_t index = 0; index < cells.size(); index++) {
    for (auto it = cells[index].begin(); it < cells[index].end();) {
      const std::size_t should_be_index = dvec3ToCellIndex(it->getX());
      if (should_be_index == index) {
        ++it;
        continue;
      }

      cells[should_be_index].push_back(*it);
      it = cells[index].erase(it);
    }
  }

  // apply boundary condition
  // it is assumed that GhostParticles do not have to persist, so we dont have
  // to iterate over the halo cells of Reflective Boundaries

  for (size_t dimension = 0; dimension < 6; ++dimension) {
    switch (boundaries[dimension]) {
      case LinkedCellsConfig::BoundaryType::Outflow: {
        // clear halo
        for (const size_t cell_index : halo_direction_cells[dimension]) {
          particle_count -=
              cells[cell_index].size();  // update particle count, only place
                                         // where particles are deleted
          cells[cell_index].clear();
          cells[cell_index].shrink_to_fit();
        }
        break;
      }
      case LinkedCellsConfig::BoundaryType::Reflective: {
        applyReflectiveBoundary(dimension);
        break;
      }
      case LinkedCellsConfig::Periodic: {
        applyPeriodicBoundary(dimension);
        break;
      }
      default: {
        DEBUG_PRINT_FMT("BoundaryType {} for dimension {} unknown",
                        static_cast<int>(boundaries[dimension]), dimension);
        break;
      }
    }
  }
}

void LinkedCellsContainer::singleIterator(
    const std::function<void(Particle &)> &f) {
  for (auto &c : cells) {
    for (auto &p : c) {
      f(p);
    }
  }
}

void LinkedCellsContainer::pairIterator(
    const std::function<void(Particle &, Particle &)> &f) {
  // - as x, y, z are all increasing offsets point to all neighbours in
  // positive directions
  // - for better cache usage (in flattened version) minimize max distance
  //   between values
  // - direction of traversal: z, y, x
  const std::array<ivec3, 13> offsets = {{
      // 9 x facing
      {{1, -1, -1}},
      {{1, -1, 0}},
      {{1, -1, 1}},
      {{1, 0, -1}},
      {{1, 0, 0}},
      {{1, 0, 1}},
      {{1, 1, -1}},
      {{1, 1, 0}},
      {{1, 1, 1}},
      // 3 y
      {{0, 1, -1}},
      {{0, 1, 0}},
      {{0, 1, 1}},
      // last z
      {{0, 0, 1}},
  }};

  // go over all cell indices
  for (std::size_t cell_index = 0; cell_index < cells.size(); cell_index++) {
    std::vector<Particle> &cell_particles = cells[cell_index];

    if (cell_particles.empty()) continue;

    ivec3 cellCoordinate = cellIndexToCoord(cell_index);
    DEBUG_PRINT_FMT("cell index: {}; coord = ({}, {}, {}); halo? = {}",
                    cellIndex, cellCoordinate[0], cellCoordinate[1],
                    cellCoordinate[2], isHalo(cellIndex));

    // iterate over particles inside cell
    for (std::size_t i = 0; i < cell_particles.size(); ++i) {
      for (std::size_t j = i + 1; j < cell_particles.size(); ++j) {
        const dvec3 p = cell_particles[i].getX();
        const dvec3 q = cell_particles[j].getX();
        if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
          continue;
        f(cell_particles[i], cell_particles[j]);
        DEBUG_PRINT_FMT("Intra cell pair: ({}, {})", cellParticles[i].getType(),
                        cellParticles[j].getType());
      }
    }

    // iterate over neighbouring particles
    for (auto &offset : offsets) {
      // compute neighbourIndex and check if it is valid
      const ivec3 neighbour_coord = {cellCoordinate[0] + offset[0],
                                    cellCoordinate[1] + offset[1],
                                    cellCoordinate[2] + offset[2]};

      if (!isValidCellCoordinate(neighbour_coord)) {
        DEBUG_PRINT_FMT("Invalid coord: ({}, {}, {})", neighbourCoord[0],
                        neighbourCoord[1], neighbourCoord[2])
        continue;
      }

      const size_t neighbour_index = cellCoordToIndex(neighbour_coord);
      DEBUG_PRINT_FMT(
          "Checking cell i={}; c=({}, {}, {}) for pairs (offset = ({}, {}, "
          "{}))",
          neighbourIndex, neighbourCoord[0], neighbourCoord[1],
          neighbourCoord[2], offset[0], offset[1], offset[2]);

      // go over all pairs with neighbour particles
      std::vector<Particle> &neighbour_particles = cells[neighbour_index];
      if (neighbour_particles.empty()) continue;

      for (auto &cell_particle : cell_particles) {
        for (auto &neighbour_particle : neighbour_particles) {
          auto p = cell_particle.getX();
          auto q = neighbour_particle.getX();

          if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
              d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
            continue;

          f(cell_particle, neighbour_particle);
          DEBUG_PRINT_FMT("Cross cell pair: ({}, {})", cellParticle.getType(),
                          neighbourParticle.getType())
        }
      }
    }
  }
}

void LinkedCellsContainer::boundaryIterator(
    const std::function<void(Particle &)> &f) {
  for (std::size_t index = 0; index < cells.size(); index++) {
    if (!isBoundary(index)) continue;

    for (auto &p : cells[index]) {
      f(p);
    }
  }
}

void LinkedCellsContainer::haloIterator(
    const std::function<void(Particle &)> &f) {
  for (std::size_t index = 0; index < cells.size(); index++) {
    if (!isHalo(index)) continue;

    for (auto &p : cells[index]) {
      f(p);
    }
  }
}

inline std::size_t LinkedCellsContainer::dvec3ToCellIndex(
    const dvec3 &position) const {
  const std::array<int, 3> cell_coords = {
      static_cast<int>(std::floor(position[0] / cell_dim[0])),
      static_cast<int>(std::floor(position[1] / cell_dim[1])),
      static_cast<int>(std::floor(position[2] / cell_dim[2]))};

  return cellCoordToIndex(cell_coords);
}

inline std::size_t LinkedCellsContainer::cellCoordToIndex(
    const ivec3 position) const {
  return (position[0] + 1) * (cell_count[1] * cell_count[2]) +
         (position[1] + 1) * (cell_count[2]) + (position[2] + 1);
}

inline ivec3 LinkedCellsContainer::cellIndexToCoord(
    std::size_t cell_index) const {
  const int x = static_cast<int>(cell_index / (cell_count[1] * cell_count[2]));
  cell_index = cell_index - (x * cell_count[1] * cell_count[2]);

  const int y = static_cast<int>(cell_index / cell_count[2]);
  const int z = static_cast<int>(cell_index - (y * cell_count[2]));

  return {x - 1, y - 1, z - 1};
}

inline bool LinkedCellsContainer::isValidCellCoordinate(
    const ivec3 coordinate) const {
  return (-1 <= coordinate[0] && coordinate[0] <= (cell_count[0] - 2)) &&
         (-1 <= coordinate[1] && coordinate[1] <= (cell_count[1] - 2)) &&
         (-1 <= coordinate[2] && coordinate[2] <= (cell_count[2] - 2));
}

inline bool LinkedCellsContainer::isHalo(const ivec3 cell_coord) const {
  return cell_coord[0] == -1 || cell_coord[1] == -1 || cell_coord[2] == -1 ||
         cell_coord[0] == (cell_count[0] - 2) ||
         cell_coord[1] == (cell_count[1] - 2) ||
         cell_coord[2] == (cell_count[2] - 2);
}

inline bool LinkedCellsContainer::isHalo(const std::size_t cell_index) const {
  const ivec3 cell_coord = cellIndexToCoord(cell_index);
  return isHalo(cell_coord);
}

inline bool LinkedCellsContainer::isBoundary(const ivec3 cell_coord) const {
  return (cell_coord[0] == 0 || cell_coord[1] == 0 || cell_coord[2] == 0 ||
          cell_coord[0] == (cell_count[0] - 3) ||
          cell_coord[1] == (cell_count[1] - 3) ||
          cell_coord[2] == (cell_count[2] - 3)) &&
         !isHalo(cell_coord);
}

inline bool LinkedCellsContainer::isBoundary(
    const std::size_t cell_index) const {
  const ivec3 cell_coord = cellIndexToCoord(cell_index);
  return isBoundary(cell_coord);
}

std::vector<std::size_t> LinkedCellsContainer::specialCellDirection(
    const std::size_t cell_index, const std::function<bool(std::size_t)> &f,
    const int lower_magic_number, const int upper_magic_number) const {
  if (!f(cell_index)) return {};

  std::vector<std::size_t> directions = {};
  const ivec3 cell_coord = cellIndexToCoord(cell_index);

  if (cell_coord[0] == lower_magic_number) {
    directions.push_back(xlow);  // west
  }
  if (cell_coord[0] == (cell_count[0] - upper_magic_number)) {
    directions.push_back(xhigh);  // east
  }
  if (cell_coord[1] == lower_magic_number) {
    directions.push_back(ylow);  // down
  }
  if (cell_coord[1] == (cell_count[1] - upper_magic_number)) {
    directions.push_back(yhigh);  // up
  }
  if (cell_coord[2] == lower_magic_number) {
    directions.push_back(zlow);  // south
  }
  if (cell_coord[2] == (cell_count[2] - upper_magic_number)) {
    directions.push_back(zhigh);  // north
  }

  return directions;
}

bool LinkedCellsContainer::isBoundaryTesting(
    const std::size_t cell_index) const {
  return isBoundary(cell_index);
}

bool LinkedCellsContainer::isHaloTesting(const std::size_t cell_index) const {
  return isHalo(cell_index);
}

std::size_t LinkedCellsContainer::dvec3ToCellIndex_testing(
    const dvec3 &position) const {
  return dvec3ToCellIndex(position);
}

void LinkedCellsContainer::applyReflectiveBoundary(const size_t dimension) {
  const std::size_t problematic_dimension = dimension / 2;
  const std::size_t problematic_dimension_direction = dimension % 2;
  // ensure that GhostParticle only interacts with specific particle
  // assumed: epsilon and sigma are the same as of the problematic
  // Particle, the cutoff is larger than half of sigma_factor * sigma
  for (const std::size_t cell_index : boundary_direction_cells[dimension]) {
    for (Particle &p : cells[cell_index]) {
      // check if it is too close
      double pos = p.getX()[problematic_dimension];
      const double boundary_position = domain[problematic_dimension];
      const double double_dist_to_boundary =
          2 * std::min(pos, boundary_position -
                                pos);  // if both of them are so small that
      // they would trigger the boundary, the
      // simulation itself is already broken

      if (double_dist_to_boundary < sigma_factor * p.getSigma()) {
        const double force =
            LennardJones::simpleForce(p, double_dist_to_boundary);
        dvec3 p_force = p.getF();
        p_force[problematic_dimension] +=
            force *
            std::pow(-1.0,
                     1 - problematic_dimension_direction);  // invert the
        // direction for
        // boundary in
        // ascending
        // coordinate
        // direction
        p.setF(p_force);
        DEBUG_PRINT_FMT(
            "Applied Force=[{}, {}, {}] to Particle at [{}, {}, {}]",
            p.getF()[0], p.getF()[1], p.getF()[2], p.getX()[0], p.getX()[1],
            p.getX()[2]);
      }
    }
  }
}

void LinkedCellsContainer::applyPeriodicBoundary(const size_t dimension) {
  // move particles in halo to the other side
  // calculate forces only for ?_high, so all particles are until then in
  // the right place
  const std::size_t problematic_dimension = dimension / 2;
  const std::size_t problematic_dimension_direction = dimension % 2;

  for (const std::size_t cell_index : halo_direction_cells[dimension]) {
    int counter = 0;
    for (auto it = cells[cell_index].begin(); it < cells[cell_index].end();
         ++it) {
      counter++;
      dvec3 new_pos = it->getX();
      new_pos[problematic_dimension] +=
          domain[problematic_dimension] *
          (problematic_dimension_direction % 2 == 0 ? 1 : -1);
      const std::size_t shouldBeIndex = dvec3ToCellIndex(new_pos);

      it->setX(new_pos);
      cells[shouldBeIndex].push_back(*it);
    }

    cells[cell_index].clear();
    cells[cell_index].shrink_to_fit();
  }

  // skip force calculation for lower side of the axis
  if (problematic_dimension_direction == 0) {
    return;
  }

  // iterate over all 9 / 3 cells on the other end
  for (const std::size_t cell_index : boundary_direction_cells[dimension]) {
    ivec3 cell_coordinates = cellIndexToCoord(cell_index);

    // change 3 to 9 for 3D
    for (std::size_t i = 0; i < 9; ++i) {
      ivec3 offset = index_offsets[problematic_dimension][i];
      const ivec3 cell_to_check = cell_coordinates + offset;
      bool is_adjacent_cell;
      ivec3 adjacent_cell_coordinates;
      dvec3 particle_distance_offset;

      std::tie(is_adjacent_cell, adjacent_cell_coordinates,
               particle_distance_offset) =
          reflective_warp_around(cell_to_check, dimension);

      if (!is_adjacent_cell) {
        continue;
      }

      const auto adjacent_cell_index =
          cellCoordToIndex(adjacent_cell_coordinates);

      // // account for the dimension that is checked
      // particle_distance_offset[problematic_dimension] =
      //     domain[problematic_dimension];

      // iterate over all pairs and calculate force

      for (Particle &p : cells[cell_index]) {
        for (Particle &q : cells[adjacent_cell_index]) {
          // distance check
          const dvec3 accounted_particle_distance =
              q.getX() - p.getX() + particle_distance_offset;

          if (ArrayUtils::L2InnerProduct(accounted_particle_distance) >=
              cutoff * cutoff) {
            continue;
          }

          const dvec3 applied_force = LennardJones::directionalForceWithOffset(
              p, q, accounted_particle_distance);
          p.setF(p.getF() + applied_force);
          q.setF(q.getF() - applied_force);
        }
      }
    }
  }
}

std::tuple<bool, ivec3, dvec3> LinkedCellsContainer::reflective_warp_around(
    const ivec3 cell_coordinate, const std::size_t raw_dimension) const {
  dvec3 offset = {0, 0, 0};

  // if (raw_dimension == yhigh &&
  //     boundaries[xlow] == LinkedCellsConfig::Periodic &&
  //     boundaries[ylow] == LinkedCellsConfig::Periodic &&
  //     (cell_coordinate[0] == -1 || cell_coordinate[0] == cell_count[1] - 2))
  //     {
  //   // both dimensions are periodic -> make sure that corner cells are valid
  //   // only once! skip this for the y dimension, since it was already
  //   calculated
  //   // in the x dimension
  //
  //   DEBUG_PRINT("cell should not be warped");
  //   return std::make_tuple(false, cell_coordinate, offset);
  // }

  if (isDoubleCorner(cell_coordinate, raw_dimension)) {
    DEBUG_PRINT("cell should not be warped");
    return std::make_tuple(false, cell_coordinate, offset);
  }

  ivec3 new_cell_coordinate = cell_coordinate;
  for (std::size_t dimension = 0; dimension < 3; dimension++) {
    if (cell_coordinate[dimension] == -1) {
      // low wrap around to high cell
      new_cell_coordinate[dimension] = cell_count[dimension] - 3;  // top cell
      offset[dimension] = -domain[dimension];
    } else if (cell_coordinate[dimension] == cell_count[dimension] - 2) {
      // high warp around to low cell
      new_cell_coordinate[dimension] = 0;  // bottom cell
      offset[dimension] = domain[dimension];
    } else {
      // no warp around, nothing can go wrong
      continue;
    }

    // if it is wrapped around but dimension is not periodic, this cell is not
    // adjacent
    if (boundaries[2 * dimension] != LinkedCellsConfig::Periodic) {
      return std::make_tuple(false, cell_coordinate, offset);
    }
  }

  DEBUG_PRINT_FMT("[{}, {}, {}] -> [{}, {}, {}] | [{}, {}, {}]",
                  cell_coordinate[0], cell_coordinate[1], cell_coordinate[2],
                  new_cell_coordinate[0], new_cell_coordinate[1],
                  new_cell_coordinate[2], offset[0], offset[1], offset[2]);
  return std::make_tuple(true, new_cell_coordinate, offset);
}

std::tuple<bool, ivec3, dvec3>
LinkedCellsContainer::reflectiveWarpAroundTesting(
    const ivec3 cell_coordinate, const std::size_t raw_dimension) const {
  return reflectiveWarpAround(cell_coordinate, raw_dimension);
}

std::size_t LinkedCellsContainer::cellCoordToIndexTesting(
    const ivec3 position) const {
  return cellCoordToIndex(position);
}

ivec3 LinkedCellsContainer::cellIndexToCoordTesting(
    const std::size_t cell_index) const {
  return cellIndexToCoord(cell_index);
}

double LinkedCellsContainer::getKineticEnergy() {
  double e_kin = 0.0;
  singleIterator([&e_kin](const Particle &p) {
    e_kin += p.getM() * ArrayUtils::L2InnerProduct(p.getV());
  });
  return e_kin * 0.5;
}

bool LinkedCellsContainer::isDoubleCorner(
    const ivec3 cell_coordinate, const std::size_t raw_dimension) const {
  // check whether it really is a corner
  int edge_of_x_dimensions_counter = 0;
  for (std::size_t dimension = 0; dimension < 3; dimension++) {
    if (cell_coordinate[dimension] == -1 ||
        cell_coordinate[dimension] == cell_count[dimension] - 2) {
      edge_of_x_dimensions_counter++;
    }
  }

  if (edge_of_x_dimensions_counter < 2) {
    // SpdWrapper::get()->info("[{}, {}, {}] is not a corner!", cell_coordinate[0],
    //                       cell_coordinate[1], cell_coordinate[2]);
    return false;
  }
  // SpdWrapper::get()->info("[{}, {}, {}] is a corner!", cell_coordinate[0],
  //                         cell_coordinate[1], cell_coordinate[2]);

  // TODO: make beautiful
  // lookup table via some if statements
  // ugly, but best I could think of
  if (boundaries[xlow] == LinkedCellsConfig::Periodic &&
      boundaries[ylow] == LinkedCellsConfig::Periodic &&
      boundaries[zlow] != LinkedCellsConfig::Periodic) {
    // 110
    // SpdWrapper::get()->info("isDoubleCorner - xlow ylow !zlow");
    if (raw_dimension == yhigh &&
        (cell_coordinate[0] == -1 || cell_coordinate[0] == cell_count[0] - 2)) {
      return true;
    }
    return false;
  }
  if (boundaries[xlow] == LinkedCellsConfig::Periodic &&
      boundaries[ylow] != LinkedCellsConfig::Periodic &&
      boundaries[zlow] == LinkedCellsConfig::Periodic) {
    // 101
    // SpdWrapper::get()->info("isDoubleCorner - xlow !ylow zlow");
    if (raw_dimension == zhigh &&
        (cell_coordinate[0] == -1 || cell_coordinate[0] == cell_count[0] - 2)) {
      return true;
    }
    return false;
  }
  if (boundaries[xlow] != LinkedCellsConfig::Periodic &&
      boundaries[ylow] == LinkedCellsConfig::Periodic &&
      boundaries[zlow] == LinkedCellsConfig::Periodic) {
    // 011
    // SpdWrapper::get()->info("isDoubleCorner - !xlow ylow zlow");
    if (raw_dimension == zhigh &&
        (cell_coordinate[1] == -1 || cell_coordinate[1] == cell_count[1] - 2)) {
      return true;
    }
    return false;
  }
  if (boundaries[xlow] == LinkedCellsConfig::Periodic &&
      boundaries[ylow] == LinkedCellsConfig::Periodic &&
      boundaries[zlow] == LinkedCellsConfig::Periodic) {
    // 111
    // SpdWrapper::get()->info("isDoubleCorner - xlow ylow zlow");
    if (raw_dimension == yhigh &&
        (cell_coordinate[0] == -1 || cell_coordinate[0] == cell_count[0] - 2)) {
      return true;
    }
    if (raw_dimension == zhigh &&
        (cell_coordinate[0] == -1 || cell_coordinate[0] == cell_count[0] - 2)) {
      return true;
    }
    if (raw_dimension == zhigh &&
        (cell_coordinate[1] == -1 || cell_coordinate[1] == cell_count[1] - 2)) {
      return true;
    }
    return false;
  }
  return false;
}
