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
#include "utils/SpdWrapper.h"

const double LinkedCellsContainer::sigma_factor = std::pow(2.0, 1.0 / 6.0);

LinkedCellsContainer::LinkedCellsContainer(
    const LinkedCellsConfig &linked_cells_config) {
  domain = linked_cells_config.domain;

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

  // add 2 for halo
  cell_count = {cell_count[0] + 2, cell_count[1] + 2, cell_count[2] + 2};

  cells.resize(cell_count[0] * cell_count[1] * cell_count[2]);

  this->boundary_config = linked_cells_config.boundary_config;

  DEBUG_PRINT_FMT("Num Cells: {}", cells.size());

  halo_direction_cells = {};
  boundary_direction_cells = {};
  // precalculate special cells
  for (std::size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
    auto halo_directions = halo_direction(cell_index);
    auto boundary_directions = boundary_direction(cell_index);

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

  // TODO: pretty
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

  DEBUG_PRINT_FMT(
      "Removed particle with coords ({}, {}, {}) from cell index: {}",
      p.getX()[0], p.getX()[1], p.getX()[2], index)
}

std::vector<Particle *> LinkedCellsContainer::getParticles() {
  std::vector<Particle *> res;
  singleIterator([&res](Particle &p) { res.push_back(&p); });

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
      const std::size_t shouldBeIndex = dvec3ToCellIndex(it->getX());
      if (shouldBeIndex == index) {
        ++it;
        continue;
      }

      cells[shouldBeIndex].push_back(*it);
      it = cells[index].erase(it);
    }
  }

  // apply boundary condition
  // it is assumed that GhostParticles do not have to persist, so we dont have
  // to iterate over the halo cells of Reflective Boundaries
  // TODO: this is heavily inefficient in 2D if 6 is used, make dimension
  // accessible, 4 instead of 6

  // clear halo, 4 hardcoded due to 2D simulation
  // should not be needed to perform on Reflective Boundaries, but in case of an
  // error prevents the simulation from crashing
  for (size_t dimension = 0; dimension < 4; ++dimension) {
    for (const size_t cell_index : halo_direction_cells[dimension]) {
      cells[cell_index].clear();
    }
  }

  for (size_t dimension = 0; dimension < 4; ++dimension) {
    switch (boundaries[dimension]) {
      case LinkedCellsConfig::BoundaryType::Outflow: {
        break;
      }
      case LinkedCellsConfig::BoundaryType::Reflective: {
        // the slides do not state what to do if it is an edge, so we place a
        // particle for every dimension it is too close to

        // ensure that GhostParticle only interacts with specific particle
        // assumed: epsilon and sigma are the same as of the problematic
        // Particle, the cutoff is larger than half of sigma_factor * sigma
        const std::size_t problematic_dimension = dimension / 2;
        const std::size_t problematic_dimension_direction = dimension % 2;

        for (const std::size_t cell_index :
             boundary_direction_cells[dimension]) {
          for (Particle &p : cells[cell_index]) {
            // check if it is too close
            double pos = p.getX()[problematic_dimension];
            const double boundary_position = domain[problematic_dimension];
            const double double_dist_to_boundary =
                2 *
                std::min(pos, boundary_position -
                                  pos);  // if both of them are so small that
                                         // they would trigger the boundary, the
                                         // simulation itself is already broken

            if (double_dist_to_boundary < sigma_factor * p.getSigma()) {
              const double force =
                  LennardJones::simpleForce(p, double_dist_to_boundary);
              dvec3 p_force = p.getF();
              p_force[problematic_dimension] +=
                  force *
                  std::pow(
                      -1.0,
                      1 - problematic_dimension_direction);  // invert the
                                                             // direction for
                                                             // boundary in
                                                             // ascending
                                                             // coordinate
                                                             // direction
              p.setF(p_force);
              DEBUG_PRINT_FMT(
                  "Applied Force=[{}, {}, {}] to Particle at [{}, {}, {}]",
                  p.getF()[0], p.getF()[1], p.getF()[2], p.getX()[0],
                  p.getX()[1], p.getX()[2]);
            }
          }
        }
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
  for (std::size_t cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
    std::vector<Particle> &cellParticles = cells[cellIndex];

    if (cellParticles.empty()) continue;

    ivec3 cellCoordinate = cellIndexToCoord(cellIndex);
    DEBUG_PRINT_FMT("cell index: {}; coord = ({}, {}, {}); halo? = {}",
                    cellIndex, cellCoordinate[0], cellCoordinate[1],
                    cellCoordinate[2], isHalo(cellIndex));

    // iterate over particles inside cell
    for (std::size_t i = 0; i < cellParticles.size(); ++i) {
      for (std::size_t j = i + 1; j < cellParticles.size(); ++j) {
        const dvec3 p = cellParticles[i].getX();
        const dvec3 q = cellParticles[j].getX();
        if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
          continue;
        f(cellParticles[i], cellParticles[j]);
        DEBUG_PRINT_FMT("Intra cell pair: ({}, {})", cellParticles[i].getType(),
                        cellParticles[j].getType());
      }
    }

    // iterate over neighbouring particles
    for (auto &offset : offsets) {
      // compute neighbourIndex and check if it is valid
      const ivec3 neighbourCoord = {cellCoordinate[0] + offset[0],
                                    cellCoordinate[1] + offset[1],
                                    cellCoordinate[2] + offset[2]};

      if (!isValidCellCoordinate(neighbourCoord)) {
        DEBUG_PRINT_FMT("Invalid coord: ({}, {}, {})", neighbourCoord[0],
                        neighbourCoord[1], neighbourCoord[2])
        continue;
      }

      const size_t neighbourIndex = cellCoordToIndex(neighbourCoord);
      DEBUG_PRINT_FMT(
          "Checking cell i={}; c=({}, {}, {}) for pairs (offset = ({}, {}, "
          "{}))",
          neighbourIndex, neighbourCoord[0], neighbourCoord[1],
          neighbourCoord[2], offset[0], offset[1], offset[2]);

      // go over all pairs with neighbour particles
      std::vector<Particle> &neighbourParticles = cells[neighbourIndex];
      if (neighbourParticles.empty()) continue;

      for (auto &cellParticle : cellParticles) {
        for (auto &neighbourParticle : neighbourParticles) {
          auto p = cellParticle.getX();
          auto q = neighbourParticle.getX();

          if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
              d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
            continue;

          f(cellParticle, neighbourParticle);
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
  const std::array<int, 3> cellCoords = {
      static_cast<int>(std::floor(position[0] / cell_dim[0])),
      static_cast<int>(std::floor(position[1] / cell_dim[1])),
      static_cast<int>(std::floor(position[2] / cell_dim[2]))};

  return cellCoordToIndex(cellCoords);
}

inline std::size_t LinkedCellsContainer::cellCoordToIndex(
    const ivec3 position) const {
  return (position[0] + 1) * (cell_count[1] * cell_count[2]) +
         (position[1] + 1) * (cell_count[2]) + (position[2] + 1);
}

inline ivec3 LinkedCellsContainer::cellIndexToCoord(
    std::size_t cellIndex) const {
  const int x = static_cast<int>(cellIndex / (cell_count[1] * cell_count[2]));
  cellIndex = cellIndex - (x * cell_count[1] * cell_count[2]);

  const int y = static_cast<int>(cellIndex / cell_count[2]);
  const int z = static_cast<int>(cellIndex - (y * cell_count[2]));

  return {x - 1, y - 1, z - 1};
}

inline bool LinkedCellsContainer::isValidCellCoordinate(
    const ivec3 coordinate) const {
  return (-1 <= coordinate[0] && coordinate[0] <= (cell_count[0] - 2)) &&
         (-1 <= coordinate[1] && coordinate[1] <= (cell_count[1] - 2)) &&
         (-1 <= coordinate[2] && coordinate[2] <= (cell_count[2] - 2));
}

inline bool LinkedCellsContainer::isHalo(const ivec3 cellCoord) const {
  return cellCoord[0] == -1 || cellCoord[1] == -1 || cellCoord[2] == -1 ||
         cellCoord[0] == (cell_count[0] - 2) ||
         cellCoord[1] == (cell_count[1] - 2) ||
         cellCoord[2] == (cell_count[2] - 2);
}

inline bool LinkedCellsContainer::isHalo(const std::size_t cellIndex) const {
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);
  return isHalo(cellCoord);
}

inline bool LinkedCellsContainer::isBoundary(const ivec3 cellCoord) const {
  return (cellCoord[0] == 0 || cellCoord[1] == 0 || cellCoord[2] == 0 ||
          cellCoord[0] == (cell_count[0] - 3) ||
          cellCoord[1] == (cell_count[1] - 3) ||
          cellCoord[2] == (cell_count[2] - 3)) &&
         !isHalo(cellCoord);
}

inline bool LinkedCellsContainer::isBoundary(
    const std::size_t cellIndex) const {
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);
  return isBoundary(cellCoord);
}

std::vector<std::size_t> LinkedCellsContainer::halo_direction(
    const std::size_t cellIndex) const {
  if (!isHalo(cellIndex)) return {};

  std::vector<std::size_t> directions = {};
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);

  if (cellCoord[0] == -1) {
    directions.push_back(0);  // west
  }
  if (cellCoord[0] == (cell_count[0] - 2)) {
    directions.push_back(1);  // east
  }
  if (cellCoord[1] == -1) {
    directions.push_back(2);  // down
  }
  if (cellCoord[1] == (cell_count[1] - 2)) {
    directions.push_back(3);  // up
  }
  if (cellCoord[2] == -1) {
    directions.push_back(4);  // south
  }
  if (cellCoord[2] == (cell_count[2] - 2)) {
    directions.push_back(5);  // north
  }

  return directions;
}

std::vector<std::size_t> LinkedCellsContainer::boundary_direction(
    const std::size_t cellIndex) const {
  if (!isBoundary(cellIndex)) return {};

  std::vector<std::size_t> directions = {};
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);

  if (cellCoord[0] == 0) {
    directions.push_back(0);  // west
  }
  if (cellCoord[0] == (cell_count[0] - 3)) {
    directions.push_back(1);  // east
  }
  if (cellCoord[1] == 0) {
    directions.push_back(2);  // down
  }
  if (cellCoord[1] == (cell_count[1] - 3)) {
    directions.push_back(3);  // up
  }
  if (cellCoord[2] == 0) {
    directions.push_back(4);  // south
  }
  if (cellCoord[2] == (cell_count[2] - 3)) {
    directions.push_back(5);  // north
  }

  return directions;
}

bool LinkedCellsContainer::isBoundary_testing(
    const std::size_t cellIndex) const {
  return isBoundary(cellIndex);
}

bool LinkedCellsContainer::isHalo_testing(const std::size_t cellIndex) const {
  return isHalo(cellIndex);
}

std::size_t LinkedCellsContainer::dvec3ToCellIndex_testing(
    const dvec3 &position) const {
  return dvec3ToCellIndex(position);
}
