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
      SpdWrapper::get()->info(
          "cell count is too small if reflective boundaries are used");
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
    auto halo_directions = special_cell_direction(
        cell_index, [this](const std::size_t index) { return isHalo(index); },
        -1, 2);
    auto boundary_directions = special_cell_direction(
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
      const std::size_t shouldBeIndex = dvec3ToCellIndex(it->getX());
      if (shouldBeIndex == index) {
        ++it;
        continue;
      }
      if (it->getX()[0] >= 100 || it->getF()[0] >= 10) {
        SpdWrapper::get()->info("ALAAAARM!");
        SpdWrapper::get()->info("x = [{}, {}, {}], v = [{}, {}, {}]", it->getX()[0], it->getX()[1], it->getX()[2], it->getV()[0], it->getX()[1], it->getV()[2]);
      }

      ivec3 current_cell = cellIndexToCoord(index);
      ivec3 new_cell = cellIndexToCoord(shouldBeIndex);
      cells[shouldBeIndex].push_back(*it);
      it = cells[index].erase(it);
    }
  }

  // apply boundary condition
  // it is assumed that GhostParticles do not have to persist, so we dont have
  // to iterate over the halo cells of Reflective Boundaries
  // TODO: this is heavily inefficient in 2D if 6 is used, make dimension
  // accessible, 4 instead of 6

  for (size_t dimension = 0; dimension < 4; ++dimension) {
    switch (boundaries[dimension]) {
      case LinkedCellsConfig::BoundaryType::Outflow: {
        // clear halo
        for (const size_t cell_index : halo_direction_cells[dimension]) {
          cells[cell_index].clear();
          cells[cell_index].shrink_to_fit();
        }
        break;
      }
      case LinkedCellsConfig::BoundaryType::Reflective: {
        apply_reflective_boundary(dimension);
        break;
      }
      case LinkedCellsConfig::Periodic: {
        // move particles in halo to the other side
        // calculate forces only for ?_high, so all particles are until then in
        // the right place
        const std::size_t problematic_dimension = dimension / 2;
        const std::size_t problematic_dimension_direction = dimension % 2;

        for (const std::size_t cell_index : halo_direction_cells[dimension]) {
          int counter = 0;
          for (auto it = cells[cell_index].begin();
               it < cells[cell_index].end(); ++it) {
            counter++;
            ivec3 c = cellIndexToCoord(cell_index);
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
          break;
        }

        // iterate over all 9 / 3 cells on the other end
        // for now strict 2D implementation for performance
        for (const std::size_t cell_index :
             boundary_direction_cells[dimension]) {
          ivec3 cell_coordinates = cellIndexToCoord(cell_index);

          // change 3 to 9 for 3D
          for (std::size_t i = 0; i < 3; ++i) {
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

            const auto adjacent_cell_index = cellCoordToIndex(adjacent_cell_coordinates);

            // // account for the dimension that is checked
            // particle_distance_offset[problematic_dimension] =
            //     domain[problematic_dimension];

            // iterate over all pairs and calculate force
            for (Particle &p : cells[cell_index]) {
              for (Particle &q : cells[adjacent_cell_index]) {
                // distance check
                const dvec3 accounted_particle_distance =
                    q.getX() - p.getX() + particle_distance_offset;

                if (ArrayUtils::squaredL2Norm(accounted_particle_distance) >=
                    cutoff * cutoff) {
                  continue;
                }

                const dvec3 applied_force =
                    LennardJones::directionalForceWithOffset(
                        p, q, accounted_particle_distance);
                p.setF(p.getF() + applied_force);
                q.setF(q.getF() - applied_force);
              }
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

std::vector<std::size_t> LinkedCellsContainer::special_cell_direction(
    const std::size_t cellIndex, const std::function<bool(std::size_t)> &f,
    const int lowerMagicNumber, const int upperMagicNumber) const {
  if (!f(cellIndex)) return {};

  std::vector<std::size_t> directions = {};
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);

  if (cellCoord[0] == lowerMagicNumber) {
    directions.push_back(xlow);  // west
  }
  if (cellCoord[0] == (cell_count[0] - upperMagicNumber)) {
    directions.push_back(xhigh);  // east
  }
  if (cellCoord[1] == lowerMagicNumber) {
    directions.push_back(ylow);  // down
  }
  if (cellCoord[1] == (cell_count[1] - upperMagicNumber)) {
    directions.push_back(yhigh);  // up
  }
  if (cellCoord[2] == lowerMagicNumber) {
    directions.push_back(zlow);  // south
  }
  if (cellCoord[2] == (cell_count[2] - upperMagicNumber)) {
    directions.push_back(zhigh);  // north
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

void LinkedCellsContainer::apply_reflective_boundary(const size_t dimension) {
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

std::tuple<bool, ivec3, dvec3> LinkedCellsContainer::reflective_warp_around(
    const ivec3 cell_coordinate, const std::size_t raw_dimension) const {
  // TODO: this is really bad code
  dvec3 offset = {0, 0, 0};

  if (raw_dimension == yhigh &&
      boundaries[xlow] == LinkedCellsConfig::Periodic &&
      boundaries[ylow] == LinkedCellsConfig::Periodic &&
      (cell_coordinate[0] == -1 || cell_coordinate[0] == cell_count[1] - 2)) {
    // both dimensions are periodic -> make sure that corner cells are valid
    // only once! skip this for the y dimension, since it was already calculated
    // in the x dimension

    SpdWrapper::get()->info("cell should not be warped");
    return std::make_tuple(false, cell_coordinate, offset);
  }

  ivec3 new_cell_coordinate = cell_coordinate;
  for (std::size_t dimension = 0; dimension < 2; dimension++) {
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
LinkedCellsContainer::reflective_warp_around_testing(
    const ivec3 cell_coordinate, const std::size_t raw_dimension) const {
  return reflective_warp_around(cell_coordinate, raw_dimension);
}

std::size_t LinkedCellsContainer::cellCoordToIndex_testing(
    const ivec3 position) const {
  return cellCoordToIndex(position);
}

ivec3 LinkedCellsContainer::cellIndexToCoord_testing(
    const std::size_t cellIndex) const {
  return cellIndexToCoord(cellIndex);
}

double LinkedCellsContainer::getKineticEnergy() {
  double E_kin = 0.0;
  singleIterator([&E_kin](const Particle &p) {
    E_kin += 0.5 * p.getM() * ArrayUtils::L2InnerProduct(p.getV());
  });
  return E_kin;
}
