//
// created by mcarn on 11/15/24
//
#include "LinkedCellsContainer.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "debug/debug_print.h"
#include "defs/Particle.h"
#include "forces/LennardJones.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

const double LinkedCellsContainer::sigma_factor = std::pow(2.0, 1.0 / 6.0);

LinkedCellsContainer::LinkedCellsContainer(
    const LinkedCellsConfig &linked_cells_config) {
  domain_ = linked_cells_config.domain;
  particle_count_ = 0;
  special_particle_count_ = 0;

  DEBUG_PRINT("LinkedCellsContainer instantiated");
  SpdWrapper::get()->info("domain size: ({}, {}, {})", domain_[0], domain_[1],
                          domain_[2]);

  cells_ = {};
  this->cutoff_ = linked_cells_config.cutoff_radius;
  this->is_membrane = linked_cells_config.is_membrane;

  cell_count_ = {
      std::max(static_cast<int>(std::floor(domain_[0] / cutoff_)), 1),
      std::max(static_cast<int>(std::floor(domain_[1] / cutoff_)), 1),
      std::max(static_cast<int>(std::floor(domain_[2] / cutoff_)), 1)};

  cell_dim_ = {static_cast<double>(domain_[0]) / cell_count_[0],
               static_cast<double>(domain_[1]) / cell_count_[1],
               static_cast<double>(domain_[2]) / cell_count_[2]};

  // safety check that minimum cell count is satisfied so the boundaries work as
  // expected:
  // could be ok, but then periodic calculation has to be disabled, since all
  // possible pairs are already iterated over
  for (std::size_t i = 0; i < 2; i++) {
    if (cell_count_[i] < 3) {
      SpdWrapper::get()->error(
          "Cell count is too small if periodic boundaries are used! If this "
          "is not a testing instance, please exit the simulation");
    }
  }

  // add 2 for halo
  cell_count_ = {cell_count_[0] + 2, cell_count_[1] + 2, cell_count_[2] + 2};

  cells_.resize(cell_count_[0] * cell_count_[1] * cell_count_[2]);

  this->boundary_config_ = linked_cells_config.boundary_config;

  DEBUG_PRINT_FMT("Num Cells: {}", cells_.size());

  halo_direction_cells_ = {};
  boundary_direction_cells_ = {};
  // precalculate special cells
  for (std::size_t cell_index = 0; cell_index < cells_.size(); ++cell_index) {
    auto halo_directions = specialCellDirection(
        cell_index, [this](const std::size_t index) { return isHalo(index); },
        -1, 2);
    auto boundary_directions = specialCellDirection(
        cell_index,
        [this](const std::size_t index) { return isBoundary(index); }, 0, 3);

    if (!halo_directions.empty()) {
      for (const unsigned long halo_direction : halo_directions) {
        halo_direction_cells_[halo_direction].push_back(cell_index);
      }
    }
    if (!boundary_directions.empty()) {
      for (const unsigned long boundary_direction : boundary_directions) {
        boundary_direction_cells_[boundary_direction].push_back(cell_index);
      }
    }
  }

  this->boundaries_ = {
      linked_cells_config.boundary_config.x_low,
      linked_cells_config.boundary_config.x_high,
      linked_cells_config.boundary_config.y_low,
      linked_cells_config.boundary_config.y_high,
      linked_cells_config.boundary_config.z_low,
      linked_cells_config.boundary_config.z_high,
  };
  is_membrane = linked_cells_config.is_membrane;
  INFO_FMT("Using membranes is {}", is_membrane);
  // TODO
  // const auto &[ids, time, force_values, dims] =
  // linked_cells_config.index_force_config; index_force = IndexForce(ids, time,
  // force_values);

  SpdWrapper::get()->info("cell dim: {}, {}, {}; cell count: {}, {}, {}",
                          cell_dim_[0], cell_dim_[1], cell_dim_[2],
                          cell_count_[0], cell_count_[1], cell_count_[2]);

  // std::random_device rd;
  // std::mt19937 g(rd());

  // std::shuffle(cells_.begin(), cells_.end(), g);

  initializeC18Schema();
}

void LinkedCellsContainer::addParticle(Particle &p) {
  const std::size_t index = dvec3ToCellIndex(p.getX());
  if (!isValidCellCoordinate(cellIndexToCoord(index))) {
    SpdWrapper::get()->error("Tried to add particle out of bounds");
    exit(1);
  }
  particles_.push_back(p);
  cells_[index].push_back(&particles_.back());

  this->particle_count_++;

  if (p.getType() < 0) {
    this->special_particle_count_++;
  }

  // DEBUG_PRINT_FMT("Added particle with coords ({}, {}, {}) into cell index:
  // {}",
  //                 p.getX()[0], p.getX()[1], p.getX()[2], index)
}

void LinkedCellsContainer::addParticles(
    const std::vector<Particle> &particles) {
  particles_.reserve(particles.size());
  for (Particle p : particles) {
    // SpdWrapper::get()->info("Adding Particle with Id : {}", p.getId());
    addParticle(p);
  }

  if (is_membrane) {
    SpdWrapper::get()->info("Recalculated the neighbour references");
    setNeighbourReferences();
  }

  SpdWrapper::get()->info("Added new particles");
}

void LinkedCellsContainer::removeParticle(const Particle &p) {
  SpdWrapper::get()->info("Particle Id remove: {}", p.getId());
  const std::size_t index = dvec3ToCellIndex(p.getX());
  std::vector<Particle *> &particles = cells_[index];

  particles.erase(std::remove_if(particles.begin(), particles.end(),
                                 [&p](const Particle *q) {
                                   return *q == p;
                                 }),  // Compare dereferenced Particle*
                  particles.end());

  this->particle_count_--;
  if (p.getType() < 0) {
    this->special_particle_count_--;
  }

  DEBUG_PRINT_FMT(
      "Removed particle with coords ({}, {}, {}) from cell index: {}",
      p.getX()[0], p.getX()[1], p.getX()[2], index)
  // TODO : throw std::runtime_error("Not implemented remove Partice.");
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
  for (auto &c : cells_) {
    count += c.size();
  }
  return count;
}

void LinkedCellsContainer::imposeInvariant() {
  // register in corresponding cell
  for (std::size_t index = 0; index < cells_.size(); index++) {
    for (auto it = cells_[index].begin(); it < cells_[index].end();) {
      if (*it == nullptr) {
        SpdWrapper::get()->error("Nullptr found");
        continue;
      }
      const std::size_t shouldBeIndex = dvec3ToCellIndex((*it)->getX());
      if (shouldBeIndex == index) {
        ++it;
        continue;
      }
      // SpdWrapper::get()->info("Imposing on partice with id {}", it->getId());
      ivec3 should_be_cell = cellIndexToCoord(shouldBeIndex);
      // SpdWrapper::get()->info(
      //     "Trying to register p {} at pos [{}, {}, {}] in cell [{}, {}, {}]",
      //     (*it)->getId(), (*it)->getX()[0], (*it)->getX()[1],
      //     (*it)->getX()[2], should_be_cell[0], should_be_cell[1],
      //     should_be_cell[2]);

      if (!isValidCellCoordinate(should_be_cell)) {
        std::cout << (*it)->getX()[0] << ", " << (*it)->getX()[1] << ", " << (*it)->getX()[2] << std::endl;
      }
      cells_[shouldBeIndex].push_back(*it);
      it = cells_[index].erase(it);
    }
  }

  // SpdWrapper::get()->info("Registered Particles in new cell");

  // apply boundary condition
  // it is assumed that GhostParticles do not have to persist, so we dont have
  // to iterate over the halo cells of Reflective Boundaries

  for (size_t dimension = 0; dimension < 6; ++dimension) {
    switch (boundaries_[dimension]) {
      case LinkedCellsConfig::BoundaryType::Outflow: {
        // clear halo
        for (const size_t cell_index : halo_direction_cells_[dimension]) {
          particle_count_ -=
              cells_[cell_index].size();  // update particle count, only place
                                          // where particles are deleted
          cells_[cell_index].clear();
          cells_[cell_index].shrink_to_fit();
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
                        static_cast<int>(boundaries_[dimension]), dimension);
        break;
      }
    }
  }
}

void LinkedCellsContainer::singleIterator(
    const std::function<void(Particle &)> &f) {
  for (const auto cell : cells_) {
    for (const auto p : cell) {
      f(*p);
    }
  }
}

void LinkedCellsContainer::setIndexForce(const IndexForce &index_force) {
  this->index_force = index_force;
}

// TODO: this is now unused => could be removed, are we allowed to? I guess
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
  for (std::size_t cell_index = 0; cell_index < cells_.size(); cell_index++) {
    std::vector<Particle *> cell_particles = cells_[cell_index];

    if (cell_particles.empty()) continue;

    ivec3 cell_coordinate = cellIndexToCoord(cell_index);
    DEBUG_PRINT_FMT("cell index: {}; coord = ({}, {}, {}); halo? = {}",
                    cell_index, cell_coordinate[0], cell_coordinate[1],
                    cell_coordinate[2], isHalo(cell_index))

    // iterate over particles inside cell
    for (std::size_t i = 0; i < cell_particles.size(); ++i) {
      for (std::size_t j = i + 1; j < cell_particles.size(); ++j) {
        const dvec3 p = cell_particles[i]->getX();
        const dvec3 q = cell_particles[j]->getX();
        if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff_ * cutoff_) {
          continue;
        }
        f(*cell_particles[i], *cell_particles[j]);
        // SpdWrapper::get()->info("Index pair {}/{}", cellParticles[i].getId(),
        //                        cellParticles[j].getId());
        DEBUG_PRINT_FMT("Intra cell pair: ({}, {})",
                        cell_particles[i]->getId(),
                        cell_particles[j]->getId());
      }
    }

    // iterate over neighbouring particles
    for (auto &offset : offsets) {
      // compute neighbourIndex and check if it is valid
      const ivec3 neighbour_coord = {cell_coordinate[0] + offset[0],
                                     cell_coordinate[1] + offset[1],
                                     cell_coordinate[2] + offset[2]};

      if (!isValidCellCoordinate(neighbour_coord)) {
        DEBUG_PRINT_FMT("Invalid coord: ({}, {}, {})", neighbour_coord[0],
                        neighbour_coord[1], neighbour_coord[2])
        continue;
      }

      const size_t neighbour_index = cellCoordToIndex(neighbour_coord);
      DEBUG_PRINT_FMT(
          "Checking cell i={}; c=({}, {}, {}) for pairs (offset = ({}, {}, "
          "{}))",
          neighbour_index, neighbour_coord[0], neighbour_coord[1],
          neighbour_coord[2], offset[0], offset[1], offset[2]);

      // go over all pairs with neighbour particles
      std::vector<Particle *> &neighbour_particles = cells_[neighbour_index];
      if (neighbour_particles.empty()) continue;

      for (const auto cell_particle : cell_particles) {
        for (const auto neighbour_particle : neighbour_particles) {
          auto p = cell_particle->getX();
          auto q = neighbour_particle->getX();

          if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
              d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff_ * cutoff_)
            continue;

          f(*cell_particle, *neighbour_particle);
          DEBUG_PRINT_FMT("Cross cell pair: ({}, {})", cell_particle->getType(),
                          neighbour_particle->getType())
        }
      }
    }
  }
}

void LinkedCellsContainer::boundaryIterator(
    const std::function<void(Particle &)> &f) {
  for (std::size_t index = 0; index < cells_.size(); index++) {
    if (!isBoundary(index)) continue;

    for (auto &p : cells_[index]) {
      f(*p);
    }
  }
}

void LinkedCellsContainer::haloIterator(
    const std::function<void(Particle &)> &f) {
  for (std::size_t index = 0; index < cells_.size(); index++) {
    if (!isHalo(index)) continue;

    for (auto &p : cells_[index]) {
      f(*p);
    }
  }
}

void LinkedCellsContainer::computeInteractiveForces(const
    std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces) {
  constexpr std::array<ivec3, 13> offsets = {{
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

  for (auto &colour : c18_colours_) {
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (const std::size_t cell_index : colour) {
      if (cells_[cell_index].empty()) {
        continue;
      }
      ivec3 cell_coordinate = cellIndexToCoord(cell_index);
      auto &cell = cells_[cell_index];

      // iterate over particles inside cell
      for (std::size_t p1 = 0; p1 < cell.size(); p1++) {
        for (std::size_t p2 = p1 + 1; p2 < cell.size(); p2++) {
          const auto p = cell[p1];
          const auto q = cell[p2];
          const dvec3 p_X = p->getX();
          const dvec3 q_X = q->getX();
          if (dvec3 d = {p_X[0] - q_X[0], p_X[1] - q_X[1], p_X[2] - q_X[2]};
              d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff_ * cutoff_)
            continue;
          dvec3 neighbour_force = {0, 0, 0};
          for (auto &force : interactive_forces) {
            neighbour_force = neighbour_force + force->directionalForce(*p, *q);
          }
          p->addF(neighbour_force);
          q->subF(neighbour_force);
          DEBUG_PRINT_FMT("Intra cell pair: ({}, {})", p->getType(),
                          q->getType());
        }
      }

      // iterate over neighbouring newton3 cells
      for (auto &offset : offsets) {
        // compute neighbourCoord and check if it is valid
        const ivec3 neighbourCoord = {cell_coordinate[0] + offset[0],
                                      cell_coordinate[1] + offset[1],
                                      cell_coordinate[2] + offset[2]};

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
        auto &neighbour_particles = cells_[neighbourIndex];
        if (neighbour_particles.empty()) continue;

        for (auto &cell_particle : cell) {
          for (auto &neighbour_particle : neighbour_particles) {
            auto p = cell_particle->getX();
            auto q = neighbour_particle->getX();

            if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff_ * cutoff_)
              continue;
            dvec3 neighbour_force = {0, 0, 0};
            for (auto &force : interactive_forces) {
              neighbour_force =
                  neighbour_force +
                  force->directionalForce(*cell_particle, *neighbour_particle);
            }
            cell_particle->addF(neighbour_force);
            neighbour_particle->subF(neighbour_force);
            // f(cell_particle, neighbour_particle);
            DEBUG_PRINT_FMT("Cross cell pair: ({}, {})", cell_particle->getId(),
                            neighbour_particle->getId())
          }
        }
      }
    }
    #pragma omp barrier
  }
}

void LinkedCellsContainer::computeSingularForces(
    const std::vector<std::unique_ptr<SingularForce>> &singular_forces) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (auto &p : particles_) {
    dvec3 f = {0, 0, 0};
    for (auto &force : singular_forces) {
      f = f + force->applyForce(p);
    }
    p.addF(f);
  }
}

inline std::size_t LinkedCellsContainer::dvec3ToCellIndex(
    const dvec3 &position) const {
  const std::array<int, 3> cell_coords = {
      static_cast<int>(std::floor(position[0] / cell_dim_[0])),
      static_cast<int>(std::floor(position[1] / cell_dim_[1])),
      static_cast<int>(std::floor(position[2] / cell_dim_[2]))};

  return cellCoordToIndex(cell_coords);
}

inline std::size_t LinkedCellsContainer::cellCoordToIndex(
    const ivec3 position) const {
  return (position[0] + 1) * (cell_count_[1] * cell_count_[2]) +
         (position[1] + 1) * (cell_count_[2]) + (position[2] + 1);
}

inline ivec3 LinkedCellsContainer::cellIndexToCoord(
    std::size_t cell_index) const {
  const int x =
      static_cast<int>(cell_index / (cell_count_[1] * cell_count_[2]));
  cell_index = cell_index - (x * cell_count_[1] * cell_count_[2]);

  const int y = static_cast<int>(cell_index / cell_count_[2]);
  const int z = static_cast<int>(cell_index - (y * cell_count_[2]));

  return {x - 1, y - 1, z - 1};
}

inline bool LinkedCellsContainer::isValidCellCoordinate(
    const ivec3 coordinate) const {
  return (-1 <= coordinate[0] && coordinate[0] <= (cell_count_[0] - 2)) &&
         (-1 <= coordinate[1] && coordinate[1] <= (cell_count_[1] - 2)) &&
         (-1 <= coordinate[2] && coordinate[2] <= (cell_count_[2] - 2));
}

inline bool LinkedCellsContainer::isHalo(const ivec3 cell_coord) const {
  return cell_coord[0] == -1 || cell_coord[1] == -1 || cell_coord[2] == -1 ||
         cell_coord[0] == (cell_count_[0] - 2) ||
         cell_coord[1] == (cell_count_[1] - 2) ||
         cell_coord[2] == (cell_count_[2] - 2);
}

inline bool LinkedCellsContainer::isHalo(const std::size_t cell_index) const {
  const ivec3 cell_coord = cellIndexToCoord(cell_index);
  return isHalo(cell_coord);
}

inline bool LinkedCellsContainer::isBoundary(const ivec3 cell_coord) const {
  return (cell_coord[0] == 0 || cell_coord[1] == 0 || cell_coord[2] == 0 ||
          cell_coord[0] == (cell_count_[0] - 3) ||
          cell_coord[1] == (cell_count_[1] - 3) ||
          cell_coord[2] == (cell_count_[2] - 3)) &&
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
  if (cell_coord[0] == (cell_count_[0] - upper_magic_number)) {
    directions.push_back(xhigh);  // east
  }
  if (cell_coord[1] == lower_magic_number) {
    directions.push_back(ylow);  // down
  }
  if (cell_coord[1] == (cell_count_[1] - upper_magic_number)) {
    directions.push_back(yhigh);  // up
  }
  if (cell_coord[2] == lower_magic_number) {
    directions.push_back(zlow);  // south
  }
  if (cell_coord[2] == (cell_count_[2] - upper_magic_number)) {
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

  for (const std::size_t cell_index : boundary_direction_cells_[dimension]) {
    for (const auto p : cells_[cell_index]) {
      // check if it is too close
      double pos = p->getX()[problematic_dimension];
      const double boundary_position = domain_[problematic_dimension];
      const double double_dist_to_boundary =
          2 * std::min(pos, boundary_position -
                                pos);  // if both of them are so small that
      // they would trigger the boundary, the
      // simulation itself is already broken
      // SpdWrapper::get()->info(
      //     "double dist: {} at [{}, {}, {}] adn v=[{}, {}, {}]",
      //     double_dist_to_boundary, p->getX()[0], p->getX()[1], p->getX()[2],
      //     p->getV()[0], p->getV()[1], p->getV()[2]);
      if (double_dist_to_boundary < sigma_factor * p->getSigma()) {
        const double force =
            LennardJones::simpleForce(*p, double_dist_to_boundary);
        dvec3 p_force = p->getF();
        p_force[problematic_dimension] +=
            force *
            std::pow(-1.0,
                     1 - problematic_dimension_direction);  // invert the
        // direction for
        // boundary in
        // ascending
        // coordinate
        // direction
        p->setF(p_force);
        DEBUG_PRINT_FMT(
            "Applied Force=[{}, {}, {}] to Particle at [{}, {}, {}]",
            p->getF()[0], p->getF()[1], p->getF()[2], p->getX()[0],
            p->getX()[1], p->getX()[2]);

        // SpdWrapper::get()->info(
        //     "Applied Force=[{}, {}, {}] to Particle at [{}, {}, {}]",
        //     p->getF()[0], p->getF()[1], p->getF()[2], p->getX()[0],
        //     p->getX()[1], p->getX()[2]);
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

  for (const std::size_t cell_index : halo_direction_cells_[dimension]) {
    int counter = 0;
    for (auto it = cells_[cell_index].begin(); it < cells_[cell_index].end();
         ++it) {
      counter++;
      dvec3 new_pos = (*it)->getX();
      new_pos[problematic_dimension] +=
          domain_[problematic_dimension] *
          (problematic_dimension_direction % 2 == 0 ? 1 : -1);
      const std::size_t shouldBeIndex = dvec3ToCellIndex(new_pos);

      (*it)->setX(new_pos);
      cells_[shouldBeIndex].push_back(*it);
    }

    cells_[cell_index].clear();
    cells_[cell_index].shrink_to_fit();
  }

  // skip force calculation for lower side of the axis
  if (problematic_dimension_direction == 0) {
    return;
  }

  // iterate over all 9 / 3 cells on the other end
  for (const std::size_t cell_index : boundary_direction_cells_[dimension]) {
    ivec3 cell_coordinates = cellIndexToCoord(cell_index);

    // change 3 to 9 for 3D
    for (std::size_t i = 0; i < 9; ++i) {
      ivec3 offset = index_offsets_[problematic_dimension][i];
      const ivec3 cell_to_check = cell_coordinates + offset;
      bool is_adjacent_cell;
      ivec3 adjacent_cell_coordinates;
      dvec3 particle_distance_offset;

      std::tie(is_adjacent_cell, adjacent_cell_coordinates,
               particle_distance_offset) =
          reflectiveWarpAround(cell_to_check, dimension);

      if (!is_adjacent_cell) {
        continue;
      }

      const auto adjacent_cell_index =
          cellCoordToIndex(adjacent_cell_coordinates);

      // // account for the dimension that is checked
      // particle_distance_offset[problematic_dimension] =
      //     domain[problematic_dimension];

      // iterate over all pairs and calculate force

      for (auto &p : cells_[cell_index]) {
        for (auto &q : cells_[adjacent_cell_index]) {
          // distance check
          const dvec3 accounted_particle_distance =
              q->getX() - p->getX() + particle_distance_offset;

          if (ArrayUtils::L2InnerProduct(accounted_particle_distance) >=
              cutoff_ * cutoff_) {
            continue;
          }

          const dvec3 applied_force = LennardJones::directionalForceWithOffset(
              *p, *q, accounted_particle_distance);
          p->setF(p->getF() + applied_force);
          q->setF(q->getF() - applied_force);
        }
      }
    }
  }
}

std::tuple<bool, ivec3, dvec3> LinkedCellsContainer::reflectiveWarpAround(
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
      new_cell_coordinate[dimension] = cell_count_[dimension] - 3;  // top cell
      offset[dimension] = -domain_[dimension];
    } else if (cell_coordinate[dimension] == cell_count_[dimension] - 2) {
      // high warp around to low cell
      new_cell_coordinate[dimension] = 0;  // bottom cell
      offset[dimension] = domain_[dimension];
    } else {
      // no warp around, nothing can go wrong
      continue;
    }

    // if it is wrapped around but dimension is not periodic, this cell is not
    // adjacent
    if (boundaries_[2 * dimension] != LinkedCellsConfig::Periodic) {
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
        cell_coordinate[dimension] == cell_count_[dimension] - 2) {
      edge_of_x_dimensions_counter++;
    }
  }

  if (edge_of_x_dimensions_counter < 2) {
    // SpdWrapper::get()->info("[{}, {}, {}] is not a corner!",
    // cell_coordinate[0],
    //                       cell_coordinate[1], cell_coordinate[2]);
    return false;
  }
  // SpdWrapper::get()->info("[{}, {}, {}] is a corner!", cell_coordinate[0],
  //                         cell_coordinate[1], cell_coordinate[2]);

  // TODO: make beautiful
  // lookup table via some if statements
  // ugly, but best I could think of
  if (boundaries_[xlow] == LinkedCellsConfig::Periodic &&
      boundaries_[ylow] == LinkedCellsConfig::Periodic &&
      boundaries_[zlow] != LinkedCellsConfig::Periodic) {
    // 110
    // SpdWrapper::get()->info("isDoubleCorner - xlow ylow !zlow");
    if (raw_dimension == yhigh && (cell_coordinate[0] == -1 ||
                                   cell_coordinate[0] == cell_count_[0] - 2)) {
      return true;
    }
    return false;
  }
  if (boundaries_[xlow] == LinkedCellsConfig::Periodic &&
      boundaries_[ylow] != LinkedCellsConfig::Periodic &&
      boundaries_[zlow] == LinkedCellsConfig::Periodic) {
    // 101
    // SpdWrapper::get()->info("isDoubleCorner - xlow !ylow zlow");
    if (raw_dimension == zhigh && (cell_coordinate[0] == -1 ||
                                   cell_coordinate[0] == cell_count_[0] - 2)) {
      return true;
    }
    return false;
  }
  if (boundaries_[xlow] != LinkedCellsConfig::Periodic &&
      boundaries_[ylow] == LinkedCellsConfig::Periodic &&
      boundaries_[zlow] == LinkedCellsConfig::Periodic) {
    // 011
    // SpdWrapper::get()->info("isDoubleCorner - !xlow ylow zlow");
    if (raw_dimension == zhigh && (cell_coordinate[1] == -1 ||
                                   cell_coordinate[1] == cell_count_[1] - 2)) {
      return true;
    }
    return false;
  }
  if (boundaries_[xlow] == LinkedCellsConfig::Periodic &&
      boundaries_[ylow] == LinkedCellsConfig::Periodic &&
      boundaries_[zlow] == LinkedCellsConfig::Periodic) {
    // 111
    // SpdWrapper::get()->info("isDoubleCorner - xlow ylow zlow");
    if (raw_dimension == yhigh && (cell_coordinate[0] == -1 ||
                                   cell_coordinate[0] == cell_count_[0] - 2)) {
      return true;
    }
    if (raw_dimension == zhigh && (cell_coordinate[0] == -1 ||
                                   cell_coordinate[0] == cell_count_[0] - 2)) {
      return true;
    }
    if (raw_dimension == zhigh && (cell_coordinate[1] == -1 ||
                                   cell_coordinate[1] == cell_count_[1] - 2)) {
      return true;
    }
    return false;
  }
  return false;
}

void LinkedCellsContainer::setNeighbourReferences() {
  // auto p3 = &particles_[0];
  // std::cout << "in setNeighbours Particle 1 has reference location " << p3
  //           << std::endl;

  for (Particle *p : getParticles()) {
    std::vector<std::pair<bool, size_t>> new_neighbours{};

    for (Particle *p2 : getParticles()) {
      for (auto [diag, ref] : p->getNeighbours()) {
        auto *new_p = reinterpret_cast<Particle *>(ref);

        if (p2->getId() == new_p->getId()) {
          size_t *pointer = (size_t *)p2;
          new_neighbours.emplace_back(diag, (size_t)pointer);
        }

        // if (new_p->getId() == 0) {
        //   size_t *pointer = (size_t *)new_p;
        //   std::cout << "P0 is at " << new_p << " or " << pointer << " or "
        //             << ((size_t)(pointer)) << std::endl;
        //   std::cout << "P0 is at actually at " << p3 << std::endl;
        // }
      }
    }

    p->resetNeighbours();
    p->setNeighbours(new_neighbours);
  }
}
