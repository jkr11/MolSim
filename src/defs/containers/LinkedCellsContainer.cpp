//
// created by mcarn on 11/15/24
//
#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "../../utils/SpdWrapper.h"
#include "../Particle.h"
#include "LinkedCellsContainer.h"

  LinkedCellsContainer::LinkedCellsContainer(const dvec3 domain, const double cutoff /*boundary types*/) {
    cells = {};
    this->cutoff = cutoff;

    cellCount = {std::max(static_cast<int>(std::floor(domain[0] / cutoff)), 1),
                 std::max(static_cast<int>(std::floor(domain[1] / cutoff)), 1),
                 std::max(static_cast<int>(std::floor(domain[2] / cutoff)), 1)};

    cellDim = {domain[0] / cellCount[0], domain[1] / cellCount[1],
               domain[2] / cellCount[2]};

    // add 2 for halo
    cellCount = {cellCount[0] + 2, cellCount[1] + 2, cellCount[2] + 2};

    cells.resize(cellCount[0] * cellCount[1] * cellCount[2]);

    //TODO: pretty
    SpdWrapper::get()->info("cell dim: {}, {}, {}; cell count: {}, {}, {}", cellDim[0], cellDim[1], cellDim[2], cellCount[0], cellCount[1], cellCount[2]);
  }

  void LinkedCellsContainer::addParticle(const Particle &p) {
    int index = dvec3ToCellIndex(p.getX());
    cells[index].emplace_back(p);

    SpdWrapper::get()->debug("Added particle with coords ({}, {}, {}) into cell index: {}", p.getX()[0], p.getX()[1], p.getX()[2], index);
  }

  void LinkedCellsContainer::removeParticle(const Particle &p) {
    //TODO
  }

  std::vector<Particle> LinkedCellsContainer::getParticles() const {
    //TODO
  }

  [[nodiscard]] std::size_t LinkedCellsContainer::size() const {
    std::size_t count = 0;
    for (auto &c : cells) {
      count += c.size();
    }
    return count;
  }

  void LinkedCellsContainer::singleIterator(const std::function<void(Particle&)>& f) {
    for (auto &c : cells) {
      for (auto &p : c) {
        f(p);
      }
    }
  }

  void LinkedCellsContainer::pairIterator(const std::function<void(Particle &, Particle &)>& f) {
    // - extract offsets out of loop below (should be 13 neighbours in 3d)
    // - as x, y, z are all increasing offsets point to all neighbours in
    // positive directions 
    // - for better cache usage (in flattened version) minimize max distance
    //   between values
    // - direction of traversal: z, y, x
    std::array<ivec3, 13> offsets = {{
      //9 x facing
      {{1, -1, -1}},
      {{1, -1, 0}},
      {{1, -1, 1}},
      {{1, 0, -1}},
      {{1, 0, 0}},
      {{1, 0, 1}},
      {{1, 1, -1}},
      {{1, 1, 0}},
      {{1, 1, 1}},
      //3 y
      {{0, 1, -1}},
      {{0, 1, 0}},
      {{0, 1, 1}},
      //last z
      {{0, 0, 1}},
    }};

    for (std::size_t cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
      std::vector<Particle> &cellParticles = cells[cellIndex];

      ivec3 co = cellIndexToCoord(cellIndex);
      SpdWrapper::get()->debug("cell index: {}; coord = ({}, {}, {}); halo? = {}", cellIndex, co[0], co[1], co[2], isHalo(cellIndex));

      if (cellParticles.empty()) continue;
      
      //iterate over particles inside cell
      for (std::size_t i = 0; i < cellParticles.size(); ++i) {
        for (std::size_t j = i + 1; j < cellParticles.size(); ++j) {
          f(cellParticles[i], cellParticles[j]);
          SpdWrapper::get()->debug("Intra cell pair: ({}, {})", cellParticles[i].getType(), cellParticles[j].getType());
        }
      }

      //iterate over neighbouring particles
      for (auto &offset : offsets) {
        //compute neighbourIndex and check if it is valid
        ivec3 neighbourCoord = {co[0] + offset[0], co[1] + offset[1], co[2] + offset[2]};

        //only fires if halo has particle
        if (!isValidCellCoordinate(neighbourCoord)) {
          SpdWrapper::get()->info("Invalid coord: ({}, {}, {})", neighbourCoord[0], neighbourCoord[1], neighbourCoord[2]);
          continue;
        }

        size_t neighbourIndex = cellCoordToIndex(neighbourCoord);
        SpdWrapper::get()->debug("Checking cell i={}; c=({}, {}, {}) for pairs (offset = ({}, {}, {}))", neighbourIndex, neighbourCoord[0], neighbourCoord[1], neighbourCoord[2], offset[0], offset[1], offset[2]);

        //go over all pairs with neighbour particles
        std::vector<Particle> &neighbourParticles = cells[neighbourIndex];
        if (neighbourParticles.empty()) continue;

        for (std::size_t i = 0; i < cellParticles.size(); ++i) {
          for (std::size_t j = 0; j < neighbourParticles.size(); ++j) {
            auto p = cellParticles[i].getX();
            auto q = neighbourParticles[j].getX();
            dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};

            if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
              continue;

            f(cellParticles[i], neighbourParticles[j]);
            SpdWrapper::get()->debug("Cross cell pair: ({}, {})", cellParticles[i].getType(), neighbourParticles[j].getType());
          }
        }
      }
    }
  }

  void LinkedCellsContainer::boundaryIterator(std::function<void(Particle&)> f) {
    //TODO
  }

  void LinkedCellsContainer::haloIterator(std::function<void(Particle&)> f) {
    //TODO
  }

  inline std::size_t LinkedCellsContainer::dvec3ToCellIndex(dvec3 position) {
    std::array<int, 3> cellCoords = {
      static_cast<int>(std::floor(position[0] / cellDim[0])),
      static_cast<int>(std::floor(position[1] / cellDim[1])),
      static_cast<int>(std::floor(position[2] / cellDim[2]))
    };

    return cellCoordToIndex(cellCoords);
  }

  inline std::size_t LinkedCellsContainer::cellCoordToIndex(std::array<int, 3> position) {
    return (position[0] + 1) * (cellCount[1] * cellCount[2]) + (position[1] + 1) * (cellCount[2]) +
           (position[2] + 1);
  }

  inline std::array<int, 3> LinkedCellsContainer::cellIndexToCoord(std::size_t cellIndex) {
    int x = cellIndex / (cellCount[1] * cellCount[2]);
    cellIndex = cellIndex - (x * cellCount[1] * cellCount[2]);
 
    int y = cellIndex / cellCount[2];
    int z = cellIndex - (y * cellCount[2]);

    return {x - 1, y - 1, z - 1};
  }

  inline bool LinkedCellsContainer::isValidCellCoordinate(ivec3 coordinate) {
    return (-1 <= coordinate[0] && coordinate[0] <= (cellCount[0] - 2))
    && (-1 <= coordinate[1] && coordinate[1] <= (cellCount[1] - 2)) 
    && (-1 <= coordinate[2] && coordinate[2] <= (cellCount[2] - 2));
  }

  inline bool LinkedCellsContainer::isHalo(std::size_t cellIndex) {
    ivec3 cc = cellIndexToCoord(cellIndex);
    return cc[0] == -1 || cc[1] == -1 || cc[2] == -1
      || cc[0] == (cellCount[0] - 2) || cc[1] == (cellCount[1] - 2) || cc[2] == (cellCount[2] - 2); 
  }