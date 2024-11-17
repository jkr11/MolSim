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
  std::array<int, 3> cellCount;
  dvec3 cellDim;
  double cutoff;

 public:
  explicit LinkedCellsContainer(const dvec3 domain, const double cutoff);

  ~LinkedCellsContainer() override = default;

  void addParticle(const Particle& p) override;

  void removeParticle(const Particle& p) override;

  [[nodiscard]] std::vector<Particle> getParticles() const override;

  [[nodiscard]] std::size_t size() const override;

  void singleIterator(const std::function<void(Particle&)>& f) override;

  void pairIterator(
      const std::function<void(Particle&, Particle&)>& f) override;

  void boundaryIterator(std::function<void(Particle&)> f);

  void haloIterator(std::function<void(Particle&)> f);

  //debug methods
  std::vector<std::vector<Particle>>& getCells() {
    return cells;
  }

  std::array<int, 3> getCellCount() {
    return cellCount;
  }

  dvec3 getCellDim() {
    return cellDim;
  }
      
private:
  inline std::size_t dvec3ToCellIndex(dvec3 position);
  inline std::size_t cellCoordToIndex(std::array<int, 3> position);
  inline std::array<int, 3> cellIndexToCoord(std::size_t cellIndex);
  inline bool isValidCellCoordinate(ivec3 coordinate);
  inline bool isHalo(std::size_t cellIndex);
};