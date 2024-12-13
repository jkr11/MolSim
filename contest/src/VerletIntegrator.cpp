//
// Created by mcarn on 12/11/24.
//

#include <iostream>
#include <immintrin.h>

#include "VerletIntegrator.h"
#include "../../src/defs/types.h"
#include "../../src/utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& container) {
  /*for (int i = 0; i < container.ids.size(); i++) {
    const dvec3 pos = {container.px[i], container.py[i], container.pz[i]};
    const dvec3 vel = {container.vx[i], container.vy[i], container.vz[i]};
    const dvec3 f = {container.fx[i], container.fy[i], container.fz[i]};
    const int type = container.types[i];
    const double mass = container.particleTypeInfo.mass[type];

    const dvec3 newX = pos + delta_t * vel + (delta_t * delta_t / (2 * mass)) * f;

    container.px[i] = newX[0];
    container.py[i] = newX[1];
    container.pz[i] = newX[2];

    container.ofx[i] = container.fx[i];
    container.ofy[i] = container.fy[i];
    container.ofz[i] = container.fz[i];

    container.fx[i] = 0;
    container.fy[i] = 0;
    container.fz[i] = 0;
  }

  //TODO: merge loops?
  for (int i = 0; i < container.ids.size(); i++) {
    container.ofx[i] = container.fx[i];
    container.ofy[i] = container.fy[i];
    container.ofz[i] = container.fz[i];

    container.fx[i] = 0;
    container.fy[i] = 0;
    container.fz[i] = 0;
  }

*/

  const std::size_t size = container.ids.size();
  const __m256d delta_t_vec = _mm256_set1_pd(delta_t);
  const __m256d half_delta_t_squared = _mm256_set1_pd((delta_t * delta_t) / 2.0);

  for (int i = 0; i < container.ids.size(); i += 4) {
    // Load position vectors
    __m256d px = _mm256_load_pd(&container.px[i]);
    __m256d py = _mm256_load_pd(&container.py[i]);
    __m256d pz = _mm256_load_pd(&container.pz[i]);

    // Load velocity vectors
    __m256d vx = _mm256_load_pd(&container.vx[i]);
    __m256d vy = _mm256_load_pd(&container.vy[i]);
    __m256d vz = _mm256_load_pd(&container.vz[i]);

    // Load force vectors
    __m256d fx = _mm256_load_pd(&container.fx[i]);
    __m256d fy = _mm256_load_pd(&container.fy[i]);
    __m256d fz = _mm256_load_pd(&container.fz[i]);

    // Compute mass for each particle (scalar lookup -> broadcast to SIMD)
    alignas(64) double masses[4];
    for (int j = 0; j < 4; ++j) {
      masses[j] = container.particleTypeInfo.mass[container.types[i + j]];
    }
    __m256d mass = _mm256_load_pd(masses);

    // Compute new positions: pos + delta_t * vel + (delta_t^2 / 2m) * f
    __m256d new_px = _mm256_fmadd_pd(delta_t_vec, vx, _mm256_fmadd_pd(_mm256_div_pd(half_delta_t_squared, mass), fx, px));
    __m256d new_py = _mm256_fmadd_pd(delta_t_vec, vy, _mm256_fmadd_pd(_mm256_div_pd(half_delta_t_squared, mass), fy, py));
    __m256d new_pz = _mm256_fmadd_pd(delta_t_vec, vz, _mm256_fmadd_pd(_mm256_div_pd(half_delta_t_squared, mass), fz, pz));

    // Store the results back to memory
    _mm256_store_pd(&container.px[i], new_px);
    _mm256_store_pd(&container.py[i], new_py);
    _mm256_store_pd(&container.pz[i], new_pz);

    _mm256_store_pd(&container.ofx[i], fx);
    _mm256_store_pd(&container.ofy[i], fy);
    _mm256_store_pd(&container.ofz[i], fz);

    __m256 zero = _mm256_setzero_pd();
    _mm256_store_pd(&container.fx[i], zero);
    _mm256_store_pd(&container.fy[i], zero);
    _mm256_store_pd(&container.fz[i], zero);
  }

  container.imposeInvariant();

  //TODO: single force application
  /*for (int i = 0; i < container.ids.size(); i++) {
    dvec3 f = {0, 0, 0};
    for (const auto& force : singular_forces) {
      f = f + force->applyForce(p);
    }
    p.setF(f);
  }*/

  //pairwise force application
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

  //TODO debug
  /*std::vector<std::pair<int, int>> pairs = {};
  for (int i = 0; i < container.ids.size(); i++) {
    for (int j = i + 1; j < container.ids.size(); j++) {
      dvec3 d = { container.px[i] - container.px[j],
                    container.py[i] - container.py[j],
                    container.pz[i] - container.pz[j] };

      if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > container.linkedCellsConfig.cutoff_radius * container.linkedCellsConfig.cutoff_radius)
        continue;

      pairs.emplace_back(container.ids[i], container.ids[j]);
      //std::cout << "should be pair: (" << container.ids[i] << ", " << container.ids[j] << ")" << std::endl;
    }
  }

  int pairCount = 0;
*/
  //debug end

  // go over all cell indices
  for (std::size_t cellIndex = 0; cellIndex < container.partitionStart.size(); cellIndex++) {

    //TODO: check if last partition is empty
    if (container.partitionSizes[cellIndex] == 0)
      continue;

    ivec3 cellCoordinate = container.cellIndexToCoordinate(cellIndex);
    //std::cout << "cell index: " << cellIndex << "; coord = (" << cellCoordinate[0] << ", " << cellCoordinate[1] << ", " << cellCoordinate[2] << ")" << std::endl;

    // iterate over particles inside cell
    std::size_t partitionStart = container.partitionStart[cellIndex];
    std::size_t partitionEnd = partitionStart + container.partitionSizes[cellIndex];

    //TODO: how to even SIMD this?
    for (std::size_t i = partitionStart; i < partitionEnd; ++i) {
      for (std::size_t j = i + 1; j < partitionEnd; ++j) {
        dvec3 d = { container.px[i] - container.px[j],
                    container.py[i] - container.py[j],
                    container.pz[i] - container.pz[j] };

        if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > container.linkedCellsConfig.cutoff_radius * container.linkedCellsConfig.cutoff_radius)
          continue;

        //TODO debug
        /*bool found = false;
        for (int k = 0; k < pairs.size(); k++) {
            if ((pairs[k].first == container.ids[i] && pairs[k].second == container.ids[j]) || (pairs[k].first == container.ids[j] && pairs[k].second == container.ids[i])) {
              found = true;
              break;
            }
        }

        if (!found) {
           std::cout << "Something just went wrong massively." << std::endl;
        }

        std::cout << "Intra cell pair: (" << container.ids[i] << ", " << container.ids[j] << ")" << std::endl;
        pairCount++;
        */

        //function start
        dvec3 f12 = {0.0, 0.0, 0.0};
        const double sigma1 = container.particleTypeInfo.sigma[container.types[i]];
        const double sigma2 = container.particleTypeInfo.sigma[container.types[j]];
        const double eps1 = container.particleTypeInfo.epsilon[container.types[i]];
        const double eps2 = container.particleTypeInfo.epsilon[container.types[j]];
        for (const auto& force : interactive_forces) {
          f12 = f12 + force->directionalForce({container.px[i], container.py[i], container.pz[i]},
                                              {container.px[j], container.py[j], container.pz[j]},
                                              sigma1, sigma2, eps1, eps2);
        }

        container.fx[i] += f12[0];
        container.fy[i] += f12[1];
        container.fz[i] += f12[2];

        container.fx[j] -= f12[0];
        container.fy[j] -= f12[1];
        container.fz[j] -= f12[2];
        //function end
      }
    }

    // iterate over neighbouring particles
    for (auto &offset : offsets) {
      // compute neighbourIndex and check if it is valid
      const ivec3 neighbourCoord = { cellCoordinate[0] + offset[0],
                                     cellCoordinate[1] + offset[1],
                                     cellCoordinate[2] + offset[2] };

      if (!container.isValidCellCoordinate(neighbourCoord))
        continue;

      const size_t neighbourIndex = container.cellCoordinateToIndex(neighbourCoord);
      //std::cout << "Checking cell " << neighbourIndex << "(" << neighbourCoord[0] << ", " << neighbourCoord[1] << ", " << neighbourCoord[2] << ")" << std::endl;

      // go over all pairs with neighbour particles
      if (container.partitionSizes[neighbourIndex] == 0)
        continue;

      std::size_t neighbourPartitionStart = container.partitionStart[neighbourIndex];
      std::size_t neighbourPartitionEnd = neighbourPartitionStart + container.partitionSizes[neighbourIndex];

      for (std::size_t i = partitionStart; i < partitionEnd; i++) {
        for (std::size_t j = neighbourPartitionStart; j < neighbourPartitionEnd; j++) {
          dvec3 d = { container.px[i] - container.px[j],
                    container.py[i] - container.py[j],
                    container.pz[i] - container.pz[j] };

          if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > container.linkedCellsConfig.cutoff_radius * container.linkedCellsConfig.cutoff_radius)
            continue;

          //TODO debug
          /*bool found = false;
          for (int k = 0; k < pairs.size(); k++) {
              if ((pairs[k].first == container.ids[i] && pairs[k].second == container.ids[j]) || (pairs[k].first == container.ids[j] && pairs[k].second == container.ids[i])) {
                found = true;
                break;
              }
          }

          if (!found) {
             std::cout << "Something just went wrong massively." << std::endl;
          }

          std::cout << "Cross cell pair: (" << container.ids[i] << ", " << container.ids[j] << ")" << std::endl;
          pairCount++;*/

          //function start
          dvec3 f12 = {0.0, 0.0, 0.0};
          const double sigma1 = container.particleTypeInfo.sigma[container.types[i]];
          const double sigma2 = container.particleTypeInfo.sigma[container.types[j]];
          const double eps1 = container.particleTypeInfo.epsilon[container.types[i]];
          const double eps2 = container.particleTypeInfo.epsilon[container.types[j]];
          for (const auto& force : interactive_forces) {
            f12 = f12 + force->directionalForce({container.px[i], container.py[i], container.pz[i]},
                                                {container.px[j], container.py[j], container.pz[j]},
                                                sigma1, sigma2, eps1, eps2);
          }

          container.fx[i] += f12[0];
          container.fy[i] += f12[1];
          container.fz[i] += f12[2];

          container.fx[j] -= f12[0];
          container.fy[j] -= f12[1];
          container.fz[j] -= f12[2];
          //function end
        }
      }
    }
  }

  //TODO debug
  //std::cout << pairCount << ", " << pairs.size() << std::endl;

  //update velocity
  for (int i = 0; i < container.ids.size(); i++) {
    const dvec3 pos = {container.px[i], container.py[i], container.pz[i]};
    const dvec3 vel = {container.vx[i], container.vy[i], container.vz[i]};
    const dvec3 f = {container.fx[i], container.fy[i], container.fz[i]};
    const dvec3 of = {container.ofx[i], container.ofy[i], container.ofz[i]};
    const int type = container.types[i];
    const double mass = container.particleTypeInfo.mass[type];

    const dvec3 newV = vel + (delta_t / (2 * mass) * (of + f));

    container.vx[i] = newV[0];
    container.vy[i] = newV[1];
    container.vz[i] = newV[2];
  }
}