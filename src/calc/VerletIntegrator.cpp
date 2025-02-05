//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include "../utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  // position update
  particle_container.singleIterator([this](Particle& p) {
    if (p.getType() < 0) {
      return;
    }  // ignore position update for walls
    const dvec3 new_x = p.getX() + delta_t_ * p.getV() +
                        (delta_t_ * delta_t_ / (2 * p.getM())) * (p.getF());

    p.setX(new_x);
  });

  // reset
  particle_container.singleIterator([](Particle& p) { p.updateForceInTime(); });

  particle_container.imposeInvariant();

  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& index_force : index_forces_) {
      for (const int id : index_force->getIndices()) {
        if (p.getId() == id) {
          f = f + index_force->applyForce(p, current_time_);
        }
      }
    }
    p.addF(f);
  });

#ifdef _OPENMP
  if (strategy_ == ParallelStrategy::STRATEGY_2) {
    particle_container.computeInteractiveForcesC18(interactive_forces_);
  } else if (strategy_ == ParallelStrategy::STRATEGY_1) {
    particle_container.computeInteractiveForcesForceBuffer(interactive_forces_);
  } else {
    particle_container.pairIterator([this](Particle& p, Particle& q) {
      dvec3 f = {0, 0, 0};
      for (const auto& interactive_force : interactive_forces_) {
        f = f + interactive_force->directionalForce(p, q);
      }
      p.addF(f);
      q.subF(f);
    });
  }
#else
  particle_container.pairIterator([this](Particle& p, Particle& q) {
    dvec3 f = {0, 0, 0};
    for (const auto& interactive_force : interactive_forces_) {
      f = f + interactive_force->directionalForce(p, q);
    }
    p.addF(f);
    q.subF(f);
  });
#endif

  particle_container.computeSingularForces(singular_forces_);

  // Velocity Update
  particle_container.singleIterator([this](Particle& p) {
    if (p.getType() < 0) {
      return;  // ignore velocity update for walls, theoretically
    }
    const dvec3 new_v =
        p.getV() + (delta_t_ / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });

  incrementTime();
}
