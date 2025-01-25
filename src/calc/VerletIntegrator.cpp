//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include "../utils/ArrayUtils.h"
#include "debug/debug_print.h"
#include "utils/SpdWrapper.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
// position update
  particle_container.singleIterator([this](Particle& p) {
    if (p.getType() < 0) {
      SpdWrapper::get()->error("Particle {} has -1 type ", p.getId());

      return;
    }  // ignore position update for walls

    dvec3 oldx = p.getX();
    const dvec3 new_x = p.getX() + delta_t * p.getV() +
                        (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);

    if (p.getId() == 874) {
      dvec3 deltapos = new_x - oldx;
      SpdWrapper::get()->info("874 moved to [{}, {}, {}]", deltapos[0], deltapos[1], deltapos[2]);
      SpdWrapper::get()->info("   with f = [{}, {}, {}]", p.getF()[0], p.getF()[1], p.getF()[2]);
    }
  });

  // reset
  particle_container.singleIterator([](Particle& p) { p.updateForceInTime(); });

  particle_container.imposeInvariant();

  // TODO: refactor in lower iterator? maybe pass time to all? just get the
  // global var?
  // Pull Up
  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& index_force : index_forces) {
      for (const int id : index_force->getIndeces()) {
        if (p.getId() == id) {
          //SpdWrapper::get()->info("Particle {}; -> [{}, {}, {}]", p.getId(), f[0], f[1], f[2]);
          f = f + index_force->applyForce(p, current_time);
          //SpdWrapper::get()->info("Particle {} adding [{}, {}, {}]", p.getId(), f[0], f[1], f[2]);
        }
      }
    }
    p.addF(f);
  });

  SpdWrapper::get()->info("----------------- Iteration done");

  // Lennard Jones (or truncated)
  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    dvec3 f12 = {0.0, 0.0, 0.0};
    for (const auto& force : interactive_forces) {
      f12 = f12 + force->directionalForce(p1, p2);
    }
    p1.addF(f12);
    p2.subF(f12);

  });

  // TODO DELETE
  particle_container.singleIterator([this](Particle& p) {
    if (p.getId() == 874) {
      SpdWrapper::get()->info("   current 874 f = [{}, {}, {}]", p.getF()[0], p.getF()[1], p.getF()[2]);
    }
  });

  // Gravity and or Membrane
  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& force : singular_forces) {
      f = f + force->applyForce(p);
    }

    if (p.getId() == 874) {
      SpdWrapper::get()->info("874 membrane force: [{}, {}, {}]", f[0], f[1], f[2]);
    }
    p.addF(f);
  });

  // Velocity Update
  particle_container.singleIterator([this](Particle& p) {
    if (p.getType() < 0) {
      SpdWrapper::get()->error("Particle {} has -1 type ", p.getId());
      return;  // ignore velocity update for walls, theoretically
    }

    const dvec3 new_v =
        p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });

  // TODO: remove
  particle_container.incrementTime();
}
