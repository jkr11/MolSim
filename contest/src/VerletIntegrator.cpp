//
// Created by mcarn on 12/11/24.
//

#include "VerletIntegrator.h"
#include "../../src/defs/types.h"
#include "../../src/utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& container) {
  for (int i = 0; i < container.ids.size(); i++) {
    const dvec3 pos = {container.px[i], container.py[i], container.pz[i]};
    const dvec3 vel = {container.vx[i], container.vy[i], container.vz[i]};
    const dvec3 f = {container.fx[i], container.fy[i], container.fz[i]};
    const int type = container.types[i];
    const double mass = container.particleTypeInfo.mass[type];

    const dvec3 newX = pos + delta_t * vel + (delta_t * delta_t / (2 * mass)) * f;

    container.px[i] = newX[0];
    container.py[i] = newX[1];
    container.pz[i] = newX[2];

    //TODO update force here?

  }

  for (int i = 0; i < container.ids.size(); i++) {
    container.ofx[i] = container.fx[i];
    container.ofy[i] = container.fy[i];
    container.ofz[i] = container.fz[i];

    container.fx[i] = 0;
    container.fy[i] = 0;
    container.fz[i] = 0;
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

  /*
  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    dvec3 f12 = {0.0, 0.0, 0.0};
    for (const auto& force : interactive_forces) {
      f12 = f12 + force->directionalForce(p1, p2);
    }
    // p1.setF(p1.getF() + f12);  // F_i = \sum_j F_ij
    // p2.setF(p2.getF() - f12);  // g12 = -g21
    p1.addF(f12);
    p2.subF(f12);
  });

   */

  //TODO pair iter

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