//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include "../utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  // update positions
  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_x = p.getX() + delta_t * p.getV() +
                        (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);
  });

  // advance time, old F_{t] -> F_{t-1}, F_{t} = 0
  particle_container.singleIterator([](Particle& p) {
    p.setOldF(p.getF());
    p.setF({0, 0, 0});
  });

  // update force

  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    const dvec3 g12 = force.directionalForce(p1, p2);
    p1.setF(p1.getF() + g12);  // F_i = \sum_j F_ij
    p2.setF(p2.getF() - g12);  // g12 = -g21
  });

  /*
  // alternative Iterator
  std::vector<Particle> particles = particle_container.getParticles();
  for (std::size_t i = 0; i < particle_container.size(); i++) {;
    dvec3 acc = particles[i].getF();
    for (std::size_t j = i+1; j < particle_container.size(); j++) {
      const dvec3 g12 = force.directionalForce(particles[i],particles[j]);
      acc = acc + g12;
      const dvec3 p2F = particles[j].getF() - g12;
      particles[j].setF(p2F);
    }
    particles[i].setF(acc);
  }
  particle_container.setParticles(particles);
  */

  // Now we use F_t and F_{t-1} to calculate the current velocity
  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_v =
        p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });
}
