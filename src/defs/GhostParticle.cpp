//
// Created by maximilian on 27.11.24.
//

#include "GhostParticle.h"
#include "debug/debug_print.h"

GhostParticle::GhostParticle(const std::array<double, 3> &x_arg,
      double m_arg, double _epsilon, double _sigma, Particle *_ghostParticleOf) {
  x = x_arg;
  m = m_arg;
  epsilon = _epsilon;
  sigma = _sigma;
  ghostParticleOf = _ghostParticleOf;
  DEBUG_PRINT("created ghost particle");
}

 GhostParticle::~GhostParticle() { DEBUG_PRINT("deleted ghost particle"); }

double GhostParticle::getEpsilon() const { return epsilon; }

double GhostParticle::getM() const { return m; }

double GhostParticle::getSigma() const { return sigma; }

const std::array<double, 3> &GhostParticle::getX() const { return x; }

std::string GhostParticle::toString() {
  return "GhostParticle";
}




