/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "../utils/ArrayUtils.h"
#include "debug/debug_print.h"

Particle::Particle(int type_arg) {
  type = type_arg;
  DEBUG_PRINT("Particle generated!");
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(const Particle& other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  DEBUG_PRINT("Particle generated by copy!");
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(const std::array<double, 3>& x_arg, const std::array<double, 3>& v_arg,
                   const double m_arg, int _type) {
  x = x_arg;
  v = v_arg;
  m = m_arg;
  type = _type;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
  DEBUG_PRINT("Particle generated!");
}

Particle::~Particle() {
  DEBUG_PRINT("Particle destructed!");
}

const std::array<double, 3>& Particle::getX() const {
  return x;
}

const std::array<double, 3>& Particle::getV() const {
  return v;
}

const std::array<double, 3>& Particle::getF() const {
  return f;
}

const std::array<double, 3>& Particle::getOldF() const {
  return old_f;
}

double Particle::getM() const {
  return m;
}

int Particle::getType() const {
  return type;
}

void Particle::setF(const std::array<double, 3>& F) {
  f = F;
}

void Particle::setV(const std::array<double, 3>& V) {
  v = V;
}

void Particle::setX(const std::array<double, 3>& X) {
  x = X;
}

void Particle::setOldF(const dvec3& oF) {
  old_f = oF;
}


std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
    << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

bool Particle::operator==(const Particle& other) const {
  return (x == other.x) and (v == other.v) and (f == other.f) and
    (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

std::ostream& operator<<(std::ostream& stream, Particle& p) {
  stream << p.toString();
  return stream;
}
