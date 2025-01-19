/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

Particle::Particle(int type_arg) {
  type = type_arg;
  // DEBUG_PRINT("Particle generated!");
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  epsilon = other.epsilon;
  sigma = other.sigma;
  // DEBUG_PRINT("Particle generated by copy!");
}

Particle::Particle(const std::array<double, 3> &x_arg,
                   const std::array<double, 3> &v_arg, const double m_arg,
                   const double _epsilon, const double _sigma, int _type) {
  x = x_arg;
  v = v_arg;
  m = m_arg;
  type = _type;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
  sigma = _sigma;
  epsilon = _epsilon;
  // DEBUG_PRINT("Particle generated!");
}

Particle::Particle(const std::array<double, 3> &x_arg,
                   const std::array<double, 3> &v_arg,
                   const std::array<double, 3> &f_arg,
                   const std::array<double, 3> &old_f_arg, const double m_arg,
                   const int type_arg, const double epsilon_arg,
                   const double sigma_arg)
    : x(x_arg),
      v(v_arg),
      f(f_arg),
      old_f(old_f_arg),
      m(m_arg),
      type(type_arg),
      epsilon(epsilon_arg),
      sigma(sigma_arg) {}

Particle::~Particle() { DEBUG_PRINT("Particle destructed!"); }

void Particle::updateForceInTime() {
  old_f = f;
  f = {0., 0., 0.};
}

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
         << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

bool Particle::operator==(const Particle &other) const {
  return (x == other.x) and (v == other.v) and (f == other.f) and
         (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

std::ostream &operator<<(std::ostream &stream, const Particle &p) {
  stream << p.toString();
  return stream;
}
