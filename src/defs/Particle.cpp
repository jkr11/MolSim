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

int Particle::global_id_counter = 0;

Particle::Particle(int type_arg) {
  type = type_arg;
  // DEBUG_PRINT("Particle generated!");
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other) = default;

// Move constructor?
Particle::Particle(Particle &&other) noexcept
    : x(other.x),
      v(other.v),
      f(other.f),
      old_f(other.old_f),
      m(other.m),
      type(other.type),
      epsilon(other.epsilon),
      sigma(other.sigma),
      id(other.id),
      neighbours(
          std::move(other.neighbours)) {  // Transfer ownership of shared_ptr
}

Particle::Particle(const std::array<double, 3> &x_arg,
                   const std::array<double, 3> &v_arg, const double m_arg,
                   const double _epsilon, const double _sigma, const int _type)
    : x(x_arg),
      v(v_arg),
      f({0, 0, 0}),
      old_f({0.0, 0.0, 0.0}),
      m(m_arg),
      type(_type),
      epsilon(_epsilon),
      sigma(_sigma),
      id(global_id_counter++) {}

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
      sigma(sigma_arg),
      id(global_id_counter++) {}

Particle::~Particle() = default;


const std::vector<std::pair<bool, std::shared_ptr<Particle>>> &
Particle::getNeighbours() const {
  return neighbours;
}

void Particle::pushBackNeighbour(bool diag,
                                 const std::shared_ptr<Particle> &particle) {
  neighbours.push_back({diag, particle});
}


void Particle::updateForceInTime() {
  old_f = f;
  f = {0., 0., 0.};
}

int Particle::getId() const { return id; }

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
