/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

int Particle::global_id_counter_ = 0;

Particle::Particle(const int type) {
  type_ = type;
  f_ = {0., 0., 0.};
  old_f_ = {0., 0., 0.};
}

Particle::Particle(const Particle &other) {
  x_ = other.x_;
  v_ = other.v_;
  f_ = other.f_;
  old_f_ = other.old_f_;
  m_ = other.m_;
  type_ = other.type_;
  epsilon_ = other.epsilon_;
  sigma_ = other.sigma_;
  id_ = other.id_;
  neighbours_ = other.neighbours_;
}

Particle::Particle(const std::array<double, 3> &x_arg,
                   const std::array<double, 3> &v_arg, const double m_arg,
                   const double epsilon, const double sigma, const int type) {
  x_ = x_arg;
  v_ = v_arg;
  m_ = m_arg;
  type_ = type;
  f_ = {0., 0., 0.};
  old_f_ = {0., 0., 0.};
  sigma_ = sigma;
  epsilon_ = epsilon;
  id_ = global_id_counter_++;
}

Particle::Particle(const std::array<double, 3> &x_arg,
                   const std::array<double, 3> &v_arg,
                   const std::array<double, 3> &f_arg,
                   const std::array<double, 3> &old_f_arg, const double m_arg,
                   const int type_arg, const double epsilon_arg,
                   const double sigma_arg)
    : x_(x_arg),
      v_(v_arg),
      f_(f_arg),
      old_f_(old_f_arg),
      m_(m_arg),
      type_(type_arg),
      epsilon_(epsilon_arg),
      sigma_(sigma_arg),
      id_(global_id_counter_++) {}

Particle::~Particle() = default;

const std::vector<std::pair<bool, size_t>> &Particle::getNeighbours() const {
  return neighbours_;
}

void Particle::pushBackNeighbour(bool diag, long particle) {
  neighbours_.emplace_back(diag, particle);
}

void Particle::updateForceInTime() {
  old_f_ = f_;
  f_ = {0., 0., 0.};
}

int Particle::getId() const { return id_; }

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x_ << " v: " << v_ << " f: " << f_
         << " old_f: " << old_f_ << " type: " << type_;
  return stream.str();
}

bool Particle::operator==(const Particle &other) const {
  return (x_ == other.x_) and (v_ == other.v_) and (f_ == other.f_) and
         (type_ == other.type_) and (m_ == other.m_) and
         (old_f_ == other.old_f_);
}

std::ostream &operator<<(std::ostream &stream, const Particle &p) {
  stream << p.toString();
  return stream;
}
