/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "defs/types.h"
#include "utils/ArrayUtils.h"

class Particle final {
 private:
  static int global_id_counter;
  /**
   * Position of the particle
   */
  std::array<double, 3> x_{};

  /**
   * Velocity of the particle
   */
  std::array<double, 3> v_{};

  /**
   * Force effective on this particle
   */
  std::array<double, 3> f_{};

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_f_{};

  /**
   * Mass of this particle
   */
  double m_{};

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   * negative types are immovable and are being ignored by some calculations
   */
  int type_{};

  /**
   * depth of the potential well
   * used in the caluclation of lennard-jones force
   */
  double epsilon_{};

  /**
   * distance from the particle at which the potential is zero \n
   * used in the calculation of lennard-jones force
   */
  double sigma_{};

  /**
   * unique identifier for every particle
   */
  int id_{};

  /**
   * neighbouring cells for the membranes
   */
  std::vector<std::pair<bool, size_t>> neighbours;

 public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg,
      double m_arg, double epsilon, double sigma, int type = 0);


  explicit Particle(const std::array<double, 3> &x_arg,
                    const std::array<double, 3> &v_arg,
                    const std::array<double, 3> &f_arg,
                    const std::array<double, 3> &old_f_arg, double m_arg,
                    int type_arg, double epsilon_arg, double sigma_arg);
  // ambiguous overload?
  // Particle()
  //    : neighbours(std::vector<std::pair<bool, std::shared_ptr<Particle>>>())
  //    {}

  ~Particle();

  [[nodiscard]] const std::array<double, 3> &getX() const { return x_; }

  [[nodiscard]] const std::array<double, 3> &getV() const { return v_; }

  [[nodiscard]] const std::array<double, 3> &getF() const { return f_; }

  [[nodiscard]] const std::array<double, 3> &getOldF() const { return old_f_; }

  [[nodiscard]] double getM() const { return m_; }

  [[nodiscard]] int getType() const { return type_; }

  [[nodiscard]] double getEpsilon() const { return epsilon_; }

  [[nodiscard]] double getSigma() const { return sigma_; }

  [[nodiscard]] const std::vector<std::pair<bool,size_t>> &getNeighbours()
      const;

  void setF(const std::array<double, 3> &f) { f_ = f; }

  void setX(const std::array<double, 3> &x) { x_ = x; }

  void setV(const std::array<double, 3> &v) { v_ = v; }

  void setOldF(const dvec3 &old_f) { old_f_ = old_f; }

  void setEpsilon(const double &epsilon) { epsilon_ = epsilon; }

  void setSigma(const double &sigma) { sigma_ = sigma; }

  void pushBackNeighbour(bool diag, long particle);

  void updateForceInTime();

  void subV(const dvec3 &v) { v_ = v_ - v; }

  void addV(const dvec3 &v) { v_ = v_ + v; }

  void mulV(const double &scalar) { v_ = scalar * v_; }

  void addF(const dvec3 &f) { f_ = f_ + f; }

  void subF(const dvec3 &f) { f_ = f_ - f; }

  void setNeighbours(
      const std::vector<std::pair<bool, size_t>> &new_neighbours) {
    neighbours = new_neighbours;
  }

  void resetNeighbours() { neighbours = {}; }

  [[nodiscard]] int getId() const;

  bool operator==(const Particle &other) const;

  [[nodiscard]] std::string toString() const;

  Particle &operator=(const Particle &other) {
    if (this != &other) {
      x_ = other.x_;
      v_ = other.v_;
      f_ = other.f_;
      old_f_ = other.old_f_;
      m_ = other.m_;
      type_ = other.type_;
      id_ = other.id_;
      epsilon_ = other.epsilon_;
      sigma_ = other.sigma_;
      neighbours = other.neighbours;  // Shallow copy of shared_ptr
    }
    return *this;
  }

  // Move assignment operator
  Particle &operator=(Particle &&other) noexcept {
    if (this != &other) {
      x_ = other.x_;
      v_ = other.v_;
      f_ = other.f_;
      old_f_ = other.old_f_;
      m_ = other.m_;
      type_ = other.type_;
      id_ = other.id_;
      epsilon_ = other.epsilon_;
      sigma_ = other.sigma_;
      neighbours =
          std::move(other.neighbours);  // Transfer ownership of shared_ptr
    }
    return *this;
  }
};

std::ostream &operator<<(std::ostream &stream, const Particle &p);
