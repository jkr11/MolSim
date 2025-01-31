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
  static int global_id_counter_;
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
   * depth of the potential-well
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
   * neighbouring cells for the membranes, saved as [diagonal, address]
   */
  std::vector<std::pair<bool, size_t>> neighbours_;

 public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(const std::array<double, 3> &x_arg,
           const std::array<double, 3> &v_arg, double m_arg, double epsilon,
           double sigma, int type = 0);

  explicit Particle(const std::array<double, 3> &x_arg,
                    const std::array<double, 3> &v_arg,
                    const std::array<double, 3> &f_arg,
                    const std::array<double, 3> &old_f_arg, double m_arg,
                    int type_arg, double epsilon_arg, double sigma_arg);

  ~Particle();

  [[nodiscard]] const std::array<double, 3> &getX() const { return x_; }

  [[nodiscard]] const std::array<double, 3> &getV() const { return v_; }

  [[nodiscard]] const std::array<double, 3> &getF() const { return f_; }

  [[nodiscard]] const std::array<double, 3> &getOldF() const { return old_f_; }

  [[nodiscard]] double getM() const { return m_; }

  [[nodiscard]] int getType() const { return type_; }

  [[nodiscard]] double getEpsilon() const { return epsilon_; }

  [[nodiscard]] double getSigma() const { return sigma_; }

  [[nodiscard]] const std::vector<std::pair<bool, size_t>> &getNeighbours()
      const;

  void setF(const std::array<double, 3> &f) { f_ = f; }

  void setX(const std::array<double, 3> &x) { x_ = x; }

  void setV(const std::array<double, 3> &v) { v_ = v; }

  void setOldF(const dvec3 &old_f) { old_f_ = old_f; }

  void setEpsilon(const double &epsilon) { epsilon_ = epsilon; }

  void setSigma(const double &sigma) { sigma_ = sigma; }

  /**
   * @brief Push Neighbour reference to a particle
   * @param diag Whether neighbour is diagonal
   * @param particle Particle reference in form of an long
   */
  void pushBackNeighbour(bool diag, long particle);

  /**
   * @brief OldF = F; F = 0
   */
  void updateForceInTime();

  void subV(const dvec3 &v) { v_ = v_ - v; }

  void addV(const dvec3 &v) { v_ = v_ + v; }

  void mulV(const double &scalar) { v_ = scalar * v_; }

  void addF(const dvec3 &f) { f_ = f_ + f; }

  void subF(const dvec3 &f) { f_ = f_ - f; }

  /**
   * @brief Directly set neighbours of a particle
   * @param new_neighbours Vector of neighbours (bool == diagonal, size_t ==
   * address)
   */
  void setNeighbours(
      const std::vector<std::pair<bool, size_t>> &new_neighbours) {
    neighbours_ = new_neighbours;
  }

  /**
   * @brief Clear neighbour references
   */
  void resetNeighbours() { neighbours_ = {}; }

  /**
   * @brief Retrieve particle id
   * @return particle id
   */
  [[nodiscard]] int getId() const;

  /**
   * @brief Comparison operator
   * @param other Other particle to be compare against
   * @return Whether particle is seen as equivalent
   */
  bool operator==(const Particle &other) const;

  /**
   * @brief Get debug string from particle
   * @return String
   */
  [[nodiscard]] std::string toString() const;

  /**
   * TODO
   */
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
      neighbours_ = other.neighbours_;  // Shallow copy of shared_ptr
    }
    return *this;
  }

  /**
   * TODO
   */
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
      neighbours_ =
          std::move(other.neighbours_);  // Transfer ownership of shared_ptr
    }
    return *this;
  }
};

/**
 * @brief Write Particle to stream
 * @param stream Output stream
 * @param p particle to be written
 * @return passed in stream
 */
std::ostream &operator<<(std::ostream &stream, const Particle &p);
