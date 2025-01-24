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

class Particle final {
 private:
  static int global_id_counter;
  /**
   * Position of the particle
   */
  std::array<double, 3> x{};

  /**
   * Velocity of the particle
   */
  std::array<double, 3> v{};

  /**
   * Force effective on this particle
   */
  std::array<double, 3> f{};

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_f{};

  /**
   * Mass of this particle
   */
  double m{};

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type{};

  /**
   * depth of the potential well
   * used in the caluclation of lennard-jones force
   */
  double epsilon{};

  /**
   * distance from the particle at which the potential is zero \n
   * used in the calculation of lennard-jones force
   */
  double sigma{};

  /**
   * unique identifier for every particle
   */
  int id{};

  /**
   * neighbouring cells for the membranes
   */
  std::vector<std::pair<bool, std::shared_ptr<Particle>>> neighbours;

 public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg,
      double m_arg, double _epsilon, double _sigma, int type = 0);

  Particle(Particle &&other) noexcept;  // Move constructor? cancer

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

  [[nodiscard]] const std::array<double, 3> &getX() const;

  [[nodiscard]] const std::array<double, 3> &getV() const;

  [[nodiscard]] const std::array<double, 3> &getF() const;

  [[nodiscard]] const std::array<double, 3> &getOldF() const;

  [[nodiscard]] double getM() const;

  [[nodiscard]] int getType() const;

  [[nodiscard]] double getEpsilon() const;

  [[nodiscard]] double getSigma() const;

  [[nodiscard]] const std::vector<std::pair<bool, std::shared_ptr<Particle>>> &
  getNeighbours() const;

  void setF(const std::array<double, 3> &F);

  void setX(const std::array<double, 3> &X);

  void setV(const std::array<double, 3> &V);

  void setOldF(const dvec3 &oF);

  void setEpsilon(const double &epsilon);

  void setSigma(const double &sigma);

  void pushBackNeighbour(bool diag, const std::shared_ptr<Particle> &particle);

  void updateForceInTime();

  void subV(const dvec3 &V);

  void addV(const dvec3 &V);

  void mulV(const double &scalar);

  void addF(const dvec3 &F);

  void subF(const dvec3 &F);

  [[nodiscard]] int getId() const;

  bool operator==(const Particle &other) const;

  [[nodiscard]] std::string toString() const;

  Particle &operator=(const Particle &other) {
    if (this != &other) {
      x = other.x;
      v = other.v;
      f = other.f;
      old_f = other.old_f;
      m = other.m;
      type = other.type;
      id = other.id;
      epsilon = other.epsilon;
      sigma = other.sigma;
      neighbours = other.neighbours;  // Shallow copy of shared_ptr
    }
    return *this;
  }

  // Move assignment operator
  Particle &operator=(Particle &&other) noexcept {
    if (this != &other) {
      x = other.x;
      v = other.v;
      f = other.f;
      old_f = other.old_f;
      m = other.m;
      type = other.type;
      id = other.id;
      epsilon = other.epsilon;
      sigma = other.sigma;
      neighbours =
          std::move(other.neighbours);  // Transfer ownership of shared_ptr
    }
    return *this;
  }
};

std::ostream &operator<<(std::ostream &stream, const Particle &p);
