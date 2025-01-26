/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>
#include <mutex>

#include "defs/types.h"

class Particle final {
 private:
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
   * Mutex to ensure proper force calculation
   */
  std::mutex mutex;

 public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg,
      double m_arg, double _epsilon, double _sigma, int type = 0);

  explicit Particle(const std::array<double, 3> &x_arg,
                    const std::array<double, 3> &v_arg,
                    const std::array<double, 3> &f_arg,
                    const std::array<double, 3> &old_f_arg, double m_arg,
                    int type_arg, double epsilon_arg, double sigma_arg);

  ~Particle();

  [[nodiscard]] const std::array<double, 3> &getX() const;

  [[nodiscard]] const std::array<double, 3> &getV() const;

  [[nodiscard]] const std::array<double, 3> &getF() const;

  [[nodiscard]] const std::array<double, 3> &getOldF() const;

  [[nodiscard]] double getM() const;

  [[nodiscard]] int getType() const;

  [[nodiscard]] double getEpsilon() const;

  [[nodiscard]] double getSigma() const;

  void setF(const std::array<double, 3> &F);

  void setX(const std::array<double, 3> &X);

  void setV(const std::array<double, 3> &V);

  void setOldF(const dvec3 &oF);

  void setEpsilon(const double &epsilon);

  void setSigma(const double &sigma);

  void updateForceInTime();

  void subV(const dvec3 &V);

  void addV(const dvec3 &V);

  void mulV(const double &scalar);

  void addF(const dvec3 &F);

  void subF(const dvec3 &F);

  bool operator==(const Particle &other) const;


  [[nodiscard]] std::string toString() const;

  Particle &operator=(const Particle &o) {
    x = o.x;
    v = o.v;
    f = o.f;
    old_f = o.old_f;
    m = o.m;
    type = o.type;
    epsilon = o.epsilon;
    sigma = o.sigma;
    //skip mutex
    return *this;
  }

};

std::ostream &operator<<(std::ostream &stream, const Particle &p);
