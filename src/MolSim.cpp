<<<<<<< HEAD

=======
>>>>>>> MolSim/test1
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

<<<<<<< HEAD
#include <iostream>
#include <list>
=======
#include <defs/Particle.h>
#include <array>
#include <iostream>
#include <list>
#include <vector>

#include "calc/Verlet.h"
#include "forces/gravity.h"
#include "outputWriter/VTKWriter.h"
>>>>>>> MolSim/test1

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
<<<<<<< HEAD
=======


>>>>>>> MolSim/test1
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
<<<<<<< HEAD
void plotParticles(int iteration);
=======
void plotParticles(int iteration, outputWriter::VTKWriter& writer, ParticleContainer& particle_container);
>>>>>>> MolSim/test1

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::list<Particle> particles;
<<<<<<< HEAD

int main(int argc, char *argsv[]) {

=======
Gravity gravity;

int main(int argc, char* argsv[]) {
>>>>>>> MolSim/test1
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);
<<<<<<< HEAD
=======
  ParticleContainer particle_container(particles);
  VerletIntegrator verlet_integrator(gravity, delta_t);
  outputWriter::VTKWriter writer;
>>>>>>> MolSim/test1

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
<<<<<<< HEAD
=======
    /*
>>>>>>> MolSim/test1
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();
<<<<<<< HEAD

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
=======
    */

    verlet_integrator.step(particle_container);
    iteration++;
    if (iteration % 50 == 0) {
      plotParticles(iteration, writer, particle_container);
>>>>>>> MolSim/test1
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

<<<<<<< HEAD
void calculateF() {
  std::list<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    for (auto &p2 : particles) {
      // @TODO: insert calculation of forces here!
=======
/*
void calculateF() {
  std::list<Particle>::iterator iterator;
  iterator = particles.begin();
  // double G = 6.674 / 1e11; // m^3 * kg^-1 * s^-2
  // actually, g = 1 for this task
  for (auto& p1 : particles) {
    for (auto& p2 : particles) {
      // @TODO: insert calculation of forces here!
      if (&p1 == &p2) { continue; }

      const auto& x2 = p2.getX();
      const double m2 = p2.getM();

      std::array<double, 3> r{};
      double dist = 0.0;
      for (int i = 0; i < 3; ++i) {
        r[i] = p2.getX()[i] - p1.getX()[i];
        dist += r[i] * r[i];
      }
      //dist = std::sqrt(dist);
      //what is faster? pow(1.5) of sqrt^3
      double F = (m1 * m2) / std::pow(dist, 1.5);
      std::array<double, 3> f1 = {};
      for (int i = 0; i < 3; ++i) {
        f1[i] = F * r[i];
      }

      auto g1 = gravity.directionalForce(p1, p2);
      std::cout << "Grav: " << g1 << std::endl;
      p1.setF(p1.getF() + g1);
>>>>>>> MolSim/test1
    }
  }
}

void calculateX() {
<<<<<<< HEAD
  for (auto &p : particles) {
    // @TODO: insert calculation of position updates here!
=======
  constexpr double dt = delta_t;
  for (auto& p : particles) {
    // @TODO: insert calculation of position updates here!
    auto& x = p.getX();
    const auto& v = p.getV();
    const auto& f = p.getF();
    double m = p.getM();
    std::array<double, 3> new_x{};

    new_x = x + dt * v + (dt * dt * 0.5 / m) * f;
    p.setX(new_x);
    p.setOldF(p.getF());
    p.setF({0, 0, 0});
>>>>>>> MolSim/test1
  }
}

void calculateV() {
<<<<<<< HEAD
  for (auto &p : particles) {
    // @TODO: insert calculation of veclocity updates here!
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
=======
  constexpr double dt = delta_t;
  for (auto& p : particles) {
    // @TODO: insert calculation of veclocity updates here!
    auto& v = p.getV();
    auto& F = p.getF();
    auto& oF = p.getOldF();
    double m = p.getM();
    std::array<double, 3> nv{};
    nv = v + (dt / (2 * m)) * (F * oF);
    p.setV(nv);
  }
}
*/
void plotParticles(int iteration, outputWriter::VTKWriter& vtkWriter, ParticleContainer& particle_container) {
  vtkWriter.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto& p : particle_container.getParticles()) {
    vtkWriter.plotParticle(p);
  }

  vtkWriter.writeFile("MD_vtk", iteration);
>>>>>>> MolSim/test1
}
