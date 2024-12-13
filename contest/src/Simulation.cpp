//
// Created by mcarn on 12/11/24.
//

#include <iostream>
#include <chrono>

#include "Simulation.h"
#include "CuboidGenerator.h"
#include "ParticleContainer.h"
#include "LennardJones.h"
#include "VerletIntegrator.h"
#include "VTKWriter.h"
#include "../../src/defs/types.h"

void Simulation::run(std::string filepath) {
  //TODO: get info from xml -> defered

  ParticleTypeInfo typeinfo = {
    .mass = {1, 2},
    .sigma = {1.2, 1.1},
    .epsilon = {1, 1}
  };

  LinkedCellsConfig config = {
                        .domain = {300, 54, 1},
                        .cutoff_radius = 3.0,
                        .boundary_config = {
                            .x_high = LinkedCellsConfig::Outflow,
                            .x_low = LinkedCellsConfig::Outflow,
                            .y_high = LinkedCellsConfig::Outflow,
                            .y_low = LinkedCellsConfig::Outflow,
                            .z_high = LinkedCellsConfig::Outflow,
                            .z_low = LinkedCellsConfig::Outflow,
  }};


  //init container
  ParticleContainer container(typeinfo, config, 10000);


  //TODO: get info from xml -> defered
  CuboidGenerator gen1({0.6, 2, 0}, {250, 20, 1}, 1.2, {0, 0, 0}, 0.1, 0, true);
  gen1.generate(container);

  CuboidGenerator gen2({0.6, 27, 0}, {250, 20, 1}, 1.2, {0, 0, 0}, 0.1, 1, true);
  gen2.generate(container);

  //impose invariant to get correct setup after adding particles
  container.imposeInvariant();

  //init integrator
  std::vector<std::unique_ptr<InteractiveForce>> interactiveForces;
  interactiveForces.push_back(std::make_unique<LennardJones>());

  std::vector<std::unique_ptr<SingularForce>> singularForces;

  double delta_t = 0.0005;

  VerletIntegrator integrator(interactiveForces, singularForces, delta_t);


  //init output writer
  VTKWriter outputWriter{};


  //run sim in while loop
  double end_time = 1; //TODO
  double current_time = 0;
  int iteration = 0;
  const auto startTime = std::chrono::high_resolution_clock::now();

  while (current_time <= end_time) {
    if (iteration == 1000) {
      const auto endTime = std::chrono::high_resolution_clock::now();

      const std::chrono::duration<double> elapsed = endTime - startTime;

      std::cout << elapsed.count() << " seconds" << std::endl;
      exit(0);
    }

    if (iteration % 100 == 0) {
      //std::cout << "Iteration: " << iteration << std::endl;

      /*outputWriter.initializeOutput(container.ids.size());

      for (int i = 0; i < container.ids.size(); i++) {
        const dvec3 pos = {container.px[i], container.py[i], container.pz[i]};
        const dvec3 vel = {container.vx[i], container.vy[i], container.vz[i]};
        const dvec3 of = {container.ofx[i], container.ofy[i], container.ofz[i]};
        const int type = container.types[i];

        outputWriter.plotParticle(pos, vel, of, type);
      }

      outputWriter.writeFile("output/gore", iteration);*/
    }

    integrator.step(container);

    iteration++;
    current_time = delta_t * iteration;
  }

  std::cout << "And thus we are done." << std::endl;
}

