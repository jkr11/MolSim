//
// Created by jkr on 12/19/24.
//
#include <gtest/gtest.h>

#include "../src/defs/Thermostat.h"
#include "calc/VerletIntegrator.h"
#include "defs/containers/DirectSumContainer.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "forces/SingularGravity.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/xml/XmlReader.h"

/**
 *  @brief tests holding the temperature at 1.0
 */
TEST(Thermostat, holding) {
  char arg0[] = "MolSim";
  char arg1[] = "-f";
  char arg2[] = "../../tests/test_cuboid.xml";
  char* argv[] = {arg0, arg1, arg2};
  auto [name, step, checkpoint] = CLArgumentParser::parse(3, argv);

  Arguments arguments;
  std::vector<Particle> particles;

  XmlReader::read(particles, name, arguments);

  std::unique_ptr<ParticleContainer> container;
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    container = std::make_unique<LinkedCellsContainer>(linked_cells_data);
    container->addParticles(particles);
    container->imposeInvariant();
  } else if (std::holds_alternative<DirectSumConfig>(
                 arguments.container_data)) {
    container = std::make_unique<DirectSumContainer>();
    container->addParticles(particles);
  } else {
    SpdWrapper::get()->error("Unrecognized container type");
    throw std::runtime_error("Unrecognized container type");
  }
  std::cout << particles.size() << " particles" << std::endl;
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  for (auto& config : arguments.interactive_force_types) {
    if (std::holds_alternative<LennardJonesConfig>(config)) {
      interactive_forces.push_back(std::make_unique<LennardJones>());
    } else if (std::holds_alternative<GravityConfig>(config)) {
      interactive_forces.push_back(std::make_unique<Gravity>());
    } else {
      SpdWrapper::get()->error("Unrecognized interactive_force_type");
    }
  }

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  for (auto config : arguments.singular_force_types) {
    if (std::holds_alternative<SingularGravityConfig>(config)) {
      const auto& [g] = std::get<SingularGravityConfig>(config);
      singular_forces.push_back(
          std::move(std::make_unique<SingularGravity>(g)));
    } else {
      SpdWrapper::get()->error("Unrecognized singular force");
    }
  }

  VerletIntegrator verlet_integrator(interactive_forces, singular_forces,
                                     arguments.delta_t);

  std::unique_ptr<Thermostat> thermostat;
  if (arguments.use_thermostat) {
    thermostat = std::make_unique<Thermostat>(arguments.thermostat_config);
  }
  double cur_temp = Thermostat::getTemperature(*container);
  EXPECT_NEAR(cur_temp, 1, 1);
  for (std::size_t i = 0; i < 1000; i++) {
    verlet_integrator.step(*container);
  }

  double temp = Thermostat::getTemperature(*container);
  EXPECT_NE(temp, 0.0);

  thermostat->setTemperature(*container);

  EXPECT_NEAR(thermostat->getTemperature(*container), 1.0, 1e-6);

  for (std::size_t i = 0; i < 10000; i++) {
    if (i != 0 && i % thermostat->n_thermostat) {
      thermostat->setTemperature(*container);
      EXPECT_NEAR(thermostat->getTemperature(*container), thermostat->t_target,
                  1e-6);
    }
  }
}

/**
 *  @brief test cooling from 20 to 1.0 with absolute values
 */
TEST(Thermostat, cooling) {
  char arg0[] = "MolSim";
  char arg1[] = "-f";
  char arg2[] = "../../tests/test_cuboid_cooling.xml";
  char* argv[] = {arg0, arg1, arg2};
  auto [name, step, checkpoint] = CLArgumentParser::parse(3, argv);

  Arguments arguments;
  std::vector<Particle> particles;

  XmlReader::read(particles, name, arguments);

  std::unique_ptr<ParticleContainer> container;
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    container = std::make_unique<LinkedCellsContainer>(linked_cells_data);
    container->addParticles(particles);
    container->imposeInvariant();
  } else if (std::holds_alternative<DirectSumConfig>(
                 arguments.container_data)) {
    container = std::make_unique<DirectSumContainer>();
    container->addParticles(particles);
  } else {
    SpdWrapper::get()->error("Unrecognized container type");
    throw std::runtime_error("Unrecognized container type");
  }
  std::cout << particles.size() << " particles" << std::endl;
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  for (auto& config : arguments.interactive_force_types) {
    if (std::holds_alternative<LennardJonesConfig>(config)) {
      interactive_forces.push_back(std::make_unique<LennardJones>());
    } else if (std::holds_alternative<GravityConfig>(config)) {
      interactive_forces.push_back(std::make_unique<Gravity>());
    } else {
      SpdWrapper::get()->error("Unrecognized interactive_force_type");
    }
  }

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  for (auto config : arguments.singular_force_types) {
    if (std::holds_alternative<SingularGravityConfig>(config)) {
      const auto& [g] = std::get<SingularGravityConfig>(config);
      singular_forces.push_back(
          std::move(std::make_unique<SingularGravity>(g)));
    } else {
      SpdWrapper::get()->error("Unrecognized singular force");
    }
  }

  VerletIntegrator verlet_integrator(interactive_forces, singular_forces,
                                     arguments.delta_t);

  std::unique_ptr<Thermostat> thermostat;
  if (arguments.use_thermostat) {
    thermostat = std::make_unique<Thermostat>(arguments.thermostat_config);
  }
  double cur_temp = Thermostat::getTemperature(*container);
  EXPECT_NEAR(cur_temp, 22, 4);
  for (std::size_t i = 0; i < 1000; i++) {
    verlet_integrator.step(*container);
  }

  double temp = Thermostat::getTemperature(*container);
  EXPECT_NE(temp, 0.0);

  thermostat->setTemperature(*container);

  EXPECT_NEAR(thermostat->getTemperature(*container), 1.0, 1e-6);
}

/**
 *  @brief test heating from 1.0 to 10.0 with absolute values
 */
TEST(Thermostat, heating) {
  char arg0[] = "MolSim";
  char arg1[] = "-f";
  char arg2[] = "../../tests/test_cuboid_heating.xml";
  char* argv[] = {arg0, arg1, arg2};
  auto [name, step, checkpoint] = CLArgumentParser::parse(3, argv);

  Arguments arguments;
  std::vector<Particle> particles;

  XmlReader::read(particles, name, arguments);

  std::unique_ptr<ParticleContainer> container;
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    container = std::make_unique<LinkedCellsContainer>(linked_cells_data);
    container->addParticles(particles);
    container->imposeInvariant();
  } else if (std::holds_alternative<DirectSumConfig>(
                 arguments.container_data)) {
    container = std::make_unique<DirectSumContainer>();
    container->addParticles(particles);
  } else {
    SpdWrapper::get()->error("Unrecognized container type");
    throw std::runtime_error("Unrecognized container type");
  }
  std::cout << particles.size() << " particles" << std::endl;
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  for (auto& config : arguments.interactive_force_types) {
    if (std::holds_alternative<LennardJonesConfig>(config)) {
      interactive_forces.push_back(std::make_unique<LennardJones>());
    } else if (std::holds_alternative<GravityConfig>(config)) {
      interactive_forces.push_back(std::make_unique<Gravity>());
    } else {
      SpdWrapper::get()->error("Unrecognized interactive_force_type");
    }
  }

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  for (auto config : arguments.singular_force_types) {
    if (std::holds_alternative<SingularGravityConfig>(config)) {
      const auto& [g] = std::get<SingularGravityConfig>(config);
      singular_forces.push_back(
          std::move(std::make_unique<SingularGravity>(g)));
    } else {
      SpdWrapper::get()->error("Unrecognized singular force");
    }
  }

  VerletIntegrator verlet_integrator(interactive_forces, singular_forces,
                                     arguments.delta_t);

  std::unique_ptr<Thermostat> thermostat;
  if (arguments.use_thermostat) {
    thermostat = std::make_unique<Thermostat>(arguments.thermostat_config);
  }
  double cur_temp = Thermostat::getTemperature(*container);
  EXPECT_NEAR(cur_temp, 1.0, 1);
  for (std::size_t i = 0; i < 1000; i++) {
    verlet_integrator.step(*container);
  }

  double temp = Thermostat::getTemperature(*container);
  EXPECT_NE(temp, 0.0);

  thermostat->setTemperature(*container);

  EXPECT_NEAR(thermostat->getTemperature(*container), 10, 1e-6);
}

/*
 * this test the gradual thermostat, so when deltaT is set. In this case its 2.0
 * as specified, so we would expect that it converges after at 3k iterations (a
 * lot earlier actually but this is good enough) if n_thermostat is 100
 */
TEST(Thermostat, gradual) {
  char arg0[] = "MolSim";
  char arg1[] = "-f";
  char arg2[] = "../../tests/test_cuboid_gradual.xml";
  char* argv[] = {arg0, arg1, arg2};
  auto [name, step, checkpoint] = CLArgumentParser::parse(3, argv);

  Arguments arguments;
  std::vector<Particle> particles;

  XmlReader::read(particles, name, arguments);

  std::unique_ptr<ParticleContainer> container;
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    container = std::make_unique<LinkedCellsContainer>(linked_cells_data);
    container->addParticles(particles);
    container->imposeInvariant();
  } else if (std::holds_alternative<DirectSumConfig>(
                 arguments.container_data)) {
    container = std::make_unique<DirectSumContainer>();
    container->addParticles(particles);
  } else {
    SpdWrapper::get()->error("Unrecognized container type");
    throw std::runtime_error("Unrecognized container type");
  }
  std::cout << particles.size() << " particles" << std::endl;
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  for (auto& config : arguments.interactive_force_types) {
    if (std::holds_alternative<LennardJonesConfig>(config)) {
      interactive_forces.push_back(std::make_unique<LennardJones>());
    } else if (std::holds_alternative<GravityConfig>(config)) {
      interactive_forces.push_back(std::make_unique<Gravity>());
    } else {
      SpdWrapper::get()->error("Unrecognized interactive_force_type");
    }
  }

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  for (auto config : arguments.singular_force_types) {
    if (std::holds_alternative<SingularGravityConfig>(config)) {
      const auto& [g] = std::get<SingularGravityConfig>(config);
      singular_forces.push_back(
          std::move(std::make_unique<SingularGravity>(g)));
    } else {
      SpdWrapper::get()->error("Unrecognized singular force");
    }
  }

  VerletIntegrator verlet_integrator(interactive_forces, singular_forces,
                                     arguments.delta_t);

  std::unique_ptr<Thermostat> thermostat;
  if (arguments.use_thermostat) {
    thermostat = std::make_unique<Thermostat>(arguments.thermostat_config);
  }
  double cur_temp = Thermostat::getTemperature(*container);
  EXPECT_NEAR(cur_temp, 1.0, 1);
  int diffs_added = 0;
  for (std::size_t i = 0; i <= 3000; i++) {
    verlet_integrator.step(*container);
    if (i > 0 && i % thermostat->n_thermostat == 0) {
      std::cout << i << "Iterations " << std::endl;
      diffs_added++;
      double temp = Thermostat::getTemperature(*container);
      thermostat->setTemperature(*container);
      EXPECT_NEAR(thermostat->getTemperature(*container),
                  temp + thermostat->delta_temp, thermostat->delta_temp);
    }
  }

  EXPECT_NEAR(thermostat->getTemperature(*container), thermostat->t_target, 1);
}
