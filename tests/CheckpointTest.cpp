//
// Created by jkr on 1/22/25.
//
#include <gtest/gtest.h>

#include <vector>

#include "calc/VerletIntegrator.h"
#include "defs/Simulation.h"
#include "defs/containers/DirectSumContainer.h"
//#include "defs/containers/LinkedCellsContainer.cpp"
#include "defs/containers/LinkedCellsContainer.h"
#include "defs/containers/ParticleContainer.h"
#include "defs/types.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "io/file/in/xml/XmlReader.h"
#include "io/file/out/XmlWriter.h"
#include "testUtil.h"

/**
 *  @brief tests holding the temperature at 1.0
 */
TEST(Checkpoint, cuboid) {
  char arg0[] = "./MolSim";
  char arg1[] = "-f";
  char arg2[] = "../../../tests/checkpoint_input_test.xml";
  char* argv[] = {arg0, arg1, arg2};
   auto [name, step, write_checkpoint] = CLArgumentParser::parse(3, argv);
  //const char * name = arg2;

  Arguments arguments;
  std::vector<Particle> particles;

  XmlReader::read(particles, name, arguments);

  LinkedCellsConfig config = {.domain = {300, 300, 1},
                              .cutoff_radius = 3.0,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  const std::unique_ptr<ParticleContainer> container =
      std::make_unique<LinkedCellsContainer>(config);
  container->addParticles(particles);
  container->imposeInvariant();

  std::cout << particles.size() << " particles" << std::endl;
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  interactive_forces.push_back(std::make_unique<LennardJones>());


  std::vector<std::unique_ptr<SingularForce>> singular_forces;

  VerletIntegrator verlet_integrator(interactive_forces, singular_forces,
                                     0.1);

  if constexpr (true) {
    XmlWriter::writeFile(*container, "../../../input/checkpoint_test.xml");
  }
  std::cout << "New file trying .... " << std::endl;
  char arg01[] = "MolSim";
  char arg11[] = "-f";
  char arg21[] = "../../../tests/checkpoint_output_test.xml";
  char arg31[] ="-c";
  char* argv1[] = {arg01, arg11, arg21, arg31};
  auto [name1, step1, write_checkpoint1] = CLArgumentParser::parse(4, argv1);

  Arguments arguments1;
  std::vector<Particle> particles1;
  XmlReader::read(particles1, name1, arguments1);

  ASSERT_EQ(particles.size(), particles1.size());

  std::sort(particles.begin(), particles.end(), [](const Particle& a, const Particle& b) {
        return a.getX() < b.getX();
    });

  std::sort(particles1.begin(), particles1.end(), [](const Particle& a, const Particle& b) {
        return a.getX() < b.getX();
    });

  for (size_t i = 0; i < particles.size(); ++i) {
    ASSERT_EQ_VEC3(particles[i].getX(), particles1[i].getX(), "Vectors not equal at index " + std::to_string(i));
  }

}
