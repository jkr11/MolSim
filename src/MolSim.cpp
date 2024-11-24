#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "defs/containers/DirectSumContainer.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "forces/SingularGravity.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/CuboidReader.h"
#include "io/file/in/xml/XmlReader.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "spdlog/fmt/bundled/chrono.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

int main(int argc, char* argv[]) {
  SpdWrapper::get()->info("Application started");

  Arguments arguments = {
      //.input_file = "",
      .t_end = 5,
      .delta_t = 0.0002,
      //.output_time_step_size = 1,
      .log_level = "info",
      .force_type = Arguments::LennardJones,
      .singular_force_type = Arguments::SingularGravity,
      .container_type = Arguments::LinkedCells,
      .singular_force_data = SingularGravityConfig{.g = 11.4},
      .container_data =
          LinkedCellsConfig{.domain = {100, 100, 100},
                            .cutoff_radius = 3.0,
                            .boundary_config =
                                {
                                    .north = LinkedCellsConfig::Outflow,
                                    .south = LinkedCellsConfig::Outflow,
                                    .east = LinkedCellsConfig::Outflow,
                                    .west = LinkedCellsConfig::Outflow,
                                    .up = LinkedCellsConfig::Outflow,
                                    .down = LinkedCellsConfig::Outflow,
                                }},
  };
  auto [input_file, step_size] = CLArgumentParser::parse(argc, argv);
  SpdWrapper::get()->info("Step size: {}", step_size);
  //  TODO: Should we change this so it doesnt get read here but the reader
  //  instantiates the container and then writes the shapes to the container?
  std::vector<Particle> particles;

  XmlReader::read(particles, input_file, arguments);

  arguments.printConfiguration();

  std::unique_ptr<ParticleContainer> container;
  // So this exists pretty cool we can save another field in the struct
  // this is also much safer I guess
  // TODO: i really hope lrz is c++ 17 and not 14 else we can revert all of this
  // yikes
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    // TODO: make it possible to just pass the entire struct into this or
    // write a translation unit
    container = std::make_unique<LinkedCellsContainer>(
        linked_cells_data.domain, linked_cells_data.cutoff_radius);
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

  std::unique_ptr<BidirectionalForce> force;
  if (arguments.force_type == Arguments::Gravity) {
    force = std::make_unique<Gravity>();
  } else if (arguments.force_type == Arguments::LennardJones) {
    force = std::make_unique<LennardJones>();
  }

  std::unique_ptr<SingularForce> singular_force;
  if (std::holds_alternative<SingularGravityConfig>(
          arguments.singular_force_data)) {
    const auto& [g] =
        std::get<SingularGravityConfig>(arguments.singular_force_data);
    singular_force = std::make_unique<SingularGravity>(g);
  } else {
    SpdWrapper::get()->error("Unrecognized singular force");
  }

  VerletIntegrator verlet_integrator(*force, *singular_force,
                                     arguments.delta_t);
  outputWriter::VTKWriter writer;

  const std::string outputDirectory =
      createOutputDirectory("./output/", argc, argv);

  double current_time = 0;
  int iteration = 0;
  int writes = 0;
  int percentage = 0;
  double next_output_time = 0;

  while (current_time <= arguments.t_end) {
    verlet_integrator.step(*container);

    if (current_time >= next_output_time) {
      plotParticles(outputDirectory, iteration, writer, *container);
      writes++;
      next_output_time = writes * step_size;
      // check if next percentage complete
      if (const double t = 100 * current_time / arguments.t_end;
          t >= percentage) {
        percentage++;
        SpdWrapper::get()->info("[{:.0f} %]: Iteration {}",
                                100 * current_time / arguments.t_end,
                                iteration);
      }
    }

    iteration++;
    current_time = arguments.delta_t * iteration;  // + start_time
  }
  std::cout << std::endl;
  SpdWrapper::get()->info("Output written. Terminating...");

  return 0;
}
