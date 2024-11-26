#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "defs/containers/DirectSumContainer.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/CuboidReader.h"
#include "io/file/in/xml/XmlReader.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "spdlog/fmt/bundled/chrono.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

int main(const int argc, char* argv[]) {
  SpdWrapper::get()->info("Application started");

  Arguments arguments = {
      .t_end = 5,
      .delta_t = 0.0002,
      .force_type = Arguments::LennardJones,
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

  //  TODO: Should we change this so the reader only reads configs? probably.
  std::vector<Particle> particles;

  XmlReader::read(particles, input_file, arguments);

  printConfiguration(arguments);
  SpdWrapper::get()->info("Step size: {}", step_size);
  // TODO: we should pass step size to reader but right now its useful for
  // testing

  std::unique_ptr<ParticleContainer> container;
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
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

  std::unique_ptr<Force> force;
  if (arguments.force_type == Arguments::Gravity) {
    force = std::make_unique<Gravity>();
  } else if (arguments.force_type == Arguments::LennardJones) {
    force = std::make_unique<LennardJones>();
  }

  VerletIntegrator verlet_integrator(*force, arguments.delta_t);
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
