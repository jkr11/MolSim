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
#include "spdlog/stopwatch.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"
int main(const int argc, char* argv[]) {
#ifndef BENCHMARK
  SpdWrapper::get()->info("Application started");
#endif
  Arguments arguments = {
      .t_end = 5,
      .delta_t = 0.0002,
      .force_type = Arguments::LennardJones,
      .container_data =
          LinkedCellsConfig{.domain = {100, 100, 100},
                            .cutoff_radius = 3.0,
                            .boundary_config =
                                {
                                    .x_high = LinkedCellsConfig::Outflow,
                                    .x_low = LinkedCellsConfig::Outflow,
                                    .y_high = LinkedCellsConfig::Outflow,
                                    .y_low = LinkedCellsConfig::Outflow,
                                    .z_high = LinkedCellsConfig::Outflow,
                                    .z_low = LinkedCellsConfig::Outflow,
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
  spdlog::stopwatch stopwatch;
  const auto start_time = std::chrono::high_resolution_clock::now();
  while (current_time <= arguments.t_end) {
    verlet_integrator.step(*container);
#ifndef BENCHMARK
    if (current_time >= next_output_time) {
      plotParticles(outputDirectory, iteration, writer, *container);
      writes++;
      next_output_time = writes * step_size;
      // check if next percentage complete
      if (const double t = 100 * current_time / arguments.t_end;
          t >= percentage) {
        percentage++;
        auto elapsed = stopwatch.elapsed();
        auto eta = (elapsed / percentage) * 100 - elapsed;
        auto h = std::chrono::duration_cast<std::chrono::hours>(eta).count();
        eta -= std::chrono::hours(h);
        auto m = std::chrono::duration_cast<std::chrono::minutes>(eta).count();
        eta -= std::chrono::minutes(m);
        auto s = std::chrono::duration_cast<std::chrono::seconds>(eta).count();

        SpdWrapper::get()->info(
            "[{:.0f} %]: Iteration {:<12} | [ETA: {}:{:02}:{:02}]",
            100 * current_time / arguments.t_end, iteration, h, m, s);
      }
    }
#endif
    iteration++;
    current_time = arguments.delta_t * iteration;
  }
  const auto end_time = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed_time = end_time - start_time;
#ifdef BENCHMARK
  std::cout << "Simulation Time: " << elapsed_time.count() << " seconds"
            << std::endl;
#endif
  SpdWrapper::get()->info("Output written. Terminating...");

  return 0;
}
