#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "defs/Simulation.h"
#include "defs/containers/DirectSumContainer.h"
#include "forces/LennardJones.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/CuboidReader.h"
#include "io/file/in/xml/XmlReader.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "spdlog/stopwatch.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

int main(int argc, char *argv[]) {
  SpdWrapper::get()->info("Application started");
  // ok so this is really ugly but i think we can merge both structs into one
  // we dont need additional file readers anymore so we can just use XMLReader
  Simulation simulation_params = {
      .delta_t = 0.014,
      .t_end = 10,
      .cutoff_radius = 0.05,
      .domain = {10, 10, 10},
  };
  // This is looking a lot better now, we can probably merge these now
  Arguments arguments = {
      "",                                // file
      10,                                // t_end
      0.014,                             // delta_t
      1,                                 // output_time_step_size
      "info",                            // logLevel
      std::make_unique<LennardJones>(),  // force
  };

  if (CLArgumentParser::parse(argc, argv, arguments) != 0) {
    exit(EXIT_FAILURE);
  }
  const auto reader = std::make_unique<XmlReader>(simulation_params);
  // SpdWrapper::get()->info("t_end: {}, delta_t: {}, output_time_step_size:
  // {}",
  //                         arguments.t_end, arguments.delta_t,
  //                        arguments.output_time_step_size);

  std::vector<Particle> particles;
  reader->read(particles, arguments.inputFile);

  auto [delta_t, t_end, cutoff_radius, domain] = reader.get()->pass();

  SpdWrapper::get()->info("t_end: {}, delta_t: {}, cutoff_radius: {}", t_end,
                          delta_t, cutoff_radius);

  DirectSumContainer container(particles);
  SpdWrapper::get()->info("particles.size: {}", particles.size());
  VerletIntegrator verlet_integrator(*arguments.force, delta_t);
  outputWriter::VTKWriter writer;

  const std::string outputDirectory =
      createOutputDirectory("./output/", argc, argv);

  const spdlog::stopwatch sw;
  double current_time = 0;  // start time is always 0
  int iteration = 0;
  int writes = 0;
  int percentage = 0;
  double next_output_time = 0;

  while (current_time <= t_end) {
    verlet_integrator.step(container);

    if (current_time >= next_output_time) {
      plotParticles(outputDirectory, iteration, writer, container);
      writes++;
      next_output_time = writes * arguments.output_time_step_size;

      // check if next percentage complete
      if (const double t = 100 * current_time / t_end; t >= percentage) {
        percentage++;
        SpdWrapper::get()->info("[{:.0f} %]: Iteration {}",
                                100 * current_time / t_end, iteration);
      }
    }

    iteration++;
    current_time = delta_t * iteration;  // + start_time
  }
  std::cout << std::endl;
  SpdWrapper::get()->info("Output written. Terminating...");

  return 0;
}
