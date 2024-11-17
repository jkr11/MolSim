#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "defs/containers/DirectSumContainer.h"
#include "forces/LennardJones.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/CuboidReader.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "spdlog/stopwatch.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

int main(int argc, char *argv[]) {
  SpdWrapper::get()->info("Application started");

  Arguments arguments = {
      "",                                // file
      10,                                // t_end
      0.014,                             // delta_t
      1,                                 // output_time_step_size
      "info",                            // logLevel
      std::make_unique<LennardJones>(),  // force
      std::make_unique<CuboidReader>(),  // Reader
  };

  if (CLArgumentParser::parse(argc, argv, arguments) != 0) {
    exit(EXIT_FAILURE);
  }

  SpdWrapper::get()->info("t_end: {}, delta_t: {}, output_time_step_size: {}",
                          arguments.t_end, arguments.delta_t,
                          arguments.output_time_step_size);

  // setup Simulation
  // FileReader::readFile(particles, input_file);
  DirectSumContainer container;
  /*TODO: arguments.reader->read(particleContainer.getParticlesReference(),
                         arguments.inputFile);*/

  VerletIntegrator verlet_integrator(*arguments.force, arguments.delta_t);
  outputWriter::VTKWriter writer;

  const std::string outputDirectory =
      createOutputDirectory("./output/", argc, argv);

  const spdlog::stopwatch sw;
  double current_time = 0;  // start time is always 0
  int iteration = 0;
  int writes = 0;
  int percentage = 0;
  double next_output_time = 0;

  while (current_time <= arguments.t_end) {
    verlet_integrator.step(container);

    if (current_time >= next_output_time) {
      plotParticles(outputDirectory, iteration, writer, container);
      writes++;
      next_output_time = writes * arguments.output_time_step_size;

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
