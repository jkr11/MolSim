#include <chrono>
#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "forces/Force.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/CuboidReader.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

// TODO: move
constexpr int output_interval = 32;

Gravity gravity;
LennardJones lennardjones;

int main(int argc, char *argv[]) {
  SpdWrapper::get()->info("Application started");

  Arguments arguments = {
      "",      // file
      100,     // t_end
      0.014,   // delta_t
      1,       // output_time_step_size
      "info",  // logLevel
  };

  if (CLArgumentParser::parse(argc, argv, arguments) != 0) {
    exit(EXIT_FAILURE);
  }

  SpdWrapper::get()->info("t_end: {}, delta_t: {}, output_time_step_size: {}",
                          arguments.t_end, arguments.delta_t,
                          arguments.output_time_step_size);

  // setup Simulation
  // FileReader::readFile(particles, input_file);
  ParticleContainer particleContainer;
  CuboidReader::read(particleContainer.getParticlesReference(),
                     arguments.inputFile);

  VerletIntegrator verlet_integrator(lennardjones, arguments.delta_t);
  outputWriter::VTKWriter writer;

  const std::string outputDirectory = createOutputDirectory("./output/", argc, argv);

  double current_time = 0;  // start time is always 0
  int iteration = 0;
  int writes = 0;

  while (current_time <= arguments.t_end) {
    verlet_integrator.step(particleContainer);
    // if (current_time >= writes * output_time_step_size) {

    if (writes % output_interval == 0) {
      plotParticles(outputDirectory, iteration, writer, particleContainer);
#ifdef DEBUG
      SpdWrapper::get()->debug("Iteration {} finished.", iteration);
#else
      const double completion_percentage = 100 * current_time / arguments.t_end;
      const std::string output_string =
          "\r[" + std::to_string(completion_percentage) +
          " %]: " + std::to_string(iteration) +
          " iterations finished at t=" + std::to_string(current_time);
      std::cout << output_string << std::flush;
#endif
    }
    writes++;
    //}

    iteration++;
    current_time = arguments.delta_t * iteration;  // + start_time
  }

  SpdWrapper::get()->info("Output written. Terminating...");

  return 0;
}
