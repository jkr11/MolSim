#include <complex>
#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "defs/Simulation.h"
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

int main(int argc, char* argv[]) {
  SpdWrapper::get()->info("Application started");
  // ok so this is really ugly but i think we can merge both structs into one
  // we dont need additional file readers anymore so we can just use XMLReader

  Arguments arguments = {
      //.input_file = "",                  // file
      .t_end = 5,         // t_end
      .delta_t = 0.0002,  // delta_t
      //.output_time_step_size = 1,        // output_time_step_size
      .log_level = "info",                    // logLevel
      .force_type = Arguments::LennardJones,  // force
      .domain = {100, 100, 100},
      .cutoff_radius = 3.0,
      .container_type = Arguments::LinkedCells,
  };  // TODO: figure out if the . assignement in structs is valid C++17

  auto [input_file, step_size] = CLArgumentParser::parse(argc, argv);
  // const auto reader = std::make_unique<XmlReader>();
  //  TODO: Should we change this so it doesnt get read here but the reader
  //  instantiates the container and then writes the shapes to the container?
  std::vector<Particle> particles;
  // reader->read(particles, input_file);
  XmlReader::read(particles, input_file, arguments);
  SpdWrapper::get()->info("Particles size {}", particles.size());
  /*
  auto [delta_t, t_end, cutoff_radius, domain, force_type, container_type] =
      reader.get()->pass();
      */

  SpdWrapper::get()->info("t_end: {}, delta_t: {}, cutoff_radius: {}",
                          arguments.t_end, arguments.delta_t,
                          arguments.cutoff_radius);

  // maybe we can make this nicer, this is the best i can come up with right now
  std::unique_ptr<ParticleContainer> container;
  if (arguments.container_type == Arguments::LinkedCells) {
    container = std::make_unique<LinkedCellsContainer>(arguments.domain,
                                                       arguments.cutoff_radius);
    container->addParticles(particles);
    container->imposeInvariant();
  } else if (arguments.container_type == Arguments::DirectSum) {
    container = std::make_unique<DirectSumContainer>();
    container->addParticles(particles);
  }

  // set force type
  std::unique_ptr<Force> force;
  if (arguments.force_type == Arguments::Gravity) {
    force = std::make_unique<Gravity>();
  } else if (arguments.force_type == Arguments::LennardJones) {
    force = std::make_unique<LennardJones>();
  }
  SpdWrapper::get()->info("particles.size: {}", particles.size());
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
