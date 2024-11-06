#include <getopt.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <list>
#include <vector>

#include "FileReader.h"
#include "calc/VerletIntegrator.h"
#include "defs/Particle.h"
#include "forces/Force.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "io/CLArgumentParser.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include "utils/CuboidReader.h"
#include "utils/SpdWrapper.h"

/**** forward declaration of the calculation functions ****/
void plotParticles(int iteration, outputWriter::VTKWriter& vtkWriter,
                   ParticleContainer& particle_container);
void prepareOutputDirectory(int argsc, char* argv[]);

constexpr int output_interval = 32;

// parent directory, subdirectory will be appended at runtime
std::string output_directory = "./output/";

std::list<Particle> particles;

Gravity gravity;
LennardJones lennardjones;

int main(const int argc, char* argv[]) {
  SpdWrapper::get()->info("Application started");

  struct Arguments arguments = {
      "",     // file
      100,    // t_end
      0.014,  // delta_t
      1,      // output_time_step_size
      0,      // logging_level -> TODO
  };

  if (CLArgumentParser::parse(argc, argv, arguments) != 0) {
    exit(EXIT_FAILURE);
  }

  /*SpdWrapper::get()->info(arguments.inputFile);
  SpdWrapper::get()->info(arguments.t_end);
  SpdWrapper::get()->info(arguments.delta_t);
  SpdWrapper::get()->info(arguments.output_time_step_size);
  SpdWrapper::get()->info(arguments.loggingLevel);*/

  // TODO: Change spdLogLevel / parse properly

  prepareOutputDirectory(argc, argv);
  SpdWrapper::get()->info("t_end: {}, delta_t: {}, output_time_step_size: {}",
                          arguments.t_end, arguments.delta_t,
                          arguments.output_time_step_size);

  // FileReader::readFile(particles, input_file);

  // setup Simulation
  ParticleContainer particle_container(particles);
  VerletIntegrator verlet_integrator(lennardjones, arguments.delta_t);
  outputWriter::VTKWriter writer;
  CuboidReader::readCuboidFile(particle_container, arguments.inputFile);
  // FileReader::readFile(particles, input_file);

  // current_time = start_time, but start_time couldn't be changed by the user
  // anyways
  double current_time = 0;
  int iteration = 0;
  int writes = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time <= arguments.t_end) {
    verlet_integrator.step(particle_container);
    // if (current_time >= writes * output_time_step_size) {

    if (writes % output_interval == 0) {
      plotParticles(iteration, writer, particle_container);
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
    // starttime remove here also
    current_time = arguments.delta_t * iteration;
  }

  SpdWrapper::get()->info("Output written. Terminating...");

  return 0;
}

void plotParticles(const int iteration, outputWriter::VTKWriter& vtkWriter,
                   ParticleContainer& particle_container) {
  vtkWriter.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto& p : particle_container.getParticles()) {
    vtkWriter.plotParticle(p);
    // SpdWrapper::get()->info("Plotted");
  }

  vtkWriter.writeFile(output_directory + "/MD_vtk", iteration);
}

void prepareOutputDirectory(const int argsc, char* argv[]) {
  // source for getting time:
  // https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
  const auto currentTime = std::chrono::high_resolution_clock::now();
  const std::time_t now = std::chrono::system_clock::to_time_t(currentTime);
  const std::tm localTime = *std::localtime(&now);
  std::ostringstream timeString;

  timeString << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S/");
  std::string output_sub_directory = timeString.str();

  // output at outputdirectory/currentTime
  const std::filesystem::path output_directory_path =
      output_directory + output_sub_directory;
  output_directory = std::string(output_directory_path);

  if (!is_directory(output_directory_path)) {
    create_directories(output_directory_path);
    SpdWrapper::get()->info("Output at {}", output_directory_path.string());
  }

  // save configuration (input) for future use
  std::ofstream config(output_directory + "configuration.txt");
  for (int i = 0; i < argsc; i++) {
    config << argv[i] << " ";
  }
  config << std::endl;
  config.close();
}