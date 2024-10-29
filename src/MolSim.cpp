#include <defs/Particle.h>
#include <getopt.h>

#include <filesystem>
#include <iostream>
#include <list>
#include <vector>
#include <chrono>

#include "FileReader.h"
#include "calc/Verlet.h"
#include "forces/gravity.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

/**** forward declaration of the calculation functions ****/
void plotParticles(int iteration, outputWriter::VTKWriter& vtkWriter,
                   ParticleContainer& particle_container);
void printUsage(const std::string& additionalNote,
                const std::string& programName);
void prepareOutputDirectory(int argsc, char* argv[]);

int output_interval = 10; // in relation to vtk writes, seems to be decent
constexpr double start_time = 0;
double t_end = 100;
double delta_t = 0.014;
double output_time_step_size = 0.01;

// parent directory, subdirectory will be appended at runtime
std::string output_directory = "./output/";

const std::string red = "\033[31m";
const std::string reset = "\033[0m";

// TODO: what data structure to pick?
std::list<Particle> particles;

Gravity gravity;

int main(const int argc, char* argsv[]) {
  // read optional arguments
  std::string input_file;
  int opt;

  while ((opt = getopt(argc, argsv, "f:t:d:s:")) != -1) {
    switch (opt) {
      case 'h':
        printUsage("Display Help page, no execution", argsv[0]);
      case 'f':
        input_file = optarg;
        break;
      case 't':
        try {
          t_end = std::stod(optarg);
        } catch (...) {
          printUsage("Invalid argument '" + std::string(optarg) +
                     "' for [-t <double>]",
                     argsv[0]);
        }
        break;
      case 'd':
        try {
          delta_t = std::stod(optarg);
        } catch (...) {
          printUsage(
              "Invalid argument '" + std::string(optarg) +
              "' for [-d <double>]", argsv[0]);
        }
        break;
      case 's':
        try {
          output_time_step_size = std::stod(optarg);
        } catch (...) {
          printUsage("Invalid argument '" + std::string(optarg) +
                     "' for [-s <double>]",
                     argsv[0]);
        }
        break;
      default:
        printUsage("unsupported flag '-" + std::string(optarg) + "' detected",
                   argsv[0]);
    }
  }

  if (!std::filesystem::exists(input_file)) {
    printUsage("Input File '" + input_file + "' does not exist", argsv[0]);
  }
  if (input_file.empty()) {
    printUsage("input file " + input_file + " is empty!", argsv[0]);
  }

  prepareOutputDirectory(argc, argsv);
  std::cout << "t_end: " << t_end << ", delta_t: " << delta_t <<
      ", output_time_step_size: " << output_time_step_size << std::endl;

  FileReader fileReader;
  fileReader.readFile(particles, input_file);

  // setup Simulation
  ParticleContainer particle_container(particles);
  VerletIntegrator verlet_integrator(gravity, delta_t);
  outputWriter::VTKWriter writer;

  double current_time = start_time;
  double next_output_time = start_time;

  int iteration = 0;
  int writes = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < t_end) {
    verlet_integrator.step(particle_container);
    iteration++;
    if (current_time >= next_output_time) {
      plotParticles(iteration, writer, particle_container);
      next_output_time = current_time + output_time_step_size;

      if (writes % output_interval == 0) {
        // clean up print if DEBUG is enabled
#ifdef DEBUG
        std::cout << "Iteration " << iteration << " finished." << std::endl;
#else
        const double completion_percentage = 100 * current_time / t_end;
        const std::string output_string =
            "\r[" + std::to_string(completion_percentage) +
            " %]: " + std::to_string(iteration) + " iterations finished";
        std::cout << output_string << std::flush;
#endif
      }

      writes++;
    }

    current_time += delta_t;
  }

  std::cout << std::endl;
  std::cout << "output written. Terminating..." << std::endl;

  return 0;
}

void plotParticles(const int iteration, outputWriter::VTKWriter& vtkWriter,
                   ParticleContainer& particle_container) {
  vtkWriter.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto& p : particle_container.getParticles()) {
    vtkWriter.plotParticle(p);
  }

  vtkWriter.writeFile(output_directory + "/MD_vtk", iteration);
}


void printUsage(const std::string& additionalNote,
                const std::string& programName) {
  std::cerr << red << "[Error:] " << additionalNote << reset << "\n";
  std::cout
      << "Usage: " << programName << " [options]\n"
      << "Options:\n"
      << "  -h                Show this help message\n"
      << "  -f <filename>     Specify the input file\n"
      << "  [-t <double>]     Specify the simulation end time (t_end), default=100\n"
      << "  [-d <double>]     Specify the simulation delta time (t_delta), default=0.014\n"
      << "  [-s <double>]     Specify how often the output will be written (step_size), default=50\n"
      << "                    note that this is independent of the time resolution (t_delta) and dependent of the simulation time"
      << "\nExample:\n"
      << "  " << programName
      << " -f ./input/eingabe-sonne.txt -t 100 -d 0.14\n";

  exit(EXIT_FAILURE);
}


void prepareOutputDirectory(const int argsc, char* argv[]) {
  // source for getting time: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
  const auto currentTime = std::chrono::high_resolution_clock::now();
  const std::time_t now =
      std::chrono::system_clock::to_time_t(currentTime);
  const std::tm localTime = *std::localtime(&now);
  std::ostringstream timeString;

  timeString << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S/");
  std::string output_sub_directory = timeString.str();

  // output at outputdirectory/currentTime
  // if init due to code style
  const std::filesystem::path output_directory_path =
      output_directory + output_sub_directory;
  output_directory = std::string(output_directory_path);

  if (!is_directory(
      output_directory_path)) {
    create_directories(output_directory_path);
    std::cout << "Output at: " << output_directory_path
        << std::endl;
  }

  // save configuration (input) for future use
  std::ofstream config(output_directory + "configuration.txt");
  for (int i = 0; i < argsc; i++) {
    config << argv[i] << " ";
  }
  config << std::endl;
  config.close();
}