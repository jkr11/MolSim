#include <getopt.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <list>
#include <vector>

#include "FileReader.h"
#include "calc/VerletIntegrator.h"
#include "defs/Particle.h"
#include "forces/Gravity.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include "utils/CuboidReader.h"

/**** forward declaration of the calculation functions ****/
void plotParticles(int iteration, outputWriter::VTKWriter& vtkWriter,
                   ParticleContainer& particle_container);
void printUsage(const std::string& additionalNote,
                const std::string& programName);
void prepareOutputDirectory(int argsc, char* argv[]);

constexpr int output_interval = 32;
// print every x vtk writes, seems to be decent
constexpr double start_time = 0;
double t_end = 100;
double delta_t = 0.014;
double output_time_step_size = 1;

// parent directory, subdirectory will be appended at runtime
std::string output_directory = "./output/";

const std::string red = "\033[31m";
const std::string reset = "\033[0m";

std::list<Particle> particles;

Gravity gravity;

int main(const int argc, char* argsv[]) {
  // read optional arguments
  std::string input_file;
  int opt;

  while ((opt = getopt(argc, argsv, "hf:t:d:s:")) != -1) {
    try {
      if (optarg == nullptr) {
        throw std::logic_error("missing option after flag");
      }

      switch (opt) {
        case 'h':
          printUsage("Display Help page, no execution", argsv[0]);
        case 'f':
          input_file = optarg;
          break;
        case 't':
          t_end = std::stod(optarg);
          break;
        case 'd':
          delta_t = std::stod(optarg);
          break;
        case 's':
          output_time_step_size = std::stod(optarg);
          break;
        default:
          printUsage("unsupported flag '-" +
                         std::string(1, static_cast<char>(opt)) + "' detected",
                     argsv[0]);
      }
    } catch (const std::invalid_argument&) {
      printUsage("Invalid arg for option -" +
                     std::string(1, static_cast<char>(opt)) + ": '" +
                     std::string(optarg) + "'",
                 argsv[0]);
    } catch (const std::out_of_range&) {
      printUsage("Out-of-range value for option -" +
                     std::string(1, static_cast<char>(opt)) + ": '" +
                     std::string(optarg) + "'",
                 argsv[0]);
    } catch (const std::logic_error&) {
      printUsage(" ^^", argsv[0]);
    }
  }

  if (!std::filesystem::exists(input_file) ||
      std::filesystem::is_directory(input_file)) {
    printUsage("Input File '" + input_file + "' does not exist", argsv[0]);
  }
  if (input_file.empty()) {
    printUsage("input file " + input_file + " is empty!", argsv[0]);
  }

  prepareOutputDirectory(argc, argsv);
  std::cout << "t_end: " << t_end << ", delta_t: " << delta_t
            << ", output_time_step_size: " << output_time_step_size
            << std::endl;

  //FileReader::readFile(particles, input_file);

  // setup Simulation
  ParticleContainer particle_container(particles);
  VerletIntegrator verlet_integrator(gravity, delta_t);
  outputWriter::VTKWriter writer;
  CuboidReader::readCuboidFile(particle_container, input_file);
  double current_time = start_time;

  int iteration = 0;
  int writes = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time <= t_end) {
    verlet_integrator.step(particle_container);
    if (current_time >= writes * output_time_step_size) {
      plotParticles(iteration, writer, particle_container);

      if (writes % output_interval == 0) {
#ifdef DEBUG
        std::cout << "Iteration " << iteration << " finished." << std::endl;
#else
        const double completion_percentage = 100 * current_time / t_end;
        const std::string output_string =
            "\r[" + std::to_string(completion_percentage) +
            " %]: " + std::to_string(iteration) +
            " iterations finished at t=" + std::to_string(current_time);
        std::cout << output_string << std::flush;
#endif
      }

      writes++;
    }

    iteration++;
    current_time = start_time + delta_t * iteration;
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
    std::cout << "particles plotted\n" << std::endl;
  }

  vtkWriter.writeFile(output_directory + "/MD_vtk", iteration);
}

void printUsage(const std::string& additionalNote,
                const std::string& programName) {
  std::cerr << red << "[Error:] " << additionalNote << reset << "\n";
  std::cout << "Usage: " << programName << " [options]\n"
            << "Options:\n"
            << "  -h                Show this help message\n"
            << "  -f <filename>     Specify the input file\n"
            << "  [-t <double>]     Specify the simulation end time (t_end), "
               "default=100\n"
            << "  [-d <double>]     Specify the simulation delta time "
               "(t_delta), default=0.014\n"
            << "  [-s <double>]     Specify how often the output will be "
               "written (step_size), default=1\n"
            << "                    note that this is independent of the time "
               "resolution (t_delta) and dependent on the simulation time"
            << "\nExample:\n"
            << "  " << programName
            << " -f ./input/eingabe-sonne.txt -t 100 -d 0.14\n";

  exit(EXIT_FAILURE);
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
    std::cout << "Output at: " << output_directory_path << std::endl;
  }

  // save configuration (input) for future use
  std::ofstream config(output_directory + "configuration.txt");
  for (int i = 0; i < argsc; i++) {
    config << argv[i] << " ";
  }
  config << std::endl;
  config.close();
}