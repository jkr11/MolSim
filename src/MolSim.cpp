#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"


#include <defs/Particle.h>
#include <iostream>
#include <list>
#include <vector>
#include <filesystem>
#include <getopt.h>

#include "calc/Verlet.h"
#include "forces/gravity.h"
#include "outputWriter/VTKWriter.h"


/**** forward declaration of the calculation functions ****/
void plotParticles(int iteration, outputWriter::VTKWriter& vtkWriter, ParticleContainer& particle_container);
void printUsage(const std::string& additionalNote, const std::string& programName);

constexpr int output_interval = 10000; // seems to be a decent rate, can be changed to arbitrary numbers
constexpr double start_time = 0;
double t_end = 100;
double delta_t = 0.014;

const std::string output = "output/";
std::string output_directory = "./output/";

const std::string red = "\033[31m";
const std::string reset = "\033[0m";

// TODO: what data structure to pick?
std::list<Particle> particles;

Gravity gravity;

int main(const int argc, char* argsv[]) {
  //read optional arguments
  std::string input_file;
  int opt;

  while ((opt = getopt(argc, argsv, "f:t:d:")) != -1) {
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
        printUsage("Invalid argument '" + std::string(optarg) + "' for [-t <double>]", argsv[0]);
      }
      break;
    case 'd':
      try {
        delta_t = std::stod(optarg);
      } catch (...) {
        printUsage("Invalid argument '" + std::string(optarg) + "' for [-d <double>]", argsv[0]);
      }
      break;
    default:
      printUsage("unsupported flag detected", argsv[0]);
    }
  }

  if (input_file.empty()) {
    printUsage("Erroneous program call, input_file empty!", argsv[0]);
  }

  // so we cant really do that incase the input is a directory
  // output_directory += input_file;
  const std::filesystem::path output_directory_path = output_directory;

  if (!std::filesystem::exists(input_file)) {
    printUsage("Input File '" + input_file + "' does not exist", argsv[0]);
  }

  // TODO: fix make this project specific, else just pass in the output path
  if (!is_directory(output_directory_path)) {
    create_directories(output_directory_path);
    std::cout << "Created directory for output: " << output_directory_path << std::endl;
  }

  // further arguments
  // check here if defaults have been used, optional type is c++17 sadly, so double copy??
  std::cout << "t_end: " << t_end << ", delta_t: " << delta_t << std::endl;

  FileReader fileReader;
  fileReader.readFile(particles,  input_file);

  //setup Simulation
  ParticleContainer particle_container(particles);
  VerletIntegrator verlet_integrator(gravity, delta_t);
  outputWriter::VTKWriter writer;

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < t_end) {
    verlet_integrator.step(particle_container);
    iteration++;
    if (iteration % 50 == 0) {
      plotParticles(iteration, writer, particle_container);
    }

    if (iteration % output_interval == 0) {
      // clean up print if DEBUG is enabled
#ifdef DEBUG
      std::cout << "Iteration " << iteration << " finished." << std::endl;
#else
      const double completion_percentage = 100 * current_time / t_end;
      const std::string output_string = "\r[" + std::to_string(completion_percentage) + " %]: "
        + std::to_string(iteration) + " iterations finished";
      std::cout << output_string << std::flush;
#endif
    }


    current_time += delta_t;
  }

  std::cout << std::endl;
  std::cout << "output written. Terminating..." << std::endl;

  return 0;
}

void plotParticles(const int iteration, outputWriter::VTKWriter& vtkWriter, ParticleContainer& particle_container) {
  vtkWriter.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto& p : particle_container.getParticles()) {
    vtkWriter.plotParticle(p);
  }

  vtkWriter.writeFile(output_directory + "/MD_vtk", iteration);
}

void printUsage(const std::string& additionalNote, const std::string& programName) {
  std::cerr << red << "[Error:] " << additionalNote << reset << "\n";
  std::cout << "Usage: " << programName << " [options]\n"
              << "Options:\n"
              << "  -h                Show this help message\n"
              << "  -f <filename>     Specify the input file\n"
              << "  [-t <double>]     Specify the simulation end time (t_end))\n"
              << "  [-d <double>]     Specify the simulation delta time (t_delta)\n"
              << "\nExample:\n"
              << "  " << programName << " -f ./input/eingabe-sonne.txt -t 100 -d 0.14\n";

  exit(EXIT_FAILURE);
}
