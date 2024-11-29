//
// Created by mcarn on 11/6/24.
//
#pragma once
#ifndef OUTPUTHELPER_H
#define OUTPUTHELPER_H
#include <filesystem>
#include <fstream>

#include "VTKWriter.h"
#include "utils/SpdWrapper.h"

/**
 * @brief prints the current state of the system to viewable files
 * @param outputDirectory specifies the base for the output folders, on which
 * other folders named after date are saved
 * @param iteration the current iteration from the main simulation loop
 * @param vtkWriter writes a .vtu a file
 * @param particle_container contains the entirety of the current particles
 */
inline void plotParticles(const std::string &outputDirectory,
                          const int iteration,
                          outputWriter::VTKWriter &vtkWriter,
                          ParticleContainer &particle_container) {
  vtkWriter.initializeOutput(static_cast<int>(particle_container.size()));

  particle_container.singleIterator(
      [&vtkWriter](const Particle &p) { vtkWriter.plotParticle(p); });

  vtkWriter.writeFile(outputDirectory + "/MD_vtk", iteration);
}

/**
 * @brief creates a timestamped directory in ./output/ containing the files at
 * the currents timestep and the specification used to attain the result.
 * @param outputDirectory the output directory named after the date
 * @param argc argc passed from main
 * @param argv argv passed from main
 * @return name of the output directory path
 */
inline std::string createOutputDirectory(const std::string &outputDirectory,
                                         int argc, char *argv[]) {
  // source for getting time:
  // https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
  const auto currentTime = std::chrono::high_resolution_clock::now();
  const std::time_t now = std::chrono::system_clock::to_time_t(currentTime);
  const std::tm localTime = *std::localtime(&now);
  std::ostringstream timeString;

  // output into 'outputDirectory/current_time/'
  timeString << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S/");
  const std::filesystem::path output_directory_path =
      outputDirectory + timeString.str();

  if (!is_directory(output_directory_path)) {
    create_directories(output_directory_path);
    SpdWrapper::get()->info("Output at {}", output_directory_path.string());
  }

  // save configuration (input) for future use
  std::ofstream config(output_directory_path.string() + "/configuration.txt");
  for (int i = 0; i < argc; i++) {
    config << argv[i] << " ";
  }
  config << std::endl;
  config.close();

  return output_directory_path.string();
}

#endif  // OUTPUTHELPER_H
