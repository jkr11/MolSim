//
// Created by mcarn on 11/6/24.
//
#pragma once
#ifndef OUTPUTHELPER_H
#define OUTPUTHELPER_H
#include <fstream>

#include "VTKWriter.h"
#include "utils/SpdWrapper.h"

/**
 * @brief prints the current state of the system to viewable files
 * @param output_directory specifies the base for the output folders, on which
 * other folders named after date are saved
 * @param iteration the current iteration from the main simulation loop
 * @param vtk_writer writes a .vtu a file
 * @param particle_container contains the entirety of the current particles
 */
inline void plotParticles(const std::string &output_directory,
                          const int iteration,
                          outputWriter::VTKWriter &vtk_writer,
                          ParticleContainer &particle_container) {
  vtk_writer.initializeOutput(static_cast<int>(particle_container.size()));

  particle_container.singleIterator(
      [&vtk_writer](const Particle &p) { vtk_writer.plotParticle(p); });

  vtk_writer.writeFile(output_directory + "/MD_vtk", iteration);
}

/**
 * @brief creates a timestamped directory in ./output/ containing the files at
 * the currents timestep and the specification used to attain the result.
 * @param output_directory the output directory named after the date
 * @param argc argc passed from main
 * @param argv argv passed from main
 * @return name of the output directory path
 */
inline std::string createOutputDirectory(const std::string &output_directory,
                                         int argc, char *argv[]) {
  // source for getting time:
  // https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
  const auto current_time = std::chrono::high_resolution_clock::now();
  const std::time_t now = std::chrono::system_clock::to_time_t(current_time);
  const std::tm local_time = *std::localtime(&now);
  std::ostringstream time_string;
  std::string output_directory_path = "./output";

  // save configuration (input) for future use
  std::ofstream config(output_directory_path + "/configuration.txt");
  for (int i = 0; i < argc; i++) {
    config << argv[i] << " ";
  }
  config << std::endl;
  config.close();

  return output_directory_path;
}

#endif  // OUTPUTHELPER_H
