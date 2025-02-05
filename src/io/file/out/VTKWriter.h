/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "defs/Particle.h"
#include "vtk-unstructured.h"

namespace outputWriter {

/**
 * This class implements the functionality to generate vtk output from
 * particles.
 */
class VTKWriter final {
 public:
  VTKWriter();

  ~VTKWriter();

  /**
   * set up internal data structures and prepare to plot a particle.
   */
  void initializeOutput(int num_particles);

  /**
   * plot type, mass, position, velocity and force of a particle.
   *
   * @note: initializeOutput() must have been called before.
   */
  void plotParticle(const Particle &p) const;

  /**
   * writes the final output file.
   *
   * @param filename the base name of the file to be written.
   * @param iteration the number of the current iteration,
   *        which is used to generate an unique filename
   */
  void writeFile(const std::string &filename, int iteration) const;

 private:
  VTKFile_t *vtk_file_{};
};

}  // namespace outputWriter
