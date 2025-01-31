//
// Created by jkr on 12/3/24.
//

#ifndef XMLWRITER_H
#define XMLWRITER_H
#pragma once
#include "defs/containers/ParticleContainer.h"

/**
 * @brief writes particles to a checkpoint file in xml format
 */
class XmlWriter {
 public:
  /**
   * @brief constructs an XML Writer
   */
  XmlWriter();
  /**
   * @brief Destructor for XmlWrither
   */
  ~XmlWriter();

  /**
   * @brief writes the particles from the contianer ot the file
   * @param particle_container the particles to be written
   * @param filepath the file to be written to
   */
  static void writeFile(ParticleContainer& particle_container,
                        const std::string& filepath);
};
#endif  // XMLWRITER_H
