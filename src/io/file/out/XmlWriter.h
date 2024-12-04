//
// Created by jkr on 12/3/24.
//

#ifndef XMLWRITER_H
#define XMLWRITER_H
#pragma once
#include "defs/containers/ParticleContainer.h"

class XmlWriter {
 public:
  XmlWriter();
  ~XmlWriter();

  static void writeFile(ParticleContainer& particle_container,
                 const std::string& filepath);
};
#endif  // XMLWRITER_H
