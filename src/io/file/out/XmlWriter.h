//
// Created by jkr on 12/3/24.
//

#ifndef XMLWRITER_H
#define XMLWRITER_H
#pragma once
#include "checkpoint-schema.hxx"
#include "io/file/in/xml/input.hxx"
#include "defs/Particle.h"
#include "defs/Simulation.h"
class XmlWriter {
 public:
  XmlWriter();
  ~XmlWriter();

  void writeFile(const std::vector<Particle>& particles, const Arguments& args,
                 std::size_t iteration);
};
#endif  // XMLWRITER_H
