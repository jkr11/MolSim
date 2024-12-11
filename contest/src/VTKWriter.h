//
// Created by mcarn on 12/11/24.
//

#pragma once

#include <iostream>

#include "../../src/defs/types.h"
#include "vtk-unstructured.h"
#include "ParticleContainer.h"

class VTKWriter {
private:
  VTKFile_t *vtkFile{};

public:
  VTKWriter() = default;
  ~VTKWriter() = default;

  void writeFile(const std::string &filename, int iteration) const;
  void initializeOutput(int numParticles);
  void plotParticle(const dvec3 pos, const dvec3 vel, const dvec3 force, int type) const;
};
