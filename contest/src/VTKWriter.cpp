/*
 * VTKWriter.cpp
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#include "VTKWriter.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>

void VTKWriter::initializeOutput(int numParticles) {
  vtkFile = new VTKFile_t("UnstructuredGrid");

  // per point, we add type, position, velocity and force
  PointData pointData;
  DataArray_t velocity(type::Float32, "velocity", 3);
  DataArray_t forces(type::Float32, "force", 3);
  DataArray_t type(type::Int32, "type", 1);
  pointData.DataArray().push_back(velocity);
  pointData.DataArray().push_back(forces);
  pointData.DataArray().push_back(type);

  CellData cellData;  // we don't have cell data => leave it empty

  // 3 coordinates
  Points points;
  DataArray_t pointCoordinates(type::Float32, "points", 3);
  points.DataArray().push_back(pointCoordinates);

  Cells cells;  // we don't have cells, => leave it empty
  // for some reasons, we have to add a dummy entry for paraview
  DataArray_t cells_data(type::Float32, "types", 0);
  cells.DataArray().push_back(cells_data);

  PieceUnstructuredGrid_t piece(pointData, cellData, points, cells,
                                numParticles, 0);
  UnstructuredGrid_t unstructuredGrid(piece);
  vtkFile->UnstructuredGrid(unstructuredGrid);
}

void VTKWriter::writeFile(const std::string &filename, int iteration) const {
  std::stringstream strstr;
  strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration
         << ".vtu";

  std::ofstream file(strstr.str().c_str());
  VTKFile(file, *vtkFile);
  delete vtkFile;
}

void VTKWriter::plotParticle(const dvec3 pos, const dvec3 vel, const dvec3 force, int type) const {


  PointData::DataArray_sequence &pointDataSequence =
      vtkFile->UnstructuredGrid()->Piece().PointData().DataArray();
  PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

  dataIterator->push_back(vel[0]);
  dataIterator->push_back(vel[1]);
  dataIterator->push_back(vel[2]);
  // cout << "Appended velocity data in: " << dataIterator->Name();

  ++dataIterator;
  dataIterator->push_back(force[0]);
  dataIterator->push_back(force[1]);
  dataIterator->push_back(force[2]);
  // cout << "Appended force data in: " << dataIterator->Name();

  ++dataIterator;
  dataIterator->push_back(type);

  Points::DataArray_sequence &pointsSequence =
      vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
  const Points::DataArray_iterator pointsIterator = pointsSequence.begin();

  pointsIterator->push_back(pos[0]);
  pointsIterator->push_back(pos[1]);
  pointsIterator->push_back(pos[2]);
}
