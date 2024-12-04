//
// Created by jkr on 12/3/24.
//
#include "XmlWriter.h"

#include "defs/Thermostat.h"
#include "io/file/out/checkpoint-schema.hxx"

XmlWriter::XmlWriter() = default;

XmlWriter::~XmlWriter() = default;

template <typename VecType, typename XmlType>
inline XmlType wrapVec(const VecType& vec, const std::string& name) {
  XmlType xmlVec;
  xmlVec.x() = vec[0];
  xmlVec.y() = vec[1];
  xmlVec.z() = vec[2];
  return xmlVec;
}

void XmlWriter::writeFile(const std::vector<Particle>& particles,
                          const Arguments& args, std::size_t iteration) {

}