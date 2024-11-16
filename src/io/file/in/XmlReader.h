//
// Created by jkr on 11/16/24.
//

#ifndef XMLREADER_H
#define XMLREADER_H
#pragma once
#include "FileReader.h"

class XmlReader final : public FileReader {
 public:
  XmlReader() = default;
  void read(std::vector<Particle> &particles,
            const std::string &filepath) override;
};
#endif  // XMLREADER_H
