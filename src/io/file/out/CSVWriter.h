//
// Created by maximilian on 19.01.25.
//
#pragma once
#include <fstream>
#include <vector>

#include "utils/SpdWrapper.h"

class CSVWriter {
  std::ofstream file;

 public:
  explicit CSVWriter(const std::string& fileName) {
    file.open(fileName);
    if (!file.is_open()) {
      SpdWrapper::get()->error("Failed to open CSV file");
      throw std::ios_base::failure("Failed to open file");
    }
  }

  ~CSVWriter() {
    closeFile();
  }

  void writeLine(const double iteration, const std::vector<std::string>& data) {
    file << iteration << ",";
    for (size_t i = 0; i < data.size(); ++i) {
      file << data[i];
      if (i < data.size() - 1) {
        file << ",";
      }
    }
    file << "\n";
  }

  void closeFile() {
    file.close();
  }
};
