//
// Created by maximilian on 19.01.25.
//
#pragma once
#include <fstream>
#include <vector>

#include "utils/SpdWrapper.h"

class CSVWriter {
  std::ofstream file_;

 public:
  explicit CSVWriter(const std::string& file_name) {
    file_.open(file_name);
    if (!file_.is_open()) {
      SpdWrapper::get()->error("Failed to open CSV output file");
      throw std::ios_base::failure("Failed to open file");
    }
  }

  ~CSVWriter() {
    closeFile();
  }

  void writeLine(const double time, const std::vector<std::string>& data) {
    file_ << time << ",";
    for (size_t i = 0; i < data.size(); ++i) {
      file_ << data[i];
      if (i < data.size() - 1) {
        file_ << ",";
      }
    }
    file_ << "\n";
    file_.flush();
  }

  void closeFile() {
    file_.close();
  }
};
