//
// Created by maximilian on 19.01.25.
//
#pragma once
#include <fstream>
#include <vector>

#include "utils/SpdWrapper.h"

/**
 * @brief Writes to any csv file
 */
class CSVWriter {
  std::ofstream file_;

 public:
  /**
   * Checks if the given file is open
   * @param file_name the target file
   */
  explicit CSVWriter(const std::string& file_name) {
    SpdWrapper::get()->info("Opening {}...", file_name);
    file_.open(file_name);
    if (!file_.is_open()) {
      SpdWrapper::get()->error("Failed to open CSV output file");
      throw std::ios_base::failure("Failed to open file");
    }
  }

  /**
   * @brief Destructor
   * @note Closes open files
   */
  ~CSVWriter() { closeFile(); }

  /**
   * @brief writes a vector of data into a csv
   * @param time timestamp for data in the 0th column
   * @param data the data that is written
   */
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

  /**
   * @brief closes the specified file
   */
  void closeFile() { file_.close(); }
};
