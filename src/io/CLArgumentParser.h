//
// Created by mcarn on 11/6/24.
//

#ifndef CLARGUMENTPARSER_H
#define CLARGUMENTPARSER_H

#include <algorithm>
#include <memory>

#include "forces/Force.h"

/**
 * @brief struct to hold command line arguments
 */
struct Arguments {
  std::string inputFile;
  double t_end;
  double delta_t;
  double output_time_step_size;
  std::string logLevel;
  std::unique_ptr<Force> force;
};

/**
 * @brief Static class to encapsulate CLI argument parsin
 * @note Does also set the logging level
 */
class CLArgumentParser {
 public:
  /**
   * @brief Fills the passed arguments struct with the parsed arguments from the
   * CLI
   * @note prints usage if failed
   * @param argc Directly from main
   * @param argv Directly from main
   * @param arguments argument struct (should hold default arguments)
   * @return -1 Failure, 0 Success
   */
  static int parse(int argc, char *argv[], Arguments &arguments);

  /**
   * @brief print usage
   * @param additionalNote additional info to customize printout
   * @param programName name of the program (argv[0])
   */
  static void printUsage(const std::string &additionalNote,
                         const std::string &programName);
};

inline std::string toLower(const std::string &s) {
  std::string res = s;
  std::transform(res.begin(), res.end(), res.begin(),
                 [](const unsigned char c) { return std::tolower(c); });
  return res;
}

#endif  // CLARGUMENTPARSER_H
