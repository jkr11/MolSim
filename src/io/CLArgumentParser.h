//
// Created by mcarn on 11/6/24.
//

#ifndef CLARGUMENTPARSER_H
#define CLARGUMENTPARSER_H

#include <algorithm>
#include <filesystem>
#include <tuple>

#include "defs/Simulation.h"

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
   * @returns path to the input file
   */
  static std::tuple<std::filesystem::path, double> parse(int argc,
                                                         char *argv[]);

  /**
   * @brief print usage
   * @param additionalNote additional info to customize printout
   * @param programName name of the program (argv[0])
   */
  static void printUsage(const std::string &additionalNote,
                         const std::string &programName);

  /**
   * @brief Validates the existence and validity of an input file.
   * @param file_path The path to the input file to be validated.
   * @throws std::invalid_argument if the file does not exist, is a directory,
   * or is empty.
   * @note verification that this is actually an xml file is done in XmlReader.
   */
  static void validateInputFile(const std::filesystem::path &file_path);
};

/**
 * @brief return the lowercase version of any (valid) string inputted
 * @param s reference to the input string
 * @return lowercase version of s
 */
inline std::string toLower(const std::string &s) {
  std::string res = s;
  std::transform(res.begin(), res.end(), res.begin(),
                 [](const unsigned char c) { return std::tolower(c); });
  return res;
}

#endif  // CLARGUMENTPARSER_H
