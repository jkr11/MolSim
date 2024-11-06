//
// Created by mcarn on 11/6/24.
//

#ifndef CLARGUMENTPARSER_H
#define CLARGUMENTPARSER_H

#include <iostream>

struct Arguments {
  std::string inputFile;
  double t_end;
  double delta_t;
  double output_time_step_size;
  std::string logLevel;
};

class CLArgumentParser {
 public:
  static int parse(int argc, char* argv[], Arguments& arguments);
  static void printUsage(const std::string& additionalNote,
                         const std::string& programName);
};

#endif  // CLARGUMENTPARSER_H