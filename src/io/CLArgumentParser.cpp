//
// Created by mcarn on 11/6/24.
//

#include "CLArgumentParser.h"

#include <getopt.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "utils/SpdWrapper.h"

int CLArgumentParser::parse(int argc, char* argv[], Arguments& arguments) {
  struct option long_options[] = {
      {"help", no_argument, nullptr, 'h'},
      {"file", required_argument, nullptr, 'f'},
      {"t_end", required_argument, nullptr, 't'},
      {"delta_t", required_argument, nullptr, 'd'},
      {"step_size", required_argument, nullptr, 's'},
      {"logging", required_argument, nullptr, 'l'},
      {nullptr, 0, nullptr, 0}};

  char opt;
  int option_index = 0;

  while ((opt = getopt_long(argc, argv, "hf:t:d:s:l:", long_options,
                            &option_index)) != -1) {
    try {
      if ((opt == 'f' || opt == 't' || opt == 'd' || opt == 's') &&
          optarg == nullptr) {
        throw std::logic_error("missing option after flag");
      }

      switch (opt) {
        case 'h':
          printUsage("Display Help page, no execution", argv[0]);
          return -1;
        case 'f':
          arguments.inputFile = optarg;
          break;
        case 't':
          arguments.t_end = std::stod(optarg);
          break;
        case 'd':
          arguments.delta_t = std::stod(optarg);
          break;
        case 's':
          arguments.output_time_step_size = std::stod(optarg);
          break;
        case 'l':
          // convert string with level to int / enum TODO
          arguments.loggingLevel = std::atoi(optarg);
          break;
        default:
          printUsage("unsupported flag '-" +
                         std::string(1, static_cast<char>(opt)) + "' detected",
                     argv[0]);
          return -1;
      }
    } catch (const std::invalid_argument&) {
      printUsage("Invalid arg for option -" +
                     std::string(1, static_cast<char>(opt)) + ": '" +
                     std::string(optarg) + "'",
                 argv[0]);
      return -1;
    } catch (const std::out_of_range&) {
      printUsage("Out-of-range value for option -" +
                     std::string(1, static_cast<char>(opt)) + ": '" +
                     std::string(optarg) + "'",
                 argv[0]);

      return -1;
    } catch (const std::logic_error&) {
      printUsage(" ^^", argv[0]);
      return -1;
    }
  }

  // validate input file (move to IO validator?)
  if (!std::filesystem::exists(arguments.inputFile) ||
      std::filesystem::is_directory(arguments.inputFile)) {
    // this seems kinda weird, it has to do with input, but not with parsing =>
    // extract printUsage to somewhere else?
    CLArgumentParser::printUsage(
        "Input File '" + arguments.inputFile + "' does not exist", argv[0]);
    return -1;
  }

  auto iss = std::ifstream(arguments.inputFile);
  if (iss.peek() == std::ifstream::traits_type::eof()) {
    // this seems kinda weird, it has to do with input, but not with parsing =>
    // extract printUsage to somewhere else?
    CLArgumentParser::printUsage(
        "input file " + arguments.inputFile + " is empty!", argv[0]);
    return -1;
  }

  return 0;
}

void CLArgumentParser::printUsage(const std::string& additionalNote,
                                  const std::string& programName) {
  // std::cerr << red << "[Error:] " << additionalNote << reset << "\n";
  SpdWrapper::get()->error(additionalNote);
  SpdWrapper::get()->info(
      "Usage: {} [options]\n"
      "Options:\n"
      "  -h                Show this help message\n"
      "  -f <filename>     Specify the input file\n"
      "  [-t <double>]     Specify the simulation end time (t_end), "
      "default=100\n"
      "  [-d <double>]     Specify the simulation delta time (t_delta), "
      "default=0.014\n"
      "  [-s <double>]     Specify how often the output will be written "
      "(step_size), default=1\n"
      "                    note that this is independent of the time "
      "resolution (t_delta) and dependent on the simulation time\n"
      "Example:\n"
      "  {} -f ./input/eingabe-sonne.txt -t 100 -d 0.14\n",
      programName, programName);
}
