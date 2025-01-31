//
// Created by mcarn on 11/6/24.
//

#include "CLArgumentParser.h"

#include <getopt.h>

//#include <filesystem>
#include <fstream>
#include <iostream>
#include <tuple>

#include "spdlog/fmt/bundled/chrono.h"
#include "utils/SpdWrapper.h"

std::tuple<std::string, double, bool> CLArgumentParser::parse(
    const int argc, char *argv[]) {
  SpdWrapper::get()->info("Parsing arguments");
  const option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                 {"file", required_argument, nullptr, 'f'},
                                 {"step_size", required_argument, nullptr, 's'},
                                 {"loglevel", required_argument, nullptr, 'l'},
                                 {"checkpoint", no_argument, nullptr, 'c'},
                                 {nullptr, 0, nullptr, 0}};

  int opt;
  int option_index = 0;

  std::string input_file{};
  double step_size = 0.5;
  bool write_checkpoint = false;

  while ((opt = getopt_long(argc, argv, "hf:s:l:c", long_options,
                            &option_index)) != -1) {
    SpdWrapper::get()->info("Parsing options");
    try {
      if ((opt == 'f' || opt == 't' || opt == 'd' || opt == 's') &&
          optarg == nullptr) {
        throw std::invalid_argument("invalid argument for option -" +
                                    std::string(1, static_cast<char>(opt)));
      }
      SpdWrapper::get()->info("Parsing -{} opt, {} optind",
                              static_cast<char>(opt), option_index);
      switch (opt) {
        case 'h':
          printUsage("Display Help page, no execution", argv[0]);
          exit(EXIT_FAILURE);
        case 'f':
          input_file = optarg;
          break;
        case 's':
          step_size = std::stod(optarg);
          break;
        case 'l':
          if (SpdWrapper::setLogLevel(optarg) != 0) {
            throw std::invalid_argument(
                "invalid argument for option --loglevel" + std::string(optarg));
          }
          break;
        case 'c':
          write_checkpoint = true;
          break;
        default:
          throw std::invalid_argument("Unsupported option: -" +
                                      std::string(1, static_cast<char>(opt)));
      }
    } catch (const std::exception &e) {
      printUsage(e.what(), argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  SpdWrapper::get()->info("Parsing");

  try {
    validateInputFile(input_file);
  } catch (const std::exception &e) {
    printUsage(e.what(), argv[0]);
    exit(EXIT_FAILURE);
  }

  SpdWrapper::get()->info("Resetting getopt vars ... ");

  optind = 0;        // Reset position
  optarg = nullptr;  // Reset argument pointer
  optopt = 0;        // Reset last option character

  return {input_file, step_size, write_checkpoint};
}

void CLArgumentParser::validateInputFile(
    const std::string &file_path) {

  SpdWrapper::get()->info("Check if empty");
  if (std::ifstream iss(file_path);
      iss.peek() == std::ifstream::traits_type::eof()) {
    throw std::invalid_argument("Input file '" + std::string(file_path) +
                                "' is empty!");
  }
}

void CLArgumentParser::printUsage(const std::string &additional_note,
                                  const std::string &program_name) {
  SpdWrapper::get()->set_level(spdlog::level::err);
  SpdWrapper::get()->error(additional_note);
  SpdWrapper::get()->error(
      "Usage: {} [options]\n"
      "Options:\n"
      "  --help | -h                     Show this help message\n"
      "  --file | -f <filename>          Specify the input file\n"
      "  [--step_size | -s <double>]     Specify how often the output will be "
      "written wrt. time"
      "(step_size), default=1\n"
      "                                  Note that this is independent of the "
      "time "
      "resolution (t_delta) and dependent on the simulation time\n"
      "  [--loglevel | -l <level>]       Specify the log level, default=info, "
      "valid=[off, error, warn, info, debug, trace]\n"
      "  [--checkpoint | -c ]            Specify if the end state will be "
      "written to a checkpoint file "
      "Example:\n"
      "  {} -f <relative input location>.xml --loglevel info --step_size "
      "0.01\n",
      program_name, program_name);
}
