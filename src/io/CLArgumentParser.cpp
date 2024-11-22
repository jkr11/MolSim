//
// Created by mcarn on 11/6/24.
//

#include "CLArgumentParser.h"

#include <getopt.h>

#include <filesystem>
#include <fstream>

#include "defs/Simulation.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "spdlog/fmt/bundled/chrono.h"
#include "utils/SpdWrapper.h"

int CLArgumentParser::parse(int argc, char *argv[], Arguments &arguments) {
  const option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                 {"file", required_argument, nullptr, 'f'},
                                 {"t_end", required_argument, nullptr, 't'},
                                 {"delta_t", required_argument, nullptr, 'd'},
                                 {"step_size", required_argument, nullptr, 's'},
                                 {"loglevel", required_argument, nullptr, 'l'},
                                 {"force", required_argument, nullptr, 'F'},
                                 {"reader", required_argument, nullptr, 'R'},
                                 {"container", required_argument, nullptr, 'C'},
                                 {nullptr, 0, nullptr, 0}};

  int opt;
  int option_index = 0;

  while ((opt = getopt_long(argc, argv, "hf:t:d:s:l:F:R:", long_options,
                            &option_index)) != -1) {
    try {
      if ((opt == 'f' || opt == 't' || opt == 'd' || opt == 's') &&
          optarg == nullptr) {
        throw std::invalid_argument("invalid argument for option -" +
                                    std::string(1, static_cast<char>(opt)));
      }

      switch (opt) {
        case 'h':
          printUsage("Display Help page, no execution", argv[0]);
          return -1;
        case 'f':
          arguments.input_file = optarg;
          break;
        case 't':
          arguments.t_end = parseDouble(optarg, "t_end");
          break;
        case 'd':
          arguments.delta_t = parseDouble(optarg, "delta_t");
          break;
        case 's':
          arguments.output_time_step_size = parseDouble(optarg, "step_size");
          break;
        case 'l':
          if (SpdWrapper::setLogLevel(optarg) != 0) {
            throw std::invalid_argument(
                "invalid argument for option --loglevel" + std::string(optarg));
          }
          arguments.log_level = optarg;
          break;
        case 'F':
          arguments.force_type = parseForceType(optarg);
          break;
        case 'R':
          SpdWrapper::get()->info(
              "This is deprecated and ignored as we only need XMLReader now, "
              "note that only xml files are allowed.");
          break;
        case 'C':
          arguments.container_type = parseContainerType(optarg);
          break;
        default:
          throw std::invalid_argument("Unsupported option: -" +
                                      std::string(1, static_cast<char>(opt)));
      }
    } catch (const std::exception &e) {
      printUsage(e.what(), argv[0]);
      return -1;
    }
  }

  try {
    validateInputFile(arguments.input_file);
  } catch (const std::exception &e) {
    printUsage(e.what(), argv[0]);
    return -1;
  }

  return 0;
}

double CLArgumentParser::parseDouble(const char *arg,
                                     const std::string &option_name) {
  try {
    return std::stod(arg);
  } catch (const std::exception &e) {
    throw std::invalid_argument("Invalid numeric value for option '" +
                                option_name + "': " + std::string(arg));
  }
}

Arguments::ForceType CLArgumentParser::parseForceType(const std::string &arg) {
  if (const std::string f = toLower(arg); f == "lennardjones") {
    return Arguments::LennardJones;
  } else {
    if (f == "gravity") {
      return Arguments::Gravity;
    }
    throw std::invalid_argument("Unknown Force Type: '" + arg);
  }
}

Arguments::ContainerType CLArgumentParser::parseContainerType(
    const std::string &arg) {
  if (const std::string c = toLower(arg); c == "linkedcells") {
    return Arguments::LinkedCells;
  } else {
    if (c == "directsum") {
      return Arguments::DirectSum;
    }
    throw std::invalid_argument("Unknown container type: " + arg);
  }
}

void CLArgumentParser::validateInputFile(const std::string &file_path) {
  if (!std::filesystem::exists(file_path) ||
      std::filesystem::is_directory(file_path)) {
    printUsage("File does not exist", file_path);
    throw std::invalid_argument("Input file '" + file_path +
                                "' does not exist or is a directory");
  }

  if (std::ifstream iss(file_path);
      iss.peek() == std::ifstream::traits_type::eof()) {
    throw std::invalid_argument("Input file '" + file_path + "' is empty!");
  }
}

void CLArgumentParser::printUsage(const std::string &additionalNote,
                                  const std::string &programName) {
  // std::cerr << red << "[Error:] " << additionalNote << reset << "\n";

  SpdWrapper::get()->set_level(spdlog::level::err);
  SpdWrapper::get()->error(additionalNote);
  SpdWrapper::get()->error(
      "Usage: {} [options]\n"
      "Options:\n"
      "  --help | -h                     Show this help message\n"
      "  --file | -f <filename>          Specify the input file\n"
      "  [--t_end | -t <double>]         Specify the simulation end time "
      "(t_end), "
      "default=100\n"
      "  [--delta_t | -d <double>]       Specify the simulation delta time "
      "(t_delta), "
      "default=0.014\n"
      "  [--step_size | -s <double>]     Specify how often the output will be "
      "written wrt. time"
      "(step_size), default=1\n"
      "                                  Note that this is independent of the "
      "time "
      "resolution (t_delta) and dependent on the simulation time\n"
      "  [--loglevel | -l <level>]       Specify the log level, default=info, "
      "valid=[off, error, warn, info, debug, trace]\n"
      "  [--force | -F <forceType>]      Specify what force to use, "
      "default=lennardjones, "
      "forceType=[lennardjones,gravity]\n"
      "  [--reader | -R <readerType>]    Specify reader type, "
      "default=cuboidreader, "
      "readerType=[cuboidreader,defaultreader]"
      "Example:\n"
      "  {} -f ../input/test.cuboid -t 5 -d 0.0002 -s 0.01\n",
      programName, programName);
}
