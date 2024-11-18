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
                                 {nullptr, 0, nullptr, 0}};

  int opt;
  int option_index = 0;

  while ((opt = getopt_long(argc, argv, "hf:t:d:s:l:F:R:", long_options,
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
          arguments.logLevel = optarg;
          break;
        case 'F': {
          if (const std::string f = toLower(optarg); f == "lennardjones") {
            arguments.force = std::make_unique<LennardJones>();
          } else if (f == "gravity") {
            arguments.force = std::make_unique<Gravity>();
          } else {
            SpdWrapper::get()->error("Unknown Force Type: {}", f);
            exit(EXIT_FAILURE);
          }
          break;
        }
        case 'R': {
          SpdWrapper::get()->info(
              "This is deprecated as we only need XMLReader now, note that "
              "only xml files are allowed.");
        }
        default:
          printUsage("unsupported flag '-" +
                         std::string(1, static_cast<char>(opt)) + "' detected",
                     argv[0]);
          exit(EXIT_FAILURE);
      }
    } catch (const std::invalid_argument &) {
      printUsage("Invalid arg for option -" +
                     std::string(1, static_cast<char>(opt)) + ": '" +
                     std::string(optarg) + "'",
                 argv[0]);
      return -1;
    } catch (const std::out_of_range &) {
      printUsage("Out-of-range value for option -" +
                     std::string(1, static_cast<char>(opt)) + ": '" +
                     std::string(optarg) + "'",
                 argv[0]);

      return -1;
    } catch (const std::logic_error &e) {
      printUsage(e.what(), argv[0]);
      return -1;
    }
  }

  // validate input file
  if (!std::filesystem::exists(arguments.inputFile) ||
      std::filesystem::is_directory(arguments.inputFile)) {
    printUsage("Input File '" + arguments.inputFile + "' does not exist",
               argv[0]);
    return -1;
  }

  if (auto iss = std::ifstream(arguments.inputFile);
      iss.peek() == std::ifstream::traits_type::eof()) {
    printUsage("input file " + arguments.inputFile + " is empty!", argv[0]);
    return -1;
  }

  // change loglevel
  if (SpdWrapper::setLogLevel(arguments.logLevel) != 0) {
    printUsage("Log level invalid.", argv[0]);
    return -1;
  }

  return 0;
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
