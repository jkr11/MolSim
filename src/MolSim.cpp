#include <filesystem>

#include "io/CLArgumentParser.h"
#include "io/file/in/CuboidReader.h"
#include "io/file/in/xml/XmlReader.h"
#include "io/file/out/OutputHelper.h"
#include "spdlog/fmt/bundled/chrono.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

int main(const int argc, char* argv[]) {
  SpdWrapper::get()->info("Application started");

  Arguments arguments = {
      .t_end = 5,
      .delta_t = 0.0002,
      .force_type = Arguments::LennardJones,
      .container_data =
          LinkedCellsConfig{.domain = {100, 100, 100},
                            .cutoff_radius = 3.0,
                            .boundary_config =
                                {
                                    .x_high = LinkedCellsConfig::Outflow,
                                    .x_low = LinkedCellsConfig::Outflow,
                                    .y_high = LinkedCellsConfig::Outflow,
                                    .y_low = LinkedCellsConfig::Outflow,
                                    .z_high = LinkedCellsConfig::Outflow,
                                    .z_low = LinkedCellsConfig::Outflow,
                                }},
  };
  auto [input_file, step_size] = CLArgumentParser::parse(argc, argv);
  XmlReader::read(input_file, arguments);

  printConfiguration(arguments);

  SpdWrapper::get()->info("Step size: {}", step_size);

  const std::string outputDirectory =
      createOutputDirectory("./output/", argc, argv);

  Simulation simulation(std::move(arguments), outputDirectory, step_size);
  simulation.initParams();
  // simulation.initParticles();
  // simulation.run();

  return 0;
}
