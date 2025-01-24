#include <filesystem>
#include <iostream>

#include "calc/VerletIntegrator.h"
#include "defs/Thermostat.h"
#include "defs/containers/DirectSumContainer.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/Gravity.h"
#include "forces/HarmonicForce.h"
#include "forces/IndexForce.h"
#include "forces/LennardJones.h"
#include "forces/SingularGravity.h"
#include "io/CLArgumentParser.h"
#include "io/file/in/xml/XmlReader.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "io/file/out/XmlWriter.h"
#include "spdlog/stopwatch.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"
#include "utils/Statistics.h"

int main(const int argc, char* argv[]) {
#ifndef BENCHMARK
  SpdWrapper::get()->info("Application started");
#endif
  Arguments arguments = {
      .t_end = 5,
      .delta_t = 0.0002,
      .thermostat_config =
          {
              .T_init = 0.5,
              .T_target = 0.5,
              .deltaT = 0.5,
              .n_thermostat = 1000,
              .use_relative = false,
          },
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
      .statistics_config =
          StatisticsConfig{
              .calc_stats = false,
              .x_bins = 0,
              .y_bins = 0,
              .output_interval = 0,
              .velocity_output_location = "",
              .density_output_location = "",
          },
  };
  auto [input_file, step_size, write_checkpoint] =
      CLArgumentParser::parse(argc, argv);

  std::vector<Particle> particles;

  XmlReader::read(particles, input_file, arguments);

  printConfiguration(arguments);
  SpdWrapper::get()->info("Step size: {}", step_size);

  // TODO: delte in config
  std::vector<std::unique_ptr<IndexForce>> index_forces;
  for (const auto& config : arguments.index_force_configs) {
    const auto& [coords, ids, time, force_values] = config;
    SpdWrapper::get()->info("Adding index force with indices {}", ids.size());
    index_forces.push_back(
        std::make_unique<IndexForce>(ids, time, force_values));
  }

  std::unique_ptr<ParticleContainer> container;
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    // linked_cells_data.particle_count = particles.size();
    // linked_cells_data.special_particle_count =
    //     std::count_if(particles.begin(), particles.end(),
    //                   [](const Particle& p) { return p.getType() < 0; });
    container = std::make_unique<LinkedCellsContainer>(linked_cells_data);
    container->addParticles(particles);
    container->imposeInvariant();
  } else if (std::holds_alternative<DirectSumConfig>(
                 arguments.container_data)) {
    container = std::make_unique<DirectSumContainer>();
    container->addParticles(particles);
  } else {
    SpdWrapper::get()->error("Unrecognized container type");
    throw std::runtime_error("Unrecognized container type");
  }
  SpdWrapper::get()->info("{} particles", particles.size());
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  for (auto& config : arguments.interactive_force_types) {
    if (std::holds_alternative<LennardJonesConfig>(config)) {
      interactive_forces.push_back(std::make_unique<LennardJones>());
    } else if (std::holds_alternative<GravityConfig>(config)) {
      interactive_forces.push_back(std::make_unique<Gravity>());
    } else {
      SpdWrapper::get()->error("Unrecognized interactive_force_type");
    }
  }

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  for (auto config : arguments.singular_force_types) {
    if (std::holds_alternative<SingularGravityConfig>(config)) {
      const auto& [g,a] = std::get<SingularGravityConfig>(config);
      singular_forces.push_back(
          std::move(std::make_unique<SingularGravity>(g,a)));
    } else if (std::holds_alternative<HarmonicForceConfig>(config)) {
      const auto& [r, k] = std::get<HarmonicForceConfig>(config);
      singular_forces.push_back(
          std::move(std::make_unique<HarmonicForce>(r, k)));
    } else {
      SpdWrapper::get()->error("Unrecognized singular force");
    }
  }

  VerletIntegrator verlet_integrator(interactive_forces, singular_forces,
                                     index_forces, arguments.delta_t);
  outputWriter::VTKWriter writer;
  std::unique_ptr<Thermostat> thermostat;
  const std::string outputDirectory =
      createOutputDirectory("./output/", argc, argv);
  if (arguments.use_thermostat) {
    thermostat = std::make_unique<Thermostat>(arguments.thermostat_config);
  }

  double current_time = 0;
  int iteration = 0;
  auto number_of_particles = particles.size();
  const auto start_time = std::chrono::high_resolution_clock::now();
  size_t particle_updates = 0;

#ifndef BENCHMARK
  int writes = 0;
  int percentage = 0;
  double next_output_time = 0;
  spdlog::stopwatch stopwatch;
  auto time_of_last_mups = start_time;
  Statistics statistics(
      arguments.statistics_config.x_bins, arguments.statistics_config.y_bins,
      *container,
      outputDirectory + "/" + arguments.statistics_config.density_output_location,
      outputDirectory + "/" + arguments.statistics_config.velocity_output_location);

#endif

  while (current_time <= arguments.t_end) {
    verlet_integrator.step(*container);
    if (arguments.use_thermostat) {
      if (iteration % thermostat->n_thermostat == 0 && iteration > 0) {
        //SpdWrapper::get()->info("Applying thermostat in iteration [{}] / time [{:.4}/{}]", iteration, current_time, arguments.t_end);
        thermostat->setTemperature(*container);
      }
    }
#ifdef BENCHMARK  // these are the first 1000 iterations for the contest
    if (iteration == 1000) {
      const auto first_1k = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> elapsed = first_1k - start_time;
      std::cout << "First 1k iterations took: " << elapsed.count() << " seconds"
                << std::endl;
      const auto mups = static_cast<double>(number_of_particles) * 1000 *
                        (1.0 / elapsed.count());
      std::cout << "MMUPS for first 1k iterations: " << mups * (1.0 / 1e6)
                << std::endl;
    }
#endif

    particle_updates += container->getParticleCount();

#ifndef BENCHMARK
    if (current_time >= next_output_time) {
      plotParticles(outputDirectory, iteration, writer, *container);
      writes++;
      next_output_time = writes * step_size;
      // check if next percentage complete
      if (const double t = 100 * current_time / arguments.t_end;
          t >= percentage) {
        auto elapsed = stopwatch.elapsed();
        auto eta = (elapsed / percentage) * 100 - elapsed;
        auto h = std::chrono::duration_cast<std::chrono::hours>(eta).count();
        eta -= std::chrono::hours(h);
        auto m = std::chrono::duration_cast<std::chrono::minutes>(eta).count();
        eta -= std::chrono::minutes(m);
        auto s = std::chrono::duration_cast<std::chrono::seconds>(eta).count();

        if (percentage == 0) {
          h = 0;
          m = 0;
          s = 0;
        }

        auto current_time_hrc = std::chrono::high_resolution_clock::now();
        auto microseconds =
            std::chrono::duration_cast<std::chrono::microseconds>(
                current_time_hrc - time_of_last_mups)
                .count();
        double mmups = particle_updates *
                       (1.0 / static_cast<double>(microseconds));
        time_of_last_mups = current_time_hrc;
        particle_updates = 0;

        // mmups are unaccounted for write time, therefore it is always a lower
        // bound
        SpdWrapper::get()->info(
            "[{:<3.0f}%]: Iteration {:<12} | [ETA: {}:{:02}:{:02}], [average "
            "MMUPS since last log: {:02}]",
            100 * current_time / arguments.t_end, iteration, h, m, s, mmups);

        percentage++;
      }
    }

    if (arguments.statistics_config.calc_stats &&
        iteration % arguments.statistics_config.output_interval == 0) {
      statistics.writeStatistics(current_time);
    }
#endif
    iteration++;
    current_time = arguments.delta_t * iteration;
  }

  // Writes the finished simulations particle state into a checkpoint file
  if (write_checkpoint) {
    XmlWriter::writeFile(*container, "./output/test.xml");
  }

#ifdef BENCHMARK
  const auto end_time = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed_time = end_time - start_time;
  std::cout << "Simulation Time: " << elapsed_time.count() << " seconds"
            << std::endl;
#endif
  auto current_time_hrc = std::chrono::high_resolution_clock::now();
  auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
                          current_time_hrc - start_time)
                          .count();
  double mmups = iteration * static_cast<double>(number_of_particles) *
                 (1.0 / static_cast<double>(microseconds));
  std::cout << "MMUPS: " << mmups << std::endl;

#ifndef BENCHMARK
  statistics.closeFiles();
#endif
  SpdWrapper::get()->info("Output written. Terminating...");

  return 0;
}