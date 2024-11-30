//
// Created by jkr on 11/28/24.
//
#include "Simulation.h"

#include <spdlog/spdlog.h>

#include <utility>

#include "calc/VerletIntegrator.h"
#include "containers/DirectSumContainer.h"
#include "containers/LinkedCellsContainer.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "io/file/out/OutputHelper.h"
#include "io/file/out/VTKWriter.h"
#include "spdlog/stopwatch.h"

Simulation::Simulation(Arguments arguments, std::string output_directory, const double step_size)
    : arguments(std::move(arguments)),
      output_directory(std::move(output_directory)), step_size(step_size) {}

Simulation::~Simulation() = default;

void Simulation::initParams() {
  if (std::holds_alternative<LinkedCellsConfig>(arguments.container_data)) {
    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(arguments.container_data);
    container = std::make_unique<LinkedCellsContainer>(linked_cells_data);
  } else if (std::holds_alternative<DirectSumConfig>(
                 arguments.container_data)) {
    container = std::make_unique<DirectSumContainer>();
  } else {
    SpdWrapper::get()->error("Unrecognized container type");
    throw std::runtime_error("Unrecognized container type");
  }

  std::unique_ptr<Force> force;
  if (arguments.force_type == Arguments::Gravity) {
    force = std::make_unique<Gravity>();
  } else if (arguments.force_type == Arguments::LennardJones) {
    force = std::make_unique<LennardJones>();
  }

  const VerletIntegrator verlet_integrator(*force, arguments.delta_t);
  integrator = std::make_unique<VerletIntegrator>(verlet_integrator);
}

void Simulation::initParticles() {
  for (const auto& generator : arguments.generator_configs) {
    generator->generate(particles);
  }
}

void Simulation::run() const {
  outputWriter::VTKWriter writer;
  double current_time = 0;
  int iteration = 0;
  int writes = 0;
  int percentage = 0;
  double next_output_time = 0;
  spdlog::stopwatch stopwatch;
  while (current_time <= arguments.t_end) {
    integrator->step(*container);

    if (current_time >= next_output_time) {
      plotParticles(output_directory, iteration, writer, *container);
      writes++;
      next_output_time = writes * step_size;
      // check if next percentage complete
      if (const double t = 100 * current_time / arguments.t_end;
          t >= percentage) {
        percentage++;
        auto elapsed = stopwatch.elapsed();
        auto eta = (elapsed / percentage) * 100 - elapsed;
        auto h = std::chrono::duration_cast<std::chrono::hours>(eta).count();
        eta -= std::chrono::hours(h);
        auto m = std::chrono::duration_cast<std::chrono::minutes>(eta).count();
        eta -= std::chrono::minutes(m);
        auto s = std::chrono::duration_cast<std::chrono::seconds>(eta).count();

        SpdWrapper::get()->info(
            "[{:.0f} %]: Iteration {:<12} | [ETA: {}:{:02}:{:02}]",
            100 * current_time / arguments.t_end, iteration, h, m, s);
      }
    }

    iteration++;
    current_time = arguments.delta_t * iteration;
  }
  SpdWrapper::get()->info("Output written. Terminating...");
}