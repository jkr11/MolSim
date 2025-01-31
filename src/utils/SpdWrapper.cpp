//
// Created by jkr on 11/1/24.
//
#include "SpdWrapper.h"

#include <spdlog/sinks/stdout_color_sinks.h>

#include "io/CLArgumentParser.h"

std::shared_ptr<spdlog::logger> SpdWrapper::instance_ = configure();

std::shared_ptr<spdlog::logger> SpdWrapper::get() { return instance_; }

std::shared_ptr<spdlog::logger> SpdWrapper::configure() {
  auto color_console_sink =
      std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  color_console_sink->set_color(spdlog::level::info, color_console_sink->blue);
  color_console_sink->set_color(spdlog::level::warn,
                                color_console_sink->yellow);
  color_console_sink->set_color(spdlog::level::err, color_console_sink->red);
  color_console_sink->set_color(spdlog::level::critical,
                                color_console_sink->red);
  color_console_sink->set_color(spdlog::level::trace,
                                color_console_sink->green);
  spdlog::init_thread_pool(8192, 1);
  auto async_logger = std::make_shared<spdlog::async_logger>(
      "asyncLogger", spdlog::sinks_init_list{color_console_sink},
      spdlog::thread_pool(), spdlog::async_overflow_policy::block);
  async_logger->set_level(spdlog::level::info);
  spdlog::register_logger(async_logger);
  return async_logger;
}

int SpdWrapper::setLogLevel(std::string level) {
  level = toLower(level);

  if (level == "trace") {
    instance_->set_level(spdlog::level::trace);
    return 0;
  } else if (level == "debug") {
    instance_->set_level(spdlog::level::debug);
    return 0;
  } else if (level == "info") {
    instance_->set_level(spdlog::level::info);
    return 0;
  } else if (level == "warn") {
    instance_->set_level(spdlog::level::warn);
    return 0;
  } else if (level == "error") {
    instance_->set_level(spdlog::level::err);
    return 0;
  } else if (level == "off") {
    instance_->set_level(spdlog::level::off);
    return 0;
  }

  return -1;
}