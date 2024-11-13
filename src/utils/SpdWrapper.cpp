//
// Created by jkr on 11/1/24.
//
#include "SpdWrapper.h"

#include <spdlog/sinks/stdout_color_sinks.h>

#include "io/CLArgumentParser.h"

std::shared_ptr<spdlog::logger> SpdWrapper::instance = configure();

std::shared_ptr<spdlog::logger> SpdWrapper::get() { return instance; }

std::shared_ptr<spdlog::logger> SpdWrapper::configure() {
  auto colorConsoleSink =
      std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  colorConsoleSink->set_color(spdlog::level::info, colorConsoleSink->blue);
  colorConsoleSink->set_color(spdlog::level::warn, colorConsoleSink->yellow);
  colorConsoleSink->set_color(spdlog::level::err, colorConsoleSink->red);
  colorConsoleSink->set_color(spdlog::level::critical, colorConsoleSink->red);
  colorConsoleSink->set_color(spdlog::level::trace, colorConsoleSink->green);
  spdlog::init_thread_pool(8192, 1);
  auto asyncLogger = std::make_shared<spdlog::async_logger>(
      "asnycLogger", spdlog::sinks_init_list{colorConsoleSink},
      spdlog::thread_pool(), spdlog::async_overflow_policy::block);
  asyncLogger->set_level(spdlog::level::info);
  spdlog::register_logger(asyncLogger);
  return asyncLogger;
}

int SpdWrapper::setLogLevel(std::string level) {
  level = toLower(level);

  if (level == "trace") {
    instance->set_level(spdlog::level::trace);
    return 0;
  } else if (level == "debug") {
    instance->set_level(spdlog::level::debug);
    return 0;
  } else if (level == "info") {
    instance->set_level(spdlog::level::info);
    return 0;
  } else if (level == "warn") {
    instance->set_level(spdlog::level::warn);
    return 0;
  } else if (level == "error") {
    instance->set_level(spdlog::level::err);
    return 0;
  } else if (level == "off") {
    instance->set_level(spdlog::level::off);
    return 0;
  }

  return -1;
}