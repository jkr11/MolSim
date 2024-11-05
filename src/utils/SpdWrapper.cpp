//
// Created by jkr on 11/1/24.
//
#include "SpdWrapper.h"
#include <spdlog/sinks/stdout_color_sinks.h>

std::shared_ptr<spdlog::logger> SpdWrapper::instance = SpdWrapper::configure();

std::shared_ptr<spdlog::logger> SpdWrapper::get() {
  return instance;
}

std::shared_ptr<spdlog::logger> SpdWrapper::configure() {
  auto colorConsoleSink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  colorConsoleSink->set_color(spdlog::level::info, colorConsoleSink->blue);
  colorConsoleSink->set_color(spdlog::level::warn, colorConsoleSink->yellow);
  colorConsoleSink->set_color(spdlog::level::err, colorConsoleSink->red);
  colorConsoleSink->set_color(spdlog::level::critical, colorConsoleSink->red);
  colorConsoleSink->set_color(spdlog::level::trace, colorConsoleSink->green);
  spdlog::init_thread_pool(8192, 1);
  auto asyncLogger = std::make_shared<spdlog::async_logger>(
      "asnycLogger", spdlog::sinks_init_list{colorConsoleSink}, spdlog::thread_pool(), spdlog::async_overflow_policy::block);
  asyncLogger->set_level(spdlog::level::info);
  spdlog::register_logger(asyncLogger);
  return asyncLogger;
}