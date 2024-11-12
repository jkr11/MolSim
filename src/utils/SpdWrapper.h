//
// Created by jkr on 11/1/24.
//
#pragma once
#ifndef SPDWRAPPER_H
#define SPDWRAPPER_H

#include <spdlog/async.h>
#include <spdlog/spdlog.h>

#include <string>
/**
 * @brief class exposing a globally accessible spdlog::logger
 *
 * SpdWrapper manages a singleton instance of spdlog::logger that is accesible
 * from everywhere in the application
 */
class SpdWrapper {
 public:
  /**
   * @brief retrieves the singleton instance
   * @return A shared pointer to the singleton instance
   */
  static std::shared_ptr<spdlog::logger> get();
  /**
   * @brief sets the loglevel for the singleton instance
   * @param level string value indicating the desired log-level, e.g.
   * "info","debug",etc.
   * @return an integer rating the success of the operation, 0 if good
   */
  static int setLogLevel(std::string level);

 private:
  /**
   * @brief configures the singleton instance and registers it with spdlog
   * @return shared pointer to the configured asnyc_logger
   */
  static std::shared_ptr<spdlog::logger> configure();
  /**
   * @brief holds the singleton instance of the logger
   */
  static std::shared_ptr<spdlog::logger> instance;
};

#endif  // SPDWRAPPER_H
