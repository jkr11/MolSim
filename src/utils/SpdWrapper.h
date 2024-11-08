//
// Created by jkr on 11/1/24.
//
#pragma once
#ifndef SPDWRAPPER_H
#define SPDWRAPPER_H

#include <spdlog/async.h>
#include <spdlog/spdlog.h>

#include <string>

class SpdWrapper {
 public:
  static std::shared_ptr<spdlog::logger> get();
  static int setLogLevel(std::string level);

 private:
  static std::shared_ptr<spdlog::logger> configure();
  static std::shared_ptr<spdlog::logger> instance;
};

#endif  // SPDWRAPPER_H
