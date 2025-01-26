//
// Created by maximilian on 23.10.24.
//
#pragma once

#ifndef DEBUG_H
#define DEBUG_H

#ifdef DEBUG
#include <iostream>  // Also careful here, when clion is not in debug mode it tells you to remove this
#endif
#include "utils/SpdWrapper.h"

#ifdef DEBUG
#define DEBUG_PRINT(msg) SpdWrapper::get()->debug(msg);
#define DEBUG_PRINT_FMT(msg, ...) SpdWrapper::get()->debug(msg, __VA_ARGS__);
#else
#define DEBUG_PRINT(msg)
#define DEBUG_PRINT_FMT(msg, ...)
#endif

#define INFO(msg) SpdWrapper::get()->info(msg);
#define INFO_FMT(msg, ...) SpdWrapper::get()->info(msg, __VA_ARGS__);

// TODO: move or rename
template <typename T>
inline void InfoVec(std::string msg, std::array<T, 3> vec) {
  INFO_FMT("{} -- [{},{},{}]", msg, vec[0], vec[1], vec[2]);
}
#endif  // DEBUG_H
