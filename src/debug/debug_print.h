//
// Created by maximilian on 23.10.24.
//
#pragma once

#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>

#include "utils/SpdWrapper.h"

#ifdef DEBUG
#define DEBUG_PRINT(msg) SpdWrapper::get()->debug(msg);
#define DEBUG_PRINT_FMT(msg, ...) SpdWrapper::get()->debug(msg, __VA_ARGS__);
#else
#define DEBUG_PRINT(msg)
#define DEBUG_PRINT_FMT(msg, ...)
#endif

#endif  // DEBUG_H
