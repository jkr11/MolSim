//
// Created by maximilian on 23.10.24.
//
#pragma once

#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>

#ifdef DEBUG
#define inline     \
  DEBUG_PRINT(msg) \
  std::cout << "\t\033[38;2;255;165;0m " << msg << "\033[0m" << std::endl
#else
#define DEBUG_PRINT(msg)
#endif

#endif  // DEBUG_H
