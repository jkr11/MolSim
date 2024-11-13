//
// Created by maximilian on 23.10.24.
//
#pragma once

#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>

#ifdef DEBUG
#define DEBUG_PRINT(msg) SpdWrapper::get()->debug(msg);
#else
#define DEBUG_PRINT(msg)
#endif

#endif  // DEBUG_H
