// Copyright 2024 James Squires

#pragma once

// #include "SC_PlugIn.hpp"
// #include "SC_PlugIn.h"

// static InterfaceTable* ft;

namespace Constants {
constexpr double TWO_PI = 6.28318530718;
constexpr double PI = 3.141592653589793238462643;
constexpr double PI_SQRD = 9.86960440109;
}  // namespace Constants

/**
namespace JSCDSP {


template <typename NumType>
struct SC_Audio_Buffer
{
  SC_Audio_Buffer<NumType>(unsigned int size, SCUnit* unit_gen) :
belongs_to(unit_gen), data((NumType*)RTAlloc(unit_gen->mWorld, size *
sizeof(NumType))), buffer_size(size) {}; ~SC_Audio_Buffer<NumType>() {
RTFree(belongs_to->mWorld, data); }; unsigned int buffer_size; NumType* data;
  SCUnit* belongs_to;

};

} // JSCDSP

**/
