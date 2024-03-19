// Copyright 2024 James Squires
/**
#pragma once
constexpr double square(double a) {
  return a * a;
}

static double filter_gain() {

  double temp[FILTER_TAP_NUM] = {};

  for (int i = 0; i < FILTER_TAP_NUM; ++i) {
    temp[i] = square(filter_taps[i]);
  }

  double sum = 0.0;
  for (int i = 0; i < FILTER_TAP_NUM; ++i) {
    sum += temp[i];
  }

  return 1 / sum;
}
**/
