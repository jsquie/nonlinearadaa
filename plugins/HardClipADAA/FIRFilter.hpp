// Copyright 2024 James Squires
#pragma once

#define FILTER_TAP_NUM 48

constexpr int up_delay = (FILTER_TAP_NUM / 2);
constexpr int down_delay = (FILTER_TAP_NUM / 2) + 1;
constexpr double up_taps[FILTER_TAP_NUM] = {
    -0.0064715474097890545, 0.006788724784527351,  -0.007134125572070907,
    0.007511871271766723,   -0.007926929217098087, 0.00838534118242672,
    -0.00889453036904902,   0.009463720022395613,  -0.010104514094437885,
    0.010831718180021,      -0.011664525313602769, 0.012628270948224513,
    -0.013757103575462731,  0.015098181413680897,  -0.01671851963595936,
    0.01871667093508393,    -0.021243750540180146, 0.024543868940610197,
    -0.0290386730354654,    0.035524608815134716,  -0.045708348639099484,
    0.06402724397938601,    -0.10675158913607562,  0.32031404953367254,
    0.32031404953367254,    -0.10675158913607562,  0.06402724397938601,
    -0.045708348639099484,  0.035524608815134716,  -0.0290386730354654,
    0.024543868940610197,   -0.021243750540180146, 0.01871667093508393,
    -0.01671851963595936,   0.015098181413680897,  -0.013757103575462731,
    0.012628270948224513,   -0.011664525313602769, 0.010831718180021,
    -0.010104514094437885,  0.009463720022395613,  -0.00889453036904902,
    0.00838534118242672,    -0.007926929217098087, 0.007511871271766723,
    -0.007134125572070907,  0.006788724784527351,  -0.0064715474097890545,
}; 
constexpr double foldScaleCoef = 0.5031597730627207;
constexpr double upScaleCoef = 1 / foldScaleCoef;
