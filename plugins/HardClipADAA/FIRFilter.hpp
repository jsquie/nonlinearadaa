// Copyright 2024 James Squires
#pragma once

/*

FIR filter designed with
http://t-filter.appspot.com

sampling frequency: 88200 Hz

* 0 Hz - 22250 Hz
  gain = 1
  desired ripple = 1 dB
  actual ripple = 0.6841297988737765 dB

* 32000 Hz - 44100 Hz
  gain = 0
  desired attenuation = -120 dB
  actual attenuation = -120.82101778440254 dB

*/

#define UP_FILTER_TAP_NUM 33
#define UP_FILTER_TAP_NUM_F 33.0f

constexpr int M_up = UP_FILTER_TAP_NUM;
constexpr float fM_up = UP_FILTER_TAP_NUM_F;

constexpr double up_filter_taps[UP_FILTER_TAP_NUM] = {
    -0.0012246091154708878, -0.006616426278879418,  -0.01273110279401792,
    -0.00605363630895313,   0.010022443831607273,   0.005056295420376784,
    -0.015811649653449722,  -0.0018617274019090756, 0.0255058609629138,
    -0.0077493311209674304, -0.03771969766161676,   0.030327526731206006,
    0.05006541556063978,    -0.08292061503807464,   -0.05936057831504514,
    0.3099774184144736,     0.5628274952286104,     0.3099774184144736,
    -0.05936057831504514,   -0.08292061503807464,   0.05006541556063978,
    0.030327526731206006,   -0.03771969766161676,   -0.0077493311209674304,
    0.0255058609629138,     -0.0018617274019090756, -0.015811649653449722,
    0.005056295420376784,   0.010022443831607273,   -0.00605363630895313,
    -0.01273110279401792,   -0.006616426278879418,  -0.0012246091154708878};

#define DOWN_FILTER_TAP_NUM 209
#define DOWN_FILTER_TAP_NUM_F 209.0f

constexpr int M_down = DOWN_FILTER_TAP_NUM;
constexpr float fM_down = DOWN_FILTER_TAP_NUM_F;


constexpr double down_filter_taps[DOWN_FILTER_TAP_NUM] = {
    -0.0016550762341430308,  -0.0041665939310652385,  -0.003594191447486389,
    0.00251719066089494,     0.009204232105743072,    0.00891685978767867,
    0.0025929990340545586,   -0.0009241688686624757,  0.001780417290846837,
    0.0040583237055417075,   0.0012718880881058892,   -0.0015563530079133385,
    0.0002721535795504213,   0.0025035930967126315,   0.0005429619868192588,
    -0.0017453123023537607,  0.000022984138603901914, 0.0020537985664696805,
    0.00013860481511378278,  -0.0018438080070986866,  0.00013467879098560712,
    0.001974941764859404,    -0.00018853242904265353, -0.0019542465345495637,
    0.0003762055310051086,   0.0020296451360248486,   -0.0005169276482732894,
    -0.0020684792621014857,  0.0006994817682676137,   0.0021265551428543772,
    -0.0008822257189953905,  -0.0021674203444939926,  0.0010935813680486932,
    0.0022115890943483034,   -0.001314750098015205,   -0.0022412051839012007,
    0.0015560490952100908,   0.0022569517475881764,   -0.0018209974880386698,
    -0.002266301757028872,   0.0020937340035129225,   0.0022477597114187653,
    -0.0023927521827856475,  -0.0022146568586112986,  0.0027050436381634005,
    0.002156314932451038,    -0.003031065300689767,   -0.002064650858536239,
    0.0033792643077078716,   0.0019474808728353638,   -0.0037376568321684455,
    -0.0017910129887243696,  0.004112825480457824,    0.0015960522407615228,
    -0.004502958711563328,   -0.0013612537988873558,  0.004903329017739637,
    0.001079361495260796,    -0.005313440162938722,   -0.000745420496662418,
    0.00573623700664673,     0.00036334614253456674,  -0.006156309094600619,
    0.00008808856155454664,  0.006589753951165301,    -0.0005956547623096374,
    -0.007014724209099923,   0.0011845526388683388,   0.007445042007534717,
    -0.0018497468084161031,  -0.007868705069367566,   0.0026066792019242466,
    0.008285357331304153,    -0.0034688911024595,     -0.008698582991922315,
    0.0044411237414403925,   0.009092986565733343,    -0.005557414337636592,
    -0.009483553419766278,   0.006823202063313745,    0.009848519291251277,
    -0.008289679120818303,   -0.010201859443027823,   0.009988471477688032,
    0.01053012360277123,     -0.011992197126080751,   -0.010834696922933549,
    0.014396832948754698,    0.011115959799250405,    -0.017345990204051556,
    -0.01136588090176962,    0.021087101188298568,    0.01158952167903816,
    -0.026030444783605375,   -0.011780011049371056,   0.03296919343261422,
    0.011939078466010828,    -0.04359427254454018,    -0.012063942524891857,
    0.0623135104924757,      0.012154411619644348,    -0.10529117509190798,
    -0.01220839945864449,    0.3180391210647029,      0.5122272478248557,
    0.3180391210647029,      -0.01220839945864449,    -0.10529117509190798,
    0.012154411619644348,    0.0623135104924757,      -0.012063942524891857,
    -0.04359427254454018,    0.011939078466010828,    0.03296919343261422,
    -0.011780011049371056,   -0.026030444783605375,   0.01158952167903816,
    0.021087101188298568,    -0.01136588090176962,    -0.017345990204051556,
    0.011115959799250405,    0.014396832948754698,    -0.010834696922933549,
    -0.011992197126080751,   0.01053012360277123,     0.009988471477688032,
    -0.010201859443027823,   -0.008289679120818303,   0.009848519291251277,
    0.006823202063313745,    -0.009483553419766278,   -0.005557414337636592,
    0.009092986565733343,    0.0044411237414403925,   -0.008698582991922315,
    -0.0034688911024595,     0.008285357331304153,    0.0026066792019242466,
    -0.007868705069367566,   -0.0018497468084161031,  0.007445042007534717,
    0.0011845526388683388,   -0.007014724209099923,   -0.0005956547623096374,
    0.006589753951165301,    0.00008808856155454664,  -0.006156309094600619,
    0.00036334614253456674,  0.00573623700664673,     -0.000745420496662418,
    -0.005313440162938722,   0.001079361495260796,    0.004903329017739637,
    -0.0013612537988873558,  -0.004502958711563328,   0.0015960522407615228,
    0.004112825480457824,    -0.0017910129887243696,  -0.0037376568321684455,
    0.0019474808728353638,   0.0033792643077078716,   -0.002064650858536239,
    -0.003031065300689767,   0.002156314932451038,    0.0027050436381634005,
    -0.0022146568586112986,  -0.0023927521827856475,  0.0022477597114187653,
    0.0020937340035129225,   -0.002266301757028872,   -0.0018209974880386698,
    0.0022569517475881764,   0.0015560490952100908,   -0.0022412051839012007,
    -0.001314750098015205,   0.0022115890943483034,   0.0010935813680486932,
    -0.0021674203444939926,  -0.0008822257189953905,  0.0021265551428543772,
    0.0006994817682676137,   -0.0020684792621014857,  -0.0005169276482732894,
    0.0020296451360248486,   0.0003762055310051086,   -0.0019542465345495637,
    -0.00018853242904265353, 0.001974941764859404,    0.00013467879098560712,
    -0.0018438080070986866,  0.00013860481511378278,  0.0020537985664696805,
    0.000022984138603901914, -0.0017453123023537607,  0.0005429619868192588,
    0.0025035930967126315,   0.0002721535795504213,   -0.0015563530079133385,
    0.0012718880881058892,   0.0040583237055417075,   0.001780417290846837,
    -0.0009241688686624757,  0.0025929990340545586,   0.00891685978767867,
    0.009204232105743072,    0.00251719066089494,     -0.003594191447486389,
    -0.0041665939310652385,  -0.0016550762341430308};

constexpr double square(double a) { return a * a; }

constexpr double up_filter_gain() {
  double temp[UP_FILTER_TAP_NUM] = {};

  for (int i = 0; i < UP_FILTER_TAP_NUM; ++i) {
    temp[i] = square(up_filter_taps[i]);
  }

  double sum = 0.0;
  for (int i = 0; i < UP_FILTER_TAP_NUM; ++i) {
    sum += temp[i];
  }

  return 1 / sum;
}

constexpr double down_filter_gain() {
  double temp[DOWN_FILTER_TAP_NUM] = {};

  for (int i = 0; i < DOWN_FILTER_TAP_NUM; ++i) {
    temp[i] = square(down_filter_taps[i]);
  }

  double sum = 0.0;
  for (int i = 0; i < DOWN_FILTER_TAP_NUM; ++i) {
    sum += temp[i];
  }

  return 1 / sum;
}
