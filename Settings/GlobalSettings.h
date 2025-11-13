#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

// Bin edge control
const double GLOBAL_epsilon = 0.0000000000000000001;
const std::array<std::array<float, 2>, 2> drawnWindowAuto = {{{-999, -999}, {-999, -999}}}; // {{{xmin, xmax}, {ymin, ymax}}}
const std::array<std::array<float, 2>, 3> drawnWindow2DAuto = {{{-999, -999}, {-999, -999}, {-999, -999}}}; // {{{xmin, xmax}, {ymin, ymax}, {zmin, zmax}}}
const std::array<std::array<float, 2>, 2> legendPlacementAuto = {{{-999, -999}, {-999, -999}}}; // {{{x1, y1}, {x2, y2}}}
const std::array<float, 2> contextPlacementAuto = {{-999, -999}}; // {{x_topleft, y_topleft}}
const TString dummyName[1] = {""};
const int N_SigmaBarlow = 2;

int hmcCollHistIsObsolete = false;

double th2ContoursNone[1] = {-999.}; int contourNumberNone = -99;

#endif
