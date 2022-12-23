#pragma once

#include <cosmology/odesolver/odesolver.h>
#include <cosmology/utility/math.h>
#include "parameters.h"

struct HomogenousHistory {
    std::vector<float> eta; // conformal time
    std::vector<float> a; // scale factor
};

HomogenousHistory getHomogenousHistory(float aStart, float etaEnd, int points);

std::vector<float> getTransferFunctions(float kMode) {
    OdeSolver<float, float> solver;


}
