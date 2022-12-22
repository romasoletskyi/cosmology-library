#pragma once

#include "parameters.h"
#include "odesolver.h"
#include "utility.h"

struct HomogenousHistory {
    std::vector<float> eta; // conformal time
    std::vector<float> a; // scale factor
};

HomogenousHistory getHomogenousHistory(float aStart, float etaEnd, int points);

std::vector<float> getTransferFunctions(float kMode) {
    OdeSolver<float, float> solver;


}
