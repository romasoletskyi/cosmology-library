#pragma once

#include <vector>
#include <cosmology/utility/spline.h>

struct HomogenousHistory {
    Spline a; // scale factor
    Spline Xe; // free electron rate
};

HomogenousHistory getHomogenousHistory(int pointsNumber);

std::vector<double> getTransferFunctions(double kMode);
