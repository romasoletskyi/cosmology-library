#pragma once

#include <vector>
#include <cosmology/utility/spline.h>

struct HomogenousHistory {
    Spline a; // scale factor
    Spline Xe; // free electron rate
};

void writeHomogenousHistoryTo(std::ostream& stream, const HomogenousHistory& history);

HomogenousHistory readHomogenousHistoryFrom(std::istream& stream);

HomogenousHistory getHomogenousHistory(int pointsNumber);

std::vector<double> getTransferFunctions(double kMode);
