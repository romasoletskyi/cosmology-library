#pragma once

#include <vector>
#include <complex>
#include <cosmology/utility/spline.h>

struct HomogenousHistory {
    Spline a; // scale factor
    Spline Xe; // free electron rate
    Spline tau; // optical depth
};

void writeHomogenousHistoryTo(std::ostream& stream, const HomogenousHistory& history);

HomogenousHistory readHomogenousHistoryFrom(std::istream& stream);

HomogenousHistory getHomogenousHistory(int pointsNumber);

std::vector<std::complex<double>> getGeneratingFunction(const HomogenousHistory& history, double wavenumber);
