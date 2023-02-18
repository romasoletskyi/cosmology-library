#pragma once

#include <vector>
#include <complex>
#include <cosmology/parameters.h>
#include <cosmology/odesolver/odesolver.h>
#include <cosmology/utility/spline.h>

struct HomogenousHistory {
    Spline a; // scale factor
    Spline Xe; // free electron rate
    Spline tau; // optical depth
};

void writeHomogenousHistoryTo(std::ostream &stream, const HomogenousHistory &history);

HomogenousHistory readHomogenousHistoryFrom(std::istream &stream);

HomogenousHistory getHomogenousHistory(int pointsNumber, const CosmoParameters &cosmo);

Solution<std::complex<double>>
getPerturbations(const CosmoParameters &cosmo, const HomogenousHistory &history, double wavenumber);

Spline getGeneratingFunction(const HomogenousHistory &history, const Solution<std::complex<double>> &solution);
