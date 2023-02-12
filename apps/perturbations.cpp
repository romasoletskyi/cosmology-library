#include <iostream>
#include <fstream>
#include <cmath>

#include <cosmology/perturbations.h>
#include <cosmology/utility/math.h>

// source for spherical bessel functions - https://github.com/joeydumont/complex_bessel
#include "complex_bessel.h"

int main() {
    /*auto history = getHomogenousHistory(1000);
    std::ofstream file(std::string(PROJECT_PATH) + "/data/homogenous-history");
    writeHomogenousHistoryTo(file, history);*/

    std::ifstream file(std::string(PROJECT_PATH) + "/data/homogenous-history");
    auto history = readHomogenousHistoryFrom(file);
    auto grid = history.a.getKnots();

    int lMax = 1500;
    int kPoints = 80;
    int etaFactor = 2;
    int etaPoints = static_cast<int>(grid.size() - 1) / etaFactor + 1;

    double etaEnd = grid.back();
    double kMax = 2 * lMax / etaEnd;

    double k[kPoints];
    double source[kPoints][etaPoints];
    std::vector<double> factoredGrid;
    factoredGrid.reserve(etaPoints);

    k[0] = 0.1;
    for (int i = 1; i < kPoints; ++i) {
        k[i] = i * kMax / kPoints;
    }

    for (int i = 0; i < kPoints; ++i) {
        auto S = getGeneratingFunction(history, k[i]);
        for(int j = 0; j < etaPoints; ++j) {
            source[i][j] = S[j * etaFactor].real();
            factoredGrid.push_back(grid[j * etaFactor]);
        }
    }

    int kFinePoints = 2 * lMax;
    double kFine[kFinePoints];
    auto sourceFine = new double[kFinePoints * etaPoints];

    for (int i = 0; i < kFinePoints; ++i) {
        kFine[i] = (i + 1) * kMax / kFinePoints;
    }

    for (int j = 0; j < etaPoints; ++j) {
        std::vector<std::pair<double, double>> points;
        points.reserve(kPoints);

        for (int i = 0; i < kPoints; ++i) {
            points.emplace_back(k[i], source[i][j]);
        }

        auto spline = SplineBuilder().buildFromPoints(points);
        for (int i = 0; i < kFinePoints; ++i) {
            sourceFine[i * etaPoints + j] = spline.evaluate(kFine[i]);
        }
    }

    std::vector<std::pair<int, double>> spectrum;

    for (int l = 20; l <= lMax; l += 20) {
        double buffer[etaPoints];
        double sourceK[kFinePoints];

        for (int i = 0; i < kFinePoints; ++i) {
            for (int j = 0; j < etaPoints; ++j) {
                buffer[j] = sourceFine[i * etaPoints + j] * sp_bessel::sph_besselJ(l, kFine[i] * (etaEnd -  factoredGrid[j])).real();
            }
            sourceK[i] = std::pow(integrate(factoredGrid.data(), buffer, etaPoints), 2) / kFine[i];
        }

        spectrum.emplace_back(l, integrate(kFine, sourceK, kFinePoints));
        auto [L, CL] = spectrum.back();
        std::cout << L <<  " " << CL << " " << CL * L * (L + 1) << std::endl;
    }
    delete[] sourceFine;
    return 0;
}