#include <iostream>
#include <fstream>
#include <cmath>

#include <cosmology/perturbations.h>
#include <cosmology/monte_carlo.h>
#include <cosmology/utility/math.h>

// source for spherical bessel functions - https://github.com/joeydumont/complex_bessel
#include "complex_bessel.h"

std::tuple<std::vector<int>, std::vector<double>, std::vector<double>>
readExperimentalCMBSpectrum(const std::string &path) {
    auto file = std::ifstream(path);
    int size;
    file >> size;

    std::vector<int> multipole;
    std::vector<double> spectrum;
    std::vector<double> error;

    for (int i = 0; i < size; ++i) {
        int l;
        double cl, err;
        file >> l >> cl >> err;

        multipole.push_back(l);
        spectrum.push_back(cl);
        error.push_back(err);
    }
    file.close();

    return {multipole, spectrum, error};
}

void printParameters(std::ostream &stream, const CosmoParameters &cosmo) {
    stream << cosmo.omegaBaryon << " " << cosmo.omegaCold << " " << cosmo.omegaLambda << " " << cosmo.h << " "
           << cosmo.ns << " " << cosmo.cmbScale << "\n";
}

void fitCMBSpectrum(const std::string &experimentalPath, std::ostream &logStream) {
    CosmoParameters start;
    auto grid = getHomogenousHistory(1000, start).a.getKnots();

    int lMax = 3000;
    int kPoints = 40;
    int etaFactor = 10;
    int etaPoints = static_cast<int>(grid.size() - 1) / etaFactor + 1;

    double etaEnd = grid.back();
    double kMax = 2 * lMax / etaEnd;

    double k[kPoints];
    double source[kPoints * etaPoints];

    k[0] = 0.1;
    for (int i = 1; i < kPoints; ++i) {
        k[i] = i * kMax / kPoints;
    }

    std::vector<double> factoredGrid;
    factoredGrid.reserve(etaPoints);

    for (int j = 0; j < etaPoints; ++j) {
        factoredGrid.push_back(grid[j * etaFactor]);
    }

    int kFinePoints = 2 * lMax;
    double kFine[kFinePoints];
    auto sourceFine = new double[kFinePoints * etaPoints];

    for (int i = 0; i < kFinePoints; ++i) {
        kFine[i] = (i + 1) * kMax / kFinePoints;
    }

    auto [multipole, spectrum, error] = readExperimentalCMBSpectrum(experimentalPath);
    auto bessel = new double[multipole.size() * kFinePoints * etaPoints];

    for (int l = 0; l < multipole.size(); ++l) {
        for (int i = 0; i < kFinePoints; ++i) {
            for (int j = 0; j < etaPoints; ++j) {
                bessel[l * kFinePoints * etaPoints + i * etaPoints + j] =
                        sp_bessel::sph_besselJ(multipole[l], kFine[i] * (etaEnd - factoredGrid[j])).real();
            }
        }
    }

    logStream << "omega_b omega_c omega_lambda H ns\n";

    Sampler<double, CosmoParameters> sampler(spectrum, error, [&](const CosmoParameters &cosmo) {
        auto history = getHomogenousHistory(1000, cosmo);

        for (int i = 0; i < kPoints; ++i) {
            auto solution = getPerturbations(cosmo, history, k[i]);
            auto S = getGeneratingFunction(history, solution);
            for (int j = 0; j < etaPoints; ++j) {
                if (factoredGrid[j] > history.a.getKnots().back()) {
                    source[i * etaPoints + j] = 0;
                } else {
                    source[i * etaPoints + j] = S.evaluate(factoredGrid[j]);
                }
            }
        }

        for (int j = 0; j < etaPoints; ++j) {
            std::vector<std::pair<double, double>> points;
            points.reserve(kPoints);

            for (int i = 0; i < kPoints; ++i) {
                points.emplace_back(k[i], source[i * etaPoints + j]);
            }

            auto spline = SplineBuilder().buildFromPoints(points);
            for (int i = 0; i < kFinePoints; ++i) {
                sourceFine[i * etaPoints + j] = spline.evaluate(kFine[i]);
            }
        }

        std::vector<double> prediction;
        for (int l = 0; l < multipole.size(); ++l) {
            double buffer[etaPoints];
            double sourceK[kFinePoints];

            for (int i = 0; i < kFinePoints; ++i) {
                for (int j = 0; j < etaPoints; ++j) {
                    buffer[j] = sourceFine[i * etaPoints + j] * bessel[l * kFinePoints * etaPoints + i * etaPoints + j];
                }
                sourceK[i] = std::pow(integrate(factoredGrid.data(), buffer, etaPoints), 2) /
                             std::pow(kFine[i], 2 - cosmo.ns);
            }

            prediction.push_back(
                    integrate(kFine, sourceK, kFinePoints) * multipole[l] * (multipole[l] + 1) * cosmo.cmbScale);
        }

        logStream << "================================STEP================================\n";
        printParameters(logStream, cosmo);
        for (int l = 0; l < multipole.size(); ++l) {
            logStream << multipole[l] << " " << prediction[l] << "\n";
        }

        return prediction;
    }, [](const CosmoParameters &cosmo) {
        return 1.0;
    }, [](const CosmoParameters &cosmo, std::mt19937 &gen) {
        std::uniform_real_distribution<> uni(-1, 1);
        auto copy = cosmo;
        copy.h *= (1 + 0.05 * uni(gen));
        copy.cmbScale *= (1 + 0.1 * uni(gen));
        copy.ns *= (1 + 0.01 * uni(gen));
        copy.omegaCold *= (1 + 0.1 * uni(gen));
        copy.omegaBaryon *= (1 + 0.1 * uni(gen));
        copy.omegaLambda *= (1 + 0.1 * uni(gen));
        return copy;
    }, start);

    for (int i = 0; i < 1000; ++i) {
        auto [parameters, loglikelihood] = sampler.next();
        logStream << "================================NEXT================================\n";
        printParameters(logStream, parameters);
        logStream << loglikelihood << "\n";
        logStream.flush();
    }

    delete[] sourceFine;
    delete[] bessel;
}

int main() {
    std::ofstream file(std::string(PROJECT_PATH) + "/data/cmb-log");
    fitCMBSpectrum(std::string(PROJECT_PATH) + "/data/cmb-planck", file);
    file.close();
    return 0;
}