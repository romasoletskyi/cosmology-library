#include <optional>
#include <iostream>

#include <cosmology/perturbations.h>
#include <cosmology/odesolver/odesolver.h>
#include <cosmology/utility/math.h>
#include <cosmology/parameters.h>

constexpr double pi = 3.14159265358979323846;

double computeRelativeDifference(const Solution<double> &lhs, const Solution<double> &rhs) {
    double difference = 0;
    for (int i = 0; i < lhs.gridLength * lhs.variableNumber; ++i) {
        if (!(i & 0x1) && lhs.data[i] > 1) {
            break;
        }
        double relative = std::abs(lhs.data[i] - rhs.data[i]) / (1e-9 + std::abs(lhs.data[i]) + std::abs(rhs.data[i]));
        if (relative > difference) {
            difference = relative;
        }
    }
    return difference;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, bool>
progressEvolution(OdeSolver<double, double> &solver, double etaStart, double aStart, double XeStart, double &etaDelta) {
    int length = 10000;
    std::vector<double> grid;
    Solution<double> solution;
    std::optional<bool> downscaling;

    while (true) {
        auto largeGrid = linSpace(etaStart, etaStart + 2 * etaDelta, length);
        auto smallGrid = linSpace(etaStart, etaStart + etaDelta, length);

        solver.setGrid(smallGrid);
        auto smallGridSolution = solver.solve({aStart, XeStart});
        solver.setGrid(largeGrid);
        auto largeGridSolution = solver.solve({aStart, XeStart});

        if (smallGridSolution(smallGridSolution.gridLength - 1, 0) > 1) {
            grid = smallGrid;
            solution = std::move(smallGridSolution);
            break;
        }

        if (computeRelativeDifference(smallGridSolution, largeGridSolution) < 1e-2) {
            etaDelta *= 2;
            if (downscaling.has_value()) {
                if (!downscaling.value()) {
                    etaDelta /= 2;
                    grid = smallGrid;
                    solution = std::move(smallGridSolution);
                    break;
                }
            } else {
                downscaling = true;
            }
        } else {
            etaDelta /= 2;
            if (downscaling.has_value()) {
                if (downscaling.value()) {
                    grid = largeGrid;
                    solution = std::move(largeGridSolution);
                    break;
                }
            } else {
                downscaling = false;
            }
        }
    }

    bool finished = false;
    std::vector<double> etaStep, aStep, XeStep;
    etaStep.reserve(grid.size());
    aStep.reserve(grid.size());
    XeStep.reserve(grid.size());

    for (int i = 0; i < grid.size(); ++i) {
        if (solution(i, 0) > 1) {
            finished = true;
            break;
        }
        etaStep.push_back(grid[i]);
        aStep.push_back(solution(i, 0));
        XeStep.push_back(solution(i, 1));
    }

    return {etaStep, aStep, XeStep, finished};
}

std::vector<std::pair<double, double>>
samplePoints(const std::vector<double> &time, const std::vector<double> &value, int pointsNumber) {
    std::vector<std::pair<double, double>> points;
    points.reserve(pointsNumber + 1);

    size_t step = (time.size() - 1) / pointsNumber;
    for (size_t i = 0; i < time.size() - 1; i += step) {
        points.emplace_back(time[i], value[i]);
    }

    points.emplace_back(time.back(), value.back());
    return points;
}

HomogenousHistory getHomogenousHistory(int pointsNumber) {
    CosmoParameters cosmo;
    PhysParameters phys;
    OdeSolver<double, double> solver;

    solver.addVariable("a");
    solver.addExplicitFunction("da", [cosmo](double time, std::vector<double> data) {
        double a = data[0];
        return std::sqrt(cosmo.omegaLambda * std::pow(a, 4) + (cosmo.omegaBaryon + cosmo.omegaCold) * a +
                         (cosmo.omegaPhoton + cosmo.omegaNeutrino));
    }, {"a"});

    solver.addVariable("Xe");
    solver.addExplicitFunction("dXe", [cosmo, phys](double time, std::vector<double> data) {
        double a = data[0];
        double Xe = data[1];

        double T = cosmo.cmb / a / phys.evInKelvin;
        double nb = phys.critInEv * cosmo.omegaBaryon / std::pow(a, 3);
        double cross_section = 9.78 * std::pow(phys.alpha, 2) / std::pow(phys.me, 2) * std::sqrt(phys.delta / T) *
                               std::log(phys.delta / T);

        return a * phys.evInCosmoFrequency * cross_section *
               (std::pow(phys.me * T / 2 / pi, 1.5) * std::exp(-phys.delta / T) * (1 - Xe) - nb * Xe * Xe);
    }, {"a", "Xe"});

    solver.addEquation("a' = da");
    solver.addEquation("Xe' = dXe");
    solver.compile<walker::RungeKuttaWalker<double, double>>();

    double etaStart = 0;
    double etaDelta = 2e-4;
    double aStart = 6.6e-4;
    double TStart = cosmo.cmb / aStart / phys.evInKelvin;
    double nbStart = phys.critInEv * cosmo.omegaBaryon / std::pow(aStart, 3);
    double XeStart = 1 - nbStart * std::exp(phys.delta / TStart) * std::pow(phys.me * TStart / 2 / pi, -1.5);
    std::vector<double> eta, a, Xe;

    while (true) {
        auto [etaStep, aStep, XeStep, finished] = progressEvolution(solver, etaStart, aStart, XeStart, etaDelta);

        eta.insert(eta.end(), etaStep.begin() + 1, etaStep.end());
        a.insert(a.end(), aStep.begin() + 1, aStep.end());
        Xe.insert(Xe.end(), XeStep.begin() + 1, XeStep.end());

        etaStart = eta.back();
        aStart = a.back();
        XeStart = Xe.back();

        if (finished) {
            break;
        }
    }

    return HomogenousHistory{SplineBuilder().buildFromPoints(samplePoints(eta, a, pointsNumber)),
                             SplineBuilder().buildFromPoints(samplePoints(eta, Xe, pointsNumber))};
}
