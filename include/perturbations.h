#pragma once

#include "parameters.h"
#include "odesolver.h"
#include "utility.h"

struct HomogenousHistory {
    std::vector<float> eta; // conformal time
    std::vector<float> a; // scale factor
};

HomogenousHistory getHomogenousHistory(float aStart, float etaEnd, int points) {
    OdeSolver<float, float> solver;
    CosmoParameters cosmo;

    solver.addVariable("a");
    solver.addExplicitFunction("da", [cosmo](float time, std::vector<float> data) {
        float a = data[0];
        return a * a * std::sqrt(cosmo.omegaLambda + (cosmo.omegaBaryon + cosmo.omegaCold) / std::pow(a, 3) +
                                 (cosmo.omegaPhoton + cosmo.omegaNeutrino) / std::pow(a, 4));
    }, {"a"});

    solver.addEquation("a' = da");
    auto grid = linSpace<float>(0, etaEnd, points);
    solver.setGrid(grid);
    auto solution = solver.solve<walker::RungeKuttaWalker<float, float>>({aStart});

    HomogenousHistory history;
    history.a.reserve(grid.size());
    for (int i = 0; i < grid.size(); ++i) {
        history.a.push_back(solution(i, 0));
    }
    history.eta = std::move(grid);

    return history;
}

std::vector<float> getTransferFunctions(float kMode) {
    OdeSolver<float, float> solver;


}
