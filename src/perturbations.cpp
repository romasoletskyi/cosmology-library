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

std::vector<double> makeFineGrid(const std::vector<double>& grid, int factor) {
    std::vector<double> fine;
    fine.reserve(grid.size() * factor);

    for (int i = 0; i < grid.size() - 1; ++i) {
        for (int j = 0; j < factor; ++j) {
            fine.push_back(grid[i] + (grid[i + 1] - grid[i]) * j / factor);
        }
    }
    fine.push_back(grid.back());

    return fine;
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

void writeHomogenousHistoryTo(std::ostream &stream, const HomogenousHistory &history) {
    stream << history.a.getKnots().size() - 1 << "\n";

    for (size_t i = 0; i < history.a.getKnots().size() - 1; ++i) {
        stream << history.a.getKnots()[i] << " " << history.a.getPolynomials()[i] << " "
               << history.Xe.getPolynomials()[i] << "\n";
    }
    stream << history.a.getKnots().back();
}

HomogenousHistory readHomogenousHistoryFrom(std::istream &stream) {
    size_t size;
    stream >> size;

    std::vector<double> eta(size + 1);
    std::vector<CubicPolynomial> aPoly(size);
    std::vector<CubicPolynomial> XePoly(size);

    for (size_t i = 0; i < size; ++i) {
        stream >> eta[i] >> aPoly[i] >> XePoly[i];
    }
    stream >> eta[size];

    return HomogenousHistory{Spline(eta, aPoly), Spline(eta, XePoly)};
}

void getTransferFunctions(const HomogenousHistory &history) {
    CosmoParameters cosmo;
    PhysParameters phys;
    OdeSolver<std::complex<double>, double> solver;

    double k = 0.01;
    solver.addCoefficient("k", [k](double time) { return k; });
    solver.addCoefficient("i", [](double time) { return std::complex<double>(0, 1); });

    solver.addCoefficient("omega_gamma", [omegaGamma = cosmo.omegaPhoton](double time) { return omegaGamma; });
    solver.addCoefficient("omega_nu", [omegaNu = cosmo.omegaNeutrino](double time) { return omegaNu; });
    solver.addCoefficient("omega_b", [omegaB = cosmo.omegaBaryon](double time) { return omegaB; });
    solver.addCoefficient("omega_c", [omegaC = cosmo.omegaCold](double time) { return omegaC; });

    solver.addCoefficient("a", [&history](double time) { return history.a.evaluate(time); });
    solver.addCoefficient("da", [cosmo, &history](double time) {
        double a = history.a.evaluate(time);
        return std::sqrt(cosmo.omegaLambda * std::pow(a, 4) + (cosmo.omegaBaryon + cosmo.omegaCold) * a +
                         (cosmo.omegaPhoton + cosmo.omegaNeutrino));
    });
    solver.addCoefficient("dtau", [phys, cosmo, &history](double time) {
        double a = history.a.evaluate(time);
        double Xe = history.Xe.evaluate(time);
        return -phys.critInEv * phys.thomsonInEv * phys.evInCosmoFrequency * cosmo.omegaBaryon * Xe / (a * a);
    });

    double aStart = history.a.evaluate(0);
    double etaStart = 2 * (std::sqrt(
            (cosmo.omegaPhoton + cosmo.omegaNeutrino) + (cosmo.omegaBaryon + cosmo.omegaCold) * aStart) -
                           std::sqrt(cosmo.omegaPhoton + cosmo.omegaNeutrino)) / (cosmo.omegaBaryon + cosmo.omegaCold);
    solver.addCoefficient("eta", [etaStart](double time) {
        return etaStart + time;
    });

    int lExpansion = 7;
    solver.addVariable("Theta_0");
    solver.addEquation("Theta_0' = -k Theta_1 - dPhi");
    solver.addVariable("N_0");
    solver.addEquation("N_0' = -k N_1 - dPhi");

    solver.addVariable("Theta_1");
    solver.addEquation("Theta_1' = k (Theta_0 - 2 Theta_2 + Psi) / 3 + dtau (Theta_1 - i u_b / 3)");
    solver.addVariable("N_1");
    solver.addEquation("N_1' = k (N_0 - 2 N_2 + Psi) / 3");

    for (int l = 2; l < lExpansion; ++l) {
        auto index = std::to_string(l);

        solver.addVariable("Theta_" + index);
        {
            std::ostringstream stream;
            stream << "Theta_" << index << "' = dtau Theta_" << index << " + ";
            stream << "k (" << index << " Theta_" << std::to_string(l - 1) << " - ";
            stream << std::to_string(l + 1) << " Theta_" << std::to_string(l + 1) << ") / "
                   << std::to_string(2 * l + 1);
            solver.addEquation(stream.str());
        }

        solver.addVariable("N_" + index);
        {
            std::ostringstream stream;
            stream << "N_" << index << "' = ";
            stream << "k (" << index << " N_" << std::to_string(l - 1) << " - ";
            stream << std::to_string(l + 1) << " N_" << std::to_string(l + 1) << ") / " << std::to_string(2 * l + 1);
            solver.addEquation(stream.str());
        }
    }

    solver.addVariable("Theta_" + std::to_string(lExpansion));
    {
        std::ostringstream stream;
        stream << "Theta_" << std::to_string(lExpansion) << "' = k Theta_" + std::to_string(lExpansion - 1) << " + ";
        stream << "Theta_" << std::to_string(lExpansion) << " (dtau - " << std::to_string(lExpansion + 1) << " / eta)";
        solver.addEquation(stream.str());
    }

    solver.addVariable("N_" + std::to_string(lExpansion));
    {
        std::ostringstream stream;
        stream << "N_" << std::to_string(lExpansion) << "' = k N_" + std::to_string(lExpansion - 1) << " - ";
        stream << std::to_string(lExpansion + 1) << " N_" << std::to_string(lExpansion) << " / eta";
        solver.addEquation(stream.str());
    }

    solver.addVariable("delta_c");
    solver.addEquation("delta_c' = -i k u_c - 3 dPhi");
    solver.addVariable("u_c");
    solver.addEquation("u_c' = -u_c da / a - i k Psi");

    solver.addVariable("delta_b");
    solver.addEquation("delta_b' = -i k u_b - 3 dPhi");
    solver.addVariable("u_b");
    solver.addEquation("u_b' = -u_b da / a - i k Psi + 4 dtau (i Theta_1 + u_b / 3) omega_gamma / (a omega_b)");

    solver.addVariable("Phi");
    solver.addEquation("Phi' = dPhi");
    solver.addVariable("Psi");
    solver.addEquation("Psi = -12 (omega_gamma Theta_2 + omega_nu N_2) / (a * a * k * k) - Phi");
    solver.addVariable("dPhi");
    solver.addEquation(
            "dPhi = Psi da / a + a / (3 da) * (-k * k * Phi + 1.5 (omega_c delta_c + omega_b delta_b) / a + 6 (omega_gamma Theta_0 + omega_nu N_0) / (a * a))");

    solver.setGrid(history.a.getKnots());
    solver.compile<walker::BackwardEulerWalker<std::complex<double>, double>>();

    std::vector<std::complex<double>> state = {1.0 / 3, 1.0 / 3};
    state.insert(state.end(), 2 * lExpansion, 0);
    state.insert(state.end(), {1, 0, 1, 0, 2.0 / 3, -2.0 / 3, 0});

    auto solution = solver.solve(state);
    for (int j = 0; j < solution.gridLength; j += 50) {
        for (int i = 0; i < solution.variableNumber; ++i) {
            std::cout << solver.getVariables()[i] << " " << solution(j, i) << " ";
        }
        std::cout << std::endl;
    }
}
