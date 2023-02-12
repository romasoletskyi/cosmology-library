#include <optional>
#include <iostream>

#include <cosmology/odesolver/odesolver.h>
#include <cosmology/utility/math.h>
#include <cosmology/utility/utility.h>
#include <cosmology/perturbations.h>
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

std::vector<double> makeFineGrid(const std::vector<double> &grid, int factor) {
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
        double cross_section = 9.78 * phys.evInCosmoFrequency * std::pow(phys.alpha, 2) / std::pow(phys.me, 2) *
                               std::sqrt(phys.delta / T) *
                               std::log(phys.delta / T);

        double da = std::sqrt(cosmo.omegaLambda * std::pow(a, 4) + (cosmo.omegaBaryon + cosmo.omegaCold) * a +
                              (cosmo.omegaPhoton + cosmo.omegaNeutrino));
        double wavelength_in_inverse_ev = 2 * pi / (0.75 * phys.delta); // 2s -> 1s transition
        double lambda_alpha = phys.cosmoFrequencyInHz * 8 * pi * a * da /
                              ((1 - Xe) * phys.critInEv * cosmo.omegaBaryon *
                               std::pow(wavelength_in_inverse_ev, 3)); // cosmological redshifting of Lyman photons rate

        double lambda_beta = phys.cosmoFrequencyInHz * cross_section * std::pow(phys.me * T / 2 / pi, 1.5) *
                             std::exp(-phys.delta / 4 / T); // ionization rate
        double lambda_two_photons = 8.227; // 2s -> 1s 2-photon transition rate
        double cr = (lambda_two_photons + lambda_alpha) / (lambda_alpha + lambda_two_photons + lambda_beta);

        return a * cross_section * cr *
               (std::pow(phys.me * T / 2 / pi, 1.5) * std::exp(-phys.delta / T) * (1 - Xe) - nb * Xe * Xe);
    }, {"a", "Xe"});

    solver.addEquation("a' = da");
    solver.addEquation("Xe' = dXe");
    solver.compile<walker::RungeKuttaWalker<double, double>>();

    double etaDelta = 2e-4;
    double aStart = 6.6e-4;
    double etaStart = 2 * (std::sqrt(aStart * (cosmo.omegaBaryon + cosmo.omegaCold) + (cosmo.omegaPhoton + cosmo.omegaNeutrino)) - std::sqrt(cosmo.omegaPhoton + cosmo.omegaNeutrino)) / (cosmo.omegaBaryon + cosmo.omegaCold);
    double TStart = cosmo.cmb / aStart / phys.evInKelvin;
    double nbStart = phys.critInEv * cosmo.omegaBaryon / std::pow(aStart, 3);
    double XeStart = 1 - nbStart * std::exp(phys.delta / TStart) * std::pow(phys.me * TStart / 2 / pi, -1.5);

    double etaCalculatedStart = 1e-5 / std::sqrt(cosmo.omegaPhoton + cosmo.omegaNeutrino);
    auto etaCalculated = linSpace(etaCalculatedStart, etaStart, pointsNumber / 10);
    std::vector<std::pair<double, double>> aCalculated, XeCalculated;

    for (size_t i = 0; i < etaCalculated.size() - 1; ++i) {
        double eta = etaCalculated[i];
        double a = std::sqrt(cosmo.omegaPhoton + cosmo.omegaNeutrino) * eta + (cosmo.omegaBaryon + cosmo.omegaCold) * std::pow(eta, 2) / 4;
        double nb = phys.critInEv * cosmo.omegaBaryon / std::pow(a, 3);
        double T = cosmo.cmb / a / phys.evInKelvin;

        aCalculated.emplace_back(eta, a);
        XeCalculated.emplace_back(eta, 1 - nb * std::pow(phys.me * T / 2 / pi, -1.5) * std::exp(phys.delta / T));
    }

    std::vector<double> etaVec, aVec, XeVec;
    while (true) {
        auto [etaStep, aStep, XeStep, finished] = progressEvolution(solver, etaStart, aStart, XeStart, etaDelta);

        etaVec.insert(etaVec.end(), etaStep.begin() + 1, etaStep.end());
        aVec.insert(aVec.end(), aStep.begin() + 1, aStep.end());
        XeVec.insert(XeVec.end(), XeStep.begin() + 1, XeStep.end());

        etaStart = etaVec.back();
        aStart = aVec.back();
        XeStart = XeVec.back();

        if (finished) {
            break;
        }
    }

    auto aSampled = samplePoints(etaVec, aVec, pointsNumber);
    auto XeSampled = samplePoints(etaVec, XeVec, pointsNumber);

    aCalculated.insert(aCalculated.end(), aSampled.begin(), aSampled.end());
    XeCalculated.insert(XeCalculated.end(), XeSampled.begin(), XeSampled.end());

    auto aSpline = SplineBuilder().buildFromPoints(aCalculated);
    auto XeSpline = SplineBuilder().buildFromPoints(XeCalculated);

    OdeSolver<double, double> opticDepthSolver;
    opticDepthSolver.addVariable("tau");
    opticDepthSolver.addEquation("tau' = dtau");
    opticDepthSolver.addCoefficient("dtau", [phys, cosmo, &aSpline, &XeSpline](double time) {
        double a = aSpline.evaluate(time);
        double Xe = XeSpline.evaluate(time);
        return -phys.critInEv * phys.thomsonInEv * phys.evInCosmoFrequency * cosmo.omegaBaryon * Xe / (a * a);
    });

    auto grid = aSpline.getKnots();
    std::reverse(grid.begin(), grid.end());
    opticDepthSolver.setGrid(grid);
    opticDepthSolver.compile<walker::RungeKuttaWalker<double, double>>();
    auto tau = opticDepthSolver.solve({0});

    std::vector<std::pair<double, double>> points;
    points.reserve(grid.size());
    for (int i = 0; i < grid.size(); ++i) {
        points.emplace_back(grid[i], tau(i, 0));
    }

    auto tauSpline = SplineBuilder().buildFromPoints(points);
    return HomogenousHistory{std::move(aSpline), std::move(XeSpline), std::move(tauSpline)};
}

void writeHomogenousHistoryTo(std::ostream &stream, const HomogenousHistory &history) {
    stream << history.a.getKnots().size() - 1 << "\n";

    for (size_t i = 0; i < history.a.getKnots().size() - 1; ++i) {
        stream << history.a.getKnots()[i] << " " << history.a.getPolynomials()[i] << " "
               << history.Xe.getPolynomials()[i] << " " << history.tau.getPolynomials()[i] << "\n";
    }
    stream << history.a.getKnots().back();
}

HomogenousHistory readHomogenousHistoryFrom(std::istream &stream) {
    size_t size;
    stream >> size;

    std::vector<double> eta(size + 1);
    std::vector<CubicPolynomial> aPoly(size);
    std::vector<CubicPolynomial> XePoly(size);
    std::vector<CubicPolynomial> tauPoly(size);

    for (size_t i = 0; i < size; ++i) {
        stream >> eta[i] >> aPoly[i] >> XePoly[i] >> tauPoly[i];
    }
    stream >> eta[size];

    return HomogenousHistory{Spline(eta, aPoly), Spline(eta, XePoly), Spline(eta, tauPoly)};
}

std::vector<std::complex<double>> getGeneratingFunction(const HomogenousHistory &history, double wavenumber) {
    CosmoParameters cosmo;
    PhysParameters phys;
    OdeSolver<std::complex<double>, double> solver;

    solver.addCoefficient("k", [wavenumber](double time) { return wavenumber; });
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

    solver.addCoefficient("eta", [](double time) {
        return time;
    });

    int lExpansion = 8;
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
    solver.addEquation("u_b' = -u_b da / a - i k Psi + 4 dtau omega_gamma (i Theta_1 + u_b / 3) / (a omega_b)");

    solver.addVariable("Phi");
    solver.addEquation("Phi' = dPhi");
    solver.addVariable("Psi");
    solver.addEquation("Psi = -12 (omega_gamma Theta_2 + omega_nu N_2) / (a * a * k * k) - Phi");
    solver.addVariable("dPhi");
    solver.addEquation(
            "dPhi = Psi da / a + a / (3 da) * (-k * k * Phi + 1.5 (omega_c delta_c + omega_b delta_b) / a + 6 (omega_gamma Theta_0 + omega_nu N_0) / (a * a))");

    solver.addVariable("S1");
    solver.addEquation("S1 = Psi - i dtau u_b / k");
    solver.addVariable("S2");
    solver.addEquation("S2 = dtau Theta_0 + dPhi");

    auto grid = history.a.getKnots();
    solver.setGrid(grid);
    solver.compile<walker::GaussLegendreWalker<std::complex<double>, double>>();

    std::vector<std::complex<double>> state = {1.0 / 3, 1.0 / 3};
    state.insert(state.end(), 2 * lExpansion, 0);
    state.insert(state.end(), {1, 0, 1, 0, 2.0 / 3, -2.0 / 3, 0, 0, 0});

    auto solution = solver.solve(state);
    std::complex<double> S1[solution.gridLength];
    std::complex<double> S2[solution.gridLength];

    std::cout << wavenumber * 4.67 << " " << solution(270, 0) << " " << solution(270, 2 * lExpansion + 6) << " "
              << solution(270, 0) + solution(270, 2 * lExpansion + 6) << std::endl;

    for (int i = 0; i < solution.gridLength; ++i) {
        double visibility = std::exp(-history.tau.evaluate(grid[i]));
        S1[i] = visibility * solution(i, solution.variableNumber - 2);
        S2[i] = visibility * solution(i, solution.variableNumber - 1);
    }

    std::complex<double> S1derivative[solution.gridLength];
    differentiate(grid.data(), S1, S1derivative, solution.gridLength);

    std::vector<std::complex<double>> S;
    S.reserve(solution.gridLength);

    for (int i = 0; i < solution.gridLength; ++i) {
        S.push_back(S1derivative[i] - S2[i]);
    }

    return S;
}
