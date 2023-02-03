#include <iostream>
#include <fstream>
#include <cosmology/perturbations.h>

int main() {
    auto history = getHomogenousHistory(1000);
    std::ofstream file(std::string(PROJECT_PATH) + "/data/HomogenousHistory");

    for (double eta: history.a.getKnots()) {
        file << eta << " " << history.a.evaluate(eta) << " " << history.Xe.evaluate(eta) << "\n";
    }
    return 0;
}