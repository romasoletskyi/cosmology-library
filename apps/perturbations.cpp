#include <iostream>
#include <fstream>
#include <cosmology/perturbations.h>

int main() {
    auto history = getHomogenousHistory(1e-3, 4.7, 1e4);
    std::ofstream file("../data/HomogenousHistory");

    for (int i = 0; i < history.eta.size(); i += 100) {
        file << history.eta[i] << " " << history.a[i] << "\n";
    }
    return 0;
}