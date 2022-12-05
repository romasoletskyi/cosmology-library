#pragma once

#include <vector>
#include <string>
#include <functional>

class OdeSolver {
public:
    void addVariable(std::string variable);
    void addCoefficient(std::pair<std::string, std::function<double(double)>>);
};
