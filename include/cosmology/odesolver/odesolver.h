#pragma once

#include <vector>
#include <functional>
#include <memory>
#include <Eigen/Dense>

#include <cosmology/odesolver/odesystem.h>
#include <cosmology/odesolver/odewalker.h>

template<class Value>
struct Solution {
    Solution(int gridLength_, int variableNumber_) : data(new Value[gridLength_ * variableNumber_]),
                                                     gridLength(gridLength_), variableNumber(variableNumber_) {
    }

    Solution() : data(nullptr), gridLength(0), variableNumber(0) {
    }

    ~Solution() {
        delete[] data;
    }

    Solution(const Solution &solution) = delete;

    Solution &operator=(const Solution &solution) = delete;

    Solution(Solution &&solution) : data(std::exchange(solution.data, nullptr)), gridLength(solution.gridLength),
                                    variableNumber(solution.variableNumber) {
    }

    Solution &operator=(Solution &&solution) noexcept {
        std::swap(data, solution.data);
        gridLength = solution.gridLength;
        variableNumber = solution.variableNumber;
        return *this;
    }

    const Value *getState(int i) const {
        return data + i * variableNumber;
    }

    Value *getState(int i) {
        return data + i * variableNumber;
    }

    Value operator()(int i, int j) const {
        return data[i * variableNumber + j];
    }

    Value &operator()(int i, int j) {
        return data[i * variableNumber + j];
    }

    Value *data;
    int gridLength, variableNumber;
};

/*
 * ODE syntax:
 * 1. Different variables/numbers have to be divided by math operators (-, +, *, /) or by space.
 *    Space between two variables is interpreted as multiplication
 *    Example: x' = 2 x + 7
 * 2. Left part contains time derivative of one variable (x' = ...) or definition of variable itself (y = 2 * x)
 * 3. addCoefficient adds function which depends on time only
 * 4. addExplicitFunction adds function on time and state at time t_n, which is going to be used to compute state
 *    at t_{n+1}. Full support is provided to functions which depend only on dynamical variables
 *    (meaning variables which have x' = ... like equation). Otherwise, order of a numerical scheme may be lower.
 *
 * Solver usage example for the equation x'(t) = x(t); 0<=t<=1
 *
 * auto solver = OdeSolver<float, float>();
 * solver.addVariable("x");
 * solver.addEquation("x' = x");
 * solver.setGrid(linSpace<float>(0, 1, 10));
 * solver.compile<walker::RungeKuttaWalker<double, double>>();
 * auto solution = solver.solve({1});
 */
template<class Value, class Time>
class OdeSolver {
public:
    using State = std::vector<Value>;
    using CoeffFunc = std::function<Value(Time)>;
    using ExplicitFunc = std::function<Value(Time, std::vector<Value>)>;

public:
    void addVariable(std::string variable) {
        checkName(variable);
        variables_.emplace_back(variable);
    }

    void addCoefficient(std::string coefficient, CoeffFunc value) {
        checkName(coefficient);
        coefficients_.emplace_back(coefficient);
        coefficientValues_.emplace_back(value);
    }

    void addExplicitFunction(std::string function, ExplicitFunc value, std::vector<std::string> variables) {
        checkName(function);
        functions_.emplace_back(function);
        functionValues_.emplace_back(value, variables);
    }

    void addEquation(std::string equation) {
        equations_.emplace_back(equation);
    }

    void setGrid(std::vector<Time> grid) {
        grid_ = std::move(grid);
    }

    template<class Walker>
    void compile() {
        walker_ = std::unique_ptr<walker::IWalker<Value, Time>>(
                new Walker(formula::SystemBuilder<Value, Time>().build(*this)));
    }

    Solution<Value> solve(const State &initialState) {
        if (variables_.size() != equations_.size()) {
            throw std::logic_error("variables and equations numbers are not equal");
        }

        Solution<Value> solution(grid_.size(), variables_.size());
        for (int i = 0; i < solution.variableNumber; ++i) {
            solution(0, i) = initialState[i];
        }

        for (size_t i = 0; i < grid_.size() - 1; ++i) {
            walker_->next(solution, grid_, i);
        }

        return solution;
    }

    const std::vector<std::string> &getVariables() const {
        return variables_;
    }

private:
    void checkName(const std::string &name) {
        if (name[0] >= '0' && name[0] <= '9') {
            throw std::logic_error("name can't start from a number");
        }
        if (namespace_.count(name)) {
            throw std::logic_error("this name was already added");
        }
        namespace_.insert(name);
    }

private:
    friend class formula::SystemBuilder<Value, Time>;

    std::unordered_set<std::string> namespace_;
    std::vector<std::string> variables_;
    std::vector<std::string> coefficients_;
    std::vector<std::string> functions_;
    std::vector<std::string> equations_;

    std::vector<CoeffFunc> coefficientValues_;
    std::vector<std::pair<ExplicitFunc, std::vector<std::string>>> functionValues_;

    std::vector<Time> grid_;
    std::unique_ptr<walker::IWalker<Value, Time>> walker_;
};
