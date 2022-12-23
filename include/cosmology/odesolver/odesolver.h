#pragma once

#include <vector>
#include <functional>
#include <Eigen/Dense>

#include <cosmology/odesolver/odeparser.h>
#include <cosmology/odesolver/odesystem.h>
#include <cosmology/odesolver/odewalker.h>

#include <cosmology/utility/graph.h>

namespace matrix {
    template<class Value, class Time>
    class SystemBuilder;
}

/*
 * ODE syntax requirements:
 * 1. All variables are space delimited (x' = 5 x)
 * 2. Left part contains only time derivative of one variable (x' = ...) or definition of variable itself (y = 2 * x)
 * 3. addCoefficient adds function which depends on time
 * 4. addExplicitFunction adds function on time and state at time t_n,
 *    which is going to be used to compute state at t_{n+1}
 * 5. Coefficients and explicit functions MUST be written before the variable (x' = k x)
 *
 * Solver usage example for the equation x'(t) = x(t); 0<=t<=1
 *
 * auto solver = OdeSolver<float, float>();
 * solver.addVariable("x");
 * solver.addEquation("x' = x");
 * solver.setGrid(linSpace<float>(0, 1, 10));
 * solver.solve<walker::RungeKuttaWalker<float, float>>({1});
 */
template<class Value, class Time>
class OdeSolver {
public:
    using State = std::vector<Value>;
    using CoeffFunc = std::function<Value(Time)>;
    using ExplicitFunc = std::function<Value(Time, std::vector<Value>)>;

    struct Solution {
        Solution(int gridLength_, int variableNumber_) : data(new Value[gridLength_ * variableNumber_]),
                                                         gridLength(gridLength_), variableNumber(variableNumber_) {
        }

        const Value *getState(int i) const {
            return data + i * variableNumber;
        }

        Value operator()(int i, int j) const {
            return data[i * variableNumber + j];
        }

        Value &operator()(int i, int j) {
            return data[i * variableNumber + j];
        }

        ~Solution() {
            delete[] data;
        }

        Value *data;
        int gridLength, variableNumber;
    };

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
    Solution solve(const State &initialState) {
        if (variables_.size() != equations_.size()) {
            throw std::logic_error("variables and equations numbers are not equal");
        }

        makeSubstitutions();
        auto walker = Walker(matrix::SystemBuilder<Value, Time>().build(*this), grid_);

        Solution solution(grid_.size(), variables_.size());
        for (int j = 0; j < solution.variableNumber; ++j) {
            solution(0, j) = initialState[j];
        }

        for (size_t i = 0; i < grid_.size() - 1; ++i) {
            walker.next(solution, i);
        }

        return solution;
    }

    /*
     * Each variable has a corresponding equation according to its left part.
     * There are two types of variables - dynamic (equation starts as x' = ...) and non-dynamic (x = ...).
     * Algorithm of substitutions is:
     * 1. Assume that there exists a topological order on dependence graph of non-dynamic variables and find it
     * 2. Substitute non-dynamic variables according to this order,
     *    now all these variables can depend only on dynamic ones
     * 3. Substitute non-dynamic variables into dynamic equations
     */
    void makeSubstitutions() {
        using namespace parser;

        Trie trie;
        for (const auto &name: namespace_) {
            trie.addString(name);
        }

        auto variableToIndex = createVectorToIndex(variables_);
        auto [dynamicVariables, nonDynamicVariables, variableToEquation] = getVariablesMapping(trie, variableToIndex);

        // edge goes from variable to its dependency
        auto dependencyGraph = buildVariableDependencyGraph(trie, variableToEquation, variableToIndex);
        auto [nonDynamicDependencyGraph, nonDynamicToIndex] = buildSubgraph(dependencyGraph, nonDynamicVariables);
        auto topologicalOrder = getTopologicalOrder(transposeGraph(nonDynamicDependencyGraph));

        for (int nonDynamicVarIndex: topologicalOrder) {
            int nonDynamicVar = nonDynamicVariables[nonDynamicVarIndex];
            for (int dependency: getVertexNeighbors(dependencyGraph, nonDynamicVar)) {
                if (nonDynamicToIndex.count(dependency)) {
                    substituteVariableTo(variableToEquation, dependency, nonDynamicVar);
                }
            }
        }

        for (int dynamicVar: dynamicVariables) {
            for (int dependency: getVertexNeighbors(dependencyGraph, dynamicVar)) {
                if (nonDynamicToIndex.count(dependency)) {
                    substituteVariableTo(variableToEquation,dependency, dynamicVar);
                }
            }
        }
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

    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
    getVariablesMapping(const parser::Trie &trie, const std::unordered_map<std::string, int> &variableToIndex) const {
        using namespace parser;

        std::vector<int> dynamicVariables;
        std::vector<int> nonDynamicVariables;
        std::vector<int> variableToEquation;
        variableToEquation.resize(variables_.size());

        for (int i = 0; i < equations_.size(); ++i) {
            std::stringstream stream(equations_[i]);
            EquationLexer parser(trie, stream);

            int variableIndex = variableToIndex.at(parser.next().second);
            variableToEquation[variableIndex] = i;

            if (parser.next().first != Lexeme::Derivative) {
                nonDynamicVariables.push_back(variableIndex);
            } else {
                dynamicVariables.push_back(variableIndex);
            }
        }

        return {dynamicVariables, nonDynamicVariables, variableToEquation};
    }

    Graph buildVariableDependencyGraph(const parser::Trie &trie,
                                       const std::vector<int> &variableToEquation,
                                       const std::unordered_map<std::string, int> &variableToIndex) const {
        using namespace parser;

        Graph graph;
        graph.edges.resize(variables_.size());

        for (int i = 0; i < variables_.size(); ++i) {
            std::stringstream stream(equations_[variableToEquation[i]]);
            EquationLexer parser(trie, stream);
            std::unordered_set<int> dependentIndices;

            auto [lexeme, variable] = parser.next();
            do {
                std::tie(lexeme, variable) = parser.next();
                if (lexeme == Lexeme::Right && variableToIndex.count(variable)) {
                    dependentIndices.insert(variableToIndex.at(variable));
                }
            } while (lexeme != Lexeme::End);

            for (int dependency: dependentIndices) {
                graph.edges[i].push_back(dependency);
            }
        }

        return graph;
    }

    void substituteVariableTo(const std::vector<int>& variableToEquation, int from, int to) {
        std::string equationFrom = equations_[variableToEquation[from]];
        std::string variable = variables_[from];
        std::string sub = equationFrom.substr(equationFrom.find('=') + 1);
        sub.erase(sub.begin(), std::find_if(sub.begin(), sub.end(), [](char c) { return !std::isspace(c); }));
        sub = "(" + sub + ")";

        std::string& equationTo = equations_[variableToEquation[to]];
        for (std::string::size_type pos{};
             std::string::npos != (pos = equationTo.find(variable.data(), pos, variable.size()));
             pos = variable.size()) {
            equationTo.replace(pos, variable.size(), sub.data(), sub.size());
        }
    }

private:
    friend class matrix::SystemBuilder<Value, Time>;

    std::unordered_set<std::string> namespace_;
    std::vector<std::string> variables_;
    std::vector<std::string> coefficients_;
    std::vector<CoeffFunc> coefficientValues_;
    std::vector<std::string> functions_;
    std::vector<std::pair<ExplicitFunc, std::vector<std::string>>> functionValues_;
    std::vector<std::string> equations_;
    std::vector<Time> grid_;
};
