#pragma once

#include <vector>
#include <map>
#include <unordered_set>
#include <string>
#include <functional>
#include <Eigen/Dense>

namespace matrix {
    template<class Value, class Time>
    class SystemBuilder;
}

/*
 * Solver usage example for the equation x'(t) = x(t); 0<=t<=1
 *
 * auto solver = OdeSolver<float, float>();
 * solver.addVariable("x");
 * solver.addEquation("x' = x");
 * solver.setGrid(linSpace<float>(0, 1, 10));
 * solver.solve({1});
 *
 * ODE syntax requirements:
 * 1. All ODEs are linear
 * 2. All variables are space delimited (x' = 5 x)
 * 3. Left part contains only time derivative of one variable (x' = ...)
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

namespace parser {
    class Trie {
    private:
        struct Node {
            std::map<char, Node *> children;

            ~Node() {
                for (const auto &[ch, child]: children) {
                    delete child;
                }
            }
        };

    public:
        ~Trie() {
            delete root_;
        }

        void addString(const std::string &string) {
            Node *node = root_;

            for (char c: string) {
                if (!node->children.count(c)) {
                    node->children[c] = new Node;
                }
                node = node->children[c];
            }
        }

        std::string consumeStream(std::istream &stream) const {
            Node *node = root_;
            std::string string;

            while (stream) {
                char c = stream.peek();

                if (node->children.count(c)) {
                    string.push_back(c);
                    node = node->children[c];
                    stream >> c;
                } else {
                    break;
                }
            }

            return string;
        }

    private:
        Node *root_ = new Node;
    };

    enum Lexeme {
        Left, Right, Plus, Minus, Derivative, Number, End
    };

    class EquationParser {
    public:
        EquationParser(const Trie &trie, std::istream &stream) : trie_(trie), stream_(stream) {
        }

        std::pair<Lexeme, std::string> next() {
            while (stream_) {
                char c = stream_.peek();
                if (c == EOF) {
                    break;
                }
                if (specialSymbols_.count(c)) {
                    stream_.get(c);
                    if (c == '+') {
                        return {Lexeme::Plus, ""};
                    } else if (c == '-') {
                        return {Lexeme::Minus, ""};
                    } else if (c == '\'') {
                        if (!left_) {
                            throw std::logic_error("can't have derivative in the right part of an equation ");
                        }
                        return {Lexeme::Derivative, ""};
                    } else if (c == '=') {
                        if (!left_) {
                            throw std::logic_error("can't have more than one '=' in an equation");
                        }
                        left_ = false;
                    }
                } else if (c >= '0' && c <= '9') {
                    return {Lexeme::Number, consumeNumber()};
                } else {
                    auto word = trie_.consumeStream(stream_);
                    if (left_) {
                        return {Lexeme::Left, word};
                    } else {
                        return {Lexeme::Right, word};
                    }
                }
            }
            return {Lexeme::End, ""};
        }

    private:
        std::string consumeNumber() {
            std::string number;

            while (stream_) {
                char c = stream_.peek();

                if (c >= '0' && c <= '9') {
                    stream_ >> c;
                    number.push_back(c);
                } else {
                    break;
                }
            }

            return number;
        }

    private:
        std::istream &stream_;
        const Trie &trie_;
        inline static const std::unordered_set<char> specialSymbols_{'+', '-', '*', '=', '\'', ' '};
        bool left_ = true;
    };
} // namespace parser

namespace matrix {
    // Stores matrix and source (A, b) of linear ode system x' = Ax + b
    template<class Value, class Time>
    class System {
    private:
        struct Matrix {
            using Mat = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;

            Mat matrix;
            std::vector<std::pair<int, int>> entries;
            std::vector<std::function<Value(Time, const Value *)>> coefficients;

            void update(Time time, const Value *state) {
                for (size_t i = 0; i < entries.size(); ++i) {
                    matrix(entries[i].first, entries[i].second) = coefficients[i](time, state);
                }
            }
        };

    public:
        using Mat = typename Matrix::Mat;

        const Mat &getMatrix(Time time, const Value *state) {
            matrix_.update(time, state);
            return matrix_.matrix;
        }

        const Mat &getSource(Time time, const Value *state) {
            source_.update(time, state);
            return source_.matrix;
        }

    private:
        friend class SystemBuilder<Value, Time>;

        Matrix matrix_;
        Matrix source_;
    };

    template<class Value, class Time>
    class SystemBuilder {
    public:
        System<Value, Time> build(const OdeSolver<Value, Time> &solver) {
            size_t varnum = solver.variables_.size();
            system.matrix_.matrix.setZero(varnum, varnum);
            system.source_.matrix.setZero(varnum, 1);

            variableToIndex = createVectorToIndex(solver.variables_);
            coefficientToIndex = createVectorToIndex(solver.coefficients_);
            functionToIndex = createVectorToIndex(solver.functions_);

            parser::Trie trie;
            for (const auto &name: solver.namespace_) {
                trie.addString(name);
            }

            for (const auto &equation: solver.equations_) {
                parseEquation(equation, trie, solver);
            }

            return system;
        }

    private:
        template<class T>
        static std::unordered_map<T, size_t> createVectorToIndex(const std::vector<T> &vector) {
            std::unordered_map<T, size_t> map;
            for (size_t i = 0; i < vector.size(); ++i) {
                map[vector[i]] = i;
            }
            return map;
        }

        void parseEquation(const std::string &equation, const parser::Trie &trie,
                           const OdeSolver<Value, Time> &solver) {
            using namespace parser;
            std::stringstream stream(equation);
            EquationParser parser(trie, stream);

            auto [lexeme, variable] = parser.next();
            size_t row = variableToIndex.at(variable);
            if (lexeme != Lexeme::Left && parser.next().first != Lexeme::Derivative) {
                throw std::logic_error("an equation has to start from variable's derivative");
            }

            int sign = 1;
            std::pair<Lexeme, std::string> word;
            std::vector<std::pair<Lexeme, std::string>> batch;

            do {
                word = parser.next();
                if (word.first == Lexeme::Right || word.first == Lexeme::Number) {
                    batch.push_back(word);
                } else {
                    if (!batch.empty()) {
                        insertBatch(sign, row, batch, solver);
                        batch.clear();
                    }
                    if (word.first == Lexeme::Minus) {
                        sign = -1;
                    } else if (word.first == Lexeme::Plus) {
                        sign = 1;
                    } else if (word.first == Lexeme::Left) {
                        throw std::logic_error("left part of an equation has to contain only one variable");
                    }
                }
            } while (word.first != Lexeme::End);
        }

        void insertBatch(int sign, size_t row,
                         const std::vector<std::pair<parser::Lexeme, std::string>> &batch,
                         const OdeSolver<Value, Time> &solver) {
            using namespace parser;

            bool containsVariable = false;
            size_t column;
            std::vector<typename OdeSolver<Value, Time>::CoeffFunc> coefficients;
            std::vector<std::pair<std::vector<int>, typename OdeSolver<Value, Time>::ExplicitFunc>> functions;
            Value constant = sign;

            for (const auto &[lex, name]: batch) {
                if (lex == Lexeme::Number) {
                    constant *= std::stof(name);
                } else {
                    if (variableToIndex.count(name)) {
                        if (containsVariable) {
                            throw std::logic_error("only linear odes are considered");
                        }
                        containsVariable = true;
                        column = variableToIndex.at(name);
                    } else if (coefficientToIndex.count(name)) {
                        coefficients.push_back(solver.coefficientValues_[coefficientToIndex.at(name)]);
                    } else if (functionToIndex.count(name)) {
                        auto [function, variableNames] = solver.functionValues_[functionToIndex.at(name)];
                        std::vector<int> indices;
                        indices.reserve(variableNames.size());
                        for (const auto &variable: variableNames) {
                            indices.push_back(variableToIndex.at(variable));
                        }
                        functions.emplace_back(indices, function);
                    } else {
                        throw std::logic_error("unknown variable/coefficient name");
                    }
                }
            }

            if (containsVariable) {
                system.matrix_.entries.emplace_back(row, column);
                system.matrix_.coefficients.emplace_back(
                        [constant, coeffs = std::move(coefficients), funcs = std::move(functions)](Time time,
                                                                                                   const Value *value) {
                            return SystemBuilder::calculateEntry(constant, coeffs, funcs, time, value);
                        });
            } else {
                system.source_.entries.emplace_back(row, 0);
                system.source_.coefficients.emplace_back(
                        [constant, coeffs = std::move(coefficients), funcs = std::move(functions)](Time time,
                                                                                                   const Value *value) {
                            return SystemBuilder::calculateEntry(constant, coeffs, funcs, time, value);
                        });
            }
        }

        static Value calculateEntry(Value multiplier,
                                    const std::vector<typename OdeSolver<Value, Time>::CoeffFunc> &coefficients,
                                    const std::vector<std::pair<std::vector<int>, typename OdeSolver<Value, Time>::ExplicitFunc>> functions,
                                    Time time, const Value *value) {
            Value entry = multiplier;
            for (const auto &coefficient: coefficients) {
                entry *= coefficient(time);
            }
            for (const auto &[indices, function]: functions) {
                std::vector<Value> values;
                values.reserve(indices.size());
                for (int index: indices) {
                    values.push_back(value[index]);
                }
                entry *= function(time, std::move(values));
            }
            return entry;
        }

    private:
        System<Value, Time> system;
        std::unordered_map<std::string, size_t> variableToIndex;
        std::unordered_map<std::string, size_t> coefficientToIndex;
        std::unordered_map<std::string, size_t> functionToIndex;
    };
}

namespace walker {
    template<class Value, class Time>
    class Walker {
    public:
        Walker(matrix::System<Value, Time> system, const std::vector<Time> &grid) : system_(system),
                                                                                    grid_(grid) {
        }

    protected:
        matrix::System<Value, Time> system_;
        const std::vector<Time> &grid_;
    };

    /*
     * Backward Euler discretizes system as:
     * (x^{n+1} - x^n) / dt = A^{n+1}x^{n+1} + b^{n+1} <=> (1 - dt A^{n+1}) x^{n+1} = x^n + dt b^{n+1}
     */
    template<class Value, class Time>
    class BackwardEulerWalker : public Walker<Value, Time> {
    public:
        using Walker<Value, Time>::Walker;

        void next(typename OdeSolver<Value, Time>::Solution &solution, size_t step) {
            using Mat = typename matrix::System<Value, Time>::Mat;

            Time time = this->grid_[step + 1];
            Time dt = time - this->grid_[step];
            Mat mat = Mat::Identity(solution.variableNumber, solution.variableNumber) -
                      dt * this->system_.getMatrix(time, solution.getState(step));
            Mat source = Eigen::Map<const Mat>(solution.getState(step), solution.variableNumber, 1) +
                         dt * this->system_.getSource(time, solution.getState(step));
            Mat nextMatState = mat.colPivHouseholderQr().solve(source);

            for (size_t j = 0; j < solution.variableNumber; ++j) {
                solution(step + 1, j) = nextMatState(j, 0);
            }
        }
    };

    // Standard 4-step Runge-Kutta method
    template<class Value, class Time>
    class RungeKuttaWalker : public Walker<Value, Time> {
    public:
        using Walker<Value, Time>::Walker;

        void next(typename OdeSolver<Value, Time>::Solution &solution, size_t step) {
            Time time = this->grid_[step];
            Time dt = this->grid_[step + 1] - time;

            Mat matState = Eigen::Map<const Mat>(solution.getState(step), solution.variableNumber, 1);
            Mat k1 = getDerivative(matState, time);
            Mat k2 = getDerivative(matState + k1 * dt / 2, time + dt / 2);
            Mat k3 = getDerivative(matState + k2 * dt / 2, time + dt / 2);
            Mat k4 = getDerivative(matState + k3 * dt, time + dt);

            Mat nextMatState = matState + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            for (size_t j = 0; j < solution.variableNumber; ++j) {
                solution(step + 1, j) = nextMatState(j, 0);
            }
        }

    private:
        using Mat = typename matrix::System<Value, Time>::Mat;

        auto getDerivative(const Mat &state, Time time) {
            return this->system_.getMatrix(time, state.data()) * state + this->system_.getSource(time, state.data());
        }
    };
} // namespace walker
