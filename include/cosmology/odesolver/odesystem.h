#pragma once

#include <Eigen/Dense>
#include <cosmology/odesolver/odeparser.h>
#include <cosmology/utility/utility.h>

template<class Value, class Time>
class OdeSolver;

namespace matrix {
    template<class Value, class Time>
    class SystemBuilder;

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
                for (auto [i, j]: entries) {
                    matrix(i, j) = 0;
                }

                for (size_t i = 0; i < entries.size(); ++i) {
                    matrix(entries[i].first, entries[i].second) += coefficients[i](time, state);
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
        void parseEquation(const std::string &equation, const parser::Trie &trie,
                           const OdeSolver<Value, Time> &solver) {
            using namespace parser;
            std::stringstream stream(equation);
            EquationLexer parser(trie, stream);

            auto [lexeme, variable] = parser.next();
            size_t row = variableToIndex.at(variable);
            if (parser.next().first != Lexeme::Derivative) {
                throw std::logic_error("an equation has to start from variable's derivative");
            }

            std::pair<Lexeme, std::string> word;
            std::vector<std::pair<Lexeme, std::string>> batch;
            std::vector<int> batchSizes = {0};
            std::vector<int> batchSigns = {1, 1};

            do {
                word = parser.next();
                if (word.first == Lexeme::Right || word.first == Lexeme::Number) {
                    batch.push_back(word);
                    ++batchSizes.back();
                } else if (word.first == Lexeme::LeftBracket) {
                    batchSizes.push_back(0);
                    batchSigns.push_back(batchSigns.back());
                } else {
                    if (batchSizes.back()) {
                        insertBatch(batchSigns.back(), row, batch, solver);
                        batch.resize(batch.size() - batchSizes.back());
                        batchSizes.back() = 0;
                        if (word.first == Lexeme::RightBracket) {
                            batchSigns.pop_back();
                            batchSizes.pop_back();
                            batch.resize(batch.size() - batchSizes.back());
                            batchSizes.back() = 0;
                        }
                    }
                    if (word.first == Lexeme::Minus) {
                        batchSigns.back() = -batchSigns[batchSigns.size() - 2];
                    } else if (word.first == Lexeme::Plus) {
                        batchSigns.back() = batchSigns[batchSigns.size() - 2];
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
        std::unordered_map<std::string, int> variableToIndex;
        std::unordered_map<std::string, int> coefficientToIndex;
        std::unordered_map<std::string, int> functionToIndex;
    };
}