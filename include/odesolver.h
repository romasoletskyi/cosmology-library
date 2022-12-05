#pragma once

#include <vector>
#include <map>
#include <unordered_set>
#include <string>
#include <functional>
#include <Eigen/Dense>

namespace walker {
    template<class Value, class Time>
    class OdeWalkerBuilder;
}

/*
 * Solver usage example for the equation x'(t) = x(t); 0<=t<=1
 *
 * auto solver = OdeSolver<float>();
 * solver.addVariable("x");
 * solver.addEquation("x' = x");
 * solver.setGrid(linSpace<OdeSolver::Value>(0, 1, 10));
 * solver.solve({1});
 *
 * ODE syntax requirements:
 * 1. All variables are space delimited (x' = 5 x)
 * 2. Left part contains only time derivative of one variable (x' = ...)
 */
template<class Value, class Time>
class OdeSolver {
public:
    using State = std::vector<Value>;
    using CoeffFunc = std::function<Value(Time)>;

public:
    void addVariable(std::string variable) {
        if (checkNumericStart(variable)) {
            throw std::logic_error("variable name can't start from a number");
        }
        variables_.emplace_back(std::move(variable));
    }

    void addCoefficient(std::string coefficient, CoeffFunc value) {
        if (checkNumericStart(coefficient)) {
            throw std::logic_error("coefficient name can't start from a number");
        }
        coefficients_.emplace_back(std::move(coefficient));
        coefficientValues_.emplace_back(std::move(value));
    }

    void addEquation(std::string equation) {
        equations_.emplace_back(std::move(equation));
    }

    void setGrid(std::vector<Time> grid) {
        grid_ = std::move(grid);
    }

    std::vector<State> solve(const State &initialState) {
        std::vector<State> states;
        states.reserve(grid_.size());

        State state = initialState;
        states.push_back(state);

        auto walker = walker::OdeWalkerBuilder<Value, Time>().build(*this);
        for (size_t i = 1; i < grid_.size(); ++i) {
            state = walker.next(state, grid_[i]);
            states.push_back(state);
        }

        return states;
    }

    const std::vector<std::string> &getVariables() const {
        return variables_;
    }

private:
    bool checkNumericStart(const std::string &string) {
        return string[0] >= '0' && string[0] <= '9';
    }

private:
    friend class walker::OdeWalkerBuilder<Value, Time>;

    std::vector<std::string> variables_;
    std::vector<std::string> coefficients_;
    std::vector<CoeffFunc> coefficientValues_;
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

namespace walker {
/*
 * Contains data of general first-order scheme for linear ODE written as
 * nextMatrix * x_{t+1} = previousMatrix * x_t + source
 */
    template<class Value, class Time>
    class OdeWalker {
    private:
        struct Matrix {
            using Mat = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;

            Mat matrix;
            std::vector<std::pair<int, int>> entries;
            std::vector<std::function<Value(Time, Time)>> coefficients;

            void update(Time previousTime, Time time) {
                for (size_t i = 0; i < entries.size(); ++i) {
                    matrix(entries[i].first, entries[i].second) = coefficients[i](previousTime, time);
                }
            }
        };

    public:
        typename OdeSolver<Value, Time>::State next(const typename OdeSolver<Value, Time>::State &state, Time time) {
            nextMatrix_.update(previousTime_, time);
            previousMatrix_.update(previousTime_, time);
            source_.update(previousTime_, time);
            previousTime_ = time;

            typename Matrix::Mat rightPart = source_.matrix + previousMatrix_.matrix * Eigen::Map<const typename Matrix::Mat>(state.data(), state.size(), 1);
            typename Matrix::Mat nextMatState = nextMatrix_.matrix.colPivHouseholderQr().solve(rightPart);

            auto nextState = state;
            for (size_t i = 0; i < nextState.size(); ++i) {
                nextState[i] = nextMatState(i, 0);
            }

            return nextState;
        }

    private:
        friend class OdeWalkerBuilder<Value, Time>;

        Matrix nextMatrix_, previousMatrix_, source_;
        Time previousTime_;
    };

// Builds backward Euler OdeWalker
    template<class Value, class Time>
    class OdeWalkerBuilder {
    public:
        OdeWalker<Value, Time> build(const OdeSolver<Value, Time> &solver) {
            size_t varnum = solver.variables_.size();
            walker.previousMatrix_.matrix.setIdentity(varnum, varnum);
            walker.nextMatrix_.matrix.setZero(varnum, varnum);
            walker.source_.matrix.setZero(varnum, 1);
            walker.previousTime_ = solver.grid_[0];

            variableToIndex = createVectorToIndex(solver.variables_);
            coefficientToIndex = createVectorToIndex(solver.coefficients_);

            if (variableToIndex.size() != solver.variables_.size() ||
                coefficientToIndex.size() != solver.coefficients_.size()) {
                throw std::logic_error("names of variables/coefficients can't repeat");
            }

            for (const auto &name: solver.variables_) {
                if (coefficientToIndex.count(name)) {
                    throw std::logic_error("variable can't be coefficient simultaneously");
                }
            }

            if (varnum != solver.equations_.size()) {
                throw std::logic_error("number of equations is not equal to number of variables");
            }

            parser::Trie trie;
            for (const auto &name: solver.variables_) {
                trie.addString(name);
            }
            for (const auto &name: solver.coefficients_) {
                trie.addString(name);
            }

            for (const auto &equation: solver.equations_) {
                parseEquation(equation, trie, solver);
            }

            return walker;
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
                    } else if (word.first == Lexeme::Left) {
                        throw std::logic_error("left part of an equation has to contain only one variable");
                    }
                }
            } while (word.first != Lexeme::End);
        }

        void insertBatch(int sign, size_t row,
                         const std::vector<std::pair<parser::Lexeme, std::string>>& batch,
                         const OdeSolver<Value, Time> &solver) {
            using namespace parser;

            bool containsVariable = false;
            size_t column;
            std::vector<typename OdeSolver<Value, Time>::CoeffFunc> coefficients;
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
                    } else {
                        throw std::logic_error("unknown variable/coefficient name");
                    }
                }
            }

            if (containsVariable) {
                walker.nextMatrix_.entries.emplace_back(row, column);
                walker.nextMatrix_.coefficients.emplace_back(
                        [row, column, constant, coeffs = std::move(coefficients)](Time previous_time,
                                                                           Time time) {
                            Value shift = OdeWalkerBuilder::calculateShift(constant, coeffs, previous_time, time);
                            return (row == column ? 1 : 0) - shift;
                        });
            } else {
                walker.source_.entries.emplace_back(row, 0);
                walker.source_.coefficients.emplace_back(
                        [constant, coeffs = std::move(coefficients)](Time previous_time, Time time) {
                            return OdeWalkerBuilder::calculateShift(constant, coeffs, previous_time, time);
                        });
            }
        }

        static Value calculateShift(Value multiplier,
                                    const std::vector<typename OdeSolver<Value, Time>::CoeffFunc> &coefficients,
                                    Time previous_time, Time time) {
            Value shift = multiplier;
            for (const auto &coefficient: coefficients) {
                shift *= coefficient(time);
            }
            return (time - previous_time) * shift;
        }

    private:
        OdeWalker<Value, Time> walker;
        std::unordered_map<std::string, size_t> variableToIndex;
        std::unordered_map<std::string, size_t> coefficientToIndex;
    };

} // namespace walker
