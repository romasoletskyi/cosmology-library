#pragma once

#include <Eigen/Dense>
#include <cosmology/odesolver/parser.h>
#include <cosmology/utility/utility.h>
#include <cosmology/utility/graph.h>

template<class Value, class Time>
class OdeSolver;

namespace formula {
    template<class Value, class Time>
    class SystemBuilder;

    // Stores matrix and source (A, b)
    template<class Value, class Time>
    class System {
    private:
        struct UpdateMatrix {
            using Matrix = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;

            Matrix matrix;
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
        using Matrix = typename UpdateMatrix::Matrix;
        using UpdateMatrix = UpdateMatrix;

        const Matrix &getMatrix(Time time, const Value *state) {
            matrix_.update(time, state);
            return matrix_.matrix;
        }

        const Matrix &getSource(Time time, const Value *state) {
            source_.update(time, state);
            return source_.matrix;
        }

    private:
        friend class SystemBuilder<Value, Time>;

        UpdateMatrix matrix_;
        UpdateMatrix source_;
    };

    template<class Value, class Time>
    class Port {
    private:
        using Matrix = typename formula::System<Value, Time>::Matrix;

    public:
        Port() = default;

        Port(std::vector<int> dynamicVariables, std::vector<int> nonDynamicVariables) : dynamicVariables_(
                std::move(dynamicVariables)), nonDynamicVariables_(std::move(nonDynamicVariables)) {
        }

        void storeTotalState(const Matrix &nextDynamicState, const Matrix &nextNonDynamicState, Value *state) const {
            for (int i = 0; i < dynamicVariables_.size(); ++i) {
                state[dynamicVariables_[i]] = nextDynamicState(i, 0);
            }

            for (int i = 0; i < nonDynamicVariables_.size(); ++i) {
                state[nonDynamicVariables_[i]] = nextNonDynamicState(i, 0);
            }
        }

        Matrix loadDynamicState(const Value *state) const {
            Matrix dynamicState(dynamicVariables_.size(), 1);

            for (int i = 0; i < dynamicVariables_.size(); ++i) {
                dynamicState(i, 0) = state[dynamicVariables_[i]];
            }

            return dynamicState;
        }

        int getDynamicVariablesNumber() const {
            return dynamicVariables_.size();
        }

    private:
        std::vector<int> dynamicVariables_;
        std::vector<int> nonDynamicVariables_;
    };

    template<class Value, class Time>
    struct GeneralSystem {
        System<Value, Time> dynamic;
        System<Value, Time> nondynamic;
        Port<Value, Time> port;
    };

    template<class Value, class Time>
    class SystemBuilder {
    public:
        GeneralSystem<Value, Time> build(OdeSolver<Value, Time> &solver) {
            initialize(solver);
            makeSubstitutions(solver);

            parseEquationsWithVariables(solver, dynamicVariables, checkDynamicEquation, system_.dynamic);
            parseEquationsWithVariables(solver, nonDynamicVariables, checkNonDynamicEquation, system_.nondynamic);
            system_.port = Port<Value, Time>(std::move(dynamicVariables), std::move(nonDynamicVariables));

            return std::move(system_);
        }

    private:
        void initialize(const OdeSolver<Value, Time> &solver) {
            using namespace parser;

            variableToIndex = createVectorToIndex(solver.variables_);
            coefficientToIndex = createVectorToIndex(solver.coefficients_);
            functionToIndex = createVectorToIndex(solver.functions_);

            for (const auto &name: solver.namespace_) {
                trie.addString(name);
            }

            variableToEquation.resize(solver.variables_.size());
            for (int i = 0; i < solver.equations_.size(); ++i) {
                std::stringstream stream(solver.equations_[i]);
                ExpressionLexer lexer(trie, stream);

                std::string variable = lexer.next().second;
                int variableIndex = variableToIndex.at(variable);
                variableToEquation[variableIndex] = i;

                if (lexer.next().first != Lexeme::Derivative) {
                    nonDynamicVariables.push_back(variableIndex);
                } else {
                    dynamicVariableToIndex[std::move(variable)] = dynamicVariables.size();
                    dynamicVariables.push_back(variableIndex);
                }
            }

            system_.dynamic.matrix_.matrix.setZero(dynamicVariables.size(), dynamicVariables.size());
            system_.dynamic.source_.matrix.setZero(dynamicVariables.size(), 1);

            system_.nondynamic.matrix_.matrix.setZero(nonDynamicVariables.size(), dynamicVariables.size());
            system_.nondynamic.source_.matrix.setZero(nonDynamicVariables.size(), 1);
        }

    private:
        /*
         * Each variable has a corresponding equation according to its left part.
         * There are two types of variables - dynamic (equation starts as x' = ...) and non-dynamic (x = ...).
         * Algorithm of substitutions is:
         * 1. Assume that there exists a topological order on dependence graph of non-dynamic variables and find it
         * 2. Substitute non-dynamic variables according to this order, now all these variables can depend only on dynamic ones
         * 3. Substitute non-dynamic variables into dynamic equations
         */
        void makeSubstitutions(OdeSolver<Value, Time> &solver) const {
            using namespace graph;

            // edge goes from variable to its dependency
            auto dependencyGraph = buildVariableDependencyGraph(solver);
            auto [nonDynamicDependencyGraph, nonDynamicToIndex] = buildSubgraph(dependencyGraph, nonDynamicVariables);
            auto topologicalOrder = getTopologicalOrder(transposeGraph(nonDynamicDependencyGraph));

            for (int nonDynamicVarIndex: topologicalOrder) {
                int nonDynamicVar = nonDynamicVariables[nonDynamicVarIndex];
                for (int dependency: getVertexNeighbors(dependencyGraph, nonDynamicVar)) {
                    if (nonDynamicToIndex.count(dependency)) {
                        substituteVariableTo(solver, dependency, nonDynamicVar);
                    }
                }
            }

            for (int dynamicVar: dynamicVariables) {
                for (int dependency: getVertexNeighbors(dependencyGraph, dynamicVar)) {
                    if (nonDynamicToIndex.count(dependency)) {
                        substituteVariableTo(solver, dependency, dynamicVar);
                    }
                }
            }
        }

        graph::Graph buildVariableDependencyGraph(const OdeSolver<Value, Time> &solver) const {
            using namespace parser;

            graph::Graph graph;
            graph.edges.resize(solver.variables_.size());

            for (int i = 0; i < solver.variables_.size(); ++i) {
                std::stringstream stream(solver.equations_[variableToEquation[i]]);
                ExpressionLexer parser(trie, stream);
                std::unordered_set<int> dependentIndices;

                auto [lexeme, variable] = parser.next();
                do {
                    std::tie(lexeme, variable) = parser.next();
                    if (lexeme == Lexeme::Variable && variableToIndex.count(variable)) {
                        dependentIndices.insert(variableToIndex.at(variable));
                    }
                } while (lexeme != Lexeme::End);

                for (int dependency: dependentIndices) {
                    graph.edges[i].push_back(dependency);
                }
            }

            return graph;
        }

        void substituteVariableTo(OdeSolver<Value, Time> &solver, int from, int to) const {
            std::string equationFrom = solver.equations_[variableToEquation[from]];
            std::string variable = solver.variables_[from];
            std::string sub = equationFrom.substr(equationFrom.find('=') + 1);
            sub.erase(sub.begin(), std::find_if(sub.begin(), sub.end(), [](char c) { return !std::isspace(c); }));
            sub = "(" + sub + ")";

            std::string &equationTo = solver.equations_[variableToEquation[to]];
            for (std::string::size_type pos{};
                 std::string::npos != (pos = equationTo.find(variable.data(), pos, variable.size()));
                 pos = variable.size()) {
                equationTo.replace(pos, variable.size(), sub.data(), sub.size());
            }
        }

    private:
        void parseEquationsWithVariables(const OdeSolver<Value, Time> &solver, const std::vector<int> &variables,
                                         std::function<void(parser::ExpressionLexer &, const std::string &)> checker,
                                         System<Value, Time> &system) {
            using namespace parser;

            for (int row = 0; row < variables.size(); ++row) {
                std::stringstream stream(solver.equations_[variableToEquation[variables[row]]]);
                ExpressionLexer lexer(trie, stream);
                checker(lexer, solver.variables_[variables[row]]);

                SyntaxTree tree;
                tree.build(lexer);

                for (const auto &batch: getBatches(tree)) {
                    updateSystem(solver, batch, row, system);
                }
            }
        }

        static void checkDynamicEquation(parser::ExpressionLexer &lexer, const std::string &variable) {
            checkEquationVariable(lexer, variable);

            parser::Token token = lexer.next();
            if (token.first != parser::Lexeme::Derivative) {
                throw std::logic_error("dynamic equation doesn't contain derivative in left part");
            }

            checkEquationEqualSign(lexer);
        }

        static void checkNonDynamicEquation(parser::ExpressionLexer &lexer, const std::string &variable) {
            checkEquationVariable(lexer, variable);
            checkEquationEqualSign(lexer);
        }

        static void checkEquationVariable(parser::ExpressionLexer &lexer, const std::string &variable) {
            parser::Token token = lexer.next();
            if (token.first != parser::Lexeme::Variable || token.second != variable) {
                throw std::logic_error("equation doesn't start from corresponding variable");
            }
        }

        static void checkEquationEqualSign(parser::ExpressionLexer &lexer) {
            parser::Token token = lexer.next();
            if (token.first != parser::Lexeme::Equal) {
                throw std::logic_error("equation doesn't have ... = ... structure");
            }
        }

        struct Batch {
            std::vector<parser::Token> direct;
            std::vector<parser::Token> inverse;
            int sign = 1;
        };

        class BatchVisitor {
        private:
            using Range = std::pair<int, int>;

        public:
            void discoverVertex(parser::SyntaxTree::Node * /*vertex*/) {
            }

            void discoverBackEdge(parser::SyntaxTree::Node */*vertex*/, parser::SyntaxTree::Node */*neighbor*/) {
            }

            void leaveVertex(parser::SyntaxTree::Node *node) {
                using namespace parser;

                switch (node->token.first) {
                    case Lexeme::Variable:
                    case Lexeme::Number: {
                        batches_.push_back({{node->token}});
                        ends_.push_back(batches_.size());
                        break;
                    }
                    case Lexeme::Plus: {
                        auto [leftRange, rightRange] = getBinaryOperatorRanges();
                        ends_.push_back(rightRange.second);
                        break;
                    }
                    case Lexeme::Minus: {
                        Range range;

                        if (node->left) {
                            range = getBinaryOperatorRanges().second;
                        } else {
                            range = getUnaryOperatorRange();
                        }

                        for (int i = range.first; i < range.second; ++i) {
                            batches_[i].sign *= -1;
                        }

                        ends_.push_back(range.second);
                        break;
                    }
                    case Lexeme::Multiplication:
                        performDistributiveOperation(multiply);
                        break;
                    case Lexeme::Division:
                        performDistributiveOperation(divide);
                        break;
                    default:
                        break;
                }
            }

            std::vector<Batch> getBatches() {
                return std::move(batches_);
            }

        private:
            Range getUnaryOperatorRange() {
                int end = ends_.back();
                ends_.pop_back();

                if (ends_.empty()) {
                    return std::make_pair(0, end);
                }

                return std::make_pair(ends_.back(), end);
            }

            std::pair<Range, Range> getBinaryOperatorRanges() {
                auto rightRange = getUnaryOperatorRange();
                auto leftRange = getUnaryOperatorRange();
                return std::make_pair(leftRange, rightRange);
            }

            std::vector<Batch> popBatches(Range range) {
                std::vector<Batch> batches;
                batches.insert(batches.end(), batches_.begin() + range.first, batches_.begin() + range.second);
                batches_.erase(batches_.begin() + range.first, batches_.begin() + range.second);
                return batches;
            }

            void performDistributiveOperation(std::function<Batch(const Batch &, const Batch &)> operation) {
                auto [leftRange, rightRange] = getBinaryOperatorRanges();
                auto rightBatches = popBatches(rightRange);
                auto leftBatches = popBatches(leftRange);

                for (const auto &leftBatch: leftBatches) {
                    for (const auto &rightBatch: rightBatches) {
                        batches_.emplace_back(operation(leftBatch, rightBatch));
                    }
                }

                ends_.push_back(batches_.size());
            }

            static Batch multiply(const Batch &lhs, const Batch &rhs) {
                auto batch = lhs;
                batch.direct.insert(batch.direct.end(), rhs.direct.begin(), rhs.direct.end());
                batch.inverse.insert(batch.inverse.end(), rhs.inverse.begin(), rhs.inverse.begin());
                batch.sign *= rhs.sign;
                return batch;
            }

            static Batch divide(const Batch &lhs, const Batch &rhs) {
                auto batch = lhs;
                batch.direct.insert(batch.direct.end(), rhs.inverse.begin(), rhs.inverse.end());
                batch.inverse.insert(batch.inverse.end(), rhs.direct.begin(), rhs.direct.end());
                batch.sign /= rhs.sign;
                return batch;
            }

        private:
            std::vector<Batch> batches_;
            std::vector<int> ends_;
        };

        std::vector<Batch> getBatches(const parser::SyntaxTree &tree) {
            BatchVisitor visitor;
            graph::DepthFirstSearch(tree, tree.getRoot(), visitor);
            return visitor.getBatches();
        }

        struct ConcreteBatch {
            Value constant;
            std::vector<typename OdeSolver<Value, Time>::CoeffFunc> directCoefficients;
            std::vector<std::pair<std::vector<int>, typename OdeSolver<Value, Time>::ExplicitFunc>> directFunctions;
            std::vector<typename OdeSolver<Value, Time>::CoeffFunc> inverseCoefficients;
            std::vector<std::pair<std::vector<int>, typename OdeSolver<Value, Time>::ExplicitFunc>> inverseFunctions;
        };

        void
        updateSystem(const OdeSolver<Value, Time> &solver, const Batch &batch, int row, System<Value, Time> &system) {
            using namespace parser;

            int column = 0;
            bool containsVariable = false;
            ConcreteBatch concreteBatch{static_cast<Value>(batch.sign)};

            for (const auto &[lexeme, name]: batch.direct) {
                if (lexeme == Lexeme::Number) {
                    concreteBatch.constant *= std::stod(name);
                } else {
                    if (dynamicVariableToIndex.count(name)) {
                        if (containsVariable) {
                            throw std::logic_error("product of two dynamic variables in the equation; "
                                                   "only quasi-linear odes are supported");
                        }
                        containsVariable = true;
                        column = dynamicVariableToIndex.at(name);
                    }
                    updateCoefficientsOrFunctions(solver, name, concreteBatch.directCoefficients,
                                                  concreteBatch.directFunctions);
                }
            }

            for (const auto &[lexeme, name]: batch.inverse) {
                if (lexeme == Lexeme::Number) {
                    concreteBatch.constant /= std::stod(name);
                } else {
                    if (dynamicVariableToIndex.count(name)) {
                        throw std::logic_error("division by dynamic variable in the equation; "
                                               "only quasi-linear odes are supported");
                    }
                    updateCoefficientsOrFunctions(solver, name, concreteBatch.inverseCoefficients,
                                                  concreteBatch.inverseFunctions);
                }
            }

            if (containsVariable) {
                updateEntry(system.matrix_, row, column, std::move(concreteBatch));
            } else {
                updateEntry(system.source_, row, column, std::move(concreteBatch));
            }
        }

        void updateCoefficientsOrFunctions(const OdeSolver<Value, Time> &solver, const std::string &name,
                                           std::vector<typename OdeSolver<Value, Time>::CoeffFunc> &coefficients,
                                           std::vector<std::pair<std::vector<int>, typename OdeSolver<Value, Time>::ExplicitFunc>> &functions) {
            if (coefficientToIndex.count(name)) {
                coefficients.push_back(solver.coefficientValues_[coefficientToIndex.at(name)]);
            } else if (functionToIndex.count(name)) {
                auto [function, variableNames] = solver.functionValues_[functionToIndex.at(name)];
                std::vector<int> indices;
                indices.reserve(variableNames.size());

                for (const auto &variable: variableNames) {
                    indices.push_back(variableToIndex.at(variable));
                }

                functions.emplace_back(indices, function);
            }
        }

        void updateEntry(typename System<Value, Time>::UpdateMatrix &matrix, int row, int column,
                         ConcreteBatch concreteBatch) {
            matrix.entries.emplace_back(row, column);
            matrix.coefficients.emplace_back([batch = std::move(concreteBatch)](Time time, const Value *value) {
                Value entry = batch.constant;

                for (const auto &coefficient: batch.directCoefficients) {
                    entry *= coefficient(time);
                }
                for (const auto &coefficient: batch.inverseCoefficients) {
                    entry /= coefficient(time);
                }

                for (const auto &[indices, function]: batch.directFunctions) {
                    entry *= evaluateFunction(indices, function, time, value);
                }
                for (const auto &[indices, function]: batch.inverseFunctions) {
                    entry /= evaluateFunction(indices, function, time, value);
                }

                return entry;
            });
        }

        static Value
        evaluateFunction(const std::vector<int> &indices, const typename OdeSolver<Value, Time>::ExplicitFunc &function,
                         Time time, const Value *value) {
            std::vector<Value> values;
            values.reserve(indices.size());

            for (int index: indices) {
                values.push_back(value[index]);
            }

            return function(time, std::move(values));
        }

    private:
        GeneralSystem<Value, Time> system_;

        std::unordered_map<std::string, int> variableToIndex;
        std::unordered_map<std::string, int> dynamicVariableToIndex;
        std::unordered_map<std::string, int> coefficientToIndex;
        std::unordered_map<std::string, int> functionToIndex;

        std::vector<int> dynamicVariables;
        std::vector<int> nonDynamicVariables;
        std::vector<int> variableToEquation;

        parser::Trie trie;
    };
} // namespace formula