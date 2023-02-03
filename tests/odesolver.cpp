#include <iostream>
#include <gtest/gtest.h>
#include <cosmology/odesolver/odesolver.h>
#include <cosmology/utility/math.h>

std::pair<parser::Lexeme, std::string> pair(parser::Lexeme lexeme, std::string string) {
    return std::make_pair(lexeme, std::move(string));
}

TEST(Lexer, Basic) {
    using namespace parser;
    std::stringstream ss("x' = x");

    Trie trie;
    trie.addString("x");

    ExpressionLexer parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Equal, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(Lexer, Numbers) {
    using namespace parser;
    std::stringstream ss("x' = 5 x");

    Trie trie;
    trie.addString("x");

    ExpressionLexer parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Equal, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "5"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Multiplication, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(Lexer, Complex) {
    using namespace parser;
    std::stringstream ss("x' = - 7 a y + 3 z");

    Trie trie;
    for (const auto &string: {"x", "y", "z", "a"}) {
        trie.addString(string);
    }

    ExpressionLexer parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Equal, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Minus, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "7"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Multiplication, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "a"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Multiplication, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "y"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Plus, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "3"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Multiplication, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Variable, "z"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(Lexer, Brackets) {
    using namespace parser;
    std::stringstream ss("x' = a (y + z) a b + b a (z + a)");

    Trie trie;
    for (const auto &string: {"x", "y", "z", "a", "b"}) {
        trie.addString(string);
    }

    ExpressionLexer parser(trie, ss);
    std::vector<Token> tokens = {{Lexeme::Variable,       "x"},
                                 {Lexeme::Derivative,     ""},
                                 {Lexeme::Equal,          ""},
                                 {Lexeme::Variable,       "a"},
                                 {Lexeme::Multiplication, ""},
                                 {Lexeme::LeftBracket,    ""},
                                 {Lexeme::Variable,       "y"},
                                 {Lexeme::Plus,           ""},
                                 {Lexeme::Variable,       "z"},
                                 {Lexeme::RightBracket,   ""},
                                 {Lexeme::Multiplication, ""},
                                 {Lexeme::Variable,       "a"},
                                 {Lexeme::Multiplication, ""},
                                 {Lexeme::Variable,       "b"},
                                 {Lexeme::Plus,           ""},
                                 {Lexeme::Variable,       "b"},
                                 {Lexeme::Multiplication, ""},
                                 {Lexeme::Variable,       "a"},
                                 {Lexeme::Multiplication, ""},
                                 {Lexeme::LeftBracket,    ""},
                                 {Lexeme::Variable,       "z"},
                                 {Lexeme::Plus,           ""},
                                 {Lexeme::Variable,       "a"},
                                 {Lexeme::RightBracket,   ""}};

    for (const auto &token: tokens) {
        ASSERT_EQ(parser.next(), token);
    }
}

TEST(SyntaxTree, Basic) {
    using namespace parser;
    std::stringstream ss("7 a y + 3 z");

    Trie trie;
    for (const auto &string: {"x", "y", "z", "a"}) {
        trie.addString(string);
    }

    ExpressionLexer lexer(trie, ss);

    SyntaxTree tree;
    tree.build(lexer);
    auto root = tree.getRoot();

    ASSERT_EQ(root->token.first, Lexeme::Plus);
    ASSERT_EQ(root->left->token.first, Lexeme::Multiplication);
    ASSERT_EQ(root->left->left->token.first, Lexeme::Multiplication);
    ASSERT_EQ(root->left->right->token, pair(Lexeme::Variable, "y"));
    ASSERT_EQ(root->left->left->left->token, pair(Lexeme::Number, "7"));
    ASSERT_EQ(root->left->left->right->token, pair(Lexeme::Variable, "a"));
    ASSERT_EQ(root->right->token.first, Lexeme::Multiplication);
    ASSERT_EQ(root->right->left->token, pair(Lexeme::Number, "3"));
    ASSERT_EQ(root->right->right->token, pair(Lexeme::Variable, "z"));
}

TEST(SyntaxTree, Brackets) {
    using namespace parser;
    std::stringstream ss("((5 + a) / -8) * (4 + b)");

    Trie trie;
    for (const auto &string: {"a", "b"}) {
        trie.addString(string);
    }

    ExpressionLexer lexer(trie, ss);

    SyntaxTree tree;
    tree.build(lexer);
    auto root = tree.getRoot();

    ASSERT_EQ(root->token.first, Lexeme::Multiplication);
    ASSERT_EQ(root->level, 0);

    ASSERT_EQ(root->left->token.first, Lexeme::Division);
    ASSERT_EQ(root->left->level, 1);

    ASSERT_EQ(root->left->left->token.first, Lexeme::Plus);
    ASSERT_EQ(root->left->left->level, 2);
    ASSERT_EQ(root->left->left->left->token, pair(Lexeme::Number, "5"));
    ASSERT_EQ(root->left->left->left->level, 2);
    ASSERT_EQ(root->left->left->right->token, pair(Lexeme::Variable, "a"));
    ASSERT_EQ(root->left->left->right->level, 2);

    ASSERT_EQ(root->left->right->token.first, Lexeme::Minus);
    ASSERT_EQ(root->left->right->level, 1);
    ASSERT_EQ(root->left->right->left, nullptr);
    ASSERT_EQ(root->left->right->right->token, pair(Lexeme::Number, "8"));
    ASSERT_EQ(root->left->right->right->level, 1);

    ASSERT_EQ(root->right->token.first, Lexeme::Plus);
    ASSERT_EQ(root->right->level, 1);
    ASSERT_EQ(root->right->left->token, pair(Lexeme::Number, "4"));
    ASSERT_EQ(root->right->left->level, 1);
    ASSERT_EQ(root->right->right->token, pair(Lexeme::Variable, "b"));
    ASSERT_EQ(root->right->right->level, 1);
}

TEST(SystemBuilder, TwoEquations) {
    OdeSolver<float, float> solver;

    solver.addVariable("x");
    solver.addVariable("y");

    solver.addCoefficient("t", [](float time) {
        return time;
    });
    solver.addCoefficient("a", [](float time) {
        return 1 - time * time;
    });

    solver.addEquation("x' = -t x + y");
    solver.addEquation("y' = a x + t y");

    auto system = formula::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto mat = system.dynamic.getMatrix(time, nullptr);
        Eigen::MatrixXf diff = mat - Eigen::MatrixXf{{-time,           1},
                                                     {1 - time * time, time}};
        ASSERT_LT(std::abs(diff(0, 0)) + std::abs(diff(0, 1)) + std::abs(diff(1, 0)) + std::abs(diff(1, 1)), 1e-5)
                                    << " " << diff;
    }

    for (auto time: linSpace<float>(0, 1, 5)) {
        auto source = system.dynamic.getSource(time, nullptr);
        ASSERT_LT(std::abs(source(0, 0)) + std::abs(source(1, 0)), 1e-5);
    }
}

TEST(SystemBuilder, NonLinear) {
    OdeSolver<float, float> solver;

    solver.addVariable("x");
    solver.addExplicitFunction("f", [](float time, std::vector<float> data) {
        return data[0] * data[0];
    }, {"x"});
    solver.addEquation("x' = f");

    auto system = formula::SystemBuilder<float, float>().build(solver);
    for (auto x: linSpace<float>(0, 1, 5)) {
        auto mat = system.dynamic.getSource(0, &x);
        ASSERT_EQ(mat(0, 0), x * x);
    }
}

TEST(SystemBuilder, Brackets) {
    OdeSolver<float, float> solver;

    solver.addVariable("x");
    solver.addVariable("y");

    solver.addCoefficient("a", [](float time) { return time * time; });
    solver.addCoefficient("b", [](float time) { return time; });

    solver.addEquation("x' = x (a - b)");
    solver.addEquation("y' = b (x + y)");

    auto system = formula::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto mat = system.dynamic.getMatrix(time, nullptr);
        ASSERT_EQ(mat(0, 0), time * time - time);
        ASSERT_EQ(mat(0, 1), 0);
        ASSERT_EQ(mat(1, 0), time);
        ASSERT_EQ(mat(1, 1), time);
    }
}

TEST(SystemBuilder, Substitution) {
    OdeSolver<float, float> solver;

    solver.addVariable("x");
    solver.addVariable("y");

    solver.addCoefficient("a", [](float time) { return time * time; });
    solver.addCoefficient("b", [](float time) { return time; });

    solver.addEquation("x' = b x + a y");
    solver.addEquation("y = b x");

    auto system = formula::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto dynamicMatrix = system.dynamic.getMatrix(time, nullptr);
        ASSERT_TRUE(dynamicMatrix.cols() == 1 && dynamicMatrix.rows() == 1);
        ASSERT_EQ(dynamicMatrix(0, 0), time + time * time * time);

        auto nonDynamicMatrix = system.nondynamic.getMatrix(time, nullptr);
        ASSERT_TRUE(nonDynamicMatrix.cols() == 1 && nonDynamicMatrix.rows() == 1);
        ASSERT_EQ(nonDynamicMatrix(0, 0), time);
    }
}

TEST(SystemBuilder, Complex) {
    OdeSolver<float, float> solver;

    for (const auto &variable: {"x", "y", "z", "dx"}) {
        solver.addVariable(variable);
    }

    solver.addCoefficient("a", [](float time) { return time; });
    solver.addCoefficient("b", [](float time) { return 2; });
    solver.addCoefficient("c", [](float time) { return -8 * time * time; });
    solver.addCoefficient("d", [](float time) { return -3 * time; });

    solver.addEquation("x' = a (dx + z)");
    solver.addEquation("z' = -3 (dx - y)");
    solver.addEquation("y = a x + b z");
    solver.addEquation("dx = c x + d z");

    // after substitutions, we have
    // x' = (a c) x + a (d + 1) z
    // z' = 3 (a - c) x + 3 (b - d) z
    auto system = formula::SystemBuilder<float, float>().build(solver);

    for (auto time: linSpace<float>(0, 1, 5)) {
        auto dynamicMatrix = system.dynamic.getMatrix(time, nullptr);
        float tol = 1e-5;

        ASSERT_TRUE(dynamicMatrix.cols() == 2 && dynamicMatrix.rows() == 2);
        ASSERT_NEAR(dynamicMatrix(0, 0), -8 * time * time * time, tol);
        ASSERT_NEAR(dynamicMatrix(0, 1), time * (-3 * time + 1), tol);
        ASSERT_NEAR(dynamicMatrix(1, 0), 3 * (time + 8 * time * time), tol);
        ASSERT_NEAR(dynamicMatrix(1, 1), 3 * (2 + 3 * time), tol);

        auto nonDynamicMatrix = system.nondynamic.getMatrix(time, nullptr);

        ASSERT_TRUE(nonDynamicMatrix.cols() == 2 && nonDynamicMatrix.rows() == 2);
        ASSERT_NEAR(nonDynamicMatrix(0, 0), time, tol);
        ASSERT_NEAR(nonDynamicMatrix(0, 1), 2, tol);
        ASSERT_NEAR(nonDynamicMatrix(1, 0), -8 * time * time, tol);
        ASSERT_NEAR(nonDynamicMatrix(1, 1), -3 * time, tol);
    }
}

template<class T>
T compareSolutions(const std::vector<T> &grid, const Solution<T> &solution,
                   const std::function<std::vector<T>(T)> &precise) {
    T maxDiff = 0;
    for (int i = 0; i < solution.gridLength; ++i) {
        T diff = 0;
        auto preciseSolution = precise(grid[i]);

        for (int j = 0; j < solution.variableNumber; ++j) {
            diff += std::abs(solution(i, j) - preciseSolution[j]);
        }

        if (diff > maxDiff) {
            maxDiff = diff;
        }
    }
    return maxDiff;
}

TEST(OdeSolver, Basic) {
    OdeSolver<float, float> solver;
    solver.addVariable("x");
    solver.addEquation("x' = x");

    size_t size = 10;
    auto grid = linSpace<float>(1, 2, size);
    float dt = grid[1] - grid[0];
    solver.setGrid(grid);
    solver.compile<walker::BackwardEulerWalker<float, float>>();

    auto solution = solver.solve({1});
    auto maxDiff = compareSolutions<float>(grid, solution, [dt](float time) {
        return std::vector<float>{std::pow(1 / (1 - dt), (time - 1) / dt)};
    });

    ASSERT_LT(maxDiff, 1e-5);
}

TEST(OdeSolver, RungeKutta) {
    OdeSolver<float, float> solver;

    solver.addVariable("x");
    solver.addVariable("y");

    solver.addCoefficient("t", [](float time) {
        return time;
    });
    solver.addCoefficient("a", [](float time) {
        return 1 - time * time;
    });

    solver.addEquation("x' = -t x + y");
    solver.addEquation("y' = a x + t y");

    size_t size = 100;
    auto grid = linSpace<float>(0, 1, size);
    solver.setGrid(grid);
    solver.compile<walker::RungeKuttaWalker<float, float>>();

    auto solution = solver.solve({1, 1});
    auto maxDiff = compareSolutions<float>(grid, solution, [](float time) {
        return std::vector({1 + time, 1 + time + time * time});
    });

    ASSERT_LT(maxDiff, 1e-6);
}

TEST(OdeSolver, RungeKuttaOrder) {
    OdeSolver<double, double> solver;

    solver.addVariable("x");
    solver.addExplicitFunction("dx", [](double time, std::vector<double> data) { return data[0]; }, {"x"});
    solver.addEquation("x' = dx");
    solver.compile<walker::RungeKuttaWalker<double, double>>();

    auto smallGrid = linSpace<double>(0, 1, 4);
    auto largeGrid = linSpace<double>(0, 1, 8);

    solver.setGrid(smallGrid);
    auto smallSolution = solver.solve({1});
    auto smallGridDiff = compareSolutions<double>(smallGrid, smallSolution,
                                                  [](double time) { return std::vector<double>{std::exp(time)}; });

    auto smallGridExactDiff = compareSolutions<double>(smallGrid, smallSolution, [](double time) {
        return std::vector<double>{std::pow(7889.0 / 6144, 4 * time)};
    });
    ASSERT_LE(smallGridExactDiff, 1e-9);

    solver.setGrid(largeGrid);
    auto largeSolution = solver.solve({1});
    auto largeGridDiff = compareSolutions<double>(largeGrid, largeSolution,
                                                  [](double time) { return std::vector<double>{std::exp(time)}; });

    ASSERT_GT(smallGridDiff / largeGridDiff, 4);
}

int main() {
    ::testing::InitGoogleTest();
    ::testing::GTEST_FLAG(filter) = "OdeSolver*RungeKuttaOrder";
    return RUN_ALL_TESTS();
}
