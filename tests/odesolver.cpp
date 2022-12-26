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

    auto system = matrix::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto mat = system.getMatrix(time, nullptr);
        Eigen::MatrixXf diff = mat - Eigen::MatrixXf{{-time,           1},
                                                     {1 - time * time, time}};
        ASSERT_LT(std::abs(diff(0, 0)) + std::abs(diff(0, 1)) + std::abs(diff(1, 0)) + std::abs(diff(1, 1)), 1e-5)
                                    << " " << diff;
    }

    for (auto time: linSpace<float>(0, 1, 5)) {
        auto source = system.getSource(time, nullptr);
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

    auto system = matrix::SystemBuilder<float, float>().build(solver);
    for (auto x: linSpace<float>(0, 1, 5)) {
        auto mat = system.getSource(0, &x);
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

    auto system = matrix::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto mat = system.getMatrix(time, nullptr);
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
    solver.makeSubstitutions();

    auto system = matrix::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto mat = system.getMatrix(time, nullptr);
        ASSERT_TRUE(mat.cols() == 1 && mat.rows() == 1);
        ASSERT_EQ(mat(0, 0), time + time * time * time);
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

    solver.makeSubstitutions();
    // after substitutions, we have
    // x' = (a c) x + a (d + 1) z
    // z' = 3 (a - c) x + 3 (b - d) z
    std::vector<std::string> leftVariables = {"x", "z"};
    ASSERT_EQ(solver.getVariables(), leftVariables);

    auto system = matrix::SystemBuilder<float, float>().build(solver);
    for (auto time: linSpace<float>(0, 1, 5)) {
        auto mat = system.getMatrix(time, nullptr);
        float tol = 1e-5;

        ASSERT_TRUE(mat.cols() == 2 && mat.rows() == 2);
        ASSERT_NEAR(mat(0, 0), -8 * time * time * time, tol);
        ASSERT_NEAR(mat(0, 1), time * (-3 * time + 1), tol);
        ASSERT_NEAR(mat(1, 0), 3 * (time + 8 * time * time), tol);
        ASSERT_NEAR(mat(1, 1), 3 * (2 + 3 * time), tol);
    }
}

float compareSolutions(const std::vector<float> &grid, const OdeSolver<float, float>::Solution &solution,
                       const std::function<std::vector<float>(float)> &precise) {
    float maxDiff = 0;
    for (int i = 0; i < solution.gridLength; ++i) {
        float diff = 0;
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

    auto solution = solver.solve<walker::BackwardEulerWalker<float, float>>({1});
    float maxDiff = compareSolutions(grid, solution, [dt](float time) {
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

    auto solution = solver.solve<walker::RungeKuttaWalker<float, float>>({1, 1});
    float maxDiff = compareSolutions(grid, solution, [](float time) {
        return std::vector({1 + time, 1 + time + time * time});
    });

    ASSERT_LT(maxDiff, 1e-6);
}

int main() {
    ::testing::InitGoogleTest();
    ::testing::GTEST_FLAG(filter) = "SyntaxTree*";
    return RUN_ALL_TESTS();
}
