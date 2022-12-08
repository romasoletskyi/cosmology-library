#include <iostream>
#include <gtest/gtest.h>
#include "../include/odesolver.h"
#include "../include/utility.h"

std::pair<parser::Lexeme, std::string> pair(parser::Lexeme lexeme, std::string string) {
    return std::make_pair(lexeme, std::move(string));
}

TEST(ParserTest, Basic) {
    using namespace parser;
    std::stringstream ss("x' = x");

    Trie trie;
    trie.addString("x");

    EquationParser parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Left, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(ParserTest, Numbers) {
    using namespace parser;
    std::stringstream ss("x' = 5 x");

    Trie trie;
    trie.addString("x");

    EquationParser parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Left, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "5"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(ParserTest, Complex) {
    using namespace parser;
    std::stringstream ss("x' = - 7 a y + 3 z");

    Trie trie;
    for (const auto &string: {"x", "y", "z", "a"}) {
        trie.addString(string);
    }

    EquationParser parser(trie, ss);
    ASSERT_EQ(parser.next(), pair(Lexeme::Left, "x"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Derivative, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Minus, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "7"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "a"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "y"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Plus, ""));
    ASSERT_EQ(parser.next(), pair(Lexeme::Number, "3"));
    ASSERT_EQ(parser.next(), pair(Lexeme::Right, "z"));
    ASSERT_EQ(parser.next(), pair(Lexeme::End, ""));
}

TEST(SystemBuilder, Complex) {
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
        ASSERT_LT(std::abs(diff(0, 0)) + std::abs(diff(0, 1)) + std::abs(diff(1, 0)) + std::abs(diff(1, 1)), 1e-5) << " " << diff;
    }

    for (auto time: linSpace<float>(0, 1, 5)) {
        auto source = system.getSource(time, nullptr);
        ASSERT_LT(std::abs(source(0, 0)) + std::abs(source(1, 0)), 1e-5);
    }
}

TEST(SystemBuilder, NonLinear) {
    OdeSolver<float, float> solver;

    solver.addVariable("x");
    solver.addExplicitFunction("f", [](float time, std::vector<float> data){
        return data[0] * data[0];
    }, {"x"});
    solver.addEquation("x' = f");

    auto system = matrix::SystemBuilder<float, float>().build(solver);
    for (auto x: linSpace<float>(0, 1, 5)) {
        auto mat = system.getSource(0, &x);
        ASSERT_EQ(mat(0, 0), x * x);
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
    // ::testing::GTEST_FLAG(filter) = "SystemBuilder*NonLinear";
    return RUN_ALL_TESTS();
}
