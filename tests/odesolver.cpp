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
    for (const auto& string: {"x", "y", "z", "a"}) {
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

TEST(OdeSolver, Basic) {
    OdeSolver<float, float> solver;
    solver.addVariable("x");
    solver.addEquation("x' = x");

    size_t size = 10;
    auto grid = linSpace<float>(0, 1, size);
    float dt = grid[1] - grid[0];
    solver.setGrid(grid);

    auto solution = solver.solve({1});

    std::vector<float> diff;
    diff.reserve(size);
    for (size_t i = 0; i < solution.size(); ++i) {
        diff.push_back(std::abs(solution[i][0] - (float) std::pow(1 / (1 - dt), i)));
    }

    ASSERT_LT(*std::max_element(diff.begin(), diff.end()), 1e-3);
}

int main() {
    ::testing::InitGoogleTest();
    ::testing::GTEST_FLAG(filter) = "OdeSolver*";
    return RUN_ALL_TESTS();
}
